//
// Copyright 2026 Gabriel Lando
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Time-domain reconstruction: each radio is resampled (polyphase), translated
// to the target baseband, and combined through complementary complex crossfade
// FIR filters (which also discard the out-of-band corners).  Radio 2 is aligned
// to radio 1 (gain + phase) using a running estimate over the shared overlap.
// All hot loops are VOLK-accelerated when available, so the pipeline sustains
// real-time at the full combined rate.
//

#include "overlap_reconstructor.hpp"
#include "bonded_resample.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>
#include <numeric>
#include <stdexcept>

#ifdef UHD_HAVE_VOLK
#    include <volk/volk.h>
#endif

namespace bonded {

namespace {

using namespace detail; // cf, cd, kPi, hamming, design_lowpass_real,
                        // dot_cc, dot_cr, poly_resampler

// --- minimal radix-2 FFT (double), used only for one-time filter design ------
void design_fft(std::vector<cd>& a, bool inverse)
{
    const size_t n = a.size();
    for (size_t i = 1, j = 0; i < n; i++) {
        size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) {
            j ^= bit;
        }
        j ^= bit;
        if (i < j) {
            std::swap(a[i], a[j]);
        }
    }
    for (size_t len = 2; len <= n; len <<= 1) {
        const double ang = (inverse ? 2.0 : -2.0) * kPi / static_cast<double>(len);
        const cd wlen(std::cos(ang), std::sin(ang));
        for (size_t i = 0; i < n; i += len) {
            cd w(1.0, 0.0);
            for (size_t k = 0; k < len / 2; k++) {
                const cd u = a[i + k];
                const cd v = a[i + k + len / 2] * w;
                a[i + k]           = u + v;
                a[i + k + len / 2] = u - v;
                w *= wlen;
            }
        }
    }
    if (inverse) {
        for (auto& x : a) {
            x /= static_cast<double>(n);
        }
    }
}

// Complex FIR designed by frequency sampling of a (complex) target response,
// windowed to L taps. Returns taps in convolution order b[0..L-1].
std::vector<cf> design_complex_fir(size_t L, double fs,
    const std::function<cd(double)>& response)
{
    size_t N = 1;
    while (N < 8 * L) {
        N <<= 1;
    }
    std::vector<cd> H(N);
    const double df = fs / static_cast<double>(N);
    for (size_t k = 0; k < N; k++) {
        const double f = (k < N / 2 ? static_cast<double>(k)
                                    : static_cast<double>(k) - static_cast<double>(N))
                         * df;
        H[k] = response(f);
    }
    design_fft(H, true); // -> impulse response (periodic)
    std::vector<cf> b(L);
    const long half = static_cast<long>(L / 2);
    for (size_t i = 0; i < L; i++) {
        const long n = static_cast<long>(i) - half; // centered, may be negative
        const size_t idx = static_cast<size_t>((n % static_cast<long>(N) + static_cast<long>(N))
                                               % static_cast<long>(N));
        const cd v = H[idx] * hamming(i, L);
        b[i] = cf(static_cast<float>(v.real()), static_cast<float>(v.imag()));
    }
    return b;
}

// Streaming FIR (decim=1), complex taps, complex data.
struct cfir
{
    std::vector<cf> taps; // window order (reversed conv taps)
    size_t K = 0;
    std::vector<cf> buf;  // persistent [tail | block] working buffer

    void init(std::vector<cf> conv_taps)
    {
        K = conv_taps.size();
        std::reverse(conv_taps.begin(), conv_taps.end());
        taps = std::move(conv_taps);
        buf.assign(K - 1, cf(0, 0)); // initial tail of zeros
    }

    void filter(const cf* x, size_t n, std::vector<cf>& out)
    {
        buf.resize(K - 1 + n);
        std::copy(x, x + n, buf.begin() + static_cast<std::ptrdiff_t>(K - 1));
        const size_t base = out.size();
        out.resize(base + n);
        for (size_t i = 0; i < n; i++) {
            out[base + i] = dot_cc(buf.data() + i, taps.data(), K);
        }
        std::copy(buf.end() - static_cast<std::ptrdiff_t>(K - 1), buf.end(),
            buf.begin()); // carry tail
        buf.resize(K - 1);
    }
};

// Streaming FIR (decim=1), real taps, complex data (for overlap analysis).
struct rfir
{
    std::vector<float> taps; // window order
    size_t K = 0;
    std::vector<cf> buf;     // persistent [tail | block] working buffer

    void init(std::vector<float> conv_taps)
    {
        K = conv_taps.size();
        std::reverse(conv_taps.begin(), conv_taps.end());
        taps = std::move(conv_taps);
        buf.assign(K - 1, cf(0, 0));
    }

    void filter(const cf* x, size_t n, std::vector<cf>& out)
    {
        buf.resize(K - 1 + n);
        std::copy(x, x + n, buf.begin() + static_cast<std::ptrdiff_t>(K - 1));
        const size_t base = out.size();
        out.resize(base + n);
        for (size_t i = 0; i < n; i++) {
            out[base + i] = dot_cr(buf.data() + i, taps.data(), K);
        }
        std::copy(buf.end() - static_cast<std::ptrdiff_t>(K - 1), buf.end(),
            buf.begin());
        buf.resize(K - 1);
    }
};

} // namespace

struct overlap_reconstructor::impl
{
    overlap_reconstructor_config cfg;
    size_t interp, decim;
    cf inc1, inc2, rot1{1, 0}, rot2{1, 0};

    poly_resampler rs1, rs2;
    cfir hlow, hhigh;
    rfir ov1, ov2; // overlap-band lowpass for alignment

    // Running alignment statistics over the overlap.
    double Saa = 0.0, Sbb = 0.0;
    cd Sxy{0, 0};
    cf alpha{1, 0};
    overlap_alignment align;
    // Alignment is static (shared clock); estimate during a warmup window then
    // freeze it and skip the overlap-analysis filters (a real-time saving).
    bool frozen = false;
    size_t processed = 0;
    static constexpr size_t kWarmup = 200000;

    // Reusable scratch buffers (retain capacity across calls; avoid per-call
    // heap churn, which otherwise dominates at small GNU Radio chunk sizes).
    std::vector<cf> r1_, r2_, o1_, o2_, hi_, out_;

    explicit impl(const overlap_reconstructor_config& c) : cfg(c)
    {
        const size_t g =
            std::gcd(static_cast<size_t>(std::llround(c.target_rate)),
                static_cast<size_t>(std::llround(c.per_radio_rate)));
        interp = static_cast<size_t>(std::llround(c.target_rate)) / g;
        decim  = static_cast<size_t>(std::llround(c.per_radio_rate)) / g;
        rs1.init(interp, decim);
        rs2.init(interp, decim);

        const double fs        = c.target_rate;
        const double half_ov   = 0.5 * c.overlap_bw;
        const double half_band = 0.5 * c.target_bw;
        // Complementary complex crossfade filters (target baseband). The
        // transition band equals the overlap, so the FIR length scales with
        // target/overlap: a wide overlap -> short, cheap filters; a narrow
        // overlap -> a sharper, longer (costlier) filter. Bounded for safety.
        const double frac = c.overlap_bw / c.target_rate; // normalized transition
        auto odd = [](long v) { return static_cast<size_t>(v | 1L); };
        const size_t L = odd(std::clamp<long>(
            std::lround(8.0 / std::max(frac, 1e-6)), 17, 257));
        hlow.init(design_complex_fir(L, fs, [=](double f) -> cd {
            if (std::abs(f) > half_band) return cd(0, 0);
            if (f <= -half_ov) return cd(1, 0);
            if (f >= half_ov)  return cd(0, 0);
            const double x = (f + half_ov) / (2.0 * half_ov);
            const double c = std::cos(0.5 * kPi * x);
            return cd(c * c, 0);
        }));
        hhigh.init(design_complex_fir(L, fs, [=](double f) -> cd {
            if (std::abs(f) > half_band) return cd(0, 0);
            if (f <= -half_ov) return cd(0, 0);
            if (f >= half_ov)  return cd(1, 0);
            const double x = (f + half_ov) / (2.0 * half_ov);
            const double s = std::sin(0.5 * kPi * x);
            return cd(s * s, 0);
        }));
        // Overlap-band lowpass (|f| < half_ov) for the alignment estimate.
        const size_t Lo = odd(std::clamp<long>(
            std::lround(6.0 / std::max(frac, 1e-6)), 17, 195));
        ov1.init(design_lowpass_real(Lo, half_ov / fs, 1.0));
        ov2.init(design_lowpass_real(Lo, half_ov / fs, 1.0));

        const double w1 = 2.0 * kPi * c.radio1_offset_hz / fs;
        const double w2 = 2.0 * kPi * c.radio2_offset_hz / fs;
        inc1 = cf(static_cast<float>(std::cos(w1)), static_cast<float>(std::sin(w1)));
        inc2 = cf(static_cast<float>(std::cos(w2)), static_cast<float>(std::sin(w2)));
    }

    void rotate(std::vector<cf>& x, cf& phase, cf inc)
    {
#ifdef UHD_HAVE_VOLK
        volk_32fc_s32fc_x2_rotator_32fc(reinterpret_cast<lv_32fc_t*>(x.data()),
            reinterpret_cast<const lv_32fc_t*>(x.data()),
            *reinterpret_cast<const lv_32fc_t*>(&inc),
            reinterpret_cast<lv_32fc_t*>(&phase),
            static_cast<unsigned>(x.size()));
#else
        for (auto& v : x) {
            v *= phase;
            phase *= inc;
        }
        const float mag = std::abs(phase);
        if (mag > 0) {
            phase /= mag;
        }
#endif
    }

    std::vector<cf> process(const cf* rx1, const cf* rx2, size_t n)
    {
        // 1) Resample 15->20 and translate to the target baseband.
        r1_.clear();
        r2_.clear();
        rs1.feed(rx1, n, r1_);
        rs2.feed(rx2, n, r2_);
        const size_t m = std::min(r1_.size(), r2_.size());
        r1_.resize(m);
        r2_.resize(m);
        rotate(r1_, rot1, inc1);
        rotate(r2_, rot2, inc2);

        // 2) Alignment estimate over the shared overlap band (warmup only).
        if (!frozen) {
            o1_.clear();
            o2_.clear();
            ov1.filter(r1_.data(), m, o1_);
            ov2.filter(r2_.data(), m, o2_);
            for (size_t i = 0; i < m; i++) {
                Saa += std::norm(o1_[i]);
                Sbb += std::norm(o2_[i]);
                Sxy += static_cast<cd>(o1_[i]) * std::conj(static_cast<cd>(o2_[i]));
            }
            if (Sbb > 0.0) {
                const double gain = std::sqrt(Saa / Sbb);
                const double ph   = std::arg(Sxy);
                alpha = cf(static_cast<float>(gain * std::cos(ph)),
                           static_cast<float>(gain * std::sin(ph)));
                align.gain        = gain;
                align.phase_rad   = ph;
                align.delay_samps = 0.0;
                align.valid       = true;
            }
            processed += m;
            if (processed >= kWarmup && Sbb > 0.0) {
                frozen = true; // hold the (static) alignment, skip overlap FIRs
            }
        }

        // 3) Align radio 2 and combine through complementary crossfade filters.
        for (size_t i = 0; i < m; i++) {
            r2_[i] *= alpha;
        }
        out_.clear();
        hi_.clear();
        hlow.filter(r1_.data(), m, out_);  // out = Hlow(r1)
        hhigh.filter(r2_.data(), m, hi_);  // hi  = Hhigh(r2_aligned)
        for (size_t i = 0; i < m; i++) {
            out_[i] += hi_[i];
        }
        return out_;
    }
};

overlap_reconstructor::overlap_reconstructor(const overlap_reconstructor_config& cfg)
    : _cfg(cfg)
{
    if (cfg.target_rate <= 0 || cfg.per_radio_rate <= 0) {
        throw std::invalid_argument("sample rates must be > 0");
    }
}

std::vector<std::complex<float>> overlap_reconstructor::reconstruct(
    const std::vector<std::complex<float>>& rx1,
    const std::vector<std::complex<float>>& rx2)
{
    reset();
    const size_t n = std::min(rx1.size(), rx2.size());
    auto out = _impl->process(rx1.data(), rx2.data(), n);
    _alignment = _impl->align;
    return out;
}

std::vector<std::complex<float>> overlap_reconstructor::process(
    const std::complex<float>* rx1, const std::complex<float>* rx2, size_t n)
{
    if (!_impl) {
        _impl = std::make_shared<impl>(_cfg);
    }
    auto out = _impl->process(rx1, rx2, n);
    _alignment = _impl->align;
    return out;
}

void overlap_reconstructor::reset()
{
    _impl = std::make_shared<impl>(_cfg);
}

// --- one-shot polyphase resampler (kept for tests / offline helpers) ---------
std::vector<std::complex<float>> rational_resample(
    const std::vector<std::complex<float>>& x, size_t interp, size_t decim)
{
    if (interp == 0 || decim == 0) {
        throw std::invalid_argument("rational_resample: interp/decim must be > 0");
    }
    if (x.empty() || interp == decim) {
        return x;
    }
    const size_t maxLM    = std::max(interp, decim);
    const size_t num_taps = 20 * maxLM + 1;
    const std::vector<float> h =
        design_lowpass_real(num_taps, 0.5 / static_cast<double>(maxLM),
            static_cast<double>(interp));
    const size_t Nx     = x.size();
    const long full_len = static_cast<long>((Nx - 1) * interp + num_taps);
    const size_t Ny     = static_cast<size_t>((full_len + decim - 1) / decim);
    std::vector<std::complex<float>> y(Ny);
    for (size_t nn = 0; nn < Ny; nn++) {
        const long m = static_cast<long>(nn) * static_cast<long>(decim);
        cd acc(0.0, 0.0);
        long t = m % static_cast<long>(interp);
        for (; t < static_cast<long>(num_taps); t += static_cast<long>(interp)) {
            const long p = (m - t) / static_cast<long>(interp);
            if (p >= 0 && p < static_cast<long>(Nx)) {
                acc += static_cast<cd>(x[static_cast<size_t>(p)])
                       * static_cast<double>(h[static_cast<size_t>(t)]);
            }
        }
        y[nn] = cf(static_cast<float>(acc.real()), static_cast<float>(acc.imag()));
    }
    return y;
}

} // namespace bonded
