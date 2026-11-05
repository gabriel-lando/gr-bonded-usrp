//
// Copyright 2026 Gabriel Lando
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Polyphase resamplers, dot-product kernels and windowed-sinc filter design
// shared by overlap_reconstructor.cpp and band_splitter.cpp (header-only, in
// bonded::detail).
//

#pragma once

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <vector>

#ifdef BONDED_HAVE_VOLK
#    include <volk/volk.h>
#endif

namespace bonded { namespace detail {

constexpr double kPi = 3.14159265358979323846;
using cf             = std::complex<float>;
using cd             = std::complex<double>;

// ---------------------------------------------------------------------------
// Windowing
// ---------------------------------------------------------------------------

inline double hamming(size_t i, size_t L)
{
    return 0.54
           - 0.46
                 * std::cos(2.0 * kPi * static_cast<double>(i)
                            / static_cast<double>(L - 1));
}

/// Real windowed-sinc lowpass; DC gain = scale.
inline std::vector<float> design_lowpass_real(size_t num_taps,
    double cutoff_norm,
    double scale)
{
    const double mid = 0.5 * static_cast<double>(num_taps - 1);
    double sum       = 0.0;
    std::vector<double> hd(num_taps);
    for (size_t i = 0; i < num_taps; i++) {
        const double m    = static_cast<double>(i) - mid;
        const double sinc = (std::abs(m) < 1e-9)
                                ? 2.0 * cutoff_norm
                                : std::sin(2.0 * kPi * cutoff_norm * m) / (kPi * m);
        hd[i] = sinc * hamming(i, num_taps);
        sum += hd[i];
    }
    std::vector<float> h(num_taps);
    for (size_t i = 0; i < num_taps; i++) {
        h[i] = static_cast<float>(hd[i] * scale / sum);
    }
    return h;
}

// ---------------------------------------------------------------------------
// VOLK dot product (scalar fallback without VOLK), for one-shot use
// ---------------------------------------------------------------------------

inline cf dot_cr(const cf* a, const float* taps, size_t K)
{
#ifdef BONDED_HAVE_VOLK
    cf r;
    volk_32fc_32f_dot_prod_32fc(reinterpret_cast<lv_32fc_t*>(&r),
        reinterpret_cast<const lv_32fc_t*>(a),
        taps,
        static_cast<unsigned>(K));
    return r;
#else
    cf r(0, 0);
    for (size_t k = 0; k < K; k++) {
        r += a[k] * taps[k];
    }
    return r;
#endif
}

// ---------------------------------------------------------------------------
// Inlined dot products for the per-output-sample hot loops.
//
// VOLK dispatches through a function pointer, which cannot be inlined; for the
// short filters used here (K ~ 13..33) that per-call overhead dominates the
// math. These inline loops are vectorized by the compiler instead.
// ---------------------------------------------------------------------------

inline cf dot_cr_hot(const cf* a, const float* taps, size_t K)
{
    const float* af = reinterpret_cast<const float*>(a); // [re,im,re,im,...]
    float acc_re = 0.0f, acc_im = 0.0f;
    for (size_t k = 0; k < K; k++) {
        acc_re += af[2 * k] * taps[k];
        acc_im += af[2 * k + 1] * taps[k];
    }
    return cf(acc_re, acc_im);
}

inline cf dot_cc_hot(const cf* a, const cf* taps, size_t K)
{
    const float* af = reinterpret_cast<const float*>(a);
    const float* tf = reinterpret_cast<const float*>(taps);
    float acc_re = 0.0f, acc_im = 0.0f;
    for (size_t k = 0; k < K; k++) {
        const float ar = af[2 * k], ai = af[2 * k + 1];
        const float tr = tf[2 * k], ti = tf[2 * k + 1];
        acc_re += ar * tr - ai * ti;
        acc_im += ar * ti + ai * tr;
    }
    return cf(acc_re, acc_im);
}

//! Same product, but over de-interleaved (split re/im) data and taps.
//
// Interleaved complex needs a shuffle per iteration to separate re from im;
// with four contiguous float streams the loop is pure FMA (about 1.9x faster
// at K=33). This is the hot loop of the fused reconstruction path.
inline cf dot_cc_split(const float* ar_p,
    const float* ai_p,
    const float* tr_p,
    const float* ti_p,
    size_t K)
{
    float acc_re = 0.0f, acc_im = 0.0f;
    for (size_t k = 0; k < K; k++) {
        acc_re += ar_p[k] * tr_p[k] - ai_p[k] * ti_p[k];
        acc_im += ar_p[k] * ti_p[k] + ai_p[k] * tr_p[k];
    }
    return cf(acc_re, acc_im);
}

// ---------------------------------------------------------------------------
// Streaming polyphase rational resampler (real prototype, complex data)
// ---------------------------------------------------------------------------

struct poly_resampler
{
    size_t interp = 1, decim = 1, kmax = 0;
    // Polyphase taps in one contiguous buffer, phase p at [p*kmax, p*kmax+kmax),
    // oldest-first (reversed), left-padded with zeros so every phase has
    // exactly kmax taps.
    std::vector<float> poly_flat; ///< L*kmax floats
    std::vector<cf> in;
    size_t head  = 0;
    long in_base = 0, out_n = 0;

    /// taps_per_phase sets the prototype length (taps_per_phase*max(L,M)+1)
    /// and so the transition width; 12 is what the RX reconstructor uses.
    /// cutoff_norm is the -6 dB point relative to the upsampled rate L*fs_in;
    /// <= 0 means the anti-alias default 0.5/max(L,M).
    void init(size_t L, size_t M, size_t taps_per_phase = 12, double cutoff_norm = 0.0)
    {
        interp             = L;
        decim              = M;
        const size_t maxLM = std::max(L, M);
        const size_t ntaps = taps_per_phase * maxLM + 1;
        if (cutoff_norm <= 0.0) {
            cutoff_norm = 0.5 / static_cast<double>(maxLM);
        }
        const std::vector<float> h =
            design_lowpass_real(ntaps, cutoff_norm, static_cast<double>(L));
        // First pass: build each phase's reversed (oldest-first) taps and kmax.
        std::vector<std::vector<float>> phases(L);
        kmax = 0;
        for (size_t phase = 0; phase < L; phase++) {
            std::vector<float> fwd;
            for (size_t t = phase; t < ntaps; t += L) {
                fwd.push_back(h[t]);
            }
            std::reverse(fwd.begin(), fwd.end()); // oldest-first for the dot
            phases[phase] = std::move(fwd);
            kmax          = std::max(kmax, phases[phase].size());
        }
        // Second pass: pack into the flat buffer, left-padded with zeros. The
        // leading zeros multiply older samples that contribute nothing, so the
        // padded dot over kmax taps is identical to the ragged dot.
        poly_flat.assign(L * kmax, 0.0f);
        for (size_t phase = 0; phase < L; phase++) {
            const std::vector<float>& fwd = phases[phase];
            float* dst = poly_flat.data() + phase * kmax + (kmax - fwd.size());
            std::copy(fwd.begin(), fwd.end(), dst);
        }
        in.clear();
        head    = 0;
        in_base = 0;
        out_n   = 0;
    }

    void feed(const cf* x, size_t n, std::vector<cf>& out)
    {
        in.insert(in.end(), x, x + n);
        out.reserve(out.size() + n * interp / decim + 2);
        const long last = in_base + static_cast<long>(in.size() - head) - 1;
        const long start_out_n = out_n;
        // Largest out_n with floor(out_n*decim/interp) <= last, in closed form:
        // floor(a/b) <= last  <=>  a <= (last+1)*interp - 1.
        long end_out_n = start_out_n;
        if (last >= 0) {
            const long L = static_cast<long>(interp);
            const long M = static_cast<long>(decim);
            const long max_out_n = ((last + 1) * L - 1) / M; // integer floor, >=0
            end_out_n = std::max(start_out_n, max_out_n + 1);
        }

        const size_t produced = static_cast<size_t>(end_out_n - start_out_n);
        const size_t base = out.size();
        out.resize(base + produced);

        const size_t K        = kmax; // constant tap count (zero-padded phases)
        const float* poly_ptr = poly_flat.data();
        const cf* in_ptr      = in.data() + head;
        const size_t in_avail = in.size() - head;
        cf* out_ptr           = out.data() + base;
        for (size_t idx = 0; idx < produced; idx++) {
            const long out_idx = start_out_n + static_cast<long>(idx);
            const long m       = out_idx * static_cast<long>(decim);
            const long pin     = m / static_cast<long>(interp);
            const float* pp =
                poly_ptr + static_cast<size_t>(m % static_cast<long>(interp)) * K;
            const long oldest = pin - static_cast<long>(K - 1);

            if (oldest >= in_base) {
                out_ptr[idx] = dot_cr_hot(
                    in_ptr + static_cast<size_t>(oldest - in_base), pp, K);
            } else {
                cf acc(0, 0);
                for (size_t k = 0; k < K; k++) {
                    const long off = oldest + static_cast<long>(k) - in_base;
                    if (off >= 0 && static_cast<size_t>(off) < in_avail) {
                        acc += pp[k] * in_ptr[static_cast<size_t>(off)];
                    }
                }
                out_ptr[idx] = acc;
            }
        }

        out_n = end_out_n;
        const long next_pin =
            (out_n * static_cast<long>(decim)) / static_cast<long>(interp);
        const long need_from = next_pin - static_cast<long>(kmax);
        while (in_base < need_from && head < in.size()) {
            head++;
            in_base++;
        }
        if (head > 65536) {
            in.erase(in.begin(),
                in.begin() + static_cast<std::ptrdiff_t>(head));
            head = 0;
        }
    }
};

// ---------------------------------------------------------------------------
// Streaming polyphase rational resampler with a COMPLEX prototype.
//
// Same scheduling as poly_resampler, but the prototype taps are complex and
// supplied by the caller. The fused reconstruction path uses one per radio,
// whose prototype embeds the resampler anti-alias, the frequency translation
// and the crossfade band-shape, so one filter pass replaces all three.
// ---------------------------------------------------------------------------
struct poly_resampler_c
{
    size_t interp = 1, decim = 1, kmax = 0;
    // Taps and history are both stored de-interleaved (all real parts, then all
    // imaginary parts) so the hot dot product runs over four contiguous float
    // streams — see dot_cc_split().
    std::vector<float> poly_r, poly_i; ///< L*kmax each, oldest-first, left-padded
    std::vector<float> in_r, in_i;
    size_t head  = 0;
    long in_base = 0, out_n = 0;
    long pending_end_out_n = 0; ///< set by prepare(), applied by commit()

    /// proto: prototype impulse response at the upsampled (interp) grid, in
    /// convolution order (same convention as poly_resampler's internal design).
    void init(size_t L, size_t M, const std::vector<cf>& proto)
    {
        interp = L;
        decim  = M;
        const size_t ntaps = proto.size();
        std::vector<std::vector<cf>> phases(L);
        kmax = 0;
        for (size_t phase = 0; phase < L; phase++) {
            std::vector<cf> fwd;
            for (size_t t = phase; t < ntaps; t += L) {
                fwd.push_back(proto[t]);
            }
            std::reverse(fwd.begin(), fwd.end());
            phases[phase] = std::move(fwd);
            kmax          = std::max(kmax, phases[phase].size());
        }
        poly_r.assign(L * kmax, 0.0f);
        poly_i.assign(L * kmax, 0.0f);
        for (size_t phase = 0; phase < L; phase++) {
            const std::vector<cf>& fwd = phases[phase];
            const size_t off = phase * kmax + (kmax - fwd.size());
            for (size_t k = 0; k < fwd.size(); k++) {
                poly_r[off + k] = fwd[k].real();
                poly_i[off + k] = fwd[k].imag();
            }
        }
        in_r.clear();
        in_i.clear();
        head    = 0;
        in_base = 0;
        out_n   = 0;
        pending_end_out_n = 0;
    }

    // Filtering in three steps, so the output loop can be split across
    // threads: prepare() appends the input and returns how many outputs this
    // batch yields; compute() writes outputs [i0, i1) of the batch (read-only,
    // so disjoint ranges may run concurrently); commit() advances the state.
    size_t prepare(const cf* x, size_t n)
    {
        const size_t prev = in_r.size();
        in_r.resize(prev + n);
        in_i.resize(prev + n);
        for (size_t i = 0; i < n; i++) {
            in_r[prev + i] = x[i].real();
            in_i[prev + i] = x[i].imag();
        }
        const long last = in_base + static_cast<long>(in_r.size() - head) - 1;
        const long start_out_n = out_n;
        long end_out_n = start_out_n;
        if (last >= 0) {
            const long L = static_cast<long>(interp);
            const long M = static_cast<long>(decim);
            const long max_out_n = ((last + 1) * L - 1) / M;
            end_out_n = std::max(start_out_n, max_out_n + 1);
        }
        pending_end_out_n = end_out_n;
        return static_cast<size_t>(end_out_n - start_out_n);
    }

    void compute(cf* out_ptr, size_t i0, size_t i1) const
    {
        const long start_out_n = out_n;
        const size_t K        = kmax;
        const float* pr_base  = poly_r.data();
        const float* pi_base  = poly_i.data();
        const float* inr      = in_r.data() + head;
        const float* ini      = in_i.data() + head;
        const size_t in_avail = in_r.size() - head;
        for (size_t idx = i0; idx < i1; idx++) {
            const long out_idx = start_out_n + static_cast<long>(idx);
            const long m       = out_idx * static_cast<long>(decim);
            const long pin     = m / static_cast<long>(interp);
            const size_t poff =
                static_cast<size_t>(m % static_cast<long>(interp)) * K;
            const long oldest = pin - static_cast<long>(K - 1);

            if (oldest >= in_base) {
                const size_t s = static_cast<size_t>(oldest - in_base);
                out_ptr[idx] =
                    dot_cc_split(inr + s, ini + s, pr_base + poff, pi_base + poff, K);
            } else {
                float acc_re = 0.0f, acc_im = 0.0f;
                for (size_t k = 0; k < K; k++) {
                    const long off = oldest + static_cast<long>(k) - in_base;
                    if (off >= 0 && static_cast<size_t>(off) < in_avail) {
                        const size_t o = static_cast<size_t>(off);
                        const float tr = pr_base[poff + k], ti = pi_base[poff + k];
                        acc_re += inr[o] * tr - ini[o] * ti;
                        acc_im += inr[o] * ti + ini[o] * tr;
                    }
                }
                out_ptr[idx] = cf(acc_re, acc_im);
            }
        }
    }

    void commit()
    {
        out_n = pending_end_out_n;
        const long next_pin =
            (out_n * static_cast<long>(decim)) / static_cast<long>(interp);
        const long need_from = next_pin - static_cast<long>(kmax);
        while (in_base < need_from && head < in_r.size()) {
            head++;
            in_base++;
        }
        if (head > 65536) {
            const std::ptrdiff_t h = static_cast<std::ptrdiff_t>(head);
            in_r.erase(in_r.begin(), in_r.begin() + h);
            in_i.erase(in_i.begin(), in_i.begin() + h);
            head = 0;
        }
    }
};

}} // namespace bonded::detail
