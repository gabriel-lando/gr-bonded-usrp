//
// Copyright 2026 Gabriel Lando
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Shared polyphase-resampler, VOLK dot-product helpers, and windowed-sinc
// filter-design primitives used by overlap_reconstructor.cpp (and vendored by
// the band_splitter dev block in the sibling test repo).
//
// All symbols live in bonded::detail and are marked inline so
// they can be defined in multiple translation units without ODR issues.
// Use `using namespace detail;` inside the calling anonymous namespace to
// get the short names (cf, kPi, poly_resampler, …).
//

#pragma once

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <vector>

#ifdef UHD_HAVE_VOLK
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
// VOLK-accelerated dot products (scalar fallback when VOLK not available)
// ---------------------------------------------------------------------------

inline cf dot_cc(const cf* a, const cf* taps, size_t K)
{
#ifdef UHD_HAVE_VOLK
    cf r;
    volk_32fc_x2_dot_prod_32fc(reinterpret_cast<lv_32fc_t*>(&r),
        reinterpret_cast<const lv_32fc_t*>(a),
        reinterpret_cast<const lv_32fc_t*>(taps),
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

inline cf dot_cr(const cf* a, const float* taps, size_t K)
{
#ifdef UHD_HAVE_VOLK
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
// Streaming polyphase rational resampler
// (real prototype filter, complex data, VOLK dot-product per output sample)
// ---------------------------------------------------------------------------

struct poly_resampler
{
    size_t interp = 1, decim = 1, kmax = 0;
    std::vector<std::vector<float>> poly; ///< poly[phase], oldest-first (reversed)
    std::vector<cf> in;
    size_t head  = 0;
    long in_base = 0, out_n = 0;

    void init(size_t L, size_t M)
    {
        interp             = L;
        decim              = M;
        const size_t maxLM = std::max(L, M);
        const size_t ntaps = 12 * maxLM + 1;
        const std::vector<float> h =
            design_lowpass_real(ntaps, 0.5 / static_cast<double>(maxLM),
                static_cast<double>(L));
        poly.assign(L, {});
        kmax = 0;
        for (size_t phase = 0; phase < L; phase++) {
            std::vector<float> fwd;
            for (size_t t = phase; t < ntaps; t += L) {
                fwd.push_back(h[t]);
            }
            std::reverse(fwd.begin(), fwd.end()); // oldest-first for the dot
            poly[phase] = std::move(fwd);
            kmax        = std::max(kmax, poly[phase].size());
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
        while (true) {
            const long m   = out_n * static_cast<long>(decim);
            const long pin = m / static_cast<long>(interp);
            if (pin > last) {
                break;
            }
            const std::vector<float>& pp =
                poly[static_cast<size_t>(m % static_cast<long>(interp))];
            const size_t K      = pp.size();
            const long   oldest = pin - static_cast<long>(K - 1);
            if (oldest >= in_base) {
                const cf* x0 =
                    in.data() + head + static_cast<size_t>(oldest - in_base);
                out.push_back(dot_cr(x0, pp.data(), K));
            } else {
                cf acc(0, 0);
                for (size_t k = 0; k < K; k++) {
                    const long off = oldest + static_cast<long>(k) - in_base;
                    if (off >= 0
                        && static_cast<size_t>(off) < in.size() - head) {
                        acc += pp[k] * in[head + static_cast<size_t>(off)];
                    }
                }
                out.push_back(acc);
            }
            out_n++;
            const long next_pin =
                (out_n * static_cast<long>(decim)) / static_cast<long>(interp);
            const long need_from = next_pin - static_cast<long>(kmax);
            while (in_base < need_from && head < in.size()) {
                head++;
                in_base++;
            }
        }
        if (head > 65536) {
            in.erase(in.begin(),
                in.begin() + static_cast<std::ptrdiff_t>(head));
            head = 0;
        }
    }
};

}} // namespace bonded::detail
