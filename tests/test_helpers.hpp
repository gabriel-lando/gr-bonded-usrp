//
// Copyright 2026 Gabriel Lando
//
// Signal helpers shared by the unit tests.
//

#pragma once

#include "bonded_resample.hpp"
#include <cmath>
#include <complex>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace bonded_test {

using cf = std::complex<float>;
using cd = std::complex<double>;
constexpr double kPi = 3.14159265358979323846;

//! One-shot polyphase resampler by interp/decim (reference implementation,
//! independent of the streaming resamplers under test).
inline std::vector<cf> rational_resample(const std::vector<cf>& x, size_t interp, size_t decim)
{
    if (interp == 0 || decim == 0) {
        throw std::invalid_argument("rational_resample: interp/decim must be > 0");
    }
    if (x.empty() || interp == decim) {
        return x;
    }
    const size_t maxLM    = std::max(interp, decim);
    const size_t num_taps = 20 * maxLM + 1;
    const std::vector<float> h = bonded::detail::design_lowpass_real(
        num_taps, 0.5 / static_cast<double>(maxLM), static_cast<double>(interp));
    const size_t Nx     = x.size();
    const long full_len = static_cast<long>((Nx - 1) * interp + num_taps);
    const size_t Ny     = static_cast<size_t>((full_len + decim - 1) / decim);
    std::vector<cf> y(Ny);
    for (size_t nn = 0; nn < Ny; nn++) {
        const long m = static_cast<long>(nn) * static_cast<long>(decim);
        cd acc(0.0, 0.0);
        for (long t = m % static_cast<long>(interp); t < static_cast<long>(num_taps);
             t += static_cast<long>(interp)) {
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

//! Complex tone exp(j(2*pi*f*n/fs + phase)) of n samples.
inline std::vector<cf> tone(double f, double fs, size_t n, double phase = 0.0, double amp = 1.0)
{
    std::vector<cf> x(n);
    for (size_t i = 0; i < n; i++) {
        const double a = 2.0 * kPi * f * static_cast<double>(i) / fs + phase;
        x[i] = cf(static_cast<float>(amp * std::cos(a)), static_cast<float>(amp * std::sin(a)));
    }
    return x;
}

//! Least-squares complex amplitude of a tone at f over x[begin, end).
inline cd tone_amplitude(const std::vector<cf>& x, double f, double fs, size_t begin, size_t end)
{
    cd acc(0, 0);
    for (size_t i = begin; i < end; i++) {
        const double a = -2.0 * kPi * f * static_cast<double>(i) / fs;
        acc += static_cast<cd>(x[i]) * cd(std::cos(a), std::sin(a));
    }
    return acc / static_cast<double>(end - begin);
}

//! Mean power of x[begin, end).
inline double mean_power(const std::vector<cf>& x, size_t begin, size_t end)
{
    double p = 0.0;
    for (size_t i = begin; i < end; i++) {
        p += std::norm(x[i]);
    }
    return p / static_cast<double>(end - begin);
}

} // namespace bonded_test
