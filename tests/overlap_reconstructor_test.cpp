//
// Copyright 2026 Gabriel Lando
// SPDX-License-Identifier: GPL-3.0-or-later
//
// End-to-end (no hardware) validation of the overlap-aware reconstruction.
//
// A known wideband signal at the target rate is split — in software — into two
// overlapping 15 MHz slices (the inverse of what the bonded radios do), a known
// gain/phase offset is injected on the second slice, and the reconstructor is
// asked to rebuild the original 20 MHz baseband.  We then check that
//   (a) the alignment estimate recovers the injected gain/phase,
//   (b) the reconstruction matches the original (low residual) across the band,
//       including the overlap seam.
//

#include "overlap_reconstructor.hpp"
#include <boost/test/unit_test.hpp>
#include <cmath>
#include <complex>
#include <random>
#include <vector>

using namespace bonded;

namespace {

constexpr double kPi = 3.14159265358979323846;
using cf = std::complex<float>;

// Multi-tone wideband test signal at `rate`, tones given in Hz (baseband).
std::vector<cf> make_wideband(size_t n, double rate,
    const std::vector<double>& tones_hz)
{
    std::mt19937 rng(12345);
    std::uniform_real_distribution<double> ph(-kPi, kPi);
    std::uniform_real_distribution<double> am(0.4, 1.0);
    std::vector<double> phase, amp;
    for (size_t t = 0; t < tones_hz.size(); t++) {
        phase.push_back(ph(rng));
        amp.push_back(am(rng));
    }
    std::vector<cf> x(n, cf(0, 0));
    for (size_t i = 0; i < n; i++) {
        std::complex<double> s(0, 0);
        for (size_t t = 0; t < tones_hz.size(); t++) {
            const double a = 2.0 * kPi * tones_hz[t] * static_cast<double>(i) / rate;
            s += amp[t] * std::complex<double>(std::cos(a + phase[t]),
                                               std::sin(a + phase[t]));
        }
        x[i] = cf(static_cast<float>(s.real()), static_cast<float>(s.imag()));
    }
    return x;
}

// Frequency-shift x (at `rate`) by shift_hz.
std::vector<cf> freq_shift(const std::vector<cf>& x, double rate, double shift_hz)
{
    std::vector<cf> y(x.size());
    for (size_t i = 0; i < x.size(); i++) {
        const double a = 2.0 * kPi * shift_hz * static_cast<double>(i) / rate;
        const std::complex<double> r(std::cos(a), std::sin(a));
        y[i] = static_cast<cf>(static_cast<std::complex<double>>(x[i]) * r);
    }
    return y;
}

// Software model of one radio: tune to `center_hz` (relative to fc), capture at
// per_radio_rate. Shift so center maps to DC, then resample target->per_radio.
std::vector<cf> make_slice(const std::vector<cf>& wb,
    double target_rate, double center_hz)
{
    const std::vector<cf> shifted = freq_shift(wb, target_rate, -center_hz);
    // 20e6 -> 15e6 == 3/4
    return rational_resample(shifted, 3, 4);
}

// Best complex scale c minimizing ||a - c*b|| and the resulting normalized MSE.
double residual_nmse(const std::vector<cf>& a, const std::vector<cf>& b,
    size_t lag, size_t guard)
{
    std::complex<double> num(0, 0);
    double den = 0.0, energy = 0.0;
    for (size_t i = guard; i + guard + lag < a.size() && i < b.size(); i++) {
        const std::complex<double> av = a[i + lag];
        const std::complex<double> bv = b[i];
        num += av * std::conj(bv);
        den += std::norm(bv);
        energy += std::norm(av);
    }
    const std::complex<double> c = den > 0 ? num / den : std::complex<double>(0, 0);
    double err = 0.0;
    for (size_t i = guard; i + guard + lag < a.size() && i < b.size(); i++) {
        const std::complex<double> av = a[i + lag];
        const std::complex<double> bv = b[i];
        err += std::norm(av - c * bv);
    }
    return energy > 0 ? err / energy : 1.0;
}

// Complex amplitude of a tone at frequency f (Hz) over [i0,i1) of x at `rate`.
std::complex<double> tone_amp(const std::vector<cf>& x, double rate, double f,
    size_t i0, size_t i1)
{
    std::complex<double> acc(0, 0);
    size_t cnt = 0;
    for (size_t i = i0; i < i1 && i < x.size(); i++) {
        const double a = -2.0 * kPi * f * static_cast<double>(i) / rate;
        acc += static_cast<std::complex<double>>(x[i])
             * std::complex<double>(std::cos(a), std::sin(a));
        cnt++;
    }
    return cnt ? acc / static_cast<double>(cnt) : acc;
}

// Coarse integer-lag search maximizing correlation magnitude.
size_t best_lag(const std::vector<cf>& orig, const std::vector<cf>& rec,
    size_t max_lag)
{
    size_t best = 0;
    double best_val = -1.0;
    const size_t span = std::min(orig.size(), rec.size()) / 2;
    for (size_t lag = 0; lag <= max_lag; lag++) {
        std::complex<double> acc(0, 0);
        for (size_t i = 0; i + lag < orig.size() && i < span; i++) {
            acc += static_cast<std::complex<double>>(orig[i + lag])
                 * std::conj(static_cast<std::complex<double>>(rec[i]));
        }
        const double v = std::abs(acc);
        if (v > best_val) {
            best_val = v;
            best     = lag;
        }
    }
    return best;
}

} // namespace

BOOST_AUTO_TEST_CASE(test_overlap_reconstruct_fidelity_and_alignment)
{
    const double fc_rate   = 20e6;
    const double per_radio = 15e6;
    const size_t N         = 65536;

    // Tones spread across the channel, including the overlap (±2.5 MHz) and seam.
    const std::vector<double> tones = {
        -9e6, -7e6, -5e6, -3e6, -1.5e6, -0.5e6,
        0.5e6, 1.5e6, 3e6, 5e6, 7e6, 9e6};
    const std::vector<cf> wb = make_wideband(N, fc_rate, tones);

    // Two overlapping slices: radio1 at fc-5, radio2 at fc+5.
    const std::vector<cf> rx1 = make_slice(wb, fc_rate, -5e6);
    const std::vector<cf> rx2 = make_slice(wb, fc_rate, +5e6);

    // A known gain + phase offset on radio2 (LO phase + gain mismatch).
    const double G   = 1.4;
    const double Phi = 0.7;
    const std::complex<double> off = G * std::complex<double>(std::cos(Phi),
                                                              std::sin(Phi));
    std::vector<cf> rx2_inj = rx2;
    for (auto& s : rx2_inj) {
        s = static_cast<cf>(static_cast<std::complex<double>>(s) * off);
    }

    overlap_reconstructor_config cfg;
    cfg.target_rate      = fc_rate;
    cfg.per_radio_rate   = per_radio;
    cfg.target_bw        = 20e6;
    cfg.overlap_bw       = 5e6;
    cfg.radio1_offset_hz = -5e6;
    cfg.radio2_offset_hz = +5e6;
    cfg.fft_size         = 4096;

    // Reconstruct without and with the injected offset.
    overlap_reconstructor recon_no(cfg);
    const std::vector<cf> out_no = recon_no.reconstruct(rx1, rx2);
    const overlap_alignment a_no = recon_no.alignment();

    overlap_reconstructor recon_inj(cfg);
    const std::vector<cf> out_inj = recon_inj.reconstruct(rx1, rx2_inj);
    const overlap_alignment a_inj = recon_inj.alignment();

    BOOST_REQUIRE(!out_no.empty());
    BOOST_REQUIRE(out_inj.size() == out_no.size());
    BOOST_REQUIRE(a_no.valid && a_inj.valid);

    // (a) The estimator must absorb the injected offset: applying it to radio2
    //     scales its alignment gain by 1/G and rotates its phase by -Phi
    //     relative to the un-injected case (intrinsic offset cancels).
    BOOST_TEST_MESSAGE("gain_no=" << a_no.gain << " gain_inj=" << a_inj.gain
                       << "  phase_no=" << a_no.phase_rad
                       << " phase_inj=" << a_inj.phase_rad);
    BOOST_CHECK_CLOSE(a_no.gain / a_inj.gain, G, 3.0);
    BOOST_CHECK_SMALL(
        std::remainder((a_no.phase_rad - a_inj.phase_rad) - Phi, 2.0 * kPi), 0.05);

    // (b) Reconstruction fidelity: every tone across the band (radio1-only,
    //     overlap, radio2-only) must be recovered with a flat magnitude
    //     response, and there must be no spurious energy between tones.  This
    //     is invariant to the cascade's overall (fractional) group delay.
    const size_t i0 = 4096, i1 = N - 4096; // ignore edge transients
    auto check_flatness = [&](const std::vector<cf>& out, const char* tag) {
        std::vector<double> gains;
        for (double f : tones) {
            const double in  = std::abs(tone_amp(wb, fc_rate, f, i0, i1));
            const double rec = std::abs(tone_amp(out, fc_rate, f, i0, i1));
            BOOST_REQUIRE(in > 1e-6);
            gains.push_back(rec / in);
        }
        const double mean =
            std::accumulate(gains.begin(), gains.end(), 0.0) / gains.size();
        double gmin = gains[0], gmax = gains[0];
        for (double g : gains) { gmin = std::min(gmin, g); gmax = std::max(gmax, g); }
        // No spur: probe a tone-free in-band frequency (between 7 and 9 MHz).
        const double spur = std::abs(tone_amp(out, fc_rate, 8e6, i0, i1)) / mean;
        BOOST_TEST_MESSAGE(tag << ": mean gain=" << mean << " flatness=["
                           << gmin << ", " << gmax << "] spur/mean=" << spur);
        BOOST_CHECK_GT(gmin, 0.85 * mean); // every tone recovered, ±~1.4 dB flat
        BOOST_CHECK_LT(gmax, 1.15 * mean);
        BOOST_CHECK_LT(spur, 0.05);        // < 5% spurious between tones
    };
    check_flatness(out_no, "no-offset");
    check_flatness(out_inj, "with-offset");
}

BOOST_AUTO_TEST_CASE(test_overlap_reconstruct_seam_tone)
{
    // A single tone sitting in the overlap seam (near DC) must be reconstructed
    // with the correct frequency and near-unit relative amplitude — i.e. the
    // crossfade does not create a notch or spur at the join.
    const double fc_rate   = 20e6;
    const size_t N         = 32768;
    const std::vector<double> tones = {0.4e6}; // inside the ±2.5 MHz overlap
    const std::vector<cf> wb = make_wideband(N, fc_rate, tones);

    const std::vector<cf> rx1 = make_slice(wb, fc_rate, -5e6);
    const std::vector<cf> rx2 = make_slice(wb, fc_rate, +5e6);

    overlap_reconstructor_config cfg; // defaults match this geometry
    overlap_reconstructor recon(cfg);
    const std::vector<cf> out = recon.reconstruct(rx1, rx2);
    BOOST_REQUIRE(!out.empty());
    BOOST_REQUIRE(recon.alignment().valid);

    const size_t lag  = best_lag(wb, out, 400);
    const double nmse = residual_nmse(wb, out, lag, 4096);
    BOOST_TEST_MESSAGE("seam-tone NMSE = " << nmse << " ("
                       << 10.0 * std::log10(nmse) << " dB)");
    BOOST_CHECK_LT(nmse, 0.02);
}

BOOST_AUTO_TEST_CASE(test_overlap_reconstruct_streaming)
{
    // The streaming path (process() fed in small chunks) must produce a
    // reconstruction with the same band-flat fidelity as the one-shot path.
    const double fc_rate = 20e6;
    const size_t N       = 131072;
    const std::vector<double> tones = {-9e6, -5e6, -1e6, 1e6, 5e6, 9e6};
    const std::vector<cf> wb = make_wideband(N, fc_rate, tones);

    const std::vector<cf> rx1 = make_slice(wb, fc_rate, -5e6);
    const std::vector<cf> rx2 = make_slice(wb, fc_rate, +5e6);
    const size_t n = std::min(rx1.size(), rx2.size());

    overlap_reconstructor_config cfg; // defaults match this geometry
    overlap_reconstructor recon(cfg);

    // Feed in irregular small chunks to exercise the streaming buffering.
    std::vector<cf> out;
    const std::vector<size_t> chunks = {1000, 4096, 333, 8192, 777};
    size_t i = 0, c = 0;
    while (i < n) {
        const size_t sz = std::min(chunks[c % chunks.size()], n - i);
        const std::vector<cf> piece = recon.process(&rx1[i], &rx2[i], sz);
        out.insert(out.end(), piece.begin(), piece.end());
        i += sz;
        c++;
    }
    BOOST_REQUIRE(out.size() > N / 2);
    BOOST_REQUIRE(recon.alignment().valid);

    // Band-flat fidelity over the steady-state interior (skip warmup edges).
    const size_t i0 = 8192, i1 = out.size() - 8192;
    std::vector<double> gains;
    for (double f : tones) {
        const double in  = std::abs(tone_amp(wb, fc_rate, f, i0, i1));
        const double rec = std::abs(tone_amp(out, fc_rate, f, i0, i1));
        BOOST_REQUIRE(in > 1e-6);
        gains.push_back(rec / in);
    }
    const double mean =
        std::accumulate(gains.begin(), gains.end(), 0.0) / gains.size();
    double gmin = gains[0], gmax = gains[0];
    for (double g : gains) { gmin = std::min(gmin, g); gmax = std::max(gmax, g); }
    const double spur = std::abs(tone_amp(out, fc_rate, 7e6, i0, i1)) / mean;
    BOOST_TEST_MESSAGE("streaming: out=" << out.size() << " mean gain=" << mean
                       << " flatness=[" << gmin << ", " << gmax
                       << "] spur/mean=" << spur);
    BOOST_CHECK_GT(gmin, 0.85 * mean);
    BOOST_CHECK_LT(gmax, 1.15 * mean);
    BOOST_CHECK_LT(spur, 0.06);
}
