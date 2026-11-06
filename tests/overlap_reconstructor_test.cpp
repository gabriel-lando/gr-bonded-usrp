//
// Copyright 2026 Gabriel Lando
//
// End-to-end (no hardware) validation of the overlap-aware reconstruction.
//
// A known wideband signal at the target rate is split in software into two
// overlapping slices (what the bonded radios capture), a known gain/phase
// offset is injected on the second slice, and the reconstructor rebuilds the
// original baseband. The tests check the alignment estimate, band flatness
// (including the overlap seam), the warmup and frozen paths, streaming,
// threading and the direct-output API.
//

#include "overlap_reconstructor.hpp"
#include "test_helpers.hpp"
#include <boost/test/unit_test.hpp>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
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
// target_rate*interp/decim (default 20 -> 15 MS/s).
std::vector<cf> make_slice(const std::vector<cf>& wb,
    double target_rate, double center_hz, size_t interp = 3, size_t decim = 4)
{
    const std::vector<cf> shifted = freq_shift(wb, target_rate, -center_hz);
    return bonded_test::rational_resample(shifted, interp, decim);
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

// Per-tone gain of `out` relative to `wb` over [i0, i1): {mean, min, max}.
struct flatness
{
    double mean, gmin, gmax;
};
flatness measure_flatness(const std::vector<cf>& wb, const std::vector<cf>& out,
    double rate, const std::vector<double>& tones, size_t i0, size_t i1)
{
    std::vector<double> gains;
    for (double f : tones) {
        const double in  = std::abs(tone_amp(wb, rate, f, i0, i1));
        const double rec = std::abs(tone_amp(out, rate, f, i0, i1));
        gains.push_back(in > 1e-9 ? rec / in : 0.0);
    }
    flatness r{std::accumulate(gains.begin(), gains.end(), 0.0) / gains.size(),
        gains[0], gains[0]};
    for (double g : gains) {
        r.gmin = std::min(r.gmin, g);
        r.gmax = std::max(r.gmax, g);
    }
    return r;
}

// Stream rx1/rx2 through process() in the given repeating chunk sizes.
std::vector<cf> stream_through(overlap_reconstructor& recon, const std::vector<cf>& rx1,
    const std::vector<cf>& rx2, const std::vector<size_t>& chunks)
{
    const size_t n = std::min(rx1.size(), rx2.size());
    std::vector<cf> out;
    size_t i = 0, c = 0;
    while (i < n) {
        const size_t sz = std::min(chunks[c % chunks.size()], n - i);
        const std::vector<cf> piece = recon.process(&rx1[i], &rx2[i], sz);
        out.insert(out.end(), piece.begin(), piece.end());
        i += sz;
        c++;
    }
    return out;
}

// Sets an environment variable for the lifetime of the object.
struct scoped_env
{
    std::string name;
    scoped_env(const char* n, const char* v) : name(n) { setenv(n, v, 1); }
    ~scoped_env() { unsetenv(name.c_str()); }
};

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

BOOST_AUTO_TEST_CASE(test_overlap_reconstruct_frozen_steady_state)
{
    // Exercise the frozen (steady-state) path: with a short warmup window the
    // reconstructor estimates+freezes alignment early, then runs its fast path
    // for the remainder. Fidelity over the frozen region must match the warmup
    // path's band-flat response (no seam/spur/gain regression at the switch).
    const double fc_rate = 20e6;
    const size_t N       = 131072;
    const std::vector<double> tones = {-9e6, -5e6, -1e6, 1e6, 5e6, 9e6};
    const std::vector<cf> wb = make_wideband(N, fc_rate, tones);

    const std::vector<cf> rx1 = make_slice(wb, fc_rate, -5e6);
    const std::vector<cf> rx2 = make_slice(wb, fc_rate, +5e6);
    const size_t n = std::min(rx1.size(), rx2.size());

    overlap_reconstructor_config cfg; // defaults match this geometry
    cfg.warmup_samples = 20000;       // freeze early so most output is frozen

    overlap_reconstructor recon(cfg);
    std::vector<cf> out;
    const std::vector<size_t> chunks = {4096, 8192, 2048, 16384};
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

    // Check fidelity well past the freeze point (warmup ~20k output samples).
    const size_t i0 = 40000, i1 = out.size() - 8192;
    BOOST_REQUIRE(i1 > i0);
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
    BOOST_TEST_MESSAGE("frozen: out=" << out.size() << " mean gain=" << mean
                       << " flatness=[" << gmin << ", " << gmax
                       << "] spur/mean=" << spur);
    BOOST_CHECK_GT(gmin, 0.85 * mean);
    BOOST_CHECK_LT(gmax, 1.15 * mean);
    BOOST_CHECK_LT(spur, 0.06);
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

BOOST_AUTO_TEST_CASE(test_overlap_reconstruct_edge_aligned_plan)
{
    // The plan used live: 20 MHz channel, 5 MHz overlap, each radio captures
    // exactly [fc-10, fc+2.5] or [fc-2.5, fc+10] at 12.5 MS/s, at fc -/+ 3.75.
    // The overlap content is symmetric (equal-amplitude tone pairs): the gain
    // estimate is a power ratio over the overlap, and each radio's anti-alias
    // roll-off attenuates the far side of it, which cancels for a flat
    // spectrum such as 802.11 but biases the gain for lopsided content.
    const double fc_rate = 20e6;
    const size_t N       = 131072;
    const std::vector<double> tones = {-8e6, -5e6, -1.5e6, -0.4e6, 0.4e6, 1.5e6, 5e6, 8e6};
    std::vector<cf> wb(N, cf(0, 0));
    for (size_t t = 0; t < tones.size(); t++) {
        const auto x = bonded_test::tone(tones[t], fc_rate, N, 0.9 * static_cast<double>(t), 0.5);
        for (size_t i = 0; i < N; i++) {
            wb[i] += x[i];
        }
    }
    const std::vector<cf> rx1 = make_slice(wb, fc_rate, -3.75e6, 5, 8);
    const std::vector<cf> rx2 = make_slice(wb, fc_rate, +3.75e6, 5, 8);

    overlap_reconstructor_config cfg;
    cfg.per_radio_rate   = 12.5e6;
    cfg.radio1_offset_hz = -3.75e6;
    cfg.radio2_offset_hz = +3.75e6;
    cfg.warmup_samples   = 20000;
    overlap_reconstructor recon(cfg);
    const std::vector<cf> out = stream_through(recon, rx1, rx2, {4096, 3000, 8192});
    BOOST_REQUIRE(recon.alignment().valid);
    BOOST_REQUIRE(out.size() > N / 2);

    const auto fl = measure_flatness(wb, out, fc_rate, tones, 40000, out.size() - 8192);
    const double spur = std::abs(tone_amp(out, fc_rate, 3e6, 40000, out.size() - 8192)) / fl.mean;
    BOOST_TEST_MESSAGE("edge-aligned: mean=" << fl.mean << " [" << fl.gmin << ", "
                                             << fl.gmax << "] spur=" << spur);
    BOOST_CHECK_GT(fl.gmin, 0.97 * fl.mean);
    BOOST_CHECK_LT(fl.gmax, 1.03 * fl.mean);
    BOOST_CHECK_LT(spur, 0.01);
}

BOOST_AUTO_TEST_CASE(test_overlap_reconstruct_threads_match_serial)
{
    // The frozen path splits each call's output range across BONDED_THREADS
    // threads; the result must match the single-threaded one (up to the tiny
    // difference between the recursive and the closed-form rotator phase).
    const double fc_rate = 20e6;
    const size_t N       = 65536;
    const std::vector<cf> wb = make_wideband(N, fc_rate, {-7e6, -1e6, 2e6, 6e6});
    const std::vector<cf> rx1 = make_slice(wb, fc_rate, -5e6);
    const std::vector<cf> rx2 = make_slice(wb, fc_rate, +5e6);

    overlap_reconstructor_config cfg;
    cfg.warmup_samples = 4096;
    const std::vector<size_t> chunks = {8192, 4096, 2048, 12000}; // >= 2048: parallel

    std::vector<cf> serial, threaded;
    {
        scoped_env env("BONDED_THREADS", "1");
        overlap_reconstructor recon(cfg);
        serial = stream_through(recon, rx1, rx2, chunks);
    }
    for (const char* nt : {"2", "3", "4", "7"}) {
        scoped_env env("BONDED_THREADS", nt);
        overlap_reconstructor recon(cfg);
        threaded = stream_through(recon, rx1, rx2, chunks);
        BOOST_REQUIRE_EQUAL(threaded.size(), serial.size());
        double err = 0.0, ref = 0.0;
        for (size_t i = 0; i < serial.size(); i++) {
            err += std::norm(threaded[i] - serial[i]);
            ref += std::norm(serial[i]);
        }
        BOOST_TEST_MESSAGE("threads=" << nt << " rel err=" << err / ref);
        BOOST_CHECK_LT(err / ref, 1e-6);
    }
    {
        scoped_env env("BONDED_PARALLEL", "0");
        overlap_reconstructor recon(cfg);
        const std::vector<cf> off = stream_through(recon, rx1, rx2, chunks);
        BOOST_REQUIRE_EQUAL(off.size(), serial.size());
        double err = 0.0;
        for (size_t i = 0; i < serial.size(); i++) {
            err += std::norm(off[i] - serial[i]);
        }
        BOOST_CHECK_EQUAL(err, 0.0);
    }
}

BOOST_AUTO_TEST_CASE(test_overlap_process_into_matches_process)
{
    // process_into() + pending() must produce what process() does, whatever
    // the destination capacity (including capacities smaller than a call's
    // output, which spill into pending()). The two combine loops may be
    // vectorized differently (-ffast-math), hence a 1e-6 tolerance.
    const double fc_rate = 20e6;
    const size_t N       = 65536;
    const std::vector<cf> wb = make_wideband(N, fc_rate, {-6e6, 0.5e6, 7e6});
    const std::vector<cf> rx1 = make_slice(wb, fc_rate, -5e6);
    const std::vector<cf> rx2 = make_slice(wb, fc_rate, +5e6);
    const size_t n = std::min(rx1.size(), rx2.size());

    overlap_reconstructor_config cfg;
    cfg.warmup_samples = 8192;
    overlap_reconstructor ref_recon(cfg);
    const std::vector<cf> ref = stream_through(ref_recon, rx1, rx2, {4096});

    for (size_t cap : {size_t(100000), size_t(5000), size_t(1)}) {
        overlap_reconstructor recon(cfg);
        std::vector<cf> out, dst(cap);
        size_t direct = 0;
        for (size_t i = 0; i < n; i += 4096) {
            const size_t sz = std::min<size_t>(4096, n - i);
            const size_t w = recon.process_into(&rx1[i], &rx2[i], sz, dst.data(), cap);
            BOOST_REQUIRE_LE(w, cap);
            out.insert(out.end(), dst.begin(), dst.begin() + static_cast<std::ptrdiff_t>(w));
            const auto& rest = recon.pending();
            out.insert(out.end(), rest.begin(), rest.end());
            direct += w;
        }
        BOOST_TEST_MESSAGE("cap=" << cap << " direct=" << direct << " total=" << out.size());
        BOOST_CHECK_GT(direct, 0u); // the frozen path did write directly
        BOOST_REQUIRE_EQUAL(out.size(), ref.size());
        for (size_t i = 0; i < out.size(); i++) {
            if (std::abs(out[i] - ref[i]) > 1e-6f * (1.0f + std::abs(ref[i]))) {
                BOOST_FAIL("cap " << cap << ": differs at sample " << i);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(test_overlap_reconstruct_restarts_from_scratch)
{
    // reconstruct() resets the state, so two calls give identical output, and
    // reset() makes process() start over too.
    const double fc_rate = 20e6;
    const std::vector<cf> wb = make_wideband(32768, fc_rate, {-4e6, 3e6});
    const std::vector<cf> rx1 = make_slice(wb, fc_rate, -5e6);
    const std::vector<cf> rx2 = make_slice(wb, fc_rate, +5e6);
    overlap_reconstructor recon(overlap_reconstructor_config{});
    const std::vector<cf> a = recon.reconstruct(rx1, rx2);
    const std::vector<cf> b = recon.reconstruct(rx1, rx2);
    BOOST_CHECK(a == b);
    recon.reset();
    const std::vector<cf> c = recon.process(rx1.data(), rx2.data(), std::min(rx1.size(), rx2.size()));
    BOOST_CHECK(a == c);
}

BOOST_AUTO_TEST_CASE(test_overlap_reconstruct_rejects_bad_config)
{
    overlap_reconstructor_config cfg;
    cfg.target_rate = 0;
    BOOST_CHECK_THROW(overlap_reconstructor{cfg}, std::invalid_argument);
    cfg            = overlap_reconstructor_config{};
    cfg.per_radio_rate = -1;
    BOOST_CHECK_THROW(overlap_reconstructor{cfg}, std::invalid_argument);
    cfg            = overlap_reconstructor_config{};
    cfg.overlap_bw = 0;
    BOOST_CHECK_THROW(overlap_reconstructor{cfg}, std::invalid_argument);
    BOOST_CHECK_NO_THROW(overlap_reconstructor{overlap_reconstructor_config{}});
}

BOOST_AUTO_TEST_CASE(test_overlap_frozen_path_continues_warmup_path)
{
    // After the freeze the fused path must produce the same samples the
    // warmup path would have (same timeline, same per-radio phase), for any
    // freeze point. A fused path restarting its own timeline at the freeze
    // changes each radio's translation phase by a different amount, which
    // leaves the frozen alignment off and notches the overlap.
    struct plan
    {
        double per_radio, offset;
        size_t L, M;
    };
    const double fc_rate = 20e6;
    const size_t N       = 98304;
    const std::vector<cf> wb = make_wideband(N, fc_rate, {-7e6, -2e6, -0.7e6, 0.4e6, 1.8e6, 6e6});
    for (const plan& p : {plan{15e6, 5e6, 3, 4}, plan{12.5e6, 3.75e6, 5, 8}}) {
        const std::vector<cf> rx1 = make_slice(wb, fc_rate, -p.offset, p.L, p.M);
        const std::vector<cf> rx2 = make_slice(wb, fc_rate, +p.offset, p.L, p.M);
        overlap_reconstructor_config cfg;
        cfg.per_radio_rate   = p.per_radio;
        cfg.radio1_offset_hz = -p.offset;
        cfg.radio2_offset_hz = +p.offset;
        cfg.warmup_samples   = 1u << 30; // never freezes: warmup path only
        // Phase tracking would refine the alignment after the freeze; this
        // test is about the timeline only.
        scoped_env no_track("BONDED_TRACK", "0");
        overlap_reconstructor ref_recon(cfg);
        const std::vector<size_t> chunks = {3000, 4096, 1111};
        const std::vector<cf> ref = stream_through(ref_recon, rx1, rx2, chunks);

        for (size_t warmup : {5000u, 10001u, 20000u, 33333u}) {
            cfg.warmup_samples = warmup;
            overlap_reconstructor recon(cfg);
            const std::vector<cf> out = stream_through(recon, rx1, rx2, chunks);
            BOOST_REQUIRE_EQUAL(out.size(), ref.size());
            // Compare from just after the freeze (the call that crosses
            // `warmup` still runs the warmup path) to the end.
            const size_t from = warmup + 8192;
            double err = 0.0, energy = 0.0;
            for (size_t i = from; i < out.size(); i++) {
                err += std::norm(out[i] - ref[i]);
                energy += std::norm(ref[i]);
            }
            BOOST_TEST_MESSAGE("plan " << p.per_radio / 1e6 << " MS/s, warmup " << warmup
                                       << ": frozen vs warmup rel err " << err / energy);
            BOOST_CHECK_LT(err / energy, 1e-6);
        }
    }
}

BOOST_AUTO_TEST_CASE(test_overlap_reconstruct_no_phase_drift)
{
    // Each radio's translation must keep its exact frequency for the whole
    // session. With float phase increments radio 1 drifted ~6 deg/s at
    // 20 MHz and radio 2 the same the other way, so a frozen alignment
    // walked off (the RX pair's "phase drift" seen live).
    const double T = 20e6, R = 12.5e6, c1 = -3.75e6, c2 = +3.75e6;
    const size_t n = static_cast<size_t>(0.4 * R);
    // Ideal coherent radios: a tone at 0 Hz seen by both (for the
    // alignment), one only radio 1 sees and one only radio 2 sees.
    std::vector<cf> rx1 = bonded_test::tone(0.0 - c1, R, n, 0.0, 0.3);
    std::vector<cf> rx2 = bonded_test::tone(0.0 - c2, R, n, 0.0, 0.3);
    const auto t1 = bonded_test::tone(-6e6 - c1, R, n, 0.0, 0.3);
    const auto t2 = bonded_test::tone(+6e6 - c2, R, n, 0.0, 0.3);
    for (size_t i = 0; i < n; i++) {
        rx1[i] += t1[i];
        rx2[i] += t2[i];
    }
    overlap_reconstructor_config cfg;
    cfg.per_radio_rate   = R;
    cfg.radio1_offset_hz = c1;
    cfg.radio2_offset_hz = c2;
    cfg.warmup_samples   = 20000;
    overlap_reconstructor recon(cfg);
    const std::vector<cf> out = stream_through(recon, rx1, rx2, {32768, 20000});
    BOOST_REQUIRE(recon.alignment().valid);
    const size_t w = 20000, early = 100000, late = out.size() - 2 * w;
    for (double f : {-6e6, 6e6}) {
        const double d = std::arg(tone_amp(out, T, f, late, late + w)
                                  / tone_amp(out, T, f, early, early + w));
        BOOST_TEST_MESSAGE("tone " << f / 1e6 << " MHz: phase moved "
                                   << d * 180 / kPi << " deg over "
                                   << (late - early) / T << " s");
        BOOST_CHECK_LT(std::abs(d) * 180 / kPi, 0.2);
    }
}

namespace {
// Radio inputs for the edge-aligned plan: per radio, `noise_n` wideband
// samples of noise of its own (no common content), then `sig_n` samples of a
// common symmetric multitone; radio 2 is rotated by -theta, so the right
// alignment phase is theta.
std::pair<std::vector<cf>, std::vector<cf>> noise_then_common(
    size_t noise_n, size_t sig_n, double theta, bool noise_after = false)
{
    const double T = 20e6;
    std::vector<cf> sig(sig_n, cf(0, 0));
    const std::vector<double> tones = {-8e6, -5e6, -1.5e6, -0.4e6, 0.4e6, 1.5e6, 5e6, 8e6};
    for (size_t t = 0; t < tones.size(); t++) {
        const auto x = bonded_test::tone(tones[t], T, sig_n, 0.9 * static_cast<double>(t), 0.3);
        for (size_t i = 0; i < sig_n; i++) sig[i] += x[i];
    }
    std::mt19937 rng(77);
    std::normal_distribution<float> g(0.0f, 0.3f);
    auto noisy = [&](const std::vector<cf>& common) {
        std::vector<cf> noise(noise_n);
        for (auto& v : noise) v = cf(g(rng), g(rng));
        std::vector<cf> w = noise_after ? common : noise;
        const auto& tail  = noise_after ? noise : common;
        w.insert(w.end(), tail.begin(), tail.end());
        return w;
    };
    std::vector<cf> rx1 = make_slice(noisy(sig), T, -3.75e6, 5, 8);
    std::vector<cf> rx2 = make_slice(noisy(sig), T, +3.75e6, 5, 8);
    const cf rot(static_cast<float>(std::cos(theta)), static_cast<float>(-std::sin(theta)));
    for (auto& v : rx2) v *= rot;
    return {rx1, rx2};
}

overlap_reconstructor_config edge_cfg()
{
    overlap_reconstructor_config cfg;
    cfg.per_radio_rate   = 12.5e6;
    cfg.radio1_offset_hz = -3.75e6;
    cfg.radio2_offset_hz = +3.75e6;
    cfg.warmup_samples   = 20000;
    return cfg;
}

//! Alignment phase a coherent warmup finds for the same content: the target
//! (the test radios add a fixed phase of their own to theta).
double coherent_phase(double theta)
{
    scoped_env no_track("BONDED_TRACK", "0");
    const auto rx = noise_then_common(0, 200000, theta);
    overlap_reconstructor recon(edge_cfg());
    (void)stream_through(recon, rx.first, rx.second, {8192});
    return recon.alignment().phase_rad;
}
} // namespace

BOOST_AUTO_TEST_CASE(test_overlap_tracks_phase_after_noise_warmup)
{
    // The warmup sees only uncorrelated noise in the overlap, so it freezes
    // on an arbitrary phase; once a common signal appears, the frozen path
    // must adopt its phase and the seam must be flat again.
    const double theta = 2.0;
    const auto rx = noise_then_common(60000, 200000, theta);
    const double want = coherent_phase(theta);
    for (const char* track : {"1", "0"}) {
        scoped_env env("BONDED_TRACK", track);
        overlap_reconstructor recon(edge_cfg());
        (void)stream_through(recon, rx.first, rx.second, {8192});
        BOOST_REQUIRE(recon.alignment().valid);
        const double err = std::remainder(recon.alignment().phase_rad - want, 2 * kPi);
        BOOST_TEST_MESSAGE("tracking " << track << ": phase error " << err * 180 / kPi << " deg");
        if (track[0] == '1') {
            BOOST_CHECK_LT(std::abs(err) * 180 / kPi, 1.0);
        } else {
            BOOST_CHECK_GT(std::abs(err) * 180 / kPi, 10.0); // the noise warmup was wrong
        }
    }
}

BOOST_AUTO_TEST_CASE(test_overlap_tracking_ignores_noise)
{
    // After a coherent warmup, a long stretch with only independent noise in
    // each radio must not move the frozen phase.
    const double theta = -1.0;
    const auto rx = noise_then_common(400000, 60000, theta, /*noise_after=*/true);
    const double want = coherent_phase(theta);
    overlap_reconstructor recon(edge_cfg());
    (void)stream_through(recon, rx.first, rx.second, {8192});
    const double err = std::remainder(recon.alignment().phase_rad - want, 2 * kPi);
    BOOST_TEST_MESSAGE("phase after noise " << err * 180 / kPi << " deg off");
    BOOST_CHECK_LT(std::abs(err) * 180 / kPi, 1.0);
}

BOOST_AUTO_TEST_CASE(test_overlap_tracking_follows_the_strongest_signal)
{
    // Each radio has its own antenna, so the right phase differs per emitter.
    // A weak continuous source (phase theta_w, seen during the warmup) and
    // strong bursts (theta_s): the alignment must end up on the bursts, the
    // signal that needs it.
    const double T = 20e6, tw = -0.5, ts = 1.8;
    const size_t n = 360000;
    const std::vector<double> tones = {-8e6, -5e6, -1.5e6, -0.4e6, 0.4e6, 1.5e6, 5e6, 8e6};
    std::vector<cf> shape(n, cf(0, 0));
    for (size_t t = 0; t < tones.size(); t++) {
        const auto x = bonded_test::tone(tones[t], T, n, 0.9 * static_cast<double>(t), 1.0);
        for (size_t i = 0; i < n; i++) shape[i] += x[i];
    }
    auto burst_on = [](size_t i) { return i > 80000 && (i % 70000) < 8000; };
    std::vector<cf> w1(n), w2(n);
    using cd = std::complex<double>;
    const cd rw(std::cos(tw), -std::sin(tw)), rs(std::cos(ts), -std::sin(ts));
    for (size_t i = 0; i < n; i++) {
        const cd base  = static_cast<cd>(shape[i]);
        const double s = burst_on(i) ? 0.3 : 0.0;
        w1[i] = static_cast<cf>(base * (0.01 + s));
        w2[i] = static_cast<cf>(base * (0.01 * rw + s * rs));
    }
    // Radio 2's rotation is applied before its slice, like a channel.
    const std::vector<cf> rx1 = make_slice(w1, T, -3.75e6, 5, 8);
    const std::vector<cf> rx2 = make_slice(w2, T, +3.75e6, 5, 8);
    overlap_reconstructor recon(edge_cfg());
    (void)stream_through(recon, rx1, rx2, {8192});
    const double want = coherent_phase(ts);
    const double err  = std::remainder(recon.alignment().phase_rad - want, 2 * kPi);
    BOOST_TEST_MESSAGE("strongest-signal tracking error " << err * 180 / kPi << " deg (weak source "
                       << std::remainder(coherent_phase(tw) - want, 2 * kPi) * 180 / kPi
                       << " deg away)");
    BOOST_CHECK_LT(std::abs(err) * 180 / kPi, 3.0);
}

BOOST_AUTO_TEST_CASE(test_overlap_skip_input_keeps_timeline)
{
    // A stretch both radios lost (a block dropped after an overflow) must be
    // stepped over, not closed up: closing it up moves each radio's
    // translation phase by a different amount, and with it the alignment.
    // After the gap the output must match an uninterrupted run, shifted by
    // the skipped outputs, once the filters' history has flushed.
    scoped_env no_track("BONDED_TRACK", "0");
    const double T = 20e6;
    const std::vector<cf> wb = make_wideband(400000, T, {-8e6, -4e6, -1.2e6, 0.3e6, 1.4e6, 5e6, 8.5e6});
    const std::vector<cf> rx1 = make_slice(wb, T, -3.75e6, 5, 8);
    const std::vector<cf> rx2 = make_slice(wb, T, +3.75e6, 5, 8);
    const size_t n = std::min(rx1.size(), rx2.size());
    overlap_reconstructor ref_recon(edge_cfg());
    const std::vector<cf> ref = stream_through(ref_recon, rx1, rx2, {8192});

    auto feed = [&](overlap_reconstructor& r, size_t from, size_t to, std::vector<cf>& out) {
        for (size_t i = from; i < to; i += 8192) {
            const size_t sz = std::min<size_t>(8192, to - i);
            const std::vector<cf> piece = r.process(&rx1[i], &rx2[i], sz);
            out.insert(out.end(), piece.begin(), piece.end());
        }
    };
    const size_t at = 60000; // input sample; the warmup is long over
    for (size_t gap : {131072u, 130000u, 7u}) {
        for (bool skip : {true, false}) {
            overlap_reconstructor recon(edge_cfg());
            std::vector<cf> out;
            feed(recon, 0, at, out);
            const size_t before = out.size();
            if (skip) {
                const std::vector<cf>& z = recon.skip_input(gap);
                out.insert(out.end(), z.begin(), z.end());
            }
            feed(recon, at + gap, n, out);
            // Outputs the skipped whole periods would have produced.
            const size_t shift = skip ? (gap / 5) * 8 : ref.size() - out.size();
            BOOST_REQUIRE_EQUAL(out.size() + shift, ref.size());
            double err = 0.0, energy = 0.0;
            for (size_t i = before + 4096; i < out.size(); i++) {
                err += std::norm(out[i] - ref[i + shift]);
                energy += std::norm(ref[i + shift]);
            }
            BOOST_TEST_MESSAGE("gap " << gap << (skip ? " skipped" : " closed up")
                                      << ": rel err after the gap " << err / energy);
            if (skip) {
                BOOST_CHECK_LT(err / energy, 1e-6);
            } else if (gap % 5 != 0) {
                // Closed up, a gap of whole resampler periods only rotates
                // both radios alike in this plan (offsets of 3/16 T); any
                // other gap moves their output grid and relative phase.
                BOOST_CHECK_GT(err / energy, 1e-2);
            }
        }
    }
}
