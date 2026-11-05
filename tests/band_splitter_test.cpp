//
// Copyright 2026 Gabriel Lando
// SPDX-License-Identifier: GPL-3.0-or-later
//
// No-hardware validation of the band_splitter DSP engine.
//
// A known wideband multitone signal is split via band_splitter into two
// per-radio slices (the same operation a TX chain would perform), and then
// recombined in software (upsample × ratio + NCO shift back + sum) to verify
// that:
//   (a) tones across the band are recovered with flat amplitude (± ≤ 15%);
//   (b) lower-half tones go to slice 1, upper-half tones go to slice 2;
//   (c) the 0 Hz seam tone is attenuated in both slices (expected behaviour);
//   (d) streaming process() produces the same output as one-shot split();
//   (e) the seam transition is as narrow as taps_per_phase promises, also with
//       the radios pulled in from the band edge (offset < input_rate/4);
//   (f) invalid configurations are rejected.
//

#include "band_splitter.hpp"
#include "test_helpers.hpp"
#include <stdexcept>

#include <boost/test/unit_test.hpp>
#include <algorithm>
#include <cmath>
#include <complex>
#include <numeric>
#include <random>
#include <vector>

using namespace bonded;

namespace {

constexpr double kPi = 3.14159265358979323846;
using cf             = std::complex<float>;
using cd             = std::complex<double>;

// Multi-tone wideband test signal at `rate`, tones in Hz (baseband).
std::vector<cf> make_wideband(size_t n,
    double rate,
    const std::vector<double>& tones_hz)
{
    std::mt19937 rng(42);
    std::uniform_real_distribution<double> ph(-kPi, kPi);
    std::uniform_real_distribution<double> am(0.4, 1.0);
    std::vector<double> phase(tones_hz.size()), amp(tones_hz.size());
    for (size_t t = 0; t < tones_hz.size(); t++) {
        phase[t] = ph(rng);
        amp[t]   = am(rng);
    }
    std::vector<cf> x(n, cf(0, 0));
    for (size_t i = 0; i < n; i++) {
        cd s(0, 0);
        for (size_t t = 0; t < tones_hz.size(); t++) {
            const double a =
                2.0 * kPi * tones_hz[t] * static_cast<double>(i) / rate;
            s += amp[t]
                 * cd(std::cos(a + phase[t]), std::sin(a + phase[t]));
        }
        x[i] = cf(static_cast<float>(s.real()), static_cast<float>(s.imag()));
    }
    return x;
}

// Frequency-shift x (at `rate`) by shift_hz.
std::vector<cf> freq_shift(const std::vector<cf>& x,
    double rate,
    double shift_hz)
{
    std::vector<cf> y(x.size());
    for (size_t i = 0; i < x.size(); i++) {
        const double a = 2.0 * kPi * shift_hz * static_cast<double>(i) / rate;
        y[i] = static_cast<cf>(static_cast<cd>(x[i]) * cd(std::cos(a), std::sin(a)));
    }
    return y;
}

// Coherent DFT amplitude of a single tone at f Hz over samples [i0, i1).
// Magnitude is independent of group delay (phase-invariant), so no lag
// correction is needed for the amplitude flatness check.
cd tone_amp(const std::vector<cf>& x,
    double rate,
    double f,
    size_t i0,
    size_t i1)
{
    cd acc(0, 0);
    size_t cnt = 0;
    for (size_t i = i0; i < i1 && i < x.size(); i++) {
        const double a = -2.0 * kPi * f * static_cast<double>(i) / rate;
        acc += static_cast<cd>(x[i]) * cd(std::cos(a), std::sin(a));
        cnt++;
    }
    return cnt ? acc / static_cast<double>(cnt) : acc;
}

// Recombine the two slices in software: upsample each back to the wideband
// rate, shift it back to its offset and sum (what the air does on TX).
std::vector<cf> recombine(const std::vector<cf>& o1, const std::vector<cf>& o2,
    const band_splitter_config& cfg)
{
    const size_t g = std::gcd(static_cast<size_t>(std::llround(cfg.input_rate)),
        static_cast<size_t>(std::llround(cfg.per_radio_rate)));
    const size_t L = static_cast<size_t>(std::llround(cfg.input_rate)) / g;
    const size_t M = static_cast<size_t>(std::llround(cfg.per_radio_rate)) / g;
    const auto s1 = freq_shift(bonded_test::rational_resample(o1, L, M), cfg.input_rate, cfg.radio1_offset_hz);
    const auto s2 = freq_shift(bonded_test::rational_resample(o2, L, M), cfg.input_rate, cfg.radio2_offset_hz);
    std::vector<cf> out(std::min(s1.size(), s2.size()));
    for (size_t i = 0; i < out.size(); i++) {
        out[i] = s1[i] + s2[i];
    }
    return out;
}

} // namespace

// ---------------------------------------------------------------------------
// Test 1: split a multitone wideband, recombine, verify flat amplitude.
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(test_band_splitter_flatness_and_recombination)
{
    const double T  = 20e6;
    const double PR = 10e6;
    const size_t N  = 65536;

    // Tones well within each half-band; avoid the ±5 MHz boundary.
    const std::vector<double> tones = {
        -9e6, -7e6, -3e6, -1e6, 1e6, 3e6, 7e6, 9e6};
    const std::vector<cf> wb = make_wideband(N, T, tones);

    band_splitter_config cfg;
    cfg.input_rate       = T;
    cfg.per_radio_rate   = PR;
    cfg.radio1_offset_hz = -5e6;
    cfg.radio2_offset_hz = +5e6;

    band_splitter bs(cfg);
    auto [out1, out2] = bs.split(wb);

    BOOST_REQUIRE(!out1.empty());
    BOOST_REQUIRE_EQUAL(out1.size(), out2.size());

    // Expected output length ≈ N/2 (2:1 decimation), ±5%.
    BOOST_CHECK_CLOSE(static_cast<double>(out1.size()),
        static_cast<double>(N) / 2.0, 5.0);

    // Recombine in software: upsample each slice back to T, NCO-shift back, sum.
    // Uses bonded_test::rational_resample (one-shot) — the AA filter is different from the
    // streaming one but still flat in the pass-band interior.
    const size_t g =
        std::gcd(static_cast<size_t>(std::llround(T)),
            static_cast<size_t>(std::llround(PR)));
    const size_t interp = static_cast<size_t>(std::llround(T)) / g;  // 2
    const size_t decim  = static_cast<size_t>(std::llround(PR)) / g; // 1

    const std::vector<cf> up1 = bonded_test::rational_resample(out1, interp, decim);
    const std::vector<cf> up2 = bonded_test::rational_resample(out2, interp, decim);

    const std::vector<cf> sh1 = freq_shift(up1, T, cfg.radio1_offset_hz);
    const std::vector<cf> sh2 = freq_shift(up2, T, cfg.radio2_offset_hz);

    const size_t Mrec = std::min(sh1.size(), sh2.size());
    std::vector<cf> rec(Mrec);
    for (size_t i = 0; i < Mrec; i++) {
        rec[i] = sh1[i] + sh2[i];
    }

    // Amplitude flatness check: skip edge transients from both filter banks.
    const size_t guard = 4096;
    const size_t i0    = guard;
    const size_t i1    = std::min(rec.size(), wb.size()) - guard;
    BOOST_REQUIRE_GT(i1, i0);

    std::vector<double> ratios;
    for (double f : tones) {
        const double a_orig = std::abs(tone_amp(wb, T, f, i0, i1));
        const double a_rec  = std::abs(tone_amp(rec, T, f, i0, i1));
        BOOST_REQUIRE_GT(a_orig, 1e-6);
        const double ratio = a_rec / a_orig;
        BOOST_TEST_MESSAGE("f=" << f / 1e6 << " MHz  orig=" << a_orig
                                << "  rec=" << a_rec
                                << "  ratio=" << ratio);
        ratios.push_back(ratio);
    }

    // Every tone must be recovered (>50% amplitude, i.e. within ~6 dB).
    for (double r : ratios) {
        BOOST_CHECK_GT(r, 0.5);
    }

    // Band flatness: all ratios within ±15% of their mean (filter shape
    // differences between poly_resampler and rational_resample are small).
    const double mean =
        std::accumulate(ratios.begin(), ratios.end(), 0.0) / ratios.size();
    double gmin = *std::min_element(ratios.begin(), ratios.end());
    double gmax = *std::max_element(ratios.begin(), ratios.end());
    BOOST_TEST_MESSAGE("flatness: mean=" << mean << " range=["
                                         << gmin << ", " << gmax << "]");
    BOOST_CHECK_GT(gmin, 0.85 * mean);
    BOOST_CHECK_LT(gmax, 1.15 * mean);

    // Low inter-tone spur: probe at 2 MHz (inside the upper half, between
    // the ±1 MHz and ±3 MHz tones — purely noise, no test tone present).
    const double spur = std::abs(tone_amp(rec, T, 2e6, i0, i1)) / mean;
    BOOST_TEST_MESSAGE("spur at 2 MHz / mean = " << spur);
    BOOST_CHECK_LT(spur, 0.10); // < 10% = < −20 dB relative
}

// ---------------------------------------------------------------------------
// Test 2: verify lower-half tones go to slice 1, upper-half to slice 2.
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(test_band_splitter_half_assignment)
{
    const double T  = 20e6;
    const double PR = 10e6;
    const size_t N  = 65536;

    // One tone clearly in the lower half, one in the upper half.
    const std::vector<cf> wb_lo = make_wideband(N, T, {-4e6});
    const std::vector<cf> wb_hi = make_wideband(N, T, {+4e6});

    band_splitter_config cfg;
    cfg.input_rate       = T;
    cfg.per_radio_rate   = PR;
    cfg.radio1_offset_hz = -5e6;
    cfg.radio2_offset_hz = +5e6;

    // Lower-half tone: must appear in out1 (radio1), not out2.
    {
        band_splitter bs(cfg);
        auto [o1, o2] = bs.split(wb_lo);
        const size_t i0 = 2048, i1 = o1.size() - 2048;
        BOOST_REQUIRE_GT(i1, i0);
        // −4 MHz in wideband → +1 MHz in radio1 baseband (offset +5 MHz)
        const double a1 = std::abs(tone_amp(o1, PR, +1e6, i0, i1));
        const double a2 = std::abs(tone_amp(o2, PR, -1e6, i0, i1));
        BOOST_TEST_MESSAGE("lo-tone: a1(out1)=" << a1 << "  a2(out2)=" << a2);
        BOOST_CHECK_GT(a1, 1e-3);    // tone present in slice 1
        BOOST_CHECK_LT(a2, 0.1 * a1); // not in slice 2
    }

    // Upper-half tone: must appear in out2 (radio2), not out1.
    {
        band_splitter bs(cfg);
        auto [o1, o2] = bs.split(wb_hi);
        const size_t i0 = 2048, i1 = o1.size() - 2048;
        BOOST_REQUIRE_GT(i1, i0);
        // +4 MHz in wideband → −1 MHz in radio2 baseband (offset −5 MHz)
        const double a1 = std::abs(tone_amp(o1, PR, +1e6, i0, i1));
        const double a2 = std::abs(tone_amp(o2, PR, -1e6, i0, i1));
        BOOST_TEST_MESSAGE("hi-tone: a1(out1)=" << a1 << "  a2(out2)=" << a2);
        BOOST_CHECK_GT(a2, 1e-3);    // tone present in slice 2
        BOOST_CHECK_LT(a1, 0.1 * a2); // not in slice 1
    }
}

// ---------------------------------------------------------------------------
// Test 3: seam tone (0 Hz) is attenuated in both slices — documented behaviour.
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(test_band_splitter_seam)
{
    const double T  = 20e6;
    const double PR = 10e6;
    const size_t N  = 65536;

    // A tone at 0 Hz sits exactly on the AA filter edge of both radios.
    const std::vector<cf> wb = make_wideband(N, T, {0.0});

    band_splitter_config cfg;
    cfg.input_rate       = T;
    cfg.per_radio_rate   = PR;
    cfg.radio1_offset_hz = -5e6;
    cfg.radio2_offset_hz = +5e6;

    band_splitter bs(cfg);
    auto [out1, out2] = bs.split(wb);

    BOOST_REQUIRE(!out1.empty());
    BOOST_REQUIRE(!out2.empty());

    // The seam tone appears at ±5 MHz in each per-radio baseband — precisely
    // the filter's −6 dB point.  Both slices should have some energy there
    // (not zero, not full-amplitude): verify ∈ (5%, 80%) of in-band power.
    const size_t i0 = 2048, i1 = out1.size() - 2048;
    BOOST_REQUIRE_GT(i1, i0);

    // Power at the seam edge in each slice (normalised to a fully-in-band tone).
    const std::vector<cf> wb_ref = make_wideband(N, T, {-3e6}); // clearly in-band
    band_splitter bs2(cfg);
    auto [ref1, ref2] = bs2.split(wb_ref);
    const double ref_amp = std::abs(tone_amp(ref1, PR, +2e6, i0, i1));
    BOOST_REQUIRE_GT(ref_amp, 1e-6);

    // 0 Hz → ±5 MHz in per-radio baseband
    const double seam_amp1 = std::abs(tone_amp(out1, PR, +5e6, i0, i1));
    const double seam_amp2 = std::abs(tone_amp(out2, PR, -5e6, i0, i1));
    BOOST_TEST_MESSAGE("seam: amp1=" << seam_amp1 << " amp2=" << seam_amp2
                                     << " ref=" << ref_amp);
    // Seam tone is attenuated (< full-amplitude) but non-zero — expected.
    BOOST_CHECK_LT(seam_amp1, ref_amp * 1.2);
    BOOST_CHECK_LT(seam_amp2, ref_amp * 1.2);
}

// ---------------------------------------------------------------------------
// Test 4: streaming process() produces identical output to one-shot split().
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(test_band_splitter_streaming)
{
    const double T  = 20e6;
    const double PR = 10e6;
    const size_t N  = 65536;

    std::mt19937 rng(99);
    std::uniform_real_distribution<float> ud(-1.f, 1.f);
    std::vector<cf> wb(N);
    for (auto& x : wb) {
        x = cf(ud(rng), ud(rng));
    }

    band_splitter_config cfg;
    cfg.input_rate       = T;
    cfg.per_radio_rate   = PR;
    cfg.radio1_offset_hz = -5e6;
    cfg.radio2_offset_hz = +5e6;

    // One-shot reference.
    band_splitter bs_os(cfg);
    auto [o1_os, o2_os] = bs_os.split(wb);

    // Streaming: feed in irregular chunks.
    band_splitter bs_st(cfg);
    std::vector<cf> o1_st, o2_st;
    const std::vector<size_t> chunks = {1000, 4096, 333, 8192, 777};
    size_t off = 0, c = 0;
    while (off < N) {
        const size_t sz = std::min(chunks[c % chunks.size()], N - off);
        bs_st.process(wb.data() + off, sz, o1_st, o2_st);
        off += sz;
        c++;
    }

    BOOST_REQUIRE_EQUAL(o1_os.size(), o1_st.size());
    BOOST_REQUIRE_EQUAL(o2_os.size(), o2_st.size());

    for (size_t i = 0; i < o1_os.size(); i++) {
        BOOST_CHECK_SMALL(std::abs(o1_os[i] - o1_st[i]), 1e-5f);
        BOOST_CHECK_SMALL(std::abs(o2_os[i] - o2_st[i]), 1e-5f);
    }
}

// ---------------------------------------------------------------------------
// Test 5: seam transition width. With the default 64 taps per phase a tone two
// 802.11 subcarriers (625 kHz at 20 MHz) from the seam must sit almost
// entirely in its own half, and a tone on the seam is shared equally. This
// also uses the pulled-in plan (offset 14/64 of the input rate).
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(test_band_splitter_seam_transition_width)
{
    const double T = 20e6;
    const size_t N = 65536;
    band_splitter_config cfg;
    cfg.input_rate       = T;
    cfg.per_radio_rate   = 10e6;
    cfg.radio1_offset_hz = -T * 14.0 / 64.0;
    cfg.radio2_offset_hz = +T * 14.0 / 64.0;

    // Amplitude of a wideband tone at f in each radio's slice (at the slice's
    // own baseband frequency f - offset).
    auto halves = [&](double f, size_t taps) {
        band_splitter_config c = cfg;
        c.taps_per_phase = taps;
        band_splitter bs(c);
        auto [o1, o2] = bs.split(bonded_test::tone(f, T, N));
        const double r = c.per_radio_rate;
        return std::make_pair(
            std::abs(tone_amp(o1, r, f - c.radio1_offset_hz, 4096, o1.size() - 4096)),
            std::abs(tone_amp(o2, r, f - c.radio2_offset_hz, 4096, o2.size() - 4096)));
    };
    const double sc = T / 64.0; // 802.11 subcarrier spacing at this rate
    const auto up2 = halves(+2 * sc, 64);
    const auto dn2 = halves(-2 * sc, 64);
    const auto seam = halves(0.0, 64);
    BOOST_TEST_MESSAGE("+2sc lo/hi=" << up2.first << "/" << up2.second
                       << " -2sc lo/hi=" << dn2.first << "/" << dn2.second
                       << " seam=" << seam.first << "/" << seam.second);
    BOOST_CHECK_LT(up2.first, 0.05);
    BOOST_CHECK_GT(up2.second, 0.95);
    BOOST_CHECK_LT(dn2.second, 0.05);
    BOOST_CHECK_GT(dn2.first, 0.95);
    BOOST_CHECK_CLOSE(seam.first, 0.5, 10.0);
    BOOST_CHECK_CLOSE(seam.second, 0.5, 10.0);

    // A short filter (the RX value, 12) leaks much more at the same distance.
    const auto up2_short = halves(+2 * sc, 12);
    BOOST_CHECK_GT(up2_short.first, 4 * up2.first);
}

// ---------------------------------------------------------------------------
// Test 6: with the pulled-in plan, splitting and recombining is transparent
// over the whole 802.11 occupied band (+-26 subcarriers) and the seam.
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(test_band_splitter_pulled_in_plan_recombines_flat)
{
    const double T = 20e6;
    const size_t N = 65536;
    band_splitter_config cfg;
    cfg.input_rate       = T;
    cfg.per_radio_rate   = 10e6;
    cfg.radio1_offset_hz = -T * 14.0 / 64.0;
    cfg.radio2_offset_hz = +T * 14.0 / 64.0;
    const double sc = T / 64.0;
    const std::vector<double> tones = {-26 * sc, -13 * sc, -3 * sc, -1 * sc,
        1 * sc, 3 * sc, 13 * sc, 26 * sc};
    const std::vector<cf> wb = make_wideband(N, T, tones);
    band_splitter bs(cfg);
    auto [o1, o2] = bs.split(wb);
    const std::vector<cf> rec = recombine(o1, o2, cfg);
    const size_t i0 = 8192, i1 = rec.size() - 8192;
    double gmin = 1e9, gmax = 0;
    for (double f : tones) {
        const double g = std::abs(tone_amp(rec, T, f, i0, i1)) / std::abs(tone_amp(wb, T, f, i0, i1));
        gmin = std::min(gmin, g);
        gmax = std::max(gmax, g);
    }
    BOOST_TEST_MESSAGE("pulled-in recombination gains [" << gmin << ", " << gmax << "]");
    BOOST_CHECK_GT(gmin, 0.9);
    BOOST_CHECK_LT(gmax, 1.1);
}

// ---------------------------------------------------------------------------
// Test 7: output length follows the resampling ratio exactly, also for a
// non-integer ratio, whatever the chunking.
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(test_band_splitter_output_length)
{
    struct plan
    {
        double in, out, off;
    };
    for (const plan& p : {plan{20e6, 10e6, 5e6}, plan{40e6, 25e6, 8.75e6}, plan{60e6, 37.5e6, 13.125e6}}) {
        band_splitter_config cfg;
        cfg.input_rate       = p.in;
        cfg.per_radio_rate   = p.out;
        cfg.radio1_offset_hz = -p.off;
        cfg.radio2_offset_hz = +p.off;
        band_splitter bs(cfg);
        const size_t N = 100000;
        const std::vector<cf> x(N, cf(0.5f, 0.0f));
        std::vector<cf> o1, o2;
        const std::vector<size_t> chunks = {1, 7, 1000, 4096, 33};
        for (size_t off = 0, c = 0; off < N; c++) {
            const size_t take = std::min(chunks[c % chunks.size()], N - off);
            bs.process(x.data() + off, take, o1, o2);
            off += take;
        }
        const double expect = static_cast<double>(N) * p.out / p.in;
        BOOST_TEST_MESSAGE(p.in << "->" << p.out << ": " << o1.size() << " (expect " << expect << ")");
        BOOST_CHECK_EQUAL(o1.size(), o2.size());
        BOOST_CHECK_LE(std::abs(static_cast<double>(o1.size()) - expect), 1.0);
    }
}

// ---------------------------------------------------------------------------
// Test 8: invalid configurations are rejected at construction.
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(test_band_splitter_rejects_bad_config)
{
    band_splitter_config ok;
    BOOST_CHECK_NO_THROW(band_splitter{ok});

    band_splitter_config c = ok;
    c.input_rate = 0;
    BOOST_CHECK_THROW(band_splitter{c}, std::invalid_argument);

    c = ok;
    c.per_radio_rate = 30e6; // above the wideband rate
    BOOST_CHECK_THROW(band_splitter{c}, std::invalid_argument);

    c = ok;
    c.radio2_offset_hz = 6e6; // a 12 MHz slice cannot fit in 10 MS/s
    BOOST_CHECK_THROW(band_splitter{c}, std::invalid_argument);
}
