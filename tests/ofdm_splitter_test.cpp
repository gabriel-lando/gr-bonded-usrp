//
// Copyright 2026 Gabriel Lando
//
// No-hardware validation of the ofdm_splitter engine. Random CP-OFDM bursts
// are split, each radio's output is put back on the wideband grid (ideal
// interpolation + translation to its centre, as the radio's DUC and LO do)
// and the two are summed with an arbitrary relative phase, as two radios that
// are not phase-coherent would arrive at a receiver. Checks:
//   (a) every subcarrier is recovered cleanly at any relative phase;
//   (b) each subcarrier goes to the radio on its side, the middle one to both;
//   (c) the radio-rate timeline: silence keeps the rate, bursts land where
//       the input put them, a burst is released without waiting for input;
//   (d) chunked streaming gives the same output as one call;
//   (e) invalid configurations are rejected.
//

#include "ofdm_splitter.hpp"
#include "test_helpers.hpp"

#include <boost/test/unit_test.hpp>
#include <algorithm>
#include <cmath>
#include <complex>
#include <random>
#include <stdexcept>
#include <vector>

using namespace bonded;

namespace {

using bonded_test::cd;
using bonded_test::cf;
using bonded_test::kPi;

constexpr size_t kN = 64, kCP = 16, kS = kN + kCP;

//! Bins of FFT-ordered spectrum X for subcarrier k.
size_t bin(long k) { return static_cast<size_t>((k % 64 + 64) % 64); }

struct burst_signal
{
    std::vector<cf> x;                 // time samples, kS per symbol
    std::vector<std::vector<cd>> X;    // per-symbol spectrum (FFT order)
};

//! Back-to-back CP-OFDM symbols with random QPSK on subcarriers +-1..+-26.
burst_signal make_burst(size_t nsym, unsigned seed)
{
    std::mt19937 rng(seed);
    std::uniform_int_distribution<int> bit(0, 1);
    burst_signal b;
    for (size_t s = 0; s < nsym; s++) {
        std::vector<cd> X(kN, cd(0, 0));
        for (long k = -26; k <= 26; k++) {
            if (k != 0) {
                X[bin(k)] = cd(bit(rng) ? 1 : -1, bit(rng) ? 1 : -1);
            }
        }
        std::vector<cf> body(kN);
        for (size_t n = 0; n < kN; n++) {
            cd acc(0, 0);
            for (size_t j = 0; j < kN; j++) {
                const double a = 2.0 * kPi * static_cast<double>(j * n) / kN;
                acc += X[j] * cd(std::cos(a), std::sin(a));
            }
            body[n] = static_cast<cf>(acc / static_cast<double>(kN) * 0.3);
        }
        b.x.insert(b.x.end(), body.end() - kCP, body.end());
        b.x.insert(b.x.end(), body.begin(), body.end());
        b.X.push_back(X);
    }
    return b;
}

//! What a receiver at fc sees: each radio interpolated back to input_rate,
//! translated to its centre, radio 2 rotated by phi.
std::vector<cf> recombine(const std::vector<cf>& o1, const std::vector<cf>& o2,
    const ofdm_splitter_config& c, double phi, size_t up, size_t down)
{
    const auto u1 = bonded_test::rational_resample(o1, up, down);
    const auto u2 = bonded_test::rational_resample(o2, up, down);
    const size_t n = std::min(u1.size(), u2.size());
    std::vector<cf> y(n);
    for (size_t i = 0; i < n; i++) {
        const double t  = static_cast<double>(i) / c.input_rate;
        const double a1 = 2.0 * kPi * c.radio1_offset_hz * t;
        const double a2 = 2.0 * kPi * c.radio2_offset_hz * t + phi;
        y[i] = static_cast<cf>(static_cast<cd>(u1[i]) * cd(std::cos(a1), std::sin(a1))
                               + static_cast<cd>(u2[i]) * cd(std::cos(a2), std::sin(a2)));
    }
    return y;
}

std::vector<cd> fft64(const std::vector<cf>& y, size_t at)
{
    std::vector<cd> Y(kN, cd(0, 0));
    for (size_t j = 0; j < kN; j++) {
        for (size_t n = 0; n < kN; n++) {
            const double a = -2.0 * kPi * static_cast<double>(j * n) / kN;
            Y[j] += static_cast<cd>(y[at + n]) * cd(std::cos(a), std::sin(a));
        }
    }
    return Y;
}

//! Worst per-subcarrier EVM (dB) of symbols 1.. of a burst received at
//! `start`, equalized on symbol 0 (as an LTF would), FFT window `skip`
//! samples into each symbol.
double worst_evm_db(const std::vector<cf>& y, const burst_signal& b, size_t start, size_t skip)
{
    const std::vector<cd> Y0 = fft64(y, start + skip);
    std::vector<double> err(kN, 0.0), pw(kN, 0.0);
    for (size_t s = 1; s < b.X.size(); s++) {
        const std::vector<cd> Y = fft64(y, start + s * kS + skip);
        for (long k = -26; k <= 26; k++) {
            if (k == 0) continue;
            const size_t j = bin(k);
            const cd H     = Y0[j] / b.X[0][j];
            err[j] += std::norm(Y[j] / H - b.X[s][j]);
            pw[j] += std::norm(b.X[s][j]);
        }
    }
    double worst = -1e9;
    for (long k = -26; k <= 26; k++) {
        if (k != 0) {
            worst = std::max(worst, 10.0 * std::log10(err[bin(k)] / pw[bin(k)]));
        }
    }
    return worst;
}

struct plan
{
    double per_radio_rate, offset_hz;
    size_t up, down; // input_rate / per_radio_rate
};

const std::vector<plan> kPlans = {
    { 12.5e6, 3.75e6, 8, 5 },   // 5/8 rate, 12 subcarriers (default)
    { 10e6, 4.0625e6, 2, 1 },   // half rate, 13 subcarriers
};

ofdm_splitter_config make_cfg(const plan& p)
{
    ofdm_splitter_config c;
    c.input_rate       = 20e6;
    c.per_radio_rate   = p.per_radio_rate;
    c.radio1_offset_hz = -p.offset_hz;
    c.radio2_offset_hz = +p.offset_hz;
    return c;
}

//! One burst of `nsym` symbols at input offset `lead`, then `tail` zeros.
struct scenario
{
    burst_signal b;
    std::vector<cf> in;
    size_t lead;
};
scenario make_scenario(size_t nsym, size_t lead, size_t tail, unsigned seed)
{
    scenario sc{make_burst(nsym, seed), {}, lead};
    sc.in.assign(lead, cf(0, 0));
    sc.in.insert(sc.in.end(), sc.b.x.begin(), sc.b.x.end());
    sc.in.resize(sc.in.size() + tail, cf(0, 0));
    return sc;
}

} // namespace

BOOST_AUTO_TEST_CASE(test_ofdm_split_any_relative_phase)
{
    // The point of the block: the receiver sees the radios' relative phase as
    // a per-half channel phase, which its equalizer removes, so dense
    // constellations survive at any phase. (The band_splitter filter split
    // reaches only about -8 dB on the subcarriers next to DC at 180 deg.)
    for (const auto& p : kPlans) {
        const auto cfg = make_cfg(p);
        const scenario sc = make_scenario(12, 800, 800, 7);
        for (double phi_deg : {0.0, 90.0, 180.0, -120.0}) {
            ofdm_splitter spl(cfg);
            std::vector<cf> o1, o2;
            spl.process(sc.in.data(), sc.in.size(), {{sc.lead, sc.b.x.size()}}, o1, o2);
            const auto y = recombine(o1, o2, cfg, phi_deg * kPi / 180.0, p.up, p.down);
            // rational_resample delays by its half length: 10*max(up,down)/down.
            const size_t delay = 10 * std::max(p.up, p.down) / p.down;
            const double evm   = worst_evm_db(y, sc.b, sc.lead + delay, 12);
            BOOST_TEST_MESSAGE("plan " << p.per_radio_rate / 1e6 << " MS/s, phi " << phi_deg
                                       << ": worst EVM " << evm << " dB");
            BOOST_CHECK_LT(evm, -30.0);
        }
    }
}

BOOST_AUTO_TEST_CASE(test_ofdm_split_assigns_subcarriers)
{
    const auto cfg = make_cfg(kPlans[0]);
    auto energy = [&](long k) {
        std::vector<cf> sym(kS);
        for (size_t n = 0; n < kS; n++) {
            const double a = 2.0 * kPi * static_cast<double>(k) * (static_cast<double>(n) - kCP) / kN;
            sym[n] = cf(static_cast<float>(std::cos(a)), static_cast<float>(std::sin(a)));
        }
        std::vector<cf> in = sym;
        in.insert(in.end(), sym.begin(), sym.end());
        ofdm_splitter spl(cfg);
        std::vector<cf> o1, o2;
        spl.process(in.data(), in.size(), {{0, in.size()}}, o1, o2);
        double e1 = 0, e2 = 0;
        for (const auto& v : o1) e1 += std::norm(v);
        for (const auto& v : o2) e2 += std::norm(v);
        return std::make_pair(e1, e2);
    };
    for (long k : {-26L, -5L, -1L}) {
        const auto e = energy(k);
        BOOST_CHECK_GT(e.first, 1.0);
        BOOST_CHECK_LT(e.second, 1e-9 * e.first);
    }
    for (long k : {1L, 5L, 26L}) {
        const auto e = energy(k);
        BOOST_CHECK_GT(e.second, 1.0);
        BOOST_CHECK_LT(e.first, 1e-9 * e.second);
    }
    // The subcarrier midway between the radio centres (DC here) is shared,
    // half the amplitude each (compared with a subcarrier that is, like DC,
    // continuous across the symbol boundary).
    const auto dc = energy(0);
    BOOST_CHECK_CLOSE(dc.first, dc.second, 1e-3);
    BOOST_CHECK_CLOSE(dc.first, energy(-4).first / 4.0, 1.0);
}

BOOST_AUTO_TEST_CASE(test_ofdm_split_timeline)
{
    const auto cfg = make_cfg(kPlans[0]);
    // Silence keeps the radio-rate timeline: n inputs -> round(n * 5/8).
    {
        ofdm_splitter spl(cfg);
        std::vector<cf> zeros(1001, cf(0, 0)), o1, o2;
        spl.process(zeros.data(), zeros.size(), {}, o1, o2);
        BOOST_CHECK_EQUAL(o1.size(), 626u);
        BOOST_CHECK_EQUAL(o2.size(), o1.size());
    }
    // A burst starts on the radio timeline where the input put it, and the
    // whole burst (its last suffix included) comes out with the input that
    // completes it: nothing waits for the next burst.
    {
        ofdm_splitter spl(cfg);
        const scenario sc = make_scenario(5, 1600, 0, 3);
        std::vector<cf> o1, o2;
        // 802.11 frames carry one extra sample (the window tail).
        std::vector<cf> in = sc.in;
        in.push_back(cf(0, 0));
        spl.process(in.data(), in.size(), {{sc.lead, sc.b.x.size() + 1}}, o1, o2);
        const size_t q0  = 1000; // 1600 * 5/8
        const size_t end = q0 + 5 * spl.radio_fft_len() + 5 * spl.radio_cp_len()
                           + spl.radio_window_len();
        BOOST_REQUIRE_GE(o1.size(), end);
        size_t first = 0;
        while (first < o1.size() && std::abs(o1[first]) == 0.0f) first++;
        BOOST_CHECK_EQUAL(first, q0);
        size_t last = o1.size();
        while (last > 0 && std::abs(o1[last - 1]) == 0.0f) last--;
        BOOST_CHECK_EQUAL(last, end);
    }
}

BOOST_AUTO_TEST_CASE(test_ofdm_split_streaming_matches_one_call)
{
    const auto cfg = make_cfg(kPlans[0]);
    // Two bursts, the second right behind the first (no gap), then silence.
    const burst_signal b1 = make_burst(4, 11), b2 = make_burst(3, 12);
    std::vector<cf> in(333, cf(0, 0));
    const size_t s1 = in.size();
    in.insert(in.end(), b1.x.begin(), b1.x.end());
    const size_t s2 = in.size();
    in.insert(in.end(), b2.x.begin(), b2.x.end());
    in.resize(in.size() + 500, cf(0, 0));
    const std::vector<std::pair<size_t, size_t>> bursts = {{s1, b1.x.size()}, {s2, b2.x.size()}};

    ofdm_splitter ref(cfg);
    std::vector<cf> r1, r2;
    std::vector<ofdm_splitter::burst> all;
    for (const auto& b : bursts) all.push_back({b.first, b.second});
    ref.process(in.data(), in.size(), all, r1, r2);

    for (size_t chunk : {1u, 7u, 79u, 80u, 1000u}) {
        ofdm_splitter spl(cfg);
        std::vector<cf> o1, o2;
        for (size_t i = 0; i < in.size(); i += chunk) {
            const size_t n = std::min(chunk, in.size() - i);
            std::vector<ofdm_splitter::burst> here;
            for (const auto& b : bursts) {
                if (b.first >= i && b.first < i + n) here.push_back({b.first - i, b.second});
            }
            spl.process(in.data() + i, n, here, o1, o2);
        }
        BOOST_REQUIRE_EQUAL(o1.size(), r1.size());
        BOOST_REQUIRE_EQUAL(o2.size(), r2.size());
        for (size_t i = 0; i < r1.size(); i++) {
            BOOST_REQUIRE_LT(std::abs(o1[i] - r1[i]), 1e-6f);
            BOOST_REQUIRE_LT(std::abs(o2[i] - r2[i]), 1e-6f);
        }
    }
}

BOOST_AUTO_TEST_CASE(test_ofdm_split_rejects_bad_config)
{
    auto bad = [](auto edit) {
        ofdm_splitter_config c = make_cfg(kPlans[0]);
        edit(c);
        BOOST_CHECK_THROW(ofdm_splitter{c}, std::invalid_argument);
    };
    bad([](ofdm_splitter_config& c) { c.radio1_offset_hz = -3.7e6; });  // not a subcarrier
    bad([](ofdm_splitter_config& c) { c.per_radio_rate = 12e6; });      // 64 * 12/20 not whole
    bad([](ofdm_splitter_config& c) { c.window_len = 17; });            // > cp_len
    bad([](ofdm_splitter_config& c) { std::swap(c.radio1_offset_hz, c.radio2_offset_hz); });
    bad([](ofdm_splitter_config& c) { c.per_radio_rate = 25e6; });      // > input_rate
    bad([](ofdm_splitter_config& c) { c.cp_len = 0; });
    BOOST_CHECK_NO_THROW(ofdm_splitter{make_cfg(kPlans[1])});
}
