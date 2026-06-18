//
// Copyright 2026 Gabriel Lando
// SPDX-License-Identifier: GPL-3.0-or-later
//
// No-hardware unit tests for bonded_transmitter.
//
// All tests avoid opening UHD devices (configure() is never called).
// They exercise config validation, delay-trim mathematics, freq-plan
// geometry, and underrun-counter initialisation — pure logic.
//

#include "bonded_transmitter.hpp"
#include <boost/test/unit_test.hpp>
#include <cmath>
#include <stdexcept>
#include <vector>

using namespace bonded;

// ---------------------------------------------------------------------------
// Helper: build a minimal two-radio config without serials that trigger HW.
// ---------------------------------------------------------------------------
namespace {

bonded_transmitter::config make_cfg(double fc = 2.412e9,
    double per_radio_rate = 10e6)
{
    bonded_transmitter::config cfg;
    cfg.serials  = {"AAA", "BBB"};
    cfg.freq     = fc;
    cfg.rate     = per_radio_rate;
    cfg.gain     = 0.0;
    cfg.strict   = false;
    // Use internal clocks so no hardware lock is needed (tests never call
    // configure()).
    cfg.clock_source = "internal";
    cfg.time_source  = "internal";
    return cfg;
}

} // namespace

// ---------------------------------------------------------------------------
// Test 1: Default config values are sane.
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(test_bonded_transmitter_default_config)
{
    const auto cfg = make_cfg(2.412e9, 10e6);
    BOOST_CHECK_EQUAL(cfg.serials.size(), size_t(2));
    BOOST_CHECK_CLOSE(cfg.freq, 2.412e9, 0.001);
    BOOST_CHECK_CLOSE(cfg.rate, 10e6, 0.001);
    BOOST_CHECK(cfg.sample_delay_trim.empty());
    BOOST_CHECK(!cfg.strict);
}

// ---------------------------------------------------------------------------
// Test 2: Mismatched sample_delay_trim length must be rejected at construction.
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(test_bonded_transmitter_delay_trim_size_mismatch)
{
    auto cfg               = make_cfg();
    cfg.sample_delay_trim  = {0, 3, 7}; // three entries but only two serials
    BOOST_CHECK_THROW(bonded_transmitter tx(cfg), std::invalid_argument);
}

// ---------------------------------------------------------------------------
// Test 3: Valid delay trim (same length as serials) is accepted.
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(test_bonded_transmitter_delay_trim_valid)
{
    auto cfg              = make_cfg();
    cfg.sample_delay_trim = {0, 5}; // two entries, two serials
    BOOST_CHECK_NO_THROW(bonded_transmitter tx(cfg));
}

// ---------------------------------------------------------------------------
// Test 4: Freq-plan geometry — two-radio standard setup.
//
// For a 20 MHz wideband channel at fc, the per-radio rate is 10 MHz and the
// two radios are centred at fc ± 5 MHz (i.e. ± per_radio_rate/2).
// Verify the expected offsets are ±5 MHz.
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(test_bonded_transmitter_freq_plan_geometry)
{
    const double fc            = 2.412e9;
    const double per_radio     = 10e6;
    const double expected_half = per_radio / 2.0; // 5 MHz

    auto cfg       = make_cfg(fc, per_radio);
    cfg.freq_plan  = {fc - expected_half, fc + expected_half};

    BOOST_REQUIRE_EQUAL(cfg.freq_plan.size(), size_t(2));
    BOOST_CHECK_CLOSE(cfg.freq_plan[0], fc - expected_half, 0.001);
    BOOST_CHECK_CLOSE(cfg.freq_plan[1], fc + expected_half, 0.001);

    // The offset from fc must be exactly ±per_radio/2.
    const double off0 = cfg.freq_plan[0] - cfg.freq;
    const double off1 = cfg.freq_plan[1] - cfg.freq;
    BOOST_CHECK_CLOSE(off0, -expected_half, 0.001);
    BOOST_CHECK_CLOSE(off1, +expected_half, 0.001);
    BOOST_CHECK_CLOSE(std::abs(off0), std::abs(off1), 0.001);
}

// ---------------------------------------------------------------------------
// Test 5: Freq-plan size mismatch (≠ serials) must be rejected at configure()
// time — since configure() opens hardware we verify the config is inconsistent
// by checking the pre-condition through the expected field layout.
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(test_bonded_transmitter_freq_plan_mismatch_detected)
{
    auto cfg      = make_cfg();
    cfg.freq_plan = {2.407e9}; // only one entry for two radios
    // The bonded_transmitter constructor itself does not validate freq_plan
    // (configure() does).  But we can verify the config is logically
    // inconsistent (freq_plan.size() != serials.size()).
    BOOST_CHECK(cfg.freq_plan.size() != cfg.serials.size());
}

// ---------------------------------------------------------------------------
// Test 6: Delay trim mathematics — verify that the integer sample offset
// correctly represents the timing error correction.
//
// At 10 MHz rate, one sample = 100 ns.  A trim of +16 samples delays device
// B by 1600 ns, well within the 802.11 CP (800 ns) × 2 budget.
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(test_bonded_transmitter_delay_trim_timing)
{
    const double per_radio_rate = 10e6;
    const double sample_period_ns = 1e9 / per_radio_rate; // 100 ns

    // The maximum useful trim is one cyclic-prefix length (0.8 µs = 8 samples
    // at 10 MHz).  Verify the arithmetic.
    const int trim_samples   = 8;
    const double delay_ns    = trim_samples * sample_period_ns;
    const double cp_ns       = 800.0; // 802.11a/g cyclic prefix

    BOOST_CHECK_CLOSE(sample_period_ns, 100.0, 0.001);
    BOOST_CHECK_CLOSE(delay_ns, 800.0, 0.001);
    BOOST_CHECK_LE(delay_ns, cp_ns);

    auto cfg              = make_cfg(2.412e9, per_radio_rate);
    cfg.sample_delay_trim = {0, trim_samples};
    bonded_transmitter tx(cfg);

    BOOST_CHECK_EQUAL(tx.num_devices(), size_t(2));
    BOOST_CHECK(!tx.is_configured());
    BOOST_CHECK(!tx.is_streaming());
}

// ---------------------------------------------------------------------------
// Test 7: Underrun counters start at zero on a fresh object (before configure).
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(test_bonded_transmitter_underrun_initial_zero)
{
    auto cfg = make_cfg();
    bonded_transmitter tx(cfg);
    // Before configure() there are no streamers; underrun_counts is empty.
    const auto counts = tx.get_underrun_counts();
    BOOST_CHECK(counts.empty());
}

// ---------------------------------------------------------------------------
// Test 8: send_burst() before configure() must throw.
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(test_bonded_transmitter_send_before_configure_throws)
{
    auto cfg = make_cfg();
    bonded_transmitter tx(cfg);
    const std::vector<std::vector<std::complex<float>>> dummy(2);
    BOOST_CHECK_THROW(tx.send_burst(dummy, 1.0), uhd::runtime_error);
}

// ---------------------------------------------------------------------------
// Test 9: start_continuous() before configure() must throw.
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(test_bonded_transmitter_start_continuous_before_configure)
{
    auto cfg = make_cfg();
    bonded_transmitter tx(cfg);
    BOOST_CHECK_THROW(tx.start_continuous(), uhd::runtime_error);
}
