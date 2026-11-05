//
// Copyright 2026 Gabriel Lando
// SPDX-License-Identifier: GPL-3.0-or-later
//
// No-hardware unit tests for bonded_transmitter: configuration checks, the
// live phase-offset API, guards against use before configure(), and the
// per-device frequency/phase correction. configure() (which opens devices)
// is never called.
//

#include "bonded_transmitter.hpp"
#include <boost/test/unit_test.hpp>
#include <uhd/exception.hpp>
#include <cmath>
#include <complex>
#include <random>
#include <stdexcept>
#include <vector>

using namespace bonded;

namespace {

using cf = std::complex<float>;
constexpr double kPi = 3.14159265358979323846;

bonded_transmitter::config make_cfg()
{
    bonded_transmitter::config cfg;
    cfg.serials      = {"AAA", "BBB"};
    cfg.freq         = 2.412e9;
    cfg.rate         = 10e6;
    cfg.clock_source = "internal";
    cfg.time_source  = "internal";
    return cfg;
}

std::vector<cf> random_samples(size_t n, unsigned seed)
{
    std::mt19937 rng(seed);
    std::uniform_real_distribution<float> u(-1.f, 1.f);
    std::vector<cf> x(n);
    for (auto& v : x) {
        v = cf(u(rng), u(rng));
    }
    return x;
}

} // namespace

BOOST_AUTO_TEST_CASE(test_bonded_transmitter_rejects_bad_delay_trim)
{
    auto cfg              = make_cfg();
    cfg.sample_delay_trim = {0, 3, 7}; // three entries, two serials
    BOOST_CHECK_THROW(bonded_transmitter tx(cfg), std::invalid_argument);
    cfg.sample_delay_trim = {0, 5};
    BOOST_CHECK_NO_THROW(bonded_transmitter tx(cfg));
}

BOOST_AUTO_TEST_CASE(test_bonded_transmitter_initial_state)
{
    bonded_transmitter tx(make_cfg());
    BOOST_CHECK_EQUAL(tx.num_devices(), size_t(2));
    BOOST_CHECK(!tx.is_configured());
    BOOST_CHECK(!tx.is_streaming());
    BOOST_CHECK_EQUAL(tx.frames_sent(), 0u);
    BOOST_CHECK(tx.get_underrun_counts().empty()); // no devices opened yet
}

BOOST_AUTO_TEST_CASE(test_bonded_transmitter_phase_offset_api)
{
    auto cfg             = make_cfg();
    cfg.phase_offset_rad = {0.25};
    bonded_transmitter tx(cfg);
    BOOST_CHECK_EQUAL(tx.get_phase_offset(0), 0.25); // from the config
    BOOST_CHECK_EQUAL(tx.get_phase_offset(1), 0.0);  // missing entry -> 0
    tx.set_phase_offset(1, -1.5);
    BOOST_CHECK_EQUAL(tx.get_phase_offset(1), -1.5);
    BOOST_CHECK_THROW(tx.set_phase_offset(2, 0.1), std::out_of_range);
    BOOST_CHECK_EQUAL(tx.get_phase_offset(7), 0.0);
}

BOOST_AUTO_TEST_CASE(test_bonded_transmitter_requires_configure)
{
    bonded_transmitter tx(make_cfg());
    BOOST_CHECK_THROW(tx.start_continuous(), uhd::runtime_error);
    const cf z(0, 0);
    BOOST_CHECK_THROW(tx.queue_samples({&z, &z}, 1), uhd::runtime_error);
    BOOST_CHECK_NO_THROW(tx.stop()); // no-op when not streaming
}

BOOST_AUTO_TEST_CASE(test_phase_correction_identity_and_constant_phase)
{
    const std::vector<cf> x = random_samples(1000, 1);
    std::vector<cf> y = x;
    detail::apply_phase_correction(y, 0.0, 0.0, 10e6, 12345);
    BOOST_CHECK(y == x);

    const double ph = 0.7;
    detail::apply_phase_correction(y, 0.0, ph, 10e6, 0);
    const cf rot(static_cast<float>(std::cos(ph)), static_cast<float>(std::sin(ph)));
    for (size_t i = 0; i < x.size(); i++) {
        BOOST_REQUIRE_SMALL(std::abs(y[i] - x[i] * rot), 1e-6f);
    }
}

BOOST_AUTO_TEST_CASE(test_phase_correction_frequency_is_continuous_across_chunks)
{
    // A frequency correction applied chunk by chunk, with the running sample
    // index as reference, must equal the correction of the whole stream.
    const double rate = 10e6, f = 12.5, ph = -0.3;
    const std::vector<cf> x = random_samples(30000, 2);

    std::vector<cf> whole = x;
    detail::apply_phase_correction(whole, f, ph, rate, 5000);

    std::vector<cf> chunked;
    uint64_t index = 5000;
    const std::vector<size_t> sizes = {1, 4, 13, 40, 1000, 4096};
    for (size_t off = 0, c = 0; off < x.size(); c++) {
        const size_t n = std::min(sizes[c % sizes.size()], x.size() - off);
        std::vector<cf> part(x.begin() + static_cast<std::ptrdiff_t>(off),
            x.begin() + static_cast<std::ptrdiff_t>(off + n));
        detail::apply_phase_correction(part, f, ph, rate, index);
        chunked.insert(chunked.end(), part.begin(), part.end());
        index += n;
        off += n;
    }
    BOOST_REQUIRE_EQUAL(chunked.size(), whole.size());
    for (size_t i = 0; i < whole.size(); i++) {
        BOOST_REQUIRE_SMALL(std::abs(chunked[i] - whole[i]), 1e-6f);
    }

    // And it is the expected rotation: sample k of the stream turns by
    // 2*pi*f*k/rate + ph.
    const size_t k = 777;
    const double a = 2.0 * kPi * f * static_cast<double>(5000 + k) / rate + ph;
    const cf expect = x[k] * cf(static_cast<float>(std::cos(a)), static_cast<float>(std::sin(a)));
    BOOST_CHECK_SMALL(std::abs(whole[k] - expect), 1e-6f);
}
