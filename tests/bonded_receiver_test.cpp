//
// Copyright 2026 Gabriel Lando
// SPDX-License-Identifier: GPL-3.0-or-later
//

#include "bonded_alignment.hpp"
#include <boost/test/unit_test.hpp>
#include <complex>
#include <deque>
#include <vector>

namespace {

struct test_chunk
{
    double timestamp;
    size_t num_samples;
    size_t offset = 0;
    std::vector<std::vector<std::complex<float>>> data;
};

test_chunk make_chunk(const double ts,
    const size_t start,
    const size_t nsamps,
    const size_t nch = 1)
{
    test_chunk c;
    c.timestamp   = ts;
    c.num_samples = nsamps;
    c.data.resize(nch);
    for (size_t ch = 0; ch < nch; ch++) {
        c.data[ch].resize(nsamps);
        for (size_t i = 0; i < nsamps; i++) {
            const float v = static_cast<float>(start + i + 1000 * ch);
            c.data[ch][i] = std::complex<float>(v, -v);
        }
    }
    return c;
}

} // namespace

using namespace bonded;

BOOST_AUTO_TEST_CASE(test_burst_alignment_jitter_tolerated)
{
    const std::vector<double> first_ts = {1.0000000, 1.0000001, 0.9999999};
    const std::vector<bool> success    = {true, true, true};

    const auto r = detail::evaluate_burst_alignment(first_ts, success, 1e-3);
    BOOST_CHECK(r.aligned);
    BOOST_CHECK(r.error_message.empty());
    BOOST_CHECK_LT(r.inter_device_spread_us, 1.0);
}

BOOST_AUTO_TEST_CASE(test_burst_alignment_missing_device)
{
    const std::vector<double> first_ts = {1.0, -1.0, 1.0};
    const std::vector<bool> success    = {true, false, true};

    const auto r = detail::evaluate_burst_alignment(first_ts, success, 1e-3);
    BOOST_CHECK(!r.aligned);
    BOOST_CHECK(!r.error_message.empty());
}

BOOST_AUTO_TEST_CASE(test_find_alignment_timestamp)
{
    std::vector<std::deque<test_chunk>> rings(3);
    rings[0].push_back(make_chunk(10.0000000, 0, 8));
    rings[1].push_back(make_chunk(10.0000002, 0, 8));
    rings[2].push_back(make_chunk(10.0000001, 0, 8));

    double align_ts = -1.0;
    std::string err;
    const bool ok = detail::find_alignment_timestamp(rings, align_ts, err);
    BOOST_CHECK(ok);
    BOOST_CHECK(err.empty());
    BOOST_CHECK_CLOSE(align_ts, 10.0000002, 0.00001);
}

BOOST_AUTO_TEST_CASE(test_extract_aligned_samples_skip_and_copy)
{
    std::deque<test_chunk> ring;
    ring.push_back(make_chunk(10.0, 0, 8, 1));

    std::vector<std::vector<std::complex<float>>> out;
    // At 10 Msps, 0.35 us => skip 3 samples after truncation.
    // Avoid exact FP boundaries (e.g., 0.3 us) that can evaluate to 2.999999...
    // and truncate one sample early.
    const auto r = detail::extract_aligned_samples(
        ring, 4, 1, 10e6, 10.0 + 3.5e-7, out);

    BOOST_CHECK(r.success);
    BOOST_CHECK(!r.zero_filled);
    BOOST_CHECK_EQUAL(out.size(), size_t(1));
    BOOST_CHECK_EQUAL(out[0].size(), size_t(4));
    BOOST_CHECK_EQUAL(out[0][0].real(), 3.0f);
    BOOST_CHECK_EQUAL(out[0][1].real(), 4.0f);
    BOOST_CHECK_EQUAL(out[0][2].real(), 5.0f);
    BOOST_CHECK_EQUAL(out[0][3].real(), 6.0f);
}

BOOST_AUTO_TEST_CASE(test_extract_aligned_samples_short_zero_fill)
{
    std::deque<test_chunk> ring;
    ring.push_back(make_chunk(20.0, 0, 2, 1));

    std::vector<std::vector<std::complex<float>>> out;
    const auto r = detail::extract_aligned_samples(ring, 5, 1, 10e6, 20.0, out);

    BOOST_CHECK(!r.success);
    BOOST_CHECK(r.zero_filled);
    BOOST_CHECK_EQUAL(out.size(), size_t(1));
    BOOST_CHECK_EQUAL(out[0].size(), size_t(5));
    BOOST_CHECK_EQUAL(out[0][0].real(), 0.0f);
    BOOST_CHECK_EQUAL(out[0][1].real(), 1.0f);
    BOOST_CHECK_EQUAL(out[0][2].real(), 0.0f);
    BOOST_CHECK_EQUAL(out[0][3].real(), 0.0f);
    BOOST_CHECK_EQUAL(out[0][4].real(), 0.0f);
}
