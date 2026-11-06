//
// Copyright 2026 Gabriel Lando
//
// Ring-buffer alignment used by bonded_receiver (no hardware): chunks are
// placed on a common time grid by timestamp, holes stay holes, and the
// "removed" count keeps the caller's O(1) bookkeeping exact.
//

#include "bonded_alignment.hpp"
#include <boost/test/unit_test.hpp>
#include <complex>
#include <deque>
#include <random>
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

BOOST_AUTO_TEST_CASE(test_extract_aligned_samples_hole_keeps_time_grid)
{
    // After an overflow the ring has a hole: 10 samples at t0, then the next
    // chunk 15 samples later. The hole must stay a hole (zeros, flagged), not
    // be closed up — closing it would shift this device 5 samples early
    // relative to the other one for the rest of the stream.
    const double rate = 10e6;
    std::deque<test_chunk> ring;
    ring.push_back(make_chunk(30.0, 0, 10, 1));
    ring.push_back(make_chunk(30.0 + 15.0 / rate, 100, 10, 1));

    std::vector<std::vector<std::complex<float>>> out;
    auto r = detail::extract_aligned_samples(ring, 20, 1, rate, 30.0, out);
    BOOST_CHECK(!r.success);
    BOOST_CHECK(r.zero_filled);
    BOOST_CHECK_EQUAL(out[0][9].real(), 9.0f);
    BOOST_CHECK_EQUAL(out[0][10].real(), 0.0f);
    BOOST_CHECK_EQUAL(out[0][14].real(), 0.0f);
    BOOST_CHECK_EQUAL(out[0][15].real(), 100.0f);
    BOOST_CHECK_EQUAL(out[0][19].real(), 104.0f);

    // The next block continues exactly at t0 + 20 samples.
    BOOST_CHECK_CLOSE(detail::next_timestamp(ring.front(), rate), 30.0 + 20.0 / rate, 1e-9);
    r = detail::extract_aligned_samples(ring, 5, 1, rate, 30.0 + 20.0 / rate, out);
    BOOST_CHECK(r.success);
    BOOST_CHECK_EQUAL(out[0][0].real(), 105.0f);
    BOOST_CHECK_EQUAL(out[0][4].real(), 109.0f);
}

BOOST_AUTO_TEST_CASE(test_extract_aligned_samples_late_start)
{
    // Data that starts after the alignment time lands at its own time slot.
    const double rate = 10e6;
    std::deque<test_chunk> ring;
    ring.push_back(make_chunk(40.0 + 2.0 / rate, 0, 8, 1));

    std::vector<std::vector<std::complex<float>>> out;
    const auto r = detail::extract_aligned_samples(ring, 4, 1, rate, 40.0, out);
    BOOST_CHECK(!r.success);
    BOOST_CHECK_EQUAL(out[0][0].real(), 0.0f);
    BOOST_CHECK_EQUAL(out[0][1].real(), 0.0f);
    BOOST_CHECK_EQUAL(out[0][2].real(), 0.0f); // sample value 0 at slot 2
    BOOST_CHECK_EQUAL(out[0][3].real(), 1.0f);
}

BOOST_AUTO_TEST_CASE(test_next_timestamp_uses_consumed_offset)
{
    test_chunk c = make_chunk(50.0, 0, 8);
    BOOST_CHECK_EQUAL(detail::next_timestamp(c, 10e6), 50.0);
    c.offset = 4; // half consumed: next sample is 4 samples later
    BOOST_CHECK_CLOSE(detail::next_timestamp(c, 10e6), 50.0 + 4.0 / 10e6, 1e-9);
}

BOOST_AUTO_TEST_CASE(test_extract_aligned_samples_contiguous_chunks_all_channels)
{
    // Back-to-back chunks are stitched seamlessly, on every channel.
    const double rate = 10e6;
    std::deque<test_chunk> ring;
    ring.push_back(make_chunk(1.0, 0, 6, 2));
    ring.push_back(make_chunk(1.0 + 6.0 / rate, 6, 6, 2));
    ring.push_back(make_chunk(1.0 + 12.0 / rate, 12, 6, 2));

    std::vector<std::vector<std::complex<float>>> out;
    const auto r = detail::extract_aligned_samples(ring, 15, 2, rate, 1.0, out);
    BOOST_CHECK(r.success);
    BOOST_CHECK(!r.zero_filled);
    BOOST_CHECK_EQUAL(r.removed, size_t(15));
    BOOST_CHECK_EQUAL(r.first_timestamp, 1.0);
    BOOST_REQUIRE_EQUAL(out.size(), size_t(2));
    for (size_t i = 0; i < 15; i++) {
        BOOST_CHECK_EQUAL(out[0][i].real(), static_cast<float>(i));
        BOOST_CHECK_EQUAL(out[1][i].real(), static_cast<float>(i + 1000));
    }
    // The last chunk is partly consumed and stays at the front.
    BOOST_REQUIRE_EQUAL(ring.size(), size_t(1));
    BOOST_CHECK_EQUAL(ring.front().offset, size_t(3));
}

BOOST_AUTO_TEST_CASE(test_extract_aligned_samples_drops_stale_chunks)
{
    // Chunks that end before the alignment time are discarded (and counted).
    const double rate = 10e6;
    std::deque<test_chunk> ring;
    ring.push_back(make_chunk(2.0, 0, 10, 1));
    ring.push_back(make_chunk(2.0 + 10.0 / rate, 10, 10, 1));
    ring.push_back(make_chunk(2.0 + 20.0 / rate, 20, 10, 1));

    std::vector<std::vector<std::complex<float>>> out;
    const auto r = detail::extract_aligned_samples(ring, 4, 1, rate, 2.0 + 21.5 / rate, out);
    BOOST_CHECK(r.success);
    // 20 dropped whole + 1 skipped at the head of the third chunk + 4 copied.
    BOOST_CHECK_EQUAL(r.removed, size_t(25));
    BOOST_CHECK_EQUAL(out[0][0].real(), 21.0f);
    BOOST_CHECK_EQUAL(out[0][3].real(), 24.0f);
}

BOOST_AUTO_TEST_CASE(test_extract_aligned_samples_empty_ring)
{
    std::deque<test_chunk> ring;
    std::vector<std::vector<std::complex<float>>> out;
    const auto r = detail::extract_aligned_samples(ring, 8, 2, 10e6, 3.0, out);
    BOOST_CHECK(!r.success);
    BOOST_CHECK(r.zero_filled);
    BOOST_CHECK_EQUAL(r.removed, size_t(0));
    BOOST_REQUIRE_EQUAL(out.size(), size_t(2));
    BOOST_CHECK_EQUAL(out[1].size(), size_t(8));
    BOOST_CHECK_EQUAL(out[1][7], std::complex<float>(0, 0));
}

BOOST_AUTO_TEST_CASE(test_extract_aligned_samples_two_devices_line_up)
{
    // Two devices whose streams start 3 samples apart: extracting both at the
    // later start time yields the same time slots (identical sample index,
    // since make_chunk numbers samples by their time slot here).
    const double rate = 10e6;
    std::deque<test_chunk> a, b;
    a.push_back(make_chunk(5.0, 0, 64, 1));
    b.push_back(make_chunk(5.0 + 3.0 / rate, 3, 64, 1));
    const double align = std::max(detail::next_timestamp(a.front(), rate),
        detail::next_timestamp(b.front(), rate));
    for (int block = 0; block < 4; block++) {
        std::vector<std::vector<std::complex<float>>> oa, ob;
        const double t = align + block * 10.0 / rate;
        const auto ra = detail::extract_aligned_samples(a, 10, 1, rate, t, oa);
        const auto rb = detail::extract_aligned_samples(b, 10, 1, rate, t, ob);
        BOOST_REQUIRE(ra.success && rb.success);
        BOOST_CHECK(oa[0] == ob[0]);
        BOOST_CHECK_EQUAL(oa[0][0].real(), static_cast<float>(3 + 10 * block));
    }
}

BOOST_AUTO_TEST_CASE(test_extract_aligned_samples_removed_accounting)
{
    // Over a long random sequence of chunks with occasional holes, the sum
    // of `removed` must equal what actually left the ring, so the receiver's
    // O(1) "available" counter never drifts.
    const double rate = 1e6;
    std::mt19937 rng(7);
    std::uniform_int_distribution<int> len(1, 300), gap(0, 9);
    std::deque<test_chunk> ring;
    size_t t = 0, total = 0;
    for (int i = 0; i < 2000; i++) {
        if (gap(rng) == 0) {
            t += static_cast<size_t>(len(rng)); // hole (lost samples)
        }
        const size_t n = static_cast<size_t>(len(rng));
        ring.push_back(make_chunk(static_cast<double>(t) / rate, t, n, 1));
        t += n;
        total += n;
    }
    size_t removed = 0, pos = 0, ok = 0, holes = 0;
    while (!ring.empty()) {
        std::vector<std::vector<std::complex<float>>> out;
        const double at = std::max(static_cast<double>(pos) / rate,
            detail::next_timestamp(ring.front(), rate));
        const auto r = detail::extract_aligned_samples(ring, 256, 1, rate, at, out);
        removed += r.removed;
        (r.success ? ok : holes)++;
        if (r.success) {
            // Samples are numbered by their time slot.
            const size_t slot = static_cast<size_t>(at * rate + 0.5);
            BOOST_REQUIRE_EQUAL(out[0][0].real(), static_cast<float>(slot));
            BOOST_REQUIRE_EQUAL(out[0][255].real(), static_cast<float>(slot + 255));
        }
        pos = static_cast<size_t>(at * rate + 0.5) + 256;
    }
    BOOST_TEST_MESSAGE("blocks ok=" << ok << " with holes=" << holes);
    BOOST_CHECK_EQUAL(removed, total);
    BOOST_CHECK_GT(ok, 0u);
    BOOST_CHECK_GT(holes, 0u);
}

BOOST_AUTO_TEST_CASE(test_correct_frequency_is_continuous_across_gaps)
{
    // A radio whose LO landed e Hz high sees a signal at the requested
    // frequency as a tone at -e. Corrected chunk by chunk from each chunk's
    // own timestamp, it must come out as one constant phasor, also across
    // gaps (an overflow) and chunk sizes that are not block multiples.
    const double rate = 50e6, e = -0.834466457; // a real B210 LO error
    const double t_start = 1234.5678;
    const std::complex<double> c0(0.6, -0.3);
    size_t pos = 0;
    double max_err = 0.0;
    for (size_t chunk : {2040u, 2040u, 1500u, 3000u, 2040u}) {
        pos += (chunk == 3000u) ? 777777u : 0u; // a gap before this chunk
        const double t0 = t_start + static_cast<double>(pos) / rate;
        std::vector<std::complex<float>> x(chunk);
        for (size_t i = 0; i < chunk; i++) {
            const double t = t0 + static_cast<double>(i) / rate;
            const double a = -2.0 * 3.141592653589793 * e * t;
            x[i] = static_cast<std::complex<float>>(c0 * std::polar(1.0, a));
        }
        bonded::detail::correct_frequency(x.data(), chunk, e, t0, rate);
        for (const auto& v : x) {
            max_err = std::max(max_err, std::abs(static_cast<std::complex<double>>(v) - c0));
        }
        pos += chunk;
    }
    BOOST_TEST_MESSAGE("max deviation from the constant " << max_err);
    BOOST_CHECK_LT(max_err, 1e-4); // float phasor rounding; a phase break would be ~1
}
