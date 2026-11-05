//
// Copyright 2026 Gabriel Lando
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Filter design, dot-product kernels and the streaming polyphase resamplers
// shared by the reconstructor and the band splitter.
//

#include "bonded_resample.hpp"
#include "test_helpers.hpp"
#include <boost/test/unit_test.hpp>
#include <cmath>
#include <complex>
#include <numeric>
#include <random>
#include <vector>

using namespace bonded::detail;
using bonded_test::tone;
using bonded_test::tone_amplitude;

namespace {

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

std::vector<cf> feed_chunked(poly_resampler& rs, const std::vector<cf>& x,
    const std::vector<size_t>& chunks)
{
    std::vector<cf> out;
    for (size_t off = 0, c = 0; off < x.size(); c++) {
        const size_t n = std::min(chunks[c % chunks.size()], x.size() - off);
        rs.feed(x.data() + off, n, out);
        off += n;
    }
    return out;
}

} // namespace

BOOST_AUTO_TEST_CASE(test_design_lowpass_real_dc_gain_and_symmetry)
{
    for (double scale : {1.0, 2.0, 5.0}) {
        const auto h = design_lowpass_real(101, 0.1, scale);
        BOOST_CHECK_CLOSE(std::accumulate(h.begin(), h.end(), 0.0), scale, 1e-4);
        for (size_t i = 0; i < h.size() / 2; i++) {
            BOOST_REQUIRE_CLOSE(h[i], h[h.size() - 1 - i], 1e-3);
        }
    }
}

BOOST_AUTO_TEST_CASE(test_dot_kernels_agree)
{
    const size_t K = 37;
    const auto a = random_samples(K, 3);
    const auto t = random_samples(K, 4);
    std::vector<float> tr(K), ti(K), ar(K), ai(K), treal(K);
    cd ref_cc(0, 0), ref_cr(0, 0);
    for (size_t k = 0; k < K; k++) {
        ar[k] = a[k].real();
        ai[k] = a[k].imag();
        tr[k] = t[k].real();
        ti[k] = t[k].imag();
        treal[k] = t[k].real();
        ref_cc += static_cast<cd>(a[k]) * static_cast<cd>(t[k]);
        ref_cr += static_cast<cd>(a[k]) * static_cast<double>(treal[k]);
    }
    BOOST_CHECK_SMALL(std::abs(static_cast<cd>(dot_cc_hot(a.data(), t.data(), K)) - ref_cc), 1e-4);
    BOOST_CHECK_SMALL(
        std::abs(static_cast<cd>(dot_cc_split(ar.data(), ai.data(), tr.data(), ti.data(), K)) - ref_cc),
        1e-4);
    BOOST_CHECK_SMALL(std::abs(static_cast<cd>(dot_cr_hot(a.data(), treal.data(), K)) - ref_cr), 1e-4);
    BOOST_CHECK_SMALL(std::abs(static_cast<cd>(dot_cr(a.data(), treal.data(), K)) - ref_cr), 1e-4);
}

BOOST_AUTO_TEST_CASE(test_poly_resampler_output_count_and_chunking)
{
    // Output count follows the ratio, and chunking does not change the output.
    struct ratio
    {
        size_t L, M;
    };
    for (const ratio r : {ratio{5, 8}, ratio{8, 5}, ratio{1, 2}, ratio{3, 4}, ratio{16, 9}}) {
        const auto x = random_samples(20000, 5);
        poly_resampler whole, chunked;
        whole.init(r.L, r.M);
        chunked.init(r.L, r.M);
        std::vector<cf> a;
        whole.feed(x.data(), x.size(), a);
        const std::vector<cf> b = feed_chunked(chunked, x, {1, 2, 17, 1000, 333});
        const double expect = static_cast<double>(x.size()) * r.L / r.M;
        BOOST_CHECK_LE(std::abs(static_cast<double>(a.size()) - expect), 1.0);
        BOOST_REQUIRE_EQUAL(a.size(), b.size());
        for (size_t i = 0; i < a.size(); i++) {
            BOOST_REQUIRE_SMALL(std::abs(a[i] - b[i]), 1e-5f);
        }
    }
}

BOOST_AUTO_TEST_CASE(test_poly_resampler_passes_and_rejects)
{
    // 20 -> 12.5 MS/s: an in-band tone keeps its amplitude, a tone above the
    // new Nyquist is suppressed.
    const double fin = 20e6, fout = 12.5e6;
    poly_resampler rs;
    rs.init(5, 8);
    std::vector<cf> y;
    const auto x = tone(3e6, fin, 40000);
    rs.feed(x.data(), x.size(), y);
    BOOST_CHECK_CLOSE(std::abs(tone_amplitude(y, 3e6, fout, 1000, y.size() - 1000)), 1.0, 2.0);

    poly_resampler rs2;
    rs2.init(5, 8);
    std::vector<cf> z;
    const auto x2 = tone(8.5e6, fin, 40000); // aliases to -4 MHz if not filtered
    rs2.feed(x2.data(), x2.size(), z);
    BOOST_CHECK_LT(bonded_test::mean_power(z, 1000, z.size() - 1000), 1e-3);
}

BOOST_AUTO_TEST_CASE(test_poly_resampler_cutoff_and_length)
{
    // A custom cutoff moves the -6 dB point; more taps per phase make the
    // transition sharper (the band splitter relies on both).
    const double fin = 20e6;
    auto gain_at = [&](double f, size_t taps, double cutoff) {
        poly_resampler rs;
        rs.init(1, 2, taps, cutoff);
        std::vector<cf> y;
        const auto x = tone(f, fin, 60000);
        rs.feed(x.data(), x.size(), y);
        return std::abs(tone_amplitude(y, f, fin / 2, 2000, y.size() - 2000));
    };
    const double cutoff = 4.375e6 / fin; // half-width 4.375 MHz
    BOOST_CHECK_CLOSE(gain_at(4.375e6, 64, cutoff), 0.5, 5.0);
    BOOST_CHECK_GT(gain_at(3.75e6, 64, cutoff), 0.95);
    BOOST_CHECK_LT(gain_at(3.75e6, 6, cutoff), 0.9); // short filter: wide transition
}

BOOST_AUTO_TEST_CASE(test_poly_resampler_c_slices_match_one_pass)
{
    // compute() over disjoint slices (what the reconstructor's threads do)
    // must give exactly the single-range result, over many batches.
    const size_t L = 8, M = 5;
    const auto proto = random_samples(13 * L + 1, 6);
    poly_resampler_c one, sliced;
    one.init(L, M, proto);
    sliced.init(L, M, proto);
    const auto x = random_samples(30000, 7);
    const std::vector<size_t> chunks = {4096, 100, 1, 9000, 333};
    for (size_t off = 0, c = 0; off < x.size(); c++) {
        const size_t n = std::min(chunks[c % chunks.size()], x.size() - off);
        const size_t p1 = one.prepare(x.data() + off, n);
        const size_t p2 = sliced.prepare(x.data() + off, n);
        BOOST_REQUIRE_EQUAL(p1, p2);
        std::vector<cf> a(p1), b(p2);
        one.compute(a.data(), 0, p1);
        const size_t k = 5;
        for (size_t s = 0; s < k; s++) {
            sliced.compute(b.data(), p2 * s / k, p2 * (s + 1) / k);
        }
        BOOST_REQUIRE(a == b);
        one.commit();
        sliced.commit();
        off += n;
    }
}
