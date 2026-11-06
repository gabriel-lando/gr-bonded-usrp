//
// Copyright 2026 Gabriel Lando
//

// Timestamp alignment of the per-device receive rings (header-only, so the
// tests can drive it with plain chunk structs).
//

#pragma once

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <deque>
#include <vector>

namespace bonded { namespace detail {

//! Time of a chunk's next unconsumed sample (its start plus what was taken).
//! ChunkT needs: double timestamp; size_t num_samples, offset;
//! std::vector<std::vector<std::complex<float>>> data (per channel).
template <typename ChunkT>
inline double next_timestamp(const ChunkT& c, const double rate)
{
    return c.timestamp + static_cast<double>(c.offset) / rate;
}

struct extract_result
{
    bool success           = false; //!< nsamps contiguous samples were copied
    bool zero_filled       = false; //!< output has zeros (hole or ran short)
    double first_timestamp = -1.0;  //!< time of out[*][0]
    //! Samples this call took out of the ring (consumed + discarded), so the
    //! caller can keep its "available" counter in O(1): rescanning a ring of
    //! tens of thousands of chunks under its mutex is too slow at high rates.
    size_t removed = 0;
};

/*!
 * Copy nsamps samples per channel starting at time align_ts out of one
 * device's ring. Chunks wholly before align_ts are dropped. Every chunk is
 * placed on the time grid by its own timestamp: a hole (lost samples after an
 * overflow) stays a hole of zeros and makes the result unsuccessful, instead of
 * pulling later samples early, which would leave this device ahead of the
 * others for the rest of the session.
 */

template <typename ChunkT>
inline extract_result extract_aligned_samples(std::deque<ChunkT>& ring,
    const size_t nsamps,
    const size_t nch,
    const double rate,
    const double align_ts,
    std::vector<std::vector<std::complex<float>>>& out)
{
    extract_result result;

    out.clear();
    out.resize(nch);
    for (size_t c = 0; c < nch; c++) {
        out[c].assign(nsamps, std::complex<float>(0.0f, 0.0f));
    }

    while (!ring.empty()) {
        const auto& front = ring.front();
        const double chunk_end_ts =
            front.timestamp + static_cast<double>(front.num_samples) / rate;
        if (chunk_end_ts <= align_ts) {
            result.removed += front.num_samples - front.offset;
            ring.pop_front();
        } else {
            break;
        }
    }

    if (ring.empty()) {
        result.first_timestamp = align_ts;
        result.zero_filled     = true;
        return result;
    }

    size_t collected = 0;
    bool gapped      = false;
    result.first_timestamp = align_ts;
    while (collected < nsamps && !ring.empty()) {
        auto& front = ring.front();
        const double front_ts = next_timestamp(front, rate);
        const double expected_ts = align_ts + static_cast<double>(collected) / rate;
        const double ahead = (front_ts - expected_ts) * rate; // samples

        if (ahead > 0.5) {
            // Hole before this chunk: skip over it in the output.
            const size_t hole = std::min(
                static_cast<size_t>(ahead + 0.5), nsamps - collected);
            collected += hole;
            gapped = true;
            continue;
        }
        size_t skip = 0;
        if (ahead < 0.0) {
            // Chunk starts before where we are: drop its stale head
            // (truncating; tiny epsilon so 2.9999999 counts as 3).
            skip = static_cast<size_t>(-ahead + 1e-6);
            if (front.offset + skip >= front.num_samples) {
                result.removed += front.num_samples - front.offset;
                ring.pop_front();
                continue;
            }
        }

        const size_t start  = front.offset + skip;
        const size_t avail  = front.num_samples - start;
        const size_t to_copy = std::min(avail, nsamps - collected);

        for (size_t c = 0; c < nch; c++) {
            std::copy(front.data[c].begin() + start,
                front.data[c].begin() + start + to_copy,
                out[c].begin() + collected);
        }
        collected += to_copy;

        result.removed += skip + to_copy;
        front.offset = start + to_copy;
        if (front.offset >= front.num_samples) {
            ring.pop_front();
        }
    }

    result.success = (collected >= nsamps) && !gapped;
    if (!result.success) {
        result.zero_filled = true;
    }

    return result;
}

/*!
 * Shift n samples, the first taken at device time t0, up by freq_hz, with the
 * phase referenced to device time 0. A radio whose LO landed freq_hz above
 * the requested frequency then delivers what an exactly tuned one would, and
 * the correction stays continuous across any gap, since every chunk is
 * anchored by its own timestamp.
 */
inline void correct_frequency(std::complex<float>* x,
    const size_t n,
    const double freq_hz,
    const double t0,
    const double rate)
{
    if (freq_hz == 0.0) {
        return;
    }
    constexpr double kTwoPi = 6.283185307179586476925;
    const double c0   = freq_hz * t0;
    const double a0   = kTwoPi * (c0 - std::floor(c0));
    const double step = kTwoPi * freq_hz / rate;
    const std::complex<float> inc(
        static_cast<float>(std::cos(step)), static_cast<float>(std::sin(step)));
    // A float phasor drifts slowly; re-anchor it every block.
    constexpr size_t kBlock = 1024;
    for (size_t i = 0; i < n; i += kBlock) {
        const double a = a0 + step * static_cast<double>(i);
        std::complex<float> p(static_cast<float>(std::cos(a)), static_cast<float>(std::sin(a)));
        const size_t end = std::min(n, i + kBlock);
        for (size_t k = i; k < end; k++) {
            x[k] *= p;
            p *= inc;
        }
    }
}

}} // namespace bonded::detail
