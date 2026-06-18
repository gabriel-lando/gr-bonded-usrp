//
// Copyright 2026 Gabriel Lando
// SPDX-License-Identifier: GPL-3.0-or-later
//

#pragma once

#include <algorithm>
#include <complex>
#include <cstddef>
#include <deque>
#include <string>
#include <vector>

namespace bonded { namespace detail {

struct burst_alignment_result
{
    double inter_device_spread_us = 0.0;
    bool aligned                  = false;
    std::string error_message;
};

inline burst_alignment_result evaluate_burst_alignment(const std::vector<double>& first_timestamps,
    const std::vector<bool>& success,
    const double align_threshold_secs)
{
    burst_alignment_result result;
    const size_t num_dev = success.size();
    if (first_timestamps.size() != num_dev) {
        result.error_message = "first_timestamps and success size mismatch";
        return result;
    }

    std::vector<double> valid_ts;
    for (size_t d = 0; d < num_dev; d++) {
        if (success[d] && first_timestamps[d] >= 0.0) {
            valid_ts.push_back(first_timestamps[d]);
        }
    }

    if (valid_ts.size() == num_dev && num_dev > 1) {
        const double ts_min         = *std::min_element(valid_ts.begin(), valid_ts.end());
        const double ts_max         = *std::max_element(valid_ts.begin(), valid_ts.end());
        result.inter_device_spread_us = (ts_max - ts_min) * 1e6;
        result.aligned              = (ts_max - ts_min) <= align_threshold_secs;
        return result;
    }
    if (valid_ts.size() == num_dev && num_dev == 1) {
        result.aligned = true;
        return result;
    }

    result.error_message = "One or more devices did not receive data";
    return result;
}

template <typename ChunkT>
inline size_t samples_available(const std::deque<ChunkT>& ring)
{
    size_t available = 0;
    for (const auto& c : ring) {
        available += c.num_samples - c.offset;
    }
    return available;
}

template <typename ChunkT>
inline bool find_alignment_timestamp(const std::vector<std::deque<ChunkT>>& ring_buffers,
    double& align_ts,
    std::string& error_message)
{
    align_ts = -1.0;
    for (size_t d = 0; d < ring_buffers.size(); d++) {
        if (ring_buffers[d].empty()) {
            error_message = "Device " + std::to_string(d) + " has no data";
            return false;
        }
        const double first_ts = ring_buffers[d].front().timestamp;
        if (first_ts > align_ts) {
            align_ts = first_ts;
        }
    }
    return true;
}

struct extract_result
{
    bool success          = false;
    bool zero_filled      = false;
    double first_timestamp = -1.0;
};

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
    bool first       = true;
    while (collected < nsamps && !ring.empty()) {
        auto& front = ring.front();
        size_t skip = 0;
        const double front_ts =
            front.timestamp + static_cast<double>(front.offset) / rate;

        if (first) {
            const double offset_secs = align_ts - front_ts;
            if (offset_secs > 0.0) {
                skip = static_cast<size_t>(offset_secs * rate);
                if (front.offset + skip >= front.num_samples) {
                    ring.pop_front();
                    continue;
                }
            }
            result.first_timestamp = align_ts;
            first                  = false;
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

        front.offset = start + to_copy;
        if (front.offset >= front.num_samples) {
            ring.pop_front();
        }
    }

    result.success = (collected >= nsamps);
    if (!result.success) {
        if (result.first_timestamp < 0.0) {
            result.first_timestamp = align_ts;
        }
        result.zero_filled = true;
    }

    return result;
}

}} // namespace bonded::detail
