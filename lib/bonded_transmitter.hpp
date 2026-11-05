//
// Copyright 2026 Gabriel Lando
// SPDX-License-Identifier: GPL-3.0-or-later
//

#pragma once

#include "bonded_api.h"
#include <uhd/types/device_addr.hpp>
#include <uhd/types/time_spec.hpp>
#include <uhd/usrp/multi_usrp.hpp>
#include <atomic>
#include <chrono>
#include <cmath>
#include <complex>
#include <condition_variable>
#include <cstdint>
#include <deque>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

namespace bonded {

namespace detail {

//! Rotate samples in place by exp(j(2*pi*freq_hz*(first_index + k)/rate +
//! phase_rad)), k = 0..n-1: the per-device frequency/phase correction, with
//! time referenced to the sample count since the stream started.
inline void apply_phase_correction(std::vector<std::complex<float>>& samp,
    double freq_hz,
    double phase_rad,
    double rate,
    uint64_t first_index)
{
    if (freq_hz == 0.0) {
        if (phase_rad == 0.0) {
            return;
        }
        const std::complex<float> rot(
            static_cast<float>(std::cos(phase_rad)), static_cast<float>(std::sin(phase_rad)));
        for (auto& v : samp) {
            v *= rot;
        }
        return;
    }
    constexpr double kTwoPi = 6.283185307179586476925286766559;
    for (size_t k = 0; k < samp.size(); k++) {
        const double t  = static_cast<double>(first_index + k) / rate;
        const double ph = kTwoPi * freq_hz * t + phase_rad;
        samp[k] *= std::complex<float>(
            static_cast<float>(std::cos(ph)), static_cast<float>(std::sin(ph)));
    }
}

} // namespace detail

/*!
 * Synchronized multi-device transmitter (the engine behind bonded_sink).
 *
 * configure() opens the devices, locks them to the shared 10 MHz reference,
 * tunes them and sets their clocks to a common time on a PPS edge.
 * start_continuous() then streams for the whole session: queue_samples()
 * appends data, and whenever no data is queued the radios get zeros, so they
 * never drop to idle (a B210 needs several hundred microseconds to ramp TX
 * power back up, longer than an 802.11 preamble). A scheduler thread hands the
 * identical block sequence to one sender thread per device, all starting at
 * the same device time; if any device underruns, every device restarts on a
 * new common start time, since an underrun shifts that device's timeline.
 *
 * Per-device sample-delay trim: a positive trim T sends T extra zeros on that
 * device after each (re)start, delaying it by T / rate seconds.
 */
class BONDED_API bonded_transmitter
{
public:
    //! Configuration for the bonded transmitter.
    struct config {
        std::vector<std::string> serials; //!< Ordered list of device serials
        std::string clock_source = "external";
        std::string time_source  = "external";
        double rate = 10e6;              //!< Per-radio TX rate (Hz)
        double freq = 100e6;             //!< Uniform center frequency (Hz)
        std::vector<double> freq_plan;   //!< Optional per-device center freqs (Hz)
        double gain = 0.0;               //!< Uniform normalized TX gain [0,1]
        std::vector<double> gain_plan;   //!< Optional per-device normalized
                                         //!< gains [0,1]; overrides `gain`.
                                         //!< Empty or same length as serials.

        //! Per-device digital frequency (Hz) and phase (rad) correction, see
        //! detail::apply_phase_correction(). Empty or one entry per serial.
        std::vector<double> freq_offset_hz;
        std::vector<double> phase_offset_rad;

        //! Extra UHD device args for every device (more TX frames in flight
        //! absorb host scheduling hiccups at high rates).
        std::string device_args = "num_send_frames=256";
        std::string subdev;              //!< Optional subdev spec
        std::string stream_args;         //!< Optional TX stream args
        bool strict        = false;      //!< Throw on lock failure vs. warn
        double lock_timeout = 5.0;       //!< Seconds to wait for ref_locked

        //! Per-device integer sample-delay trim (samples at `rate`).
        //! Positive value T: prepend T zero-samples to that device's stream.
        //! Must be either empty (no trim) or the same length as serials.
        std::vector<int> sample_delay_trim;
    };

    explicit bonded_transmitter(const config& cfg);
    ~bonded_transmitter();

    // Non-copyable
    bonded_transmitter(const bonded_transmitter&) = delete;
    bonded_transmitter& operator=(const bonded_transmitter&) = delete;

    /*!
     * Full setup: open devices, sync clocks, configure HW, create TX streamers,
     * align time via PPS.
     * \throws uhd::runtime_error on failure (or warns if !strict)
     */
    void configure();

    // --- Continuous mode ---

    /*!
     * Start the scheduler and sender threads; streaming starts 0.5 s later on
     * a common device time. Call queue_samples() to send data, stop() to end.
     */
    void start_continuous();

    /*!
     * Queue n samples per device (ptrs[d] points to device d's samples).
     * Thread-safe.
     */
    void queue_samples(const std::vector<const std::complex<float>*>& ptrs,
        size_t n);

    /*!
     * Stop streaming: whatever is still queued is sent, then the bursts end
     * and the scheduler and sender threads are joined.
     */
    void stop();

    // --- Queries ---
    size_t num_devices() const { return _config.serials.size(); }
    bool is_streaming() const { return _streaming; }
    bool is_configured() const { return _configured; }
    std::vector<size_t> get_underrun_counts() const;

    //! Frames (runs of real data between idle zero-fill) sent so far.
    uint64_t frames_sent() const;

    /*!
     * Change device d's phase rotation while streaming (applied from the next
     * queued chunk on).  The radios' LOs are not phase-coherent, so the
     * relative phase of the two halves is arbitrary per session but stable
     * within it.  Near 180 deg the subcarriers next to the seam see ISI (the
     * two split filters no longer sum to an impulse); near 0 the split is
     * transparent, so a receiver can measure the phase and feed it back here.
     */
    void set_phase_offset(size_t device, double rad);
    double get_phase_offset(size_t device) const;

private:
    // The per-device samples of one queue_samples() call, kept together so
    // the devices always receive the same sequence.
    struct tx_chunk {
        std::vector<std::vector<std::complex<float>>> per_dev; // per_dev[device]
    };

    config _config;
    bool _configured = false;
    bool _streaming  = false;
    size_t _channels_per_device = 1; // TX bonding uses 1 TX channel per device

    std::vector<uhd::usrp::multi_usrp::sptr> _devices;
    std::vector<uhd::tx_streamer::sptr> _streamers;

    std::deque<tx_chunk> _tx_queue;
    std::mutex _tx_mutex;
    std::condition_variable _tx_cv;
    std::thread _send_thread; // scheduler (_send_loop_continuous)

    // Per-device feed: the scheduler pushes the same sequence of blocks to
    // every feed; that device's sender thread drains it with send(). A
    // resync item closes the current burst and opens a new one at `start`;
    // the end item (pushed last, on stop) ends the burst and the thread.
    struct feed_item {
        std::shared_ptr<const std::vector<std::complex<float>>> buf;
        bool resync = false;
        bool end    = false;
        uhd::time_spec_t start;
    };
    struct dev_feed {
        std::mutex mu;
        std::condition_variable cv;
        std::deque<feed_item> q;
        size_t queued = 0; // samples in q
    };
    std::vector<std::unique_ptr<dev_feed>> _feeds;
    std::vector<std::thread> _dev_threads;
    std::mutex _sched_mu;
    std::condition_variable _sched_cv;
    std::atomic<bool> _resync_request{false};
    std::atomic<uint64_t> _resyncs{0};
    std::chrono::steady_clock::time_point _last_resync_req{};
    std::atomic<bool> _stop_flag{false};
    std::atomic<uint64_t> _frames_sent{0};

    // Live copy of config.phase_offset_rad (settable while streaming).
    mutable std::mutex _phase_mutex;
    std::vector<double> _phase_offset;

    // Async-message threads (one per device, from configure() until
    // destruction) count underruns and request the re-alignment.
    std::vector<std::thread> _async_threads;
    std::atomic<bool> _async_stop{false};
    std::vector<std::atomic<size_t>> _underrun_counts;

    // Internal helpers.
    void _open_devices();
    void _apply_sync();
    void _configure_hardware();
    void _create_streamers();
    void _align_time();
    void _send_loop_continuous();
    void _device_send_loop(size_t device_index);
    void _async_msg_loop(size_t device_index);
};

} // namespace bonded
