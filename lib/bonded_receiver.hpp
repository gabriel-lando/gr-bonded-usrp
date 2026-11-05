//
// Copyright 2026 Gabriel Lando
// SPDX-License-Identifier: GPL-3.0-or-later
//

#pragma once

#include "bonded_api.h"
#include <uhd/types/device_addr.hpp>
#include <uhd/usrp/multi_usrp.hpp>
#include <atomic>
#include <chrono>
#include <complex>
#include <deque>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

namespace bonded {

/*!
 * Synchronized multi-device receiver (the engine behind bonded_source).
 *
 * configure() opens the devices, locks them to the shared 10 MHz reference,
 * tunes them and sets their clocks to a common time on a PPS edge.
 * start_continuous() starts one recv thread per device, each filling a
 * timestamped ring; get_aligned_samples() hands out blocks that start at the
 * same device time on every device.
 */
class BONDED_API bonded_receiver
{
public:
    //! Configuration for the bonded receiver
    struct config {
        std::vector<std::string> serials; //!< Ordered list of device serials
        std::string clock_source = "external";
        std::string time_source  = "external";
        double rate = 10e6;              //!< Sample rate (Hz)
        double freq = 100e6;             //!< Uniform center frequency (Hz)
        std::vector<double> freq_plan;   //!< Optional per-device center freqs (Hz)
        double gain = 40.0;              //!< RX gain (dB)
        std::string subdev;              //!< Optional subdev spec (empty = default)
        std::string stream_args;         //!< Optional RX stream args (e.g. num_recv_frames=...,recv_frame_size=...)
        bool strict        = false;      //!< Throw on lock failure vs. warn
        double lock_timeout = 5.0;       //!< Seconds to wait for ref_locked
        double continuous_buffer_seconds = 4.0; //!< Per-device ring retention in continuous mode
    };

    //! One aligned block (see get_aligned_samples()).
    struct rx_result {
        //! data[device_index][channel_index] = sample buffer
        std::vector<std::vector<std::vector<std::complex<float>>>> data;
        std::vector<double> first_timestamps; //!< Per device
        std::vector<bool> success;            //!< Per device
        double inter_device_spread_us = 0.0;  //!< max - min first timestamp (µs)
        bool aligned = false;                 //!< every device returned data
        std::string error_message;            //!< Non-empty on failure
    };

    explicit bonded_receiver(const config& cfg);
    ~bonded_receiver();

    // Non-copyable
    bonded_receiver(const bonded_receiver&) = delete;
    bonded_receiver& operator=(const bonded_receiver&) = delete;

    /*!
     * Full setup: open devices, sync clocks, configure HW, create streamers,
     * align time via PPS.
     * \throws uhd::runtime_error on failure (or warns if !strict)
     */
    void configure();

    // --- Continuous mode ---

    /*!
     * Start background recv threads filling ring buffers.
     * Call get_aligned_samples() to retrieve aligned data.
     */
    void start_continuous();

    /*!
     * Wait (up to timeout_sec) until every device has nsamps samples, then
     * return nsamps per channel starting at a common device time.
     * success[d] is false where device d's block has a hole or ran short; the
     * caller should drop the whole block then, on every device.
     */
    rx_result get_aligned_samples(size_t nsamps, double timeout_sec = 1.0);

    /*!
     * Stop continuous streaming and join background threads.
     */
    void stop();

    // --- Queries ---
    size_t num_devices() const { return _config.serials.size(); }
    size_t channels_per_device() const { return _channels_per_device; }
    size_t num_channels() const { return num_devices() * _channels_per_device; }
    bool is_streaming() const { return _streaming; }
    bool is_configured() const { return _configured; }

private:
    // Internal timestamped chunk for ring buffers
    struct chunk {
        double timestamp;
        size_t num_samples;
        size_t offset = 0;
        std::vector<std::vector<std::complex<float>>> data; // [channel][samples]
    };

    config _config;
    bool _configured = false;
    bool _streaming  = false;
    size_t _channels_per_device = 2; // B210 default

    std::vector<uhd::usrp::multi_usrp::sptr> _devices;
    std::vector<uhd::rx_streamer::sptr> _streamers;

    // Continuous mode state
    std::vector<std::thread> _recv_threads;
    std::vector<std::deque<chunk>> _ring_buffers;
    std::vector<std::mutex> _ring_mutexes;
    std::vector<size_t> _ring_available_samples;
    std::atomic<bool> _stop_flag{false};

    // Opt-in acquisition profiler (BONDED_RX_PROFILE=1): are the radios not
    // delivering, or are we dropping what they deliver?
    bool _rx_profile = false;
    std::vector<std::atomic<size_t>> _prof_recv_samples;    // delivered by recv()
    std::vector<std::atomic<size_t>> _prof_served_samples;  // handed to the flowgraph
    std::vector<std::atomic<size_t>> _prof_dropped_samples; // silent ring trim
    std::vector<std::atomic<size_t>> _prof_ring_depth;      // samples resident
    std::atomic<size_t> _prof_wait_ns{0};                   // consumer poll wait
    std::atomic<size_t> _prof_lock_ns{0};                   // ring mutex hold (extract)
    std::atomic<size_t> _prof_calls{0};
    std::chrono::steady_clock::time_point _prof_last_log{};
    void _maybe_log_rx_profile();

    // Internal helpers
    void _open_devices();
    void _apply_sync();
    void _configure_hardware();
    void _create_streamers();
    void _align_time();
    void _recv_loop(size_t device_index);
};

} // namespace bonded
