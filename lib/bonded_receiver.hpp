//
// Copyright 2026 Gabriel Lando
// SPDX-License-Identifier: GPL-3.0-or-later
//

#pragma once

#include <atomic>
#include "bonded_api.h"
#include <uhd/types/device_addr.hpp>
#include <uhd/usrp/multi_usrp.hpp>
#include <complex>
#include <deque>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

namespace bonded {

/*!
 * Reusable multi-device synchronized streaming engine.
 *
 * Encapsulates the full lifecycle of bonded USB USRP (B210) reception:
 * device open → sync → configure → stream → align → deliver.
 *
 * Supports two modes:
 *   - Burst: one-shot timed capture via capture_burst()
 *   - Continuous: background recv threads with get_aligned_samples()
 *
 * Designed for direct use from standalone tools, Python (via pybind11),
 * or as the engine inside a future GNURadio OOT block.
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

    //! Result of a receive operation (burst or continuous chunk)
    struct rx_result {
        //! data[device_index][channel_index] = sample buffer
        std::vector<std::vector<std::vector<std::complex<float>>>> data;
        std::vector<double> first_timestamps; //!< Per device
        std::vector<bool> success;            //!< Per device
        double inter_device_spread_us = 0.0;  //!< max - min first timestamp (µs)
        bool aligned = false;                 //!< spread < threshold
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

    // --- Burst mode ---

    /*!
     * Capture exactly nsamps from all devices with a timed stream command.
     * \param nsamps   Number of samples per channel to capture
     * \param delay_sec Seconds in the future to schedule the burst
     * \return rx_result with captured data and alignment info
     */
    rx_result capture_burst(size_t nsamps, double delay_sec = 1.5);

    // --- Continuous mode ---

    /*!
     * Start background recv threads filling ring buffers.
     * Call get_aligned_samples() to retrieve aligned data.
     */
    void start_continuous();

    /*!
     * Block until all devices have aligned data, return nsamps per channel.
     * \param nsamps      Number of samples to retrieve per channel
     * \param timeout_sec Maximum wait time for alignment
     * \return rx_result with aligned data or error
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
    std::vector<size_t> get_overflow_counts() const;
    std::vector<size_t> get_zero_fill_counts() const;

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
    std::vector<std::atomic<size_t>> _overflow_counts;
    std::vector<std::atomic<size_t>> _zero_fill_counts;
    bool _stop_flag = false;

    // Internal helpers
    void _open_devices();
    void _apply_sync();
    void _configure_hardware();
    void _create_streamers();
    void _align_time();
    void _recv_loop(size_t device_index);
};

} // namespace bonded
