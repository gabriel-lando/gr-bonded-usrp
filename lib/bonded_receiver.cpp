//
// Copyright 2026 Gabriel Lando
//

#include "bonded_receiver.hpp"

#include "bonded_alignment.hpp"
#include <uhd/exception.hpp>
#include <uhd/types/tune_request.hpp>
#include <uhd/utils/thread.hpp>
#include <boost/format.hpp>
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <sstream>

namespace bonded {

bonded_receiver::bonded_receiver(const config& cfg) : _config(cfg) {}

bonded_receiver::~bonded_receiver()
{
    if (_streaming) {
        stop();
    }
}

// =========================================================================
// configure() — full setup pipeline
// =========================================================================
void bonded_receiver::configure()
{
    if (_configured) {
        throw uhd::runtime_error("bonded_receiver::configure() called twice");
    }
    _open_devices();
    _apply_sync();
    _configure_hardware();
    _create_streamers();
    _align_time();
    _configured = true;
}

// =========================================================================
// Continuous mode
// =========================================================================
void bonded_receiver::start_continuous()
{
    if (!_configured) {
        throw uhd::runtime_error(
            "bonded_receiver::start_continuous() called before configure()");
    }
    if (_streaming) {
        throw uhd::runtime_error("bonded_receiver: already streaming");
    }

    const size_t num_dev = _devices.size();
    _stop_flag = false;
    _ring_buffers.resize(num_dev);
    _ring_mutexes = std::vector<std::mutex>(num_dev);
    if (const char* e = std::getenv("BONDED_RX_PROFILE")) {
        _rx_profile = (std::string(e) != "0");
    }
    _prof_recv_samples    = std::vector<std::atomic<size_t>>(num_dev);
    _prof_served_samples  = std::vector<std::atomic<size_t>>(num_dev);
    _prof_dropped_samples = std::vector<std::atomic<size_t>>(num_dev);
    _prof_ring_depth      = std::vector<std::atomic<size_t>>(num_dev);
    _prof_last_log        = std::chrono::steady_clock::now();
    _ring_available_samples.assign(num_dev, 0);

    // Issue a common timed continuous stream command to all devices.
    // Using stream_now=true may fail for multi-channel streamers on some
    // devices (e.g., B210), and a timed start keeps devices aligned.
    const uhd::time_spec_t start_time =
        _devices[0]->get_time_now() + uhd::time_spec_t(0.1);
    uhd::stream_cmd_t stream_cmd(uhd::stream_cmd_t::STREAM_MODE_START_CONTINUOUS);
    stream_cmd.stream_now = false;
    stream_cmd.time_spec  = start_time;
    for (auto& s : _streamers) {
        s->issue_stream_cmd(stream_cmd);
    }

    // Launch recv threads
    _recv_threads.reserve(num_dev);
    for (size_t d = 0; d < num_dev; d++) {
        _recv_threads.emplace_back(&bonded_receiver::_recv_loop, this, d);
    }

    // Prime ring buffers so the first get_aligned_samples() calls don't hit
    // startup transients while threads are still filling initial chunks.
    const auto prime_deadline =
        std::chrono::steady_clock::now() + std::chrono::milliseconds(500);
    while (std::chrono::steady_clock::now() < prime_deadline) {
        bool all_have_data = true;
        for (size_t d = 0; d < num_dev; d++) {
            std::lock_guard<std::mutex> lock(_ring_mutexes[d]);
            if (_ring_buffers[d].empty()) {
                all_have_data = false;
                break;
            }
        }
        if (all_have_data) {
            break;
        }
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }

    _streaming = true;
}

bonded_receiver::rx_result bonded_receiver::get_aligned_samples(
    size_t nsamps, double timeout_sec)
{
    if (!_streaming) {
        throw uhd::runtime_error(
            "bonded_receiver::get_aligned_samples() called while not streaming");
    }

    const size_t num_dev = _devices.size();
    rx_result result;
    result.data.resize(num_dev);
    result.first_timestamps.resize(num_dev, -1.0);
    result.success.resize(num_dev, false);
    bool degraded = false;

    const auto deadline = std::chrono::steady_clock::now()
                          + std::chrono::duration<double>(timeout_sec);
    const auto wait_start = std::chrono::steady_clock::now();

    // Wait until all ring buffers have enough data
    while (std::chrono::steady_clock::now() < deadline) {
        bool all_ready = true;
        for (size_t d = 0; d < num_dev; d++) {
            std::lock_guard<std::mutex> lock(_ring_mutexes[d]);
            if (_ring_available_samples[d] < nsamps) {
                all_ready = false;
                break;
            }
        }
        if (all_ready) break;
        std::this_thread::sleep_for(std::chrono::microseconds(100));
    }
    if (_rx_profile) {
        _prof_wait_ns.fetch_add(
            static_cast<size_t>(std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::steady_clock::now() - wait_start).count()),
            std::memory_order_relaxed);
        _prof_calls.fetch_add(1, std::memory_order_relaxed);
    }

    // Hold every ring while choosing the start time and extracting: a full
    // ring loses its oldest chunk on every push, and if a front moved between
    // the two steps, every fetch would start on a hole and fail for as long
    // as the rings stay full.
    const auto lock_start =
        _rx_profile ? std::chrono::steady_clock::now() : std::chrono::steady_clock::time_point{};
    std::vector<std::unique_lock<std::mutex>> locks;
    locks.reserve(num_dev);
    for (size_t d = 0; d < num_dev; d++) {
        locks.emplace_back(_ring_mutexes[d]);
    }

    double align_ts = -1.0;
    for (size_t d = 0; d < num_dev; d++) {
        if (_ring_buffers[d].empty()) {
            result.error_message = "Device " + std::to_string(d) + " has no data";
            return result;
        }
        // Next unconsumed sample, not the chunk start: a partly consumed
        // front chunk would otherwise give a stale (too early) alignment.
        const double first_ts = detail::next_timestamp(_ring_buffers[d].front(), _config.rate);
        if (first_ts > align_ts) {
            align_ts = first_ts;
        }
    }

    for (size_t d = 0; d < num_dev; d++) {
        auto& ring = _ring_buffers[d];

        const auto extract = detail::extract_aligned_samples(
            ring, nsamps, _channels_per_device, _config.rate, align_ts, result.data[d]);

        _ring_available_samples[d] -=
            std::min(_ring_available_samples[d], extract.removed);

        result.first_timestamps[d] = extract.first_timestamp;
        result.success[d]          = extract.success;
        if (_rx_profile && extract.success) {
            _prof_served_samples[d].fetch_add(nsamps, std::memory_order_relaxed);
        }
        if (extract.zero_filled) {
            degraded = true;
        }
    }
    locks.clear();
    if (_rx_profile) {
        _prof_lock_ns.fetch_add(
            static_cast<size_t>(std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::steady_clock::now() - lock_start).count()),
            std::memory_order_relaxed);
        _maybe_log_rx_profile();
    }

    std::vector<double> valid_ts;
    for (size_t d = 0; d < num_dev; d++) {
        if (result.first_timestamps[d] >= 0.0) {
            valid_ts.push_back(result.first_timestamps[d]);
        }
    }
    if (valid_ts.size() == num_dev) {
        const double ts_min = *std::min_element(valid_ts.begin(), valid_ts.end());
        const double ts_max = *std::max_element(valid_ts.begin(), valid_ts.end());
        result.inter_device_spread_us = (ts_max - ts_min) * 1e6;
        result.aligned = true; // aligned by construction (same align_ts)
    }

    if (degraded) {
        result.error_message =
            "One or more devices overflowed or ran short; missing samples were zero-filled";
    }

    return result;
}

void bonded_receiver::stop()
{
    if (!_streaming) return;

    _stop_flag = true;

    // Stop streaming on all devices
    uhd::stream_cmd_t stop_cmd(uhd::stream_cmd_t::STREAM_MODE_STOP_CONTINUOUS);
    for (auto& s : _streamers) {
        s->issue_stream_cmd(stop_cmd);
    }

    // Join all threads
    for (auto& t : _recv_threads) {
        if (t.joinable()) {
            t.join();
        }
    }
    _recv_threads.clear();
    _ring_buffers.clear();
    _streaming = false;
}

// =========================================================================
// Private helpers
// =========================================================================

void bonded_receiver::_open_devices()
{
    std::cout << boost::format(
                     "[bonded_receiver] Opening %zu device(s)...\n")
                     % _config.serials.size();

    for (size_t i = 0; i < _config.serials.size(); i++) {
        uhd::device_addr_t dev_args;
        dev_args["serial"] = _config.serials[i];
        std::cout << boost::format(
                         "[bonded_receiver] Opening serial=%s ...\n")
                         % _config.serials[i];
        _devices.push_back(uhd::usrp::multi_usrp::make(dev_args));
    }
}

void bonded_receiver::_apply_sync()
{
    std::cout << "[bonded_receiver] Applying sync config to all devices...\n";
    const double lock_timeout = _config.lock_timeout;
    const bool strict         = _config.strict;

    for (size_t d = 0; d < _devices.size(); d++) {
        auto& usrp = _devices[d];
        const size_t n = usrp->get_num_mboards();
        for (size_t m = 0; m < n; m++) {
            usrp->set_clock_source(_config.clock_source, m);
            usrp->set_time_source(_config.time_source, m);
        }

        // Wait for reference lock on external/gpsdo
        if (_config.clock_source == "external" || _config.clock_source == "gpsdo") {
            const auto deadline = std::chrono::steady_clock::now()
                                  + std::chrono::duration<double>(lock_timeout);
            bool locked = false;
            while (std::chrono::steady_clock::now() < deadline) {
                bool all_locked = true;
                for (size_t m = 0; m < n; m++) {
                    try {
                        if (!usrp->get_mboard_sensor("ref_locked", m).to_bool()) {
                            all_locked = false;
                            break;
                        }
                    } catch (const uhd::lookup_error&) {
                        // sensor absent — assume locked
                    }
                }
                if (all_locked) {
                    locked = true;
                    break;
                }
                std::this_thread::sleep_for(std::chrono::milliseconds(100));
            }
            if (!locked) {
                const std::string msg =
                    "Device " + std::to_string(d) + " (serial="
                    + _config.serials[d] + "): ref not locked within "
                    + std::to_string(static_cast<int>(lock_timeout)) + " s.";
                if (strict) {
                    throw uhd::runtime_error(msg);
                }
                std::cerr << "[bonded_receiver] WARNING: " << msg << "\n";
            } else {
                std::cout << boost::format(
                                 "[bonded_receiver] Device %u: reference locked.\n")
                                 % d;
            }
        }
    }
}

void bonded_receiver::_configure_hardware()
{
    const bool use_freq_plan = !_config.freq_plan.empty();
    if (use_freq_plan && _config.freq_plan.size() != _devices.size()) {
        throw uhd::runtime_error(
            "bonded_receiver: freq_plan size must match number of devices");
    }

    for (size_t d = 0; d < _devices.size(); d++) {
        auto& usrp = _devices[d];
        if (!_config.subdev.empty()) {
            usrp->set_rx_subdev_spec(_config.subdev);
        }
        usrp->set_rx_rate(_config.rate);
        const double dev_freq = use_freq_plan ? _config.freq_plan[d] : _config.freq;
        const size_t nch = usrp->get_rx_num_channels();
        for (size_t ch = 0; ch < nch; ch++) {
            usrp->set_rx_freq(uhd::tune_request_t(dev_freq), ch);
            usrp->set_rx_gain(_config.gain, ch);
        }
        std::cout << boost::format(
                         "[bonded_receiver] Device %u: %.3f Msps, %.3f MHz, %u ch%s\n")
                         % d % (usrp->get_rx_rate() / 1e6)
                         % (usrp->get_rx_freq(0) / 1e6) % nch
                         % (use_freq_plan ? " (freq_plan)" : "");

        if (d == 0) {
            _channels_per_device = nch;
        }
    }
}

void bonded_receiver::_create_streamers()
{
    for (size_t d = 0; d < _devices.size(); d++) {
        const size_t nch = _channels_per_device;
        std::vector<size_t> ch_list(nch);
        for (size_t c = 0; c < nch; c++) {
            ch_list[c] = c;
        }
        uhd::stream_args_t sa("fc32");
        if (!_config.stream_args.empty()) {
            sa.args = uhd::device_addr_t(_config.stream_args);
        }
        sa.channels  = ch_list;
        _streamers.push_back(_devices[d]->get_rx_stream(sa));
    }
}

void bonded_receiver::_align_time()
{
    if (_config.time_source == "external" || _config.time_source == "gpsdo") {
        // B210 quirk: MCR changes reset time_source to "none".
        // Re-apply after all MCR changes from set_rx_rate/get_rx_stream.
        for (auto& dev : _devices) {
            dev->set_time_source(_config.time_source);
        }

        // Wait 500ms for FPGA PPS detection circuit to stabilise
        std::cout << "[bonded_receiver] Waiting for FPGA to stabilise after MCR change...\n";
        std::this_thread::sleep_for(std::chrono::milliseconds(500));

        // Detect PPS edge on device 0
        std::cout << "[bonded_receiver] Waiting for PPS edge on device 0...\n";
        const auto pps_deadline =
            std::chrono::steady_clock::now() + std::chrono::milliseconds(2500);
        uhd::time_spec_t last_pps = _devices[0]->get_time_last_pps();
        while (_devices[0]->get_time_last_pps() == last_pps) {
            if (std::chrono::steady_clock::now() > pps_deadline) {
                throw uhd::runtime_error(
                    "No PPS detected on device 0 within 2.5 s.\n"
                    "Check PPS cable connection to all devices.");
            }
            std::this_thread::sleep_for(std::chrono::milliseconds(1));
        }

        // PPS just fired — arm ALL devices before the next edge (~999 ms away)
        std::cout << "[bonded_receiver] PPS detected — arming all devices...\n";
        for (auto& dev : _devices) {
            dev->set_time_next_pps(uhd::time_spec_t(0.0));
        }

        // Wait for the next PPS edge to latch time 0 on every device
        std::this_thread::sleep_for(std::chrono::seconds(1));
        std::cout << "[bonded_receiver] PPS alignment done on all devices.\n";
    } else {
        // Internal time: set_time_now(0) on all devices simultaneously
        std::cout << "[bonded_receiver] Internal time: set_time_now(0) on all devices.\n";
        for (auto& dev : _devices) {
            dev->set_time_now(
                uhd::time_spec_t(0.0), uhd::usrp::multi_usrp::ALL_MBOARDS);
        }
    }
}

void bonded_receiver::_recv_loop(size_t d)
{
    uhd::set_thread_priority_safe();

    auto& streamer      = _streamers[d];
    const size_t nch    = _channels_per_device;
    const size_t buf_sz = streamer->get_max_num_samps();

    std::vector<std::vector<std::complex<float>>> bufs(
        nch, std::vector<std::complex<float>>(buf_sz));
    std::vector<std::complex<float>*> ptrs(nch);
    for (size_t c = 0; c < nch; c++) {
        ptrs[c] = bufs[c].data();
    }

    uhd::rx_metadata_t md;
    size_t pending_overflow_logs = 0;
    auto last_overflow_log = std::chrono::steady_clock::now();
    while (!_stop_flag) {
        const size_t n = streamer->recv(ptrs, buf_sz, md, 0.1);

        // An overflow is reported on a recv() that returns no samples, so
        // check the metadata before skipping empty returns.
        if (md.error_code == uhd::rx_metadata_t::ERROR_CODE_OVERFLOW) {
            pending_overflow_logs++;
            const auto now = std::chrono::steady_clock::now();
            if (now - last_overflow_log >= std::chrono::seconds(1)) {
                std::cerr << boost::format(
                                 "[bonded_receiver] Device %u: %u overflow(s) in the last second\n")
                                 % d % pending_overflow_logs;
                pending_overflow_logs = 0;
                last_overflow_log = now;
            }
        } else if (md.error_code != uhd::rx_metadata_t::ERROR_CODE_NONE
                   && md.error_code != uhd::rx_metadata_t::ERROR_CODE_TIMEOUT) {
            std::cerr << boost::format(
                             "[bonded_receiver] Device %u: error %s\n")
                             % d % md.strerror();
            continue;
        }
        if (n == 0) continue;
        if (_rx_profile) {
            _prof_recv_samples[d].fetch_add(n, std::memory_order_relaxed);
        }

        // Build chunk and push to ring buffer
        chunk c;
        c.timestamp   = md.time_spec.get_real_secs();
        c.num_samples = n;
        c.data.resize(nch);
        for (size_t ch = 0; ch < nch; ch++) {
            c.data[ch].assign(bufs[ch].begin(), bufs[ch].begin() + n);
        }

        {
            std::lock_guard<std::mutex> lock(_ring_mutexes[d]);
            _ring_buffers[d].push_back(std::move(c));
            _ring_available_samples[d] += n;
            // Keep at most continuous_buffer_seconds; the oldest data goes.
            const size_t max_chunks =
                static_cast<size_t>(_config.continuous_buffer_seconds * _config.rate / buf_sz)
                + 1;
            while (_ring_buffers[d].size() > max_chunks) {
                const size_t lost = _ring_buffers[d].front().num_samples
                                    - _ring_buffers[d].front().offset;
                _ring_available_samples[d] -= lost;
                if (_rx_profile) {
                    _prof_dropped_samples[d].fetch_add(lost, std::memory_order_relaxed);
                }
                _ring_buffers[d].pop_front();
            }
            if (_rx_profile) {
                _prof_ring_depth[d].store(_ring_available_samples[d],
                    std::memory_order_relaxed);
            }
        }
    }

    if (pending_overflow_logs != 0) {
        std::cerr << boost::format(
                         "[bonded_receiver] Device %u: %u overflow(s) before stop\n")
                         % d % pending_overflow_logs;
    }
}
void bonded_receiver::_maybe_log_rx_profile()
{
    const auto now = std::chrono::steady_clock::now();
    const double el = std::chrono::duration<double>(now - _prof_last_log).count();
    if (el < 1.0) {
        return;
    }
    _prof_last_log = now;

    const size_t calls = _prof_calls.exchange(0, std::memory_order_relaxed);
    const double wait_ms =
        _prof_wait_ns.exchange(0, std::memory_order_relaxed) / 1e6;
    const double lock_ms =
        _prof_lock_ns.exchange(0, std::memory_order_relaxed) / 1e6;

    std::ostringstream ss;
    ss << "[bonded_rx_profile] calls=" << calls
       << " wait=" << boost::format("%.1f") % (100.0 * wait_ms / (el * 1000.0)) << "%"
       << " ringlock=" << boost::format("%.1f") % (100.0 * lock_ms / (el * 1000.0)) << "%";
    for (size_t d = 0; d < _prof_recv_samples.size(); d++) {
        const double recv_msps =
            _prof_recv_samples[d].exchange(0, std::memory_order_relaxed) / el / 1e6;
        const double served_msps =
            _prof_served_samples[d].exchange(0, std::memory_order_relaxed) / el / 1e6;
        const double dropped_msps =
            _prof_dropped_samples[d].exchange(0, std::memory_order_relaxed) / el / 1e6;
        const double ring_ms =
            1e3 * _prof_ring_depth[d].load(std::memory_order_relaxed) / _config.rate;
        ss << boost::format(" | dev%u recv=%.2f served=%.2f dropped=%.2f ring=%.0fms")
                  % d % recv_msps % served_msps % dropped_msps % ring_ms;
    }
    std::cerr << ss.str() << "\n";
}

} // namespace bonded
