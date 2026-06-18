//
// Copyright 2026 Gabriel Lando
// SPDX-License-Identifier: GPL-3.0-or-later
//

#include "bonded_transmitter.hpp"

#include <uhd/exception.hpp>
#include <uhd/types/tune_request.hpp>
#include <uhd/utils/thread.hpp>
#include <boost/format.hpp>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>

namespace bonded {

bonded_transmitter::bonded_transmitter(const config& cfg) : _config(cfg)
{
    if (!cfg.sample_delay_trim.empty()
        && cfg.sample_delay_trim.size() != cfg.serials.size()) {
        throw std::invalid_argument(
            "bonded_transmitter: sample_delay_trim must be empty or match "
            "the number of serials");
    }
}

bonded_transmitter::~bonded_transmitter()
{
    if (_streaming) {
        stop();
    }
}

// =========================================================================
// configure() — full setup pipeline
// =========================================================================
void bonded_transmitter::configure()
{
    if (_configured) {
        throw uhd::runtime_error("bonded_transmitter::configure() called twice");
    }
    _open_devices();
    _apply_sync();
    _configure_hardware();
    _create_streamers();
    _align_time();
    _configured = true;
}

// =========================================================================
// Burst mode
// =========================================================================
bonded_transmitter::tx_result bonded_transmitter::send_burst(
    const std::vector<std::vector<std::complex<float>>>& samples,
    double delay_sec)
{
    if (!_configured) {
        throw uhd::runtime_error(
            "bonded_transmitter::send_burst() called before configure()");
    }
    if (_streaming) {
        throw uhd::runtime_error(
            "bonded_transmitter::send_burst() called while continuous streaming");
    }
    if (samples.size() != _devices.size()) {
        throw std::invalid_argument(
            "bonded_transmitter::send_burst(): samples must have one entry "
            "per device");
    }

    const size_t num_dev = _devices.size();
    tx_result result;
    result.success.assign(num_dev, false);
    result.underrun_counts.assign(num_dev, 0);

    const uhd::time_spec_t start_time =
        _devices[0]->get_time_now() + uhd::time_spec_t(delay_sec);

    std::cout << boost::format(
                     "[bonded_transmitter] Scheduling burst at t=%.6f s\n")
                     % start_time.get_real_secs();

    // Send from each device in parallel threads.
    std::mutex print_mu;
    auto send_fn = [&](size_t d) {
        auto& streamer       = _streamers[d];
        const auto& dev_samp = samples[d];

        // Apply integer sample-delay trim: prepend zeros if trim[d] > 0.
        std::vector<std::complex<float>> padded;
        const std::complex<float>* send_ptr = dev_samp.data();
        size_t send_n                       = dev_samp.size();
        if (!_config.sample_delay_trim.empty() && _config.sample_delay_trim[d] > 0) {
            const size_t pad = static_cast<size_t>(_config.sample_delay_trim[d]);
            padded.resize(pad + dev_samp.size(), std::complex<float>(0.0f, 0.0f));
            std::copy(dev_samp.begin(), dev_samp.end(), padded.begin() + pad);
            send_ptr = padded.data();
            send_n   = padded.size();
        }

        uhd::tx_metadata_t md;
        md.start_of_burst = true;
        md.end_of_burst   = false;
        md.has_time_spec  = true;
        md.time_spec      = start_time;

        const size_t buf_sz = streamer->get_max_num_samps();
        size_t total        = 0;
        bool ok             = true;

        while (total < send_n) {
            const bool last = (total + buf_sz >= send_n);
            if (last) {
                md.end_of_burst = true;
            }
            const size_t chunk = std::min(buf_sz, send_n - total);
            const std::vector<const std::complex<float>*> ptrs = {send_ptr + total};
            const size_t n = streamer->send(ptrs, chunk, md, 1.0);
            if (n == 0) {
                std::lock_guard<std::mutex> lock(print_mu);
                std::cerr << boost::format(
                                 "[bonded_transmitter] Device %u: send() returned 0\n")
                                 % d;
                ok = false;
                break;
            }
            md.start_of_burst = false;
            md.has_time_spec  = false;
            total += n;
        }
        result.success[d] = ok;
    };

    std::vector<std::thread> threads;
    threads.reserve(num_dev);
    for (size_t d = 0; d < num_dev; d++) {
        threads.emplace_back(send_fn, d);
    }
    for (auto& t : threads) {
        t.join();
    }

    for (size_t d = 0; d < num_dev; d++) {
        result.underrun_counts[d] = _underrun_counts[d].load();
    }
    return result;
}

// =========================================================================
// Continuous mode
// =========================================================================
void bonded_transmitter::start_continuous()
{
    if (!_configured) {
        throw uhd::runtime_error(
            "bonded_transmitter::start_continuous() called before configure()");
    }
    if (_streaming) {
        throw uhd::runtime_error("bonded_transmitter: already streaming");
    }

    const size_t num_dev = _devices.size();
    _stop_flag.store(false);
    {
        std::lock_guard<std::mutex> lk(_tx_mutex);
        _tx_queue.clear();
    }
    _underrun_counts = std::vector<std::atomic<size_t>>(num_dev);
    for (size_t d = 0; d < num_dev; d++) {
        _underrun_counts[d].store(0);
    }

    // Single consumer thread: drains the coupled queue and sends each 802.11
    // frame as one timed burst on the shared PPS time base (see
    // _send_loop_burst).  The per-device sample-delay trim is applied per
    // burst inside that loop.
    _send_thread = std::thread([this]() { _send_loop_burst(); });

    _streaming = true;
}

void bonded_transmitter::queue_samples(
    const std::vector<const std::complex<float>*>& ptrs, size_t n)
{
    if (!_streaming) {
        throw uhd::runtime_error(
            "bonded_transmitter::queue_samples() called while not streaming");
    }
    if (ptrs.size() != _devices.size()) {
        throw std::invalid_argument(
            "bonded_transmitter::queue_samples(): ptrs must have one entry per "
            "device");
    }
    tx_chunk c;
    c.per_dev.resize(ptrs.size());
    for (size_t d = 0; d < ptrs.size(); d++) {
        c.per_dev[d].assign(ptrs[d], ptrs[d] + n);
    }
    {
        std::lock_guard<std::mutex> lk(_tx_mutex);
        _tx_queue.push_back(std::move(c));
    }
    _tx_cv.notify_one();
}

void bonded_transmitter::stop()
{
    if (!_streaming) return;

    _stop_flag.store(true);

    // Wake the send thread so it can exit.
    _tx_cv.notify_all();
    if (_send_thread.joinable()) _send_thread.join();

    for (auto& t : _async_threads) {
        if (t.joinable()) t.join();
    }
    _async_threads.clear();
    {
        std::lock_guard<std::mutex> lk(_tx_mutex);
        _tx_queue.clear();
    }
    _streaming = false;
}

std::vector<size_t> bonded_transmitter::get_underrun_counts() const
{
    std::vector<size_t> counts(_underrun_counts.size(), 0);
    for (size_t d = 0; d < _underrun_counts.size(); d++) {
        counts[d] = _underrun_counts[d].load();
    }
    return counts;
}

// =========================================================================
// Private helpers
// =========================================================================

void bonded_transmitter::_open_devices()
{
    std::cout << boost::format(
                     "[bonded_transmitter] Opening %zu device(s)...\n")
                     % _config.serials.size();

    for (const auto& serial : _config.serials) {
        uhd::device_addr_t dev_args;
        dev_args["serial"] = serial;
        std::cout << boost::format(
                         "[bonded_transmitter] Opening serial=%s ...\n")
                         % serial;
        _devices.push_back(uhd::usrp::multi_usrp::make(dev_args));
    }

    // Init underrun counters here (size is now known).
    _underrun_counts = std::vector<std::atomic<size_t>>(_devices.size());
    for (auto& c : _underrun_counts) c.store(0);
}

void bonded_transmitter::_apply_sync()
{
    std::cout << "[bonded_transmitter] Applying sync config to all devices...\n";

    for (size_t d = 0; d < _devices.size(); d++) {
        auto& usrp       = _devices[d];
        const size_t nmb = usrp->get_num_mboards();
        for (size_t m = 0; m < nmb; m++) {
            usrp->set_clock_source(_config.clock_source, m);
            usrp->set_time_source(_config.time_source, m);
        }

        if (_config.clock_source == "external"
            || _config.clock_source == "gpsdo") {
            const auto deadline = std::chrono::steady_clock::now()
                                  + std::chrono::duration<double>(
                                      _config.lock_timeout);
            bool locked = false;
            while (std::chrono::steady_clock::now() < deadline) {
                bool all_locked = true;
                for (size_t m = 0; m < nmb; m++) {
                    try {
                        if (!usrp->get_mboard_sensor("ref_locked", m)
                                 .to_bool()) {
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
                    + std::to_string(static_cast<int>(_config.lock_timeout))
                    + " s.";
                if (_config.strict) {
                    throw uhd::runtime_error(msg);
                }
                std::cerr << "[bonded_transmitter] WARNING: " << msg << "\n";
            } else {
                std::cout << boost::format(
                                 "[bonded_transmitter] Device %u: reference "
                                 "locked.\n")
                                 % d;
            }
        }
    }
}

void bonded_transmitter::_configure_hardware()
{
    const bool use_freq_plan = !_config.freq_plan.empty();
    if (use_freq_plan && _config.freq_plan.size() != _devices.size()) {
        throw uhd::runtime_error(
            "bonded_transmitter: freq_plan size must match number of devices");
    }
    const bool use_gain_plan = !_config.gain_plan.empty();
    if (use_gain_plan && _config.gain_plan.size() != _devices.size()) {
        throw uhd::runtime_error(
            "bonded_transmitter: gain_plan size must match number of devices");
    }

    for (size_t d = 0; d < _devices.size(); d++) {
        auto& usrp = _devices[d];
        if (!_config.subdev.empty()) {
            usrp->set_tx_subdev_spec(_config.subdev);
        }
        usrp->set_tx_rate(_config.rate);
        const double dev_freq =
            use_freq_plan ? _config.freq_plan[d] : _config.freq;
        const double dev_gain =
            use_gain_plan ? _config.gain_plan[d] : _config.gain;
        const size_t nch = usrp->get_tx_num_channels();
        for (size_t ch = 0; ch < nch; ch++) {
            usrp->set_tx_freq(uhd::tune_request_t(dev_freq), ch);
            usrp->set_normalized_tx_gain(dev_gain, ch);
        }
        std::cout << boost::format(
                         "[bonded_transmitter] Device %u: %.3f Msps, "
                         "%.3f MHz, gain=%.2f, %u ch%s\n")
                         % d % (usrp->get_tx_rate() / 1e6)
                         % (usrp->get_tx_freq(0) / 1e6) % dev_gain % nch
                         % (use_freq_plan ? " (freq_plan)" : "");
    }
    // Bonded TX uses exactly one TX chain per device (one half-band per B210).
    // Force this regardless of how many channels the device advertises.
    _channels_per_device = 1;
}

void bonded_transmitter::_create_streamers()
{
    for (size_t d = 0; d < _devices.size(); d++) {
        std::vector<size_t> ch_list(_channels_per_device);
        for (size_t c = 0; c < _channels_per_device; c++) ch_list[c] = c;

        uhd::stream_args_t sa("fc32");
        if (!_config.stream_args.empty()) {
            sa.args = uhd::device_addr_t(_config.stream_args);
        }
        sa.channels = ch_list;
        _streamers.push_back(_devices[d]->get_tx_stream(sa));
    }

    // Start async-message listener threads (one per device).
    _async_threads.reserve(_devices.size());
    for (size_t d = 0; d < _devices.size(); d++) {
        _async_threads.emplace_back(
            &bonded_transmitter::_async_msg_loop, this, d);
    }
}

void bonded_transmitter::_align_time()
{
    if (_config.time_source == "external"
        || _config.time_source == "gpsdo") {
        // Re-apply after MCR changes (same B210 quirk as on the RX side).
        for (auto& dev : _devices) {
            dev->set_time_source(_config.time_source);
        }

        std::cout << "[bonded_transmitter] Waiting for FPGA to stabilise...\n";
        std::this_thread::sleep_for(std::chrono::milliseconds(500));

        std::cout << "[bonded_transmitter] Waiting for PPS edge on device 0...\n";
        const auto pps_deadline =
            std::chrono::steady_clock::now() + std::chrono::milliseconds(2500);
        uhd::time_spec_t last_pps = _devices[0]->get_time_last_pps();
        while (_devices[0]->get_time_last_pps() == last_pps) {
            if (std::chrono::steady_clock::now() > pps_deadline) {
                throw uhd::runtime_error(
                    "No PPS detected on device 0 within 2.5 s.\n"
                    "Check PPS cable to all devices.");
            }
            std::this_thread::sleep_for(std::chrono::milliseconds(1));
        }

        std::cout << "[bonded_transmitter] PPS detected — arming all devices...\n";
        for (auto& dev : _devices) {
            dev->set_time_next_pps(uhd::time_spec_t(0.0));
        }
        std::this_thread::sleep_for(std::chrono::seconds(1));
        std::cout << "[bonded_transmitter] PPS alignment done.\n";
    } else {
        std::cout << "[bonded_transmitter] Internal time: set_time_now(0).\n";
        for (auto& dev : _devices) {
            dev->set_time_now(uhd::time_spec_t(0.0),
                uhd::usrp::multi_usrp::ALL_MBOARDS);
        }
    }
}

void bonded_transmitter::_send_loop_burst()
{
    uhd::set_thread_priority_safe();

    const size_t num_dev = _devices.size();

    // A frame boundary is a lull in the queue.  Consecutive GR work() chunks
    // for the same 802.11 frame arrive < ~200 µs apart; the inter-frame strobe
    // gap is hundreds of ms.  A few ms cleanly separates the two cases.
    constexpr auto kFrameGap = std::chrono::milliseconds(5);
    // Schedule each burst this far in the future so its samples are copied into
    // the FPGA TX FIFO before the timestamp arrives.  A frame is only a few
    // thousand samples, so 30 ms is ample and stays well inside the ~300 ms
    // inter-frame gap.
    const double kBurstDelay = 0.03; // seconds

    // Keeps consecutive bursts strictly ordered even if two frames happen to
    // arrive closer together than kBurstDelay.
    uhd::time_spec_t last_burst_end(0.0);

    // Per-device frame accumulator, reused across frames.
    std::vector<std::vector<std::complex<float>>> frame(num_dev);

    while (!_stop_flag.load()) {
        // ---- Wait for the first chunk of the next frame --------------------
        {
            std::unique_lock<std::mutex> lk(_tx_mutex);
            _tx_cv.wait(lk, [&]() {
                return !_tx_queue.empty() || _stop_flag.load();
            });
            if (_stop_flag.load() && _tx_queue.empty()) break;
        }

        // ---- Start a fresh frame; prepend per-device delay-trim zeros ------
        // Applying the trim here (per burst) means it survives every frame,
        // unlike a one-time pre-queue that drifts away.
        for (size_t d = 0; d < num_dev; d++) {
            frame[d].clear();
            if (!_config.sample_delay_trim.empty()
                && _config.sample_delay_trim[d] > 0) {
                frame[d].assign(
                    static_cast<size_t>(_config.sample_delay_trim[d]),
                    std::complex<float>(0.0f, 0.0f));
            }
        }

        // ---- Accumulate coupled chunks until a gap ends the frame ----------
        bool got_data = false;
        while (!_stop_flag.load()) {
            std::deque<tx_chunk> batch;
            {
                std::unique_lock<std::mutex> lk(_tx_mutex);
                if (_tx_queue.empty()) {
                    _tx_cv.wait_for(lk, kFrameGap, [&]() {
                        return !_tx_queue.empty() || _stop_flag.load();
                    });
                }
                if (_tx_queue.empty()) {
                    break; // lull → frame complete
                }
                batch.swap(_tx_queue);
            }
            for (auto& chunk : batch) {
                for (size_t d = 0; d < num_dev; d++) {
                    auto& src = chunk.per_dev[d];
                    frame[d].insert(frame[d].end(), src.begin(), src.end());
                }
                got_data = true;
            }
        }

        if (!got_data) continue;

        // ---- Emit the frame as ONE timed burst on every device ------------
        // All devices share the same time_spec, so the half-bands are sample-
        // aligned over the air regardless of send order, scheduling jitter, or
        // any underrun in the preceding gap.
        const uhd::time_spec_t now_ts = _devices[0]->get_time_now();
        uhd::time_spec_t start = now_ts + uhd::time_spec_t(kBurstDelay);
        if (start < last_burst_end) {
            start = last_burst_end;
        }

        std::cerr << boost::format(
                         "[bonded_transmitter] burst: %zu samp/dev  "
                         "now=%.4f  start=%.4f\n")
                         % frame[0].size() % now_ts.get_real_secs()
                         % start.get_real_secs();

        const double t0 = start.get_real_secs();
        size_t max_len = 0;
        for (size_t d = 0; d < num_dev; d++) {
            const size_t total = frame[d].size();

            // Inter-radio phase-coherence correction: rotate this device's
            // samples by a digital NCO referenced to the shared-clock burst
            // time t0, so it cancels relative LO drift across the inter-frame
            // gap as well as within the burst.
            const double fcorr =
                (d < _config.freq_offset_hz.size()) ? _config.freq_offset_hz[d] : 0.0;
            const double pcorr =
                (d < _config.phase_offset_rad.size()) ? _config.phase_offset_rad[d] : 0.0;
            if (fcorr != 0.0 || pcorr != 0.0) {
                for (size_t k = 0; k < total; k++) {
                    const double ph =
                        2.0 * M_PI * fcorr * (t0 + static_cast<double>(k) / _config.rate)
                        + pcorr;
                    frame[d][k] *= std::complex<float>(
                        static_cast<float>(std::cos(ph)), static_cast<float>(std::sin(ph)));
                }
            }

            const std::complex<float>* ptr = frame[d].data();
            max_len                        = std::max(max_len, total);
            const size_t bs                = _streamers[d]->get_max_num_samps();

            uhd::tx_metadata_t md;
            md.start_of_burst = true;
            md.end_of_burst   = false;
            md.has_time_spec  = true;
            md.time_spec      = start;

            size_t off = 0;
            while (off < total && !_stop_flag.load()) {
                const size_t n   = std::min(bs, total - off);
                md.end_of_burst  = (off + n >= total);
                const std::vector<const std::complex<float>*> ptrs = {ptr + off};
                _streamers[d]->send(ptrs, n, md, 0.1);
                md.start_of_burst = false;
                md.has_time_spec  = false;
                off += n;
            }
        }
        last_burst_end =
            start + uhd::time_spec_t(static_cast<double>(max_len) / _config.rate);
    }
}

void bonded_transmitter::_async_msg_loop(size_t d)
{
    while (!_stop_flag.load()) {
        uhd::async_metadata_t amd;
        if (_streamers[d]->recv_async_msg(amd, 0.05)) {
            switch (amd.event_code) {
                case uhd::async_metadata_t::EVENT_CODE_UNDERFLOW:
                case uhd::async_metadata_t::EVENT_CODE_UNDERFLOW_IN_PACKET:
                    _underrun_counts[d].fetch_add(1);
                    std::cerr << boost::format(
                                     "[bonded_transmitter] Device %u: underrun "
                                     "(U)\n")
                                     % d;
                    break;
                case uhd::async_metadata_t::EVENT_CODE_TIME_ERROR:
                    std::cerr << boost::format(
                                     "[bonded_transmitter] Device %u: burst "
                                     "scheduled in the past (L/time error)\n")
                                     % d;
                    break;
                case uhd::async_metadata_t::EVENT_CODE_SEQ_ERROR:
                case uhd::async_metadata_t::EVENT_CODE_SEQ_ERROR_IN_BURST:
                    std::cerr << boost::format(
                                     "[bonded_transmitter] Device %u: sequence "
                                     "error\n")
                                     % d;
                    break;
                default:
                    break;
            }
        }
    }
}

} // namespace bonded
