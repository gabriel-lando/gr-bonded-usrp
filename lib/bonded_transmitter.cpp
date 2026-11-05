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
#include <cstdint>
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
    _phase_offset.assign(cfg.serials.size(), 0.0);
    for (size_t d = 0; d < cfg.phase_offset_rad.size() && d < cfg.serials.size(); d++) {
        _phase_offset[d] = cfg.phase_offset_rad[d];
    }
}

void bonded_transmitter::set_phase_offset(size_t device, double rad)
{
    std::lock_guard<std::mutex> lk(_phase_mutex);
    if (device >= _phase_offset.size()) {
        throw std::out_of_range("bonded_transmitter::set_phase_offset: bad device");
    }
    _phase_offset[device] = rad;
}

double bonded_transmitter::get_phase_offset(size_t device) const
{
    std::lock_guard<std::mutex> lk(_phase_mutex);
    return device < _phase_offset.size() ? _phase_offset[device] : 0.0;
}

bonded_transmitter::~bonded_transmitter()
{
    if (_streaming) {
        stop();
    }
    _async_stop.store(true);
    for (auto& t : _async_threads) {
        if (t.joinable()) t.join();
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
    // The async threads are already running: reset the counters in place.
    for (size_t d = 0; d < num_dev; d++) {
        _underrun_counts[d].store(0);
    }
    _resync_request.store(false);

    // One sender thread per device, fed by one scheduler thread that keeps
    // the TX front-end continuously streaming (see _send_loop_continuous).
    _feeds.clear();
    for (size_t d = 0; d < num_dev; d++) _feeds.push_back(std::make_unique<dev_feed>());
    for (size_t d = 0; d < num_dev; d++) {
        _dev_threads.emplace_back(&bonded_transmitter::_device_send_loop, this, d);
    }
    _send_thread = std::thread([this]() { _send_loop_continuous(); });

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

    // The scheduler hands what is still queued to the senders and then an
    // end marker; each sender sends everything before it and ends its burst.
    // (GNU Radio stops the blocks as soon as a finite source is exhausted, so
    // dropping the queue here would cut the end of the input.)
    _stop_flag.store(true);
    _tx_cv.notify_all();
    _sched_cv.notify_all();
    if (_send_thread.joinable()) _send_thread.join();
    for (auto& t : _dev_threads) {
        if (t.joinable()) t.join();
    }
    _dev_threads.clear();
    _feeds.clear();
    _streaming = false;
}

uint64_t bonded_transmitter::frames_sent() const
{
    return _frames_sent.load();
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
        uhd::device_addr_t dev_args(_config.device_args);
        dev_args["serial"] = serial;
        std::cout << boost::format(
                         "[bonded_transmitter] Opening serial=%s ...\n")
                         % serial;
        _devices.push_back(uhd::usrp::multi_usrp::make(dev_args));
    }

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
    // One TX chain per device (one half-band per B210).
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

void bonded_transmitter::_send_loop_continuous()
{
    // Scheduler: decides what is sent (queued data, or zeros when there is
    // none) and hands the identical sequence to every device's feed; one
    // sender thread per device (_device_send_loop) does the blocking send().
    // Every device gets the same sample count in the same order, which with
    // a common timed start keeps the halves sample-aligned. (A single thread
    // sending to the devices in turn cannot keep up past ~2x20 MS/s.)
    uhd::set_thread_priority_safe();

    const size_t num_dev = _devices.size();
    constexpr size_t kFillChunk = 4096; // samples/device per idle-fill block
    const auto zeros = std::make_shared<const std::vector<std::complex<float>>>(
        kFillChunk, std::complex<float>(0.0f, 0.0f));
    // Samples kept queued ahead per device: absorbs the linger below and host
    // scheduling hiccups.
    const size_t target = std::max<size_t>(4 * kFillChunk,
        static_cast<size_t>(_config.rate * 0.02));

    auto push_all = [&](std::vector<std::shared_ptr<const std::vector<std::complex<float>>>> bufs,
                        bool resync_marker, const uhd::time_spec_t& t) {
        for (size_t d = 0; d < num_dev; d++) {
            auto& f = *_feeds[d];
            {
                std::lock_guard<std::mutex> lk(f.mu);
                feed_item it;
                it.resync = resync_marker;
                it.start  = t;
                if (!resync_marker) {
                    it.buf = bufs[d];
                    f.queued += it.buf->size();
                }
                f.q.push_back(std::move(it));
            }
            f.cv.notify_one();
        }
    };
    auto min_queued = [&]() {
        size_t m = SIZE_MAX;
        for (auto& f : _feeds) {
            std::lock_guard<std::mutex> lk(f->mu);
            m = std::min(m, f->queued);
        }
        return m;
    };
    // Common (re)start: every device opens a new burst at the same future
    // device time; used at session start and after any underrun.
    auto schedule_start = [&](double lead_s) {
        const uhd::time_spec_t t = _devices[0]->get_time_now() + uhd::time_spec_t(lead_s);
        push_all({}, true, t);
        std::cout << boost::format("[bonded_transmitter] %s continuous TX at t=%.6f s\n")
                         % (_resyncs.load() ? "re-aligned" : "starting") % t.get_real_secs();
    };

    // Hands a batch of queued chunks to every feed, applying the per-device
    // frequency/phase correction. A frame is a run of data between zero-fills
    // (with the linger below, all chunks of one frame land in one run; frames
    // closer together than kLinger count as one).
    bool idle           = true; // last block scheduled was zero-fill (or nothing yet)
    uint64_t total_sent = 0;    // samples scheduled per device (equal on all devices)
    auto push_batch = [&](std::deque<tx_chunk>& batch) {
        if (batch.empty()) {
            return;
        }
        if (idle) {
            _frames_sent.fetch_add(1);
            idle = false;
        }
        for (auto& chunk : batch) {
            const size_t n = chunk.per_dev[0].size();
            std::vector<std::shared_ptr<const std::vector<std::complex<float>>>> bufs(num_dev);
            for (size_t d = 0; d < num_dev; d++) {
                auto& samp = chunk.per_dev[d];
                const double fcorr = (d < _config.freq_offset_hz.size())
                                          ? _config.freq_offset_hz[d]
                                          : 0.0;
                detail::apply_phase_correction(
                    samp, fcorr, get_phase_offset(d), _config.rate, total_sent);
                bufs[d] = std::make_shared<const std::vector<std::complex<float>>>(std::move(samp));
            }
            push_all(std::move(bufs), false, {});
            total_sent += n;
        }
    };

    schedule_start(0.5);

    // One frame reaches queue_samples() as several work() chunks; zero-filling
    // as soon as the queue is momentarily empty would cut a hole into it. So
    // right after data, wait briefly for the rest of the frame before
    // resuming zero-fill (the feeds hold far more than this).
    constexpr auto kLinger = std::chrono::milliseconds(2);
    auto last_data = std::chrono::steady_clock::now() - kLinger;

    while (!_stop_flag.load()) {
        if (_resync_request.exchange(false)) {
            _resyncs.fetch_add(1);
            std::cerr << "[bonded_transmitter] underrun -> re-aligning all devices "
                         "on a common start time\n";
            // The new burst must begin after everything already queued (feeds
            // and UHD's in-flight frames) has gone out on the old one.
            schedule_start(0.5);
        }

        // Pace on the feeds: only schedule more once some device runs low.
        {
            std::unique_lock<std::mutex> lk(_sched_mu);
            _sched_cv.wait_for(lk, std::chrono::milliseconds(5), [&] {
                return min_queued() < target || _stop_flag.load() || _resync_request.load();
            });
        }
        if (_stop_flag.load()) break;
        if (min_queued() >= target) continue;

        std::deque<tx_chunk> batch;
        {
            std::unique_lock<std::mutex> lk(_tx_mutex);
            const auto linger_end = last_data + kLinger;
            if (_tx_queue.empty() && std::chrono::steady_clock::now() < linger_end) {
                _tx_cv.wait_until(lk, linger_end, [this] {
                    return !_tx_queue.empty() || _stop_flag.load();
                });
            }
            batch.swap(_tx_queue);
        }

        if (batch.empty()) {
            push_all(std::vector<std::shared_ptr<const std::vector<std::complex<float>>>>(num_dev, zeros),
                false, {});
            total_sent += kFillChunk;
            idle = true;
            continue;
        }

        push_batch(batch);
        last_data = std::chrono::steady_clock::now();
    }

    // Stopping: send what is still queued, then tell every sender to finish.
    std::deque<tx_chunk> rest;
    {
        std::lock_guard<std::mutex> lk(_tx_mutex);
        rest.swap(_tx_queue);
    }
    push_batch(rest);
    for (auto& f : _feeds) {
        {
            std::lock_guard<std::mutex> lk(f->mu);
            feed_item end;
            end.end = true;
            f->q.push_back(std::move(end));
        }
        f->cv.notify_one();
    }
}

void bonded_transmitter::_device_send_loop(size_t d)
{
    uhd::set_thread_priority_safe();
    auto& f        = *_feeds[d];
    auto& streamer = _streamers[d];
    const size_t bs = streamer->get_max_num_samps();

    uhd::tx_metadata_t md;
    bool in_burst = false;
    const size_t trim = (!_config.sample_delay_trim.empty() && _config.sample_delay_trim[d] > 0)
                            ? static_cast<size_t>(_config.sample_delay_trim[d])
                            : 0;

    // Send exactly n samples, in packets of at most bs. sent == 0 while the
    // timed first packet of a burst is pending is normal until that time
    // arrives (the FIFO is full and the burst has not started); only once the
    // start time is well in the past do we fall back to an untimed send.
    auto send_all = [&](const std::complex<float>* ptr, size_t n) -> bool {
        size_t off       = 0;
        int consec_fails = 0;
        while (off < n) {
            const size_t chunk = std::min(bs, n - off);
            const std::vector<const std::complex<float>*> ptrs = {ptr + off};
            const size_t sent = streamer->send(ptrs, chunk, md, 0.1);
            if (sent == 0) {
                if (md.has_time_spec) {
                    if (_devices[d]->get_time_now() > md.time_spec + uhd::time_spec_t(0.05)) {
                        std::cerr << boost::format("[bonded_transmitter] Device %u: start "
                                                   "time missed, sending untimed\n") % d;
                        md.has_time_spec = false;
                    }
                    continue;
                }
                if (++consec_fails >= 50) {
                    std::cerr << boost::format(
                                     "[bonded_transmitter] Device %u: send() "
                                     "returned 0 repeatedly — giving up\n")
                                     % d;
                    return false;
                }
                continue;
            }
            consec_fails       = 0;
            md.start_of_burst = false;
            md.has_time_spec  = false;
            off += sent;
        }
        return true;
    };
    auto end_burst = [&]() {
        if (!in_burst) return;
        uhd::tx_metadata_t end_md;
        end_md.end_of_burst = true;
        const std::complex<float> z(0.0f, 0.0f);
        const std::vector<const std::complex<float>*> ptrs = {&z};
        streamer->send(ptrs, 0, end_md, 0.1);
        in_burst = false;
    };

    while (true) {
        feed_item it;
        {
            std::unique_lock<std::mutex> lk(f.mu);
            f.cv.wait(lk, [&] { return !f.q.empty(); });
            it = std::move(f.q.front());
            f.q.pop_front();
        }
        if (it.end) break;
        if (it.resync) {
            // Close the current burst and open the next one at the common
            // start time; the per-device trim is relative to that start.
            end_burst();
            md                = uhd::tx_metadata_t();
            md.start_of_burst = true;
            md.has_time_spec  = true;
            md.time_spec      = it.start;
            in_burst          = true;
            if (trim > 0) {
                const std::vector<std::complex<float>> pad(trim, std::complex<float>(0.0f, 0.0f));
                send_all(pad.data(), trim);
            }
            continue;
        }
        const bool ok = send_all(it.buf->data(), it.buf->size());
        {
            std::lock_guard<std::mutex> lk(f.mu);
            f.queued -= it.buf->size();
        }
        _sched_cv.notify_one();
        if (!ok) break;
    }
    end_burst();
}

void bonded_transmitter::_async_msg_loop(size_t d)
{
    while (!_async_stop.load()) {
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
                    // After an underrun the B210 resumes as soon as samples
                    // arrive, so this device has slipped against the others.
                    // Ask the scheduler to restart every device on a common
                    // time (debounced: devices often report the same event).
                    {
                        const auto now = std::chrono::steady_clock::now();
                        std::lock_guard<std::mutex> lk(_sched_mu);
                        if (now - _last_resync_req > std::chrono::milliseconds(500)) {
                            _last_resync_req = now;
                            _resync_request.store(true);
                        }
                    }
                    _sched_cv.notify_one();
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
