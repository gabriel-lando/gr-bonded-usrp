#include "bonded_source_impl.h"

#include "bonded_csv.hpp"
#include "bonded_receiver.hpp"
#include <gnuradio/io_signature.h>
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace gr {
namespace bonded_usrp {

bonded_source::sptr bonded_source::make(const std::string& serials_csv,
                                        double rate,
                                        double freq,
                                        double gain,
                                        const std::string& clock_source,
                                        const std::string& time_source,
                                        bool strict_lock,
                                        double lock_timeout,
                                        double fetch_timeout,
                                        const std::string& subdev,
                                        const std::string& freq_plan_csv,
                                        const std::string& stream_args)
{
    return gnuradio::make_block_sptr<bonded_source_impl>(serials_csv,
                                                rate,
                                                freq,
                                                gain,
                                                clock_source,
                                                time_source,
                                                strict_lock,
                                                lock_timeout,
                                                fetch_timeout,
                                                subdev,
                                                freq_plan_csv,
                                                stream_args);
}

std::vector<std::string> bonded_source_impl::_parse_serials(const std::string& serials_csv)
{
    auto serials = bonded::detail::split_csv(serials_csv);
    if (serials.size() < 2) {
        throw std::invalid_argument("bonded_source requires at least 2 serials in serials_csv");
    }
    return serials;
}

void bonded_source_impl::_log_fetch_issue(const std::string& reason,
                                          int noutput_items,
                                          const bonded::bonded_receiver::rx_result* result)
{
    const auto now = std::chrono::steady_clock::now();
    if (_last_fetch_log.time_since_epoch().count() != 0
        && now - _last_fetch_log < std::chrono::seconds(1)) {
        return;
    }

    _last_fetch_log = now;
    std::ostringstream ss;
    ss << "[bonded_source] " << reason << " requesting " << noutput_items
       << " samples";
    if (result != nullptr) {
        if (!result->error_message.empty()) {
            ss << "; detail: " << result->error_message;
        }
        if (!result->first_timestamps.empty()) {
            ss << "; first_ts=";
            for (size_t i = 0; i < result->first_timestamps.size(); i++) {
                if (i != 0) {
                    ss << ',';
                }
                ss << result->first_timestamps[i];
            }
        }
        if (result->aligned) {
            ss << "; spread_us=" << result->inter_device_spread_us;
        }
    }
    std::cerr << ss.str() << std::endl;
}

bonded_source_impl::bonded_source_impl(const std::string& serials_csv,
                                       double rate,
                                       double freq,
                                       double gain,
                                       const std::string& clock_source,
                                       const std::string& time_source,
                                       bool strict_lock,
                                       double lock_timeout,
                                       double fetch_timeout,
                                       const std::string& subdev,
                                       const std::string& freq_plan_csv,
                                       const std::string& stream_args)
    : gr::sync_block("bonded_source",
                     gr::io_signature::make(0, 0, 0),
                     gr::io_signature::make(
                         static_cast<int>(_parse_serials(serials_csv).size() * 2),
                         static_cast<int>(_parse_serials(serials_csv).size() * 2),
                         sizeof(gr_complex))),
      _serials(_parse_serials(serials_csv)),
      _fetch_timeout(fetch_timeout)
{
    bonded::bonded_receiver::config cfg;
    cfg.serials = _serials;
    cfg.clock_source = clock_source;
    cfg.time_source = time_source;
    cfg.rate = rate;
    cfg.freq = freq;
    cfg.gain = gain;
    cfg.subdev = subdev;
    cfg.strict = strict_lock;
    cfg.lock_timeout = lock_timeout;
    cfg.freq_plan = bonded::detail::parse_csv_doubles(freq_plan_csv, "freq_plan_csv");
    cfg.stream_args = stream_args;
    if (!cfg.freq_plan.empty() && cfg.freq_plan.size() != _serials.size()) {
        throw std::invalid_argument(
            "freq_plan_csv length must match number of serials");
    }
    _rx = std::make_unique<bonded::bonded_receiver>(cfg);

    // GNU Radio's default output buffer means a few thousand samples per
    // work() call at tens of MS/s, and the scheduler handoffs then cap the
    // throughput with no thread saturated. BONDED_SRC_BUFFER overrides it.
    long out_buf = 262144;
    if (const char* env = std::getenv("BONDED_SRC_BUFFER")) {
        const long v = std::strtol(env, nullptr, 10);
        if (v > 0) {
            out_buf = v;
        }
    }
    set_min_output_buffer(out_buf);
}

bonded_source_impl::~bonded_source_impl() = default;

bool bonded_source_impl::start()
{
    std::lock_guard<std::mutex> lock(_state_mutex);
    _rx->configure();
    _rx->start_continuous();
    _running = true;
    return true;
}

bool bonded_source_impl::stop()
{
    std::lock_guard<std::mutex> lock(_state_mutex);
    if (_running) {
        _rx->stop();
        _running = false;
    }
    return true;
}

int bonded_source_impl::work(int noutput_items,
                             gr_vector_const_void_star&,
                             gr_vector_void_star& output_items)
{
    const size_t expected_devices = _serials.size();
    const auto deadline = std::chrono::steady_clock::now()
                          + std::chrono::duration<double>(_fetch_timeout);
    // A degraded block (a hole on some device) has been consumed; fetch the
    // next one right away. Returning 0 instead would leave the block idle
    // until the scheduler retries a source that produced nothing, about twice
    // a second: too rarely to drain rings that are full, which would then
    // keep every fetch degraded for the rest of the session.
    for (;;) {
        auto result = _rx->get_aligned_samples(static_cast<size_t>(noutput_items), _fetch_timeout);

        if (result.data.size() != expected_devices) {
            _log_fetch_issue("aligned fetch returned wrong device count", noutput_items, &result);
            return 0;
        }

        bool shape_ok = true;
        for (size_t d = 0; d < expected_devices; d++) {
            // A device may provide fewer channels than output ports (e.g. a
            // single-channel subdev spec); missing ports are filled with zeros.
            const size_t channels_available =
                std::min(result.data[d].size(), _channels_per_device);
            if (channels_available == 0) {
                shape_ok = false;
                break;
            }
            for (size_t ch = 0; ch < channels_available; ch++) {
                if (result.data[d][ch].size() != static_cast<size_t>(noutput_items)) {
                    shape_ok = false;
                    break;
                }
            }
            if (!shape_ok) {
                break;
            }
        }

        if (!shape_ok) {
            _log_fetch_issue("aligned fetch returned incomplete channel shape", noutput_items, &result);
            return 0;
        }

        bool degraded = false;
        for (bool success : result.success) {
            if (!success) {
                degraded = true;
                break;
            }
        }
        if (degraded) {
            _log_fetch_issue("aligned fetch degraded; not forwarding zero-filled data", noutput_items, &result);
            if (std::chrono::steady_clock::now() < deadline) {
                continue;
            }
            return 0;
        }

        for (size_t d = 0; d < expected_devices; d++) {
            const size_t channels_available =
                std::min(result.data[d].size(), _channels_per_device);
            for (size_t ch = 0; ch < channels_available; ch++) {
                const size_t port = d * _channels_per_device + ch;
                auto* out = static_cast<gr_complex*>(output_items[port]);
                std::memcpy(out,
                            result.data[d][ch].data(),
                            sizeof(gr_complex) * static_cast<size_t>(noutput_items));
            }

            for (size_t ch = channels_available; ch < _channels_per_device; ch++) {
                const size_t port = d * _channels_per_device + ch;
                auto* out = static_cast<gr_complex*>(output_items[port]);
                std::fill_n(out, static_cast<size_t>(noutput_items), gr_complex{});
            }
        }

        return noutput_items;
    }
}

} // namespace bonded_usrp
} // namespace gr
