#include "bonded_source_impl.h"

#include "bonded_receiver.hpp"
#include <gnuradio/io_signature.h>
#include <algorithm>
#include <cctype>
#include <cstring>
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
    return std::make_shared<bonded_source_impl>(serials_csv,
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
    std::vector<std::string> serials;
    std::stringstream ss(serials_csv);
    std::string item;
    while (std::getline(ss, item, ',')) {
        item.erase(std::remove_if(item.begin(),
                                  item.end(),
                                  [](unsigned char ch) { return std::isspace(ch) != 0; }),
            item.end());
        if (!item.empty()) {
            serials.push_back(item);
        }
    }
    if (serials.size() < 2) {
        throw std::invalid_argument("bonded_source requires at least 2 serials in serials_csv");
    }
    return serials;
}

std::vector<int> bonded_source_impl::_make_output_signature(size_t nports)
{
    return std::vector<int>(nports, sizeof(gr_complex));
}

std::vector<double> bonded_source_impl::_parse_freq_plan(const std::string& freq_plan_csv)
{
    std::vector<double> freqs;
    if (freq_plan_csv.empty()) {
        return freqs;
    }

    std::stringstream ss(freq_plan_csv);
    std::string item;
    while (std::getline(ss, item, ',')) {
        item.erase(std::remove_if(item.begin(),
                                  item.end(),
                                  [](unsigned char ch) { return std::isspace(ch) != 0; }),
            item.end());
        if (!item.empty()) {
            freqs.push_back(std::stod(item));
        }
    }
    return freqs;
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
                     gr::io_signature::makev(
                         static_cast<int>(_parse_serials(serials_csv).size() * 2),
                         static_cast<int>(_parse_serials(serials_csv).size() * 2),
                         _make_output_signature(_parse_serials(serials_csv).size() * 2))),
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
    cfg.freq_plan = _parse_freq_plan(freq_plan_csv);
    cfg.stream_args = stream_args;
    if (!cfg.freq_plan.empty() && cfg.freq_plan.size() != _serials.size()) {
        throw std::invalid_argument(
            "freq_plan_csv length must match number of serials");
    }
    _rx = std::make_unique<bonded::bonded_receiver>(cfg);
}

bonded_source_impl::~bonded_source_impl() {}

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
    auto result = _rx->get_aligned_samples(static_cast<size_t>(noutput_items), _fetch_timeout);
    const size_t expected_devices = _serials.size();

    auto zero_all_outputs = [&]() {
        for (auto* out : output_items) {
            std::memset(out, 0, sizeof(gr_complex) * static_cast<size_t>(noutput_items));
        }
    };

    if (result.data.size() != expected_devices) {
        zero_all_outputs();
        return noutput_items;
    }

    bool shape_ok = true;
    for (size_t d = 0; d < expected_devices; d++) {
        // bonded_receiver may provide 1 channel/device when a freq plan is
        // active. We accept 1.._channels_per_device channels and copy what is
        // available into the corresponding output ports.
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
        zero_all_outputs();
        return noutput_items;
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

        // Keep deterministic zeros on fixed output ports that do not have
        // runtime channel data (e.g., 1 channel/device with freq_plan).
        for (size_t ch = channels_available; ch < _channels_per_device; ch++) {
            const size_t port = d * _channels_per_device + ch;
            auto* out = static_cast<gr_complex*>(output_items[port]);
            std::memset(out, 0, sizeof(gr_complex) * static_cast<size_t>(noutput_items));
        }
    }

    // Per-buffer tagging is expensive at high sample rates and is currently
    // not consumed by the capture/decoder flowgraphs. Keep disabled to avoid
    // scheduler pressure and RX overflows.

    return noutput_items;
}

} // namespace bonded_usrp
} // namespace gr
