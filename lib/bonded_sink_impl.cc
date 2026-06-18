#include "bonded_sink_impl.h"

#include "bonded_transmitter.hpp"
#include <gnuradio/io_signature.h>

#include <sstream>

namespace gr {
namespace bonded_usrp {

bonded_sink::sptr bonded_sink::make(const std::string& serials_csv,
    double rate,
    double freq,
    double gain,
    const std::string& clock_source,
    const std::string& time_source,
    bool strict_lock,
    double lock_timeout,
    const std::string& subdev,
    const std::string& freq_plan_csv,
    const std::string& stream_args,
    int delay_trim_a,
    int delay_trim_b,
    const std::string& gain_plan_csv,
    const std::string& freq_offset_csv,
    const std::string& phase_offset_csv)
{
    return gnuradio::make_block_sptr<bonded_sink_impl>(serials_csv, rate, freq,
        gain, clock_source, time_source, strict_lock, lock_timeout, subdev,
        freq_plan_csv, stream_args, delay_trim_a, delay_trim_b, gain_plan_csv,
        freq_offset_csv, phase_offset_csv);
}

// ---------------------------------------------------------------------------

bonded_sink_impl::bonded_sink_impl(const std::string& serials_csv,
    double rate,
    double freq,
    double gain,
    const std::string& clock_source,
    const std::string& time_source,
    bool strict_lock,
    double lock_timeout,
    const std::string& subdev,
    const std::string& freq_plan_csv,
    const std::string& stream_args,
    int delay_trim_a,
    int delay_trim_b,
    const std::string& gain_plan_csv,
    const std::string& freq_offset_csv,
    const std::string& phase_offset_csv)
    : gr::sync_block("bonded_sink",
          gr::io_signature::make(2, 2, sizeof(gr_complex)),
          gr::io_signature::make(0, 0, 0)),
      _rate(rate)
{
    bonded::bonded_transmitter::config cfg;
    cfg.serials       = _parse_serials(serials_csv);
    cfg.clock_source  = clock_source;
    cfg.time_source   = time_source;
    cfg.rate          = rate;
    cfg.freq          = freq;
    cfg.gain          = gain;
    cfg.gain_plan     = _parse_freq_plan(gain_plan_csv); // generic double-CSV parse
    cfg.freq_offset_hz   = _parse_freq_plan(freq_offset_csv);
    cfg.phase_offset_rad = _parse_freq_plan(phase_offset_csv);
    cfg.subdev        = subdev;
    cfg.freq_plan     = _parse_freq_plan(freq_plan_csv);
    cfg.stream_args   = stream_args;
    cfg.strict        = strict_lock;
    cfg.lock_timeout  = lock_timeout;
    cfg.sample_delay_trim = {delay_trim_a, delay_trim_b};

    _tx = std::make_unique<bonded::bonded_transmitter>(cfg);
}

bonded_sink_impl::~bonded_sink_impl() = default;

bool bonded_sink_impl::start()
{
    std::lock_guard<std::mutex> lk(_state_mutex);
    if (_running) return true;
    _tx->configure();
    _tx->start_continuous();
    _running = true;
    return true;
}

bool bonded_sink_impl::stop()
{
    std::lock_guard<std::mutex> lk(_state_mutex);
    if (!_running) return true;
    _tx->stop();
    _running = false;
    return true;
}

int bonded_sink_impl::work(int noutput_items,
    gr_vector_const_void_star& input_items,
    gr_vector_void_star& /*output_items*/)
{
    const gr_complex* in0 = static_cast<const gr_complex*>(input_items[0]);
    const gr_complex* in1 = static_cast<const gr_complex*>(input_items[1]);

    const std::vector<const std::complex<float>*> ptrs = {
        reinterpret_cast<const std::complex<float>*>(in0),
        reinterpret_cast<const std::complex<float>*>(in1)};

    _tx->queue_samples(ptrs, static_cast<size_t>(noutput_items));

    return noutput_items;
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

std::vector<std::string> bonded_sink_impl::_parse_serials(
    const std::string& csv)
{
    std::vector<std::string> out;
    std::istringstream ss(csv);
    std::string tok;
    while (std::getline(ss, tok, ',')) {
        while (!tok.empty() && tok.front() == ' ') tok.erase(tok.begin());
        while (!tok.empty() && tok.back() == ' ') tok.pop_back();
        if (!tok.empty()) out.push_back(tok);
    }
    return out;
}

std::vector<double> bonded_sink_impl::_parse_freq_plan(const std::string& csv)
{
    if (csv.empty()) return {};
    std::vector<double> out;
    std::istringstream ss(csv);
    std::string tok;
    while (std::getline(ss, tok, ',')) {
        try {
            out.push_back(std::stod(tok));
        } catch (...) {
        }
    }
    return out;
}

} // namespace bonded_usrp
} // namespace gr
