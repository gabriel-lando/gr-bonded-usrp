#include "bonded_sink_impl.h"

#include "bonded_csv.hpp"
#include "bonded_transmitter.hpp"
#include <gnuradio/io_signature.h>

#include <stdexcept>

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
          gr::io_signature::make(0, 0, 0))
{
    using bonded::detail::parse_csv_doubles;
    bonded::bonded_transmitter::config cfg;
    cfg.serials = bonded::detail::split_csv(serials_csv);
    if (cfg.serials.size() != 2) {
        throw std::invalid_argument("bonded_sink requires exactly 2 serials in serials_csv");
    }
    const auto per_device = [&](const std::string& csv, const char* name) {
        auto v = parse_csv_doubles(csv, name);
        if (!v.empty() && v.size() != cfg.serials.size()) {
            throw std::invalid_argument(
                std::string("bonded_sink: ") + name + " must be empty or have 2 values");
        }
        return v;
    };
    cfg.clock_source     = clock_source;
    cfg.time_source      = time_source;
    cfg.rate             = rate;
    cfg.freq             = freq;
    cfg.gain             = gain;
    cfg.gain_plan        = per_device(gain_plan_csv, "gain_plan_csv");
    cfg.freq_offset_hz   = per_device(freq_offset_csv, "freq_offset_csv");
    cfg.phase_offset_rad = per_device(phase_offset_csv, "phase_offset_csv");
    cfg.subdev           = subdev;
    cfg.freq_plan        = per_device(freq_plan_csv, "freq_plan_csv");
    cfg.stream_args   = stream_args;
    cfg.strict        = strict_lock;
    cfg.lock_timeout  = lock_timeout;
    cfg.sample_delay_trim = {delay_trim_a, delay_trim_b};

    _tx = std::make_unique<bonded::bonded_transmitter>(cfg);

    // Device setup takes seconds, so it runs here and not in start(): GNU
    // Radio does not order start() across blocks, and upstream message
    // sources would otherwise queue seconds of frames that then go out as
    // one burst when streaming begins.
    _tx->configure();
}

bonded_sink_impl::~bonded_sink_impl() = default;

bool bonded_sink_impl::start()
{
    std::lock_guard<std::mutex> lk(_state_mutex);
    if (_running) return true;
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

void bonded_sink_impl::set_phase_offset(int device, double rad)
{
    _tx->set_phase_offset(static_cast<size_t>(device), rad);
}

double bonded_sink_impl::phase_offset(int device) const
{
    return _tx->get_phase_offset(static_cast<size_t>(device));
}

uint64_t bonded_sink_impl::frames_sent() const { return _tx->frames_sent(); }


} // namespace bonded_usrp
} // namespace gr
