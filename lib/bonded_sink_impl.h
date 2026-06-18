#pragma once

#include <gnuradio/bonded_usrp/bonded_sink.h>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

namespace bonded {
class bonded_transmitter;
}

namespace gr {
namespace bonded_usrp {

class bonded_sink_impl : public bonded_sink
{
public:
    bonded_sink_impl(const std::string& serials_csv,
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
                     const std::string& phase_offset_csv);

    ~bonded_sink_impl() override;

    bool start() override;
    bool stop() override;

    int work(int noutput_items,
             gr_vector_const_void_star& input_items,
             gr_vector_void_star& output_items) override;

private:
    static std::vector<std::string> _parse_serials(const std::string& csv);
    static std::vector<double> _parse_freq_plan(const std::string& csv);

    double _rate;
    bool _running = false;
    std::mutex _state_mutex;
    std::unique_ptr<bonded::bonded_transmitter> _tx;
};

} // namespace bonded_usrp
} // namespace gr
