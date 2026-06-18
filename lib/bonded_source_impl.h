#pragma once

#include <gnuradio/bonded_usrp/bonded_source.h>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

namespace bonded {
class bonded_receiver;
}

namespace gr {
namespace bonded_usrp {

class bonded_source_impl : public bonded_source
{
public:
    bonded_source_impl(const std::string& serials_csv,
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
                       const std::string& stream_args);

    ~bonded_source_impl() override;

    bool start() override;
    bool stop() override;

    int work(int noutput_items,
             gr_vector_const_void_star& input_items,
             gr_vector_void_star& output_items) override;

private:
    static std::vector<std::string> _parse_serials(const std::string& serials_csv);
    static std::vector<double> _parse_freq_plan(const std::string& freq_plan_csv);
    static std::vector<int> _make_output_signature(size_t nports);

    std::vector<std::string> _serials;
    size_t _channels_per_device = 2;
    double _fetch_timeout = 1.0;

    bool _running = false;
    std::mutex _state_mutex;
    std::unique_ptr<bonded::bonded_receiver> _rx;
};

} // namespace bonded_usrp
} // namespace gr
