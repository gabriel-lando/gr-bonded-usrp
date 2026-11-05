#pragma once

#include "bonded_receiver.hpp"
#include <gnuradio/bonded_usrp/bonded_source.h>
#include <chrono>
#include <cstddef>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

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

    std::vector<std::string> _serials;
    size_t _channels_per_device = 2;
    double _fetch_timeout = 1.0;

    std::chrono::steady_clock::time_point _last_fetch_log{};

    void _log_fetch_issue(const std::string& reason,
                          int noutput_items,
                          const bonded::bonded_receiver::rx_result* result = nullptr);

    bool _running = false;
    std::mutex _state_mutex;
    std::unique_ptr<bonded::bonded_receiver> _rx;
};

} // namespace bonded_usrp
} // namespace gr
