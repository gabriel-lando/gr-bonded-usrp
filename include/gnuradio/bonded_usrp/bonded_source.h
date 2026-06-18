#pragma once

#include <gnuradio/bonded_usrp/api.h>
#include <gnuradio/sync_block.h>
#include <memory>
#include <string>

namespace gr {
namespace bonded_usrp {

class BONDED_USRP_API bonded_source : virtual public gr::sync_block
{
public:
    using sptr = std::shared_ptr<bonded_source>;

    static sptr make(const std::string& serials_csv,
                     double rate,
                     double freq,
                     double gain,
                     const std::string& clock_source = "external",
                     const std::string& time_source = "external",
                     bool strict_lock = false,
                     double lock_timeout = 5.0,
                     double fetch_timeout = 1.0,
                     const std::string& subdev = "",
                     const std::string& freq_plan_csv = "",
                     const std::string& stream_args = "");
};

} // namespace bonded_usrp
} // namespace gr
