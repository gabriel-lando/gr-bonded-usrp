#pragma once

#include <gnuradio/bonded_usrp/api.h>
#include <gnuradio/sync_block.h>
#include <memory>
#include <string>

namespace gr {
namespace bonded_usrp {

/*!
 * \brief Bonded USRP sink block (TX side).
 *
 * Accepts two per-radio complex baseband streams (one per bonded B210) and
 * transmits them simultaneously via the configured TX streamers.  Internally
 * wraps uhd::usrp::bonded::bonded_transmitter, which handles 10 MHz + PPS
 * clock/time synchronisation and a common timed TX start.
 *
 *   in0: Radio A (lower half, tuned at fc + radio1_offset, at rate)
 *   in1: Radio B (upper half, tuned at fc + radio2_offset, at rate)
 *
 * Use this block together with bonded_usrp_band_splitter to form the complete
 * TX bonding path:  wifi_phy_hier → band_splitter → bonded_sink → two B210s.
 */
class BONDED_USRP_API bonded_sink : virtual public gr::sync_block
{
public:
    using sptr = std::shared_ptr<bonded_sink>;

    static sptr make(const std::string& serials_csv,
                     double rate,
                     double freq,
                     double gain,
                     const std::string& clock_source = "external",
                     const std::string& time_source  = "external",
                     bool strict_lock    = false,
                     double lock_timeout = 5.0,
                     const std::string& subdev        = "",
                     const std::string& freq_plan_csv = "",
                     const std::string& stream_args   = "",
                     int delay_trim_a = 0,
                     int delay_trim_b = 0,
                     const std::string& gain_plan_csv  = "",
                     const std::string& freq_offset_csv  = "",
                     const std::string& phase_offset_csv = "");
};

} // namespace bonded_usrp
} // namespace gr
