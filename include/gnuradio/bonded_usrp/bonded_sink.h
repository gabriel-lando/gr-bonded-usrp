#pragma once

#include <gnuradio/bonded_usrp/api.h>
#include <gnuradio/sync_block.h>
#include <cstdint>
#include <memory>
#include <string>

namespace gr {
namespace bonded_usrp {

/*!
 * \brief Bonded USRP sink block (TX side).
 *
 * Transmits two per-radio complex baseband streams through two B210s that
 * share an external 10 MHz reference and PPS, sample-aligned on a common
 * device time.  The radios stream continuously; gaps between input bursts are
 * filled with zeros so each burst starts at full power, and an underrun on
 * either radio restarts both on a common time.
 *
 *   in0: radio A (lower half, tuned at fc + radio1_offset, at rate)
 *   in1: radio B (upper half, tuned at fc + radio2_offset, at rate)
 *
 * Typical TX chain: wideband baseband -> band_splitter -> bonded_sink.
 * The device setup runs in the constructor (it takes a few seconds).
 *
 * serials_csv must list exactly two serials; the *_csv per-radio parameters
 * are empty or two comma-separated numbers.  gain is normalized [0, 1].
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

    /*!
     * Set device's (0 = A, 1 = B) phase rotation in radians while running,
     * e.g. to cancel the arbitrary relative phase of the two halves as
     * measured at a receiver.
     */
    virtual void set_phase_offset(int device, double rad) = 0;
    virtual double phase_offset(int device) const = 0;

    //! Bursts of input (runs of data between zero-fill) sent so far.
    virtual uint64_t frames_sent() const = 0;
};

} // namespace bonded_usrp
} // namespace gr
