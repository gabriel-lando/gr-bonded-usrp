#pragma once

#include <gnuradio/bonded_usrp/api.h>
#include <gnuradio/block.h>
#include <memory>

namespace gr {
namespace bonded_usrp {

/*!
 * \brief Wideband band-splitter block (TX side).
 *
 * Splits one wideband complex baseband stream (at rate_in) at its centre into
 * two per-radio streams (each at rate_out):
 *   in:   wideband at rate_in, centred at fc
 *   out0: lower half at rate_out, centred at fc + radio1_offset
 *   out1: upper half at rate_out, centred at fc + radio2_offset
 *
 * Each radio's filter is a lowpass of half-width |offset|, so the -6 dB point
 * of both is the centre seam (the unused DC subcarrier in 802.11).
 * taps_per_phase sets that filter's length (taps_per_phase*max(L,M)+1) and so
 * how many subcarriers around the seam both radios radiate.  Requires
 * 2*|offset| <= rate_out.
 */
class BONDED_USRP_API band_splitter : virtual public gr::block
{
public:
    using sptr = std::shared_ptr<band_splitter>;

    static sptr make(double rate_in,
                     double rate_out,
                     double radio1_offset,
                     double radio2_offset,
                     int taps_per_phase = 64);
};

} // namespace bonded_usrp
} // namespace gr
