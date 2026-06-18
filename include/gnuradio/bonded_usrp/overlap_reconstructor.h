#pragma once

#include <gnuradio/bonded_usrp/api.h>
#include <gnuradio/block.h>
#include <memory>

namespace gr {
namespace bonded_usrp {

/*!
 * \brief Overlap-aware wideband reconstruction block.
 *
 * Takes the two overlapping radio slices (channel 0 of each bonded B210, at
 * rate_in) and reconstructs the contiguous target channel (at rate_out) by
 * aligning and crossfading the shared overlap region and discarding the
 * out-of-band corners.  All DSP lives in uhd::usrp::bonded::overlap_reconstructor.
 *
 *   in0: radio1 slice (tuned at fc + radio1_offset)
 *   in1: radio2 slice (tuned at fc + radio2_offset)
 *   out: reconstructed baseband centered at fc, width target_bw
 */
class BONDED_USRP_API overlap_reconstructor : virtual public gr::block
{
public:
    using sptr = std::shared_ptr<overlap_reconstructor>;

    static sptr make(double rate_in,
                     double rate_out,
                     double target_bw,
                     double overlap_bw,
                     double radio1_offset,
                     double radio2_offset,
                     int fft_size = 4096);
};

} // namespace bonded_usrp
} // namespace gr
