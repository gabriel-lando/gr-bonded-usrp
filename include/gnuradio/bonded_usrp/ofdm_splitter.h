#pragma once

#include <gnuradio/bonded_usrp/api.h>
#include <gnuradio/block.h>
#include <memory>
#include <string>

namespace gr {
namespace bonded_usrp {

/*!
 * \brief TX split of CP-OFDM bursts in the subcarrier domain.
 *
 *   in:   wideband CP-OFDM bursts at rate_in, centred at fc, each starting at
 *         a `len_tag_key` length tag (gr-ieee802-11 and gr-digital OFDM
 *         transmitters emit "packet_len")
 *   out0: radio1's subcarriers (the lower half) at rate_out, centred at
 *         fc + radio1_offset
 *   out1: radio2's subcarriers (the upper half) at rate_out, centred at
 *         fc + radio2_offset
 *
 * Every symbol is re-synthesized per radio at that radio's rate with only its
 * subcarriers, so each radio transmits a complete CP-OFDM signal and the two
 * halves stay orthogonal at the receiver whatever the radios' relative phase.
 * Use it instead of band_splitter for OFDM: band_splitter's split filter is
 * longer than the cyclic prefix, which costs dense constellations unless the
 * radios happen to radiate in phase.
 *
 * Requirements: fft_len and cp_len a whole number of samples at rate_out,
 * offsets a whole number of subcarriers (rate_in/fft_len), window_len <=
 * cp_len. Input outside tagged bursts is sent as silence.
 */
class BONDED_USRP_API ofdm_splitter : virtual public gr::block
{
public:
    using sptr = std::shared_ptr<ofdm_splitter>;
    static sptr make(double rate_in,
                     double rate_out,
                     double radio1_offset,
                     double radio2_offset,
                     int fft_len                    = 64,
                     int cp_len                     = 16,
                     int window_len                 = 8,
                     const std::string& len_tag_key = "packet_len");
};

} // namespace bonded_usrp
} // namespace gr
