#pragma once

#include <gnuradio/bonded_usrp/ofdm_splitter.h>
#include <pmt/pmt.h>
#include <cstddef>
#include <memory>
#include <vector>

namespace bonded {
class ofdm_splitter;
}

namespace gr {
namespace bonded_usrp {

class ofdm_splitter_impl : public ofdm_splitter
{
public:
    ofdm_splitter_impl(double rate_in,
                       double rate_out,
                       double radio1_offset,
                       double radio2_offset,
                       int fft_len,
                       int cp_len,
                       int window_len,
                       const std::string& len_tag_key);
    ~ofdm_splitter_impl() override;

    bool start() override;
    void forecast(int noutput_items,
                  gr_vector_int& ninput_items_required) override;
    int general_work(int noutput_items,
                     gr_vector_int& ninput_items,
                     gr_vector_const_void_star& input_items,
                     gr_vector_void_star& output_items) override;

private:
    double _rate_in;
    double _rate_out;
    pmt::pmt_t _len_tag_key;
    std::unique_ptr<bonded::ofdm_splitter> _splitter;
    // Output that did not fit in the current call (both ports equal length).
    std::vector<gr_complex> _leftover0, _leftover1;
    size_t _lead = 0;
    std::vector<gr_complex> _tmp0, _tmp1;
};

} // namespace bonded_usrp
} // namespace gr
