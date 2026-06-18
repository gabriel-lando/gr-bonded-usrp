#pragma once

#include <gnuradio/bonded_usrp/overlap_reconstructor.h>
#include <cstddef>
#include <memory>
#include <vector>

namespace bonded {
class overlap_reconstructor;
}

namespace gr {
namespace bonded_usrp {

class overlap_reconstructor_impl : public overlap_reconstructor
{
public:
    overlap_reconstructor_impl(double rate_in,
                               double rate_out,
                               double target_bw,
                               double overlap_bw,
                               double radio1_offset,
                               double radio2_offset,
                               int fft_size);
    ~overlap_reconstructor_impl() override;

    bool start() override;

    void forecast(int noutput_items, gr_vector_int& ninput_items_required) override;

    int general_work(int noutput_items,
                      gr_vector_int& ninput_items,
                      gr_vector_const_void_star& input_items,
                      gr_vector_void_star& output_items) override;

private:
    double _rate_in;
    double _rate_out;
    std::unique_ptr<bonded::overlap_reconstructor> _recon;
    std::vector<gr_complex> _leftover; // reconstructed samples not yet emitted
    size_t _lead = 0;                  // read offset into _leftover
};

} // namespace bonded_usrp
} // namespace gr
