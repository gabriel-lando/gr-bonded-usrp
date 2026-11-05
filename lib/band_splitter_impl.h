#pragma once

#include <gnuradio/bonded_usrp/band_splitter.h>
#include <cstddef>
#include <memory>
#include <vector>

namespace bonded {
class band_splitter;
}

namespace gr {
namespace bonded_usrp {

class band_splitter_impl : public band_splitter
{
public:
    band_splitter_impl(double rate_in,
                       double rate_out,
                       double radio1_offset,
                       double radio2_offset,
                       int taps_per_phase);
    ~band_splitter_impl() override;

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
    std::unique_ptr<bonded::band_splitter> _splitter;

    // Leftover buffers: samples produced by process() that didn't fit in the
    // current general_work output buffer.  Both ports always have equal length.
    std::vector<gr_complex> _leftover0, _leftover1;
    size_t _lead = 0;

    // Reusable work vectors to avoid per-call heap allocation.
    std::vector<gr_complex> _tmp0, _tmp1;
};

} // namespace bonded_usrp
} // namespace gr
