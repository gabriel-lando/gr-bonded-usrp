#pragma once

#include <gnuradio/bonded_usrp/overlap_reconstructor.h>
#include <chrono>
#include <cstddef>
#include <cstdint>
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
                               double radio2_offset);
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
    bool _direct_out = true;           // BONDED_DIRECT_OUT=0 restores the copy path
    // Input timeline from rx_time tags: time of the last tag and the input
    // samples consumed since, to measure the gap at the next tag.
    bool _have_time = false;
    double _tag_time = 0.0;
    uint64_t _since_tag = 0;
    uint64_t _tag_done = ~uint64_t(0); // offset of the tag already handled
    uint64_t _skipped_unlogged = 0;    // gap samples since the last log line
    size_t _gaps_unlogged      = 0;
    std::chrono::steady_clock::time_point _last_gap_log{};

    //! Steps over the gap an rx_time tag at offset `offset` reveals; returns
    //! the output samples written to `out` (the rest goes to _leftover).
    int _handle_rx_time(const gr::tag_t& tag, gr_complex* out, int room);
};

} // namespace bonded_usrp
} // namespace gr
