#include "overlap_reconstructor_impl.h"

#include "overlap_reconstructor.hpp" // bonded::overlap_reconstructor
#include <gnuradio/io_signature.h>

#include <algorithm>
#include <cmath>
#include <cstring>

namespace gr {
namespace bonded_usrp {

overlap_reconstructor::sptr overlap_reconstructor::make(double rate_in,
    double rate_out,
    double target_bw,
    double overlap_bw,
    double radio1_offset,
    double radio2_offset,
    int fft_size)
{
    return gnuradio::make_block_sptr<overlap_reconstructor_impl>(rate_in,
        rate_out, target_bw, overlap_bw, radio1_offset, radio2_offset, fft_size);
}

overlap_reconstructor_impl::overlap_reconstructor_impl(double rate_in,
    double rate_out,
    double target_bw,
    double overlap_bw,
    double radio1_offset,
    double radio2_offset,
    int fft_size)
    : gr::block("overlap_reconstructor",
          gr::io_signature::make(2, 2, sizeof(gr_complex)),
          gr::io_signature::make(1, 1, sizeof(gr_complex))),
      _rate_in(rate_in),
      _rate_out(rate_out)
{
    bonded::overlap_reconstructor_config cfg;
    cfg.target_rate      = rate_out;
    cfg.per_radio_rate   = rate_in;
    cfg.target_bw        = target_bw;
    cfg.overlap_bw       = overlap_bw;
    cfg.radio1_offset_hz = radio1_offset;
    cfg.radio2_offset_hz = radio2_offset;
    cfg.fft_size         = static_cast<size_t>(fft_size);
    _recon = std::make_unique<bonded::overlap_reconstructor>(cfg);

    set_relative_rate(rate_out / rate_in);
}

overlap_reconstructor_impl::~overlap_reconstructor_impl() = default;

bool overlap_reconstructor_impl::start()
{
    _recon->reset();
    _leftover.clear();
    _lead = 0;
    return true;
}

void overlap_reconstructor_impl::forecast(int noutput_items,
    gr_vector_int& ninput_items_required)
{
    // Roughly rate_in/rate_out input samples per output; the leftover queue
    // absorbs the variable per-frame yield.
    int need = static_cast<int>(std::ceil(noutput_items * _rate_in / _rate_out));
    need = std::max(1, need);
    for (auto& r : ninput_items_required) {
        r = need;
    }
}

int overlap_reconstructor_impl::general_work(int noutput_items,
    gr_vector_int& ninput_items,
    gr_vector_const_void_star& input_items,
    gr_vector_void_star& output_items)
{
    const gr_complex* in0 = static_cast<const gr_complex*>(input_items[0]);
    const gr_complex* in1 = static_cast<const gr_complex*>(input_items[1]);
    gr_complex* out        = static_cast<gr_complex*>(output_items[0]);

    int produced = 0;

    // 1) Drain previously reconstructed samples (bulk copy).
    const size_t avail = _leftover.size() - _lead;
    if (avail > 0) {
        const size_t take = std::min(avail, static_cast<size_t>(noutput_items));
        std::memcpy(out, _leftover.data() + _lead, take * sizeof(gr_complex));
        produced += static_cast<int>(take);
        _lead += take;
        if (_lead == _leftover.size()) {
            _leftover.clear();
            _lead = 0;
        }
    }

    // 2) Reconstruct from the input available this call (the time-domain
    //    reconstructor streams correctly on any chunk size).
    int consumed = 0;
    if (produced < noutput_items && _leftover.empty()) {
        const int navail = std::min(ninput_items[0], ninput_items[1]);
        if (navail > 0) {
            const std::vector<gr_complex> chunk =
                _recon->process(in0, in1, static_cast<size_t>(navail));
            consumed = navail;
            const size_t room = static_cast<size_t>(noutput_items - produced);
            const size_t take = std::min(room, chunk.size());
            std::memcpy(out + produced, chunk.data(), take * sizeof(gr_complex));
            produced += static_cast<int>(take);
            if (take < chunk.size()) {
                _leftover.assign(chunk.begin() + static_cast<std::ptrdiff_t>(take),
                    chunk.end());
                _lead = 0;
            }
        }
    }
    consume_each(consumed);
    return produced;
}

} // namespace bonded_usrp
} // namespace gr
