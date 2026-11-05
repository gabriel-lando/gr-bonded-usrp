#include "overlap_reconstructor_impl.h"

#include "overlap_reconstructor.hpp" // bonded::overlap_reconstructor
#include <gnuradio/io_signature.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <string>

namespace gr {
namespace bonded_usrp {

overlap_reconstructor::sptr overlap_reconstructor::make(double rate_in,
    double rate_out,
    double target_bw,
    double overlap_bw,
    double radio1_offset,
    double radio2_offset)
{
    return gnuradio::make_block_sptr<overlap_reconstructor_impl>(rate_in,
        rate_out, target_bw, overlap_bw, radio1_offset, radio2_offset);
}

overlap_reconstructor_impl::overlap_reconstructor_impl(double rate_in,
    double rate_out,
    double target_bw,
    double overlap_bw,
    double radio1_offset,
    double radio2_offset)
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
    _recon = std::make_unique<bonded::overlap_reconstructor>(cfg);

    set_relative_rate(rate_out / rate_in);

    // Per-call overhead only amortizes when general_work() is handed many
    // output samples at once: ask for multiples of a large chunk and an output
    // buffer several chunks deep, so a slower downstream block (the 802.11 PHY)
    // does not stall us immediately.
    size_t chunk = 8192;
    if (const char* env = std::getenv("BONDED_CHUNK")) {
        const long v = std::strtol(env, nullptr, 10);
        if (v > 0) {
            chunk = static_cast<size_t>(v);
        }
    }
    if (const char* env = std::getenv("BONDED_DIRECT_OUT")) {
        _direct_out = (std::string(env) != "0");
    }
    set_output_multiple(static_cast<int>(chunk));
    set_min_output_buffer(0, static_cast<long>(chunk * 8));
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

    // 2) Reconstruct new input, consuming only about as much as fits in the
    //    remaining output room. Taking everything available (up to
    //    bonded_source's 262144-item buffer) would park most of the output in
    //    _leftover and make this block and the one downstream run in series.
    int consumed = 0;
    if (produced < noutput_items && _leftover.empty()) {
        const int navail = std::min(ninput_items[0], ninput_items[1]);
        if (navail > 0) {
            const size_t room = static_cast<size_t>(noutput_items - produced);
            const int want =
                static_cast<int>(std::ceil(room * _rate_in / _rate_out));
            consumed = std::min(navail, std::max(want, 1));
            if (_direct_out) {
                produced += static_cast<int>(_recon->process_into(
                    in0, in1, static_cast<size_t>(consumed), out + produced, room));
                const std::vector<gr_complex>& rest = _recon->pending();
                if (!rest.empty()) {
                    _leftover.assign(rest.begin(), rest.end());
                    _lead = 0;
                }
            } else {
                // Copy path (BONDED_DIRECT_OUT=0), for A/B measurements.
                const std::vector<gr_complex>& res =
                    _recon->process_view(in0, in1, static_cast<size_t>(consumed));
                const size_t take = std::min(res.size(), room);
                std::memcpy(out + produced, res.data(), take * sizeof(gr_complex));
                produced += static_cast<int>(take);
                if (take < res.size()) {
                    _leftover.assign(res.begin() + take, res.end());
                    _lead = 0;
                }
            }
        }
    }
    consume_each(consumed);
    return produced;
}

} // namespace bonded_usrp
} // namespace gr
