// GNU Radio wrapper for bonded::ofdm_splitter: turns the length tags into
// bursts and feeds the engine no more input than the output room can hold.
//

#include "ofdm_splitter_impl.h"
#include "ofdm_splitter.hpp" // bonded::ofdm_splitter
#include <gnuradio/io_signature.h>

#include <algorithm>
#include <cmath>
#include <cstring>

namespace gr {
namespace bonded_usrp {

ofdm_splitter::sptr ofdm_splitter::make(double rate_in,
    double rate_out,
    double radio1_offset,
    double radio2_offset,
    int fft_len,
    int cp_len,
    int window_len,
    const std::string& len_tag_key)
{
    return gnuradio::make_block_sptr<ofdm_splitter_impl>(rate_in, rate_out,
        radio1_offset, radio2_offset, fft_len, cp_len, window_len, len_tag_key);
}

ofdm_splitter_impl::ofdm_splitter_impl(double rate_in,
    double rate_out,
    double radio1_offset,
    double radio2_offset,
    int fft_len,
    int cp_len,
    int window_len,
    const std::string& len_tag_key)
    : gr::block("ofdm_splitter",
          gr::io_signature::make(1, 1, sizeof(gr_complex)),
          gr::io_signature::make(2, 2, sizeof(gr_complex))),
      _rate_in(rate_in),
      _rate_out(rate_out),
      _len_tag_key(pmt::intern(len_tag_key))
{
    bonded::ofdm_splitter_config cfg;
    cfg.input_rate       = rate_in;
    cfg.per_radio_rate   = rate_out;
    cfg.radio1_offset_hz = radio1_offset;
    cfg.radio2_offset_hz = radio2_offset;
    cfg.fft_len          = static_cast<size_t>(std::max(0, fft_len));
    cfg.cp_len           = static_cast<size_t>(std::max(0, cp_len));
    cfg.window_len       = static_cast<size_t>(std::max(0, window_len));
    _splitter            = std::make_unique<bonded::ofdm_splitter>(cfg);
    set_relative_rate(rate_out / rate_in);
    set_tag_propagation_policy(TPP_DONT);
}

ofdm_splitter_impl::~ofdm_splitter_impl() = default;

bool ofdm_splitter_impl::start()
{
    _splitter->reset();
    _leftover0.clear();
    _leftover1.clear();
    _lead = 0;
    return true;
}

void ofdm_splitter_impl::forecast(int noutput_items,
    gr_vector_int& ninput_items_required)
{
    // As in band_splitter: no input is required while output is pending, or
    // a bursty source would leave it parked until the next burst.
    int need = 0;
    if (_leftover0.size() == _lead) {
        need = std::max(1,
            static_cast<int>(std::ceil(noutput_items * _rate_in / _rate_out)));
    }
    for (auto& r : ninput_items_required) {
        r = need;
    }
}

int ofdm_splitter_impl::general_work(int noutput_items,
    gr_vector_int& ninput_items,
    gr_vector_const_void_star& input_items,
    gr_vector_void_star& output_items)
{
    const gr_complex* in = static_cast<const gr_complex*>(input_items[0]);
    gr_complex* out0     = static_cast<gr_complex*>(output_items[0]);
    gr_complex* out1     = static_cast<gr_complex*>(output_items[1]);
    int produced         = 0;

    const size_t avail = _leftover0.size() - _lead;
    if (avail > 0) {
        const size_t take = std::min(avail, static_cast<size_t>(noutput_items));
        std::memcpy(out0, _leftover0.data() + _lead, take * sizeof(gr_complex));
        std::memcpy(out1, _leftover1.data() + _lead, take * sizeof(gr_complex));
        produced += static_cast<int>(take);
        _lead += take;
        if (_lead == _leftover0.size()) {
            _leftover0.clear();
            _leftover1.clear();
            _lead = 0;
        }
    }

    int consumed = 0;
    if (produced < noutput_items && _leftover0.empty()) {
        const size_t room = static_cast<size_t>(noutput_items - produced);
        const size_t want = static_cast<size_t>(
            std::ceil(static_cast<double>(room) * _rate_in / _rate_out));
        const size_t navail = std::min(static_cast<size_t>(ninput_items[0]), want);
        if (navail > 0) {
            const uint64_t nread = nitems_read(0);
            std::vector<gr::tag_t> tags;
            get_tags_in_range(tags, 0, nread, nread + navail, _len_tag_key);
            std::sort(tags.begin(), tags.end(), gr::tag_t::offset_compare);
            std::vector<bonded::ofdm_splitter::burst> bursts;
            for (const auto& t : tags) {
                if (pmt::is_integer(t.value)) {
                    bursts.push_back({ static_cast<size_t>(t.offset - nread),
                        static_cast<size_t>(std::max<long>(0, pmt::to_long(t.value))) });
                }
            }
            _tmp0.clear();
            _tmp1.clear();
            _splitter->process(in, navail, bursts, _tmp0, _tmp1);
            consumed = static_cast<int>(navail);
            const size_t take = std::min(room, _tmp0.size());
            if (take > 0) {
                std::memcpy(out0 + produced, _tmp0.data(), take * sizeof(gr_complex));
                std::memcpy(out1 + produced, _tmp1.data(), take * sizeof(gr_complex));
                produced += static_cast<int>(take);
            }
            if (take < _tmp0.size()) {
                _leftover0.assign(_tmp0.begin() + static_cast<std::ptrdiff_t>(take), _tmp0.end());
                _leftover1.assign(_tmp1.begin() + static_cast<std::ptrdiff_t>(take), _tmp1.end());
                _lead = 0;
            }
        }
    }

    consume_each(consumed);
    return produced;
}

} // namespace bonded_usrp
} // namespace gr
