#include "overlap_reconstructor_impl.h"

#include "overlap_reconstructor.hpp" // bonded::overlap_reconstructor
#include <gnuradio/io_signature.h>

#include <algorithm>
#include <cmath>
#include <iostream>
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
    _lead      = 0;
    _have_time = false;
    _since_tag = 0;
    _tag_done  = ~uint64_t(0);
    return true;
}

int overlap_reconstructor_impl::_handle_rx_time(const gr::tag_t& tag,
    gr_complex* out,
    int room)
{
    _tag_done = tag.offset;
    if (!pmt::is_tuple(tag.value) || pmt::length(tag.value) != 2) {
        return 0;
    }
    const double t = static_cast<double>(pmt::to_uint64(pmt::tuple_ref(tag.value, 0)))
                     + pmt::to_double(pmt::tuple_ref(tag.value, 1));
    int produced = 0;
    if (_have_time) {
        const long long gap = std::llround((t - _tag_time) * _rate_in)
                              - static_cast<long long>(_since_tag);
        if (gap > 0) {
            _skipped_unlogged += static_cast<uint64_t>(gap);
            _gaps_unlogged++;
            const auto now = std::chrono::steady_clock::now();
            if (now - _last_gap_log >= std::chrono::seconds(1)) {
                std::cerr << "[overlap_reconstructor] stepped over " << _gaps_unlogged
                          << " gap(s) in the radio timeline, " << _skipped_unlogged
                          << " samples" << std::endl;
                _skipped_unlogged = 0;
                _gaps_unlogged    = 0;
                _last_gap_log     = now;
            }
            const std::vector<gr_complex>& z =
                _recon->skip_input(static_cast<size_t>(gap));
            const size_t take = std::min(z.size(), static_cast<size_t>(room));
            std::memcpy(out, z.data(), take * sizeof(gr_complex));
            produced = static_cast<int>(take);
            if (take < z.size()) {
                _leftover.assign(z.begin() + take, z.end());
                _lead = 0;
            }
        } else if (gap < 0) {
            std::cerr << "[overlap_reconstructor] rx_time went back " << -gap
                      << " samples; ignored" << std::endl;
        }
    }
    _have_time = true;
    _tag_time  = t;
    _since_tag = 0;
    return produced;
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
        int navail = std::min(ninput_items[0], ninput_items[1]);

        // bonded_source tags the first sample after a gap in the radio
        // timeline (a block dropped after an overflow) with rx_time. Step
        // over the gap there, so the translation phases follow the radios'
        // time; process input only up to the next such tag.
        static const pmt::pmt_t rx_time = pmt::mp("rx_time");
        const uint64_t nread = nitems_read(0);
        std::vector<gr::tag_t> tags;
        get_tags_in_range(tags, 0, nread, nread + static_cast<uint64_t>(std::max(navail, 0)), rx_time);
        std::sort(tags.begin(), tags.end(), gr::tag_t::offset_compare);
        for (const gr::tag_t& tag : tags) {
            if (tag.offset == _tag_done) {
                continue;
            }
            if (tag.offset == nread) {
                produced += _handle_rx_time(tag, out + produced, noutput_items - produced);
                continue;
            }
            navail = static_cast<int>(tag.offset - nread);
            break;
        }

        if (navail > 0 && produced < noutput_items && _leftover.empty()) {
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
    _since_tag += static_cast<uint64_t>(consumed);
    return produced;
}

} // namespace bonded_usrp
} // namespace gr
