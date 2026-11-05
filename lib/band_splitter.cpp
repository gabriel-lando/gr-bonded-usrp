//
// Copyright 2026 Gabriel Lando
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Wideband band splitter, the TX counterpart of overlap_reconstructor.  The
// combined baseband stream is, per radio:
//   1. NCO-translated by -radio_offset_hz, bringing that radio's half to DC;
//   2. polyphase-resampled to per_radio_rate through a lowpass of half-width
//      |radio_offset_hz|, whose -6 dB point is therefore the centre seam (the
//      unused DC subcarrier in 802.11 and similar OFDM signals): a hard split.
//

#include "band_splitter.hpp"
#include "bonded_resample.hpp"

#include <cmath>
#include <numeric>
#include <stdexcept>

#ifdef BONDED_HAVE_VOLK
#    include <volk/volk.h>
#endif

namespace bonded {

namespace {
using namespace detail; // cf, kPi, poly_resampler
} // namespace

struct band_splitter::impl
{
    band_splitter_config cfg;
    size_t interp, decim;

    // NCO state: unit-phasor (phase accumulator) and per-sample increment.
    cf inc1, inc2;
    cf rot1{1.0f, 0.0f};
    cf rot2{1.0f, 0.0f};

    poly_resampler rs1, rs2;

    // Reusable scratch buffers.
    std::vector<cf> tmp1_, tmp2_;

    explicit impl(const band_splitter_config& c) : cfg(c)
    {
        const size_t g =
            std::gcd(static_cast<size_t>(std::llround(c.input_rate)),
                static_cast<size_t>(std::llround(c.per_radio_rate)));
        interp = static_cast<size_t>(std::llround(c.per_radio_rate)) / g;
        decim  = static_cast<size_t>(std::llround(c.input_rate)) / g;
        // Lowpass of half-width |offset| around the radio's centre: its
        // -6 dB edge is the seam at fc whatever the per-radio rate is.
        const double up = static_cast<double>(interp) * c.input_rate;
        rs1.init(interp, decim, c.taps_per_phase, std::abs(c.radio1_offset_hz) / up);
        rs2.init(interp, decim, c.taps_per_phase, std::abs(c.radio2_offset_hz) / up);

        // Negate offset: shift the radio's slice UP to DC.
        const double w1 =
            -2.0 * kPi * c.radio1_offset_hz / c.input_rate;
        const double w2 =
            -2.0 * kPi * c.radio2_offset_hz / c.input_rate;
        inc1 = cf(static_cast<float>(std::cos(w1)),
            static_cast<float>(std::sin(w1)));
        inc2 = cf(static_cast<float>(std::cos(w2)),
            static_cast<float>(std::sin(w2)));
    }

    void rotate(std::vector<cf>& x, cf& phase, cf inc)
    {
#ifdef BONDED_HAVE_VOLK
        volk_32fc_s32fc_x2_rotator2_32fc(
            reinterpret_cast<lv_32fc_t*>(x.data()),
            reinterpret_cast<const lv_32fc_t*>(x.data()),
            reinterpret_cast<const lv_32fc_t*>(&inc),
            reinterpret_cast<lv_32fc_t*>(&phase),
            static_cast<unsigned>(x.size()));
#else
        for (auto& v : x) {
            v *= phase;
            phase *= inc;
        }
        const float mag = std::abs(phase);
        if (mag > 0.f) {
            phase /= mag;
        }
#endif
    }

    void process(const cf* in,
        size_t n,
        std::vector<cf>& out1,
        std::vector<cf>& out2)
    {
        // Copy wideband into two scratch buffers for independent NCO rotation.
        tmp1_.assign(in, in + n);
        tmp2_.assign(in, in + n);

        // 1) NCO: translate each copy so its radio's slice is at DC.
        rotate(tmp1_, rot1, inc1);
        rotate(tmp2_, rot2, inc2);

        // 2) Polyphase decimate; the built-in AA filter is the split filter.
        rs1.feed(tmp1_.data(), n, out1);
        rs2.feed(tmp2_.data(), n, out2);
    }
};

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

band_splitter::band_splitter(const band_splitter_config& cfg) : _cfg(cfg)
{
    if (cfg.input_rate <= 0 || cfg.per_radio_rate <= 0) {
        throw std::invalid_argument("band_splitter: sample rates must be > 0");
    }
    if (cfg.per_radio_rate > cfg.input_rate) {
        throw std::invalid_argument(
            "band_splitter: per_radio_rate must be <= input_rate");
    }
    // Each radio's slice is [offset-|offset|, offset+|offset|] around fc,
    // i.e. 2|offset| wide; it must fit inside that radio's Nyquist band.
    if (2.0 * std::max(std::abs(cfg.radio1_offset_hz), std::abs(cfg.radio2_offset_hz))
        > cfg.per_radio_rate * (1.0 + 1e-9)) {
        throw std::invalid_argument(
            "band_splitter: 2*|radio offset| must be <= per_radio_rate");
    }
}

std::pair<std::vector<std::complex<float>>, std::vector<std::complex<float>>>
band_splitter::split(const std::vector<std::complex<float>>& wideband)
{
    reset();
    std::vector<std::complex<float>> out1, out2;
    _impl->process(wideband.data(), wideband.size(), out1, out2);
    return {std::move(out1), std::move(out2)};
}

void band_splitter::process(const std::complex<float>* in,
    size_t n,
    std::vector<std::complex<float>>& out1,
    std::vector<std::complex<float>>& out2)
{
    if (!_impl) {
        _impl = std::make_shared<impl>(_cfg);
    }
    _impl->process(in, n, out1, out2);
}

void band_splitter::reset()
{
    _impl = std::make_shared<impl>(_cfg);
}

} // namespace bonded
