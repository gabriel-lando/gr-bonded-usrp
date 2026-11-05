//
// Copyright 2026 Gabriel Lando
//

#pragma once

#include "bonded_api.h"
#include <complex>
#include <cstddef>
#include <memory>
#include <utility>
#include <vector>

namespace bonded {

/*!
 * Configuration for band_splitter (the TX inverse of overlap_reconstructor:
 * the wideband stream is the input, the per-radio streams the outputs).
 *
 * Each radio's half is translated to DC and filtered by a lowpass of
 * half-width |radio offset|, so the -6 dB point of both filters is the centre
 * seam (the unused DC subcarrier of 802.11 and similar OFDM signals).
 * per_radio_rate may exceed 2|offset| (oversampled radios), and |offset| may be
 * less than input_rate/4 when the signal does not fill the channel (802.11
 * uses 26 of 32 subcarriers each side), which moves the seam away from the
 * radio's own band edge.
 *
 * Classic 20 MHz plan: input_rate 20e6, per_radio_rate 10e6, offsets -/+5e6.
 */
struct band_splitter_config
{
    double input_rate       = 20e6; //!< Wideband input rate (Hz)
    double per_radio_rate   = 10e6; //!< Per-radio output rate (Hz)
    double radio1_offset_hz = -5e6; //!< radio1 center − fc (Hz), lower slice
    double radio2_offset_hz = +5e6; //!< radio2 center − fc (Hz), upper slice

    //! Split-filter length: taps_per_phase*max(L,M)+1 prototype taps.  Every
    //! subcarrier inside the filter's transition is radiated by both radios
    //! with an arbitrary relative phase (they are not phase-coherent), so it
    //! can fade; 64 keeps the transition within +-1 subcarrier of 802.11.
    size_t taps_per_phase = 64;
};

class BONDED_API band_splitter
{
public:
    explicit band_splitter(const band_splitter_config& cfg);

    /*!
     * One-shot split of a finite wideband capture (offline / unit-test).
     * Returns {slice1, slice2} at per_radio_rate.
     */
    std::pair<std::vector<std::complex<float>>, std::vector<std::complex<float>>>
    split(const std::vector<std::complex<float>>& wideband);

    /*!
     * Streaming split (real-time / GNU Radio).
     * Consumes n samples from the wideband input (at input_rate) and appends
     * the resulting per-radio samples (at per_radio_rate) to out1 and out2.
     * Call reset() to start a new stream.
     */
    void process(const std::complex<float>* in,
        size_t n,
        std::vector<std::complex<float>>& out1,
        std::vector<std::complex<float>>& out2);

    //! Clear all streaming state (starts a fresh filter history).
    void reset();

    const band_splitter_config& config() const { return _cfg; }

private:
    struct impl;

    band_splitter_config _cfg;
    std::shared_ptr<impl> _impl;
};

} // namespace bonded
