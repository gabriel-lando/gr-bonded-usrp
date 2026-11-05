//
// Copyright 2026 Gabriel Lando
// SPDX-License-Identifier: GPL-3.0-or-later
//

#pragma once

#include "bonded_api.h"
#include <complex>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

namespace bonded {

/*!
 * Configuration for ofdm_splitter, the TX split for CP-OFDM bursts.
 *
 * band_splitter cuts the wideband signal with a filter far longer than the
 * cyclic prefix. When the two radios do not radiate in phase, the receiver
 * sees that filter's two halves summed with the relative phase, and its long
 * response spreads every symbol into the next around the seam (EVM about
 * -8 dB on the subcarriers next to DC at 180 degrees, which 64-QAM cannot
 * survive). ofdm_splitter splits each OFDM symbol in the subcarrier domain
 * instead and synthesizes every radio's subcarriers directly at its own rate,
 * so each radio transmits a proper CP-OFDM signal of its own subcarriers. The
 * receiver's FFT then separates the two radios exactly, whatever their
 * relative phase and delay (within the cyclic prefix).
 *
 * The rate ratio must make fft_len and cp_len whole numbers of radio samples,
 * and each radio offset must be a whole number of subcarriers. The split is
 * at the subcarrier midway between the two radio centres (shared half/half
 * when it falls on a subcarrier, the unused DC subcarrier in 802.11).
 *
 * Default (20 MHz 802.11, same per-radio plan as the bonded RX): 12.5 MS/s
 * per radio at fc -/+ 3.75 MHz (12 subcarriers).
 */
struct ofdm_splitter_config
{
    double input_rate       = 20e6;    //!< Wideband input rate (Hz)
    double per_radio_rate   = 12.5e6;  //!< Per-radio output rate (Hz)
    double radio1_offset_hz = -3.75e6; //!< radio1 centre - fc (Hz), lower half
    double radio2_offset_hz = +3.75e6; //!< radio2 centre - fc (Hz), upper half
    size_t fft_len          = 64;      //!< Input FFT size (subcarriers)
    size_t cp_len           = 16;      //!< Input cyclic-prefix length
    //! Raised-cosine overlap between consecutive symbols, in input samples
    //! (rounded down to radio samples). It confines each radio's spectrum,
    //! whose edge at the seam is otherwise as sharp as the symbol boundary.
    //! The receiver's FFT window must start at least this far into the
    //! cyclic prefix; 8 of 16 leaves half of it for the channel and timing.
    size_t window_len = 8;
};

class BONDED_API ofdm_splitter
{
public:
    //! A burst of back-to-back CP-OFDM symbols, `length` input samples from
    //! `offset` (relative to the process() call's input, < n). A trailing partial
    //! symbol is ignored.
    struct burst
    {
        size_t offset = 0;
        size_t length = 0;
    };

    explicit ofdm_splitter(const ofdm_splitter_config& cfg);
    ~ofdm_splitter();

    /*!
     * Streaming split. Consumes n input samples (any chunking) and appends
     * the per-radio samples that are complete to out1 and out2 (equal
     * counts). Input outside the bursts becomes silence on the same timeline
     * (round(n * per_radio_rate / input_rate) samples). A burst's output is
     * released as soon as its last whole symbol has been consumed. `bursts`
     * must be sorted by offset; a burst that starts before the previous one
     * ended truncates it.
     */
    void process(const std::complex<float>* in,
        size_t n,
        const std::vector<burst>& bursts,
        std::vector<std::complex<float>>& out1,
        std::vector<std::complex<float>>& out2);

    //! Clear all streaming state.
    void reset();

    size_t radio_fft_len() const;    //!< fft_len at the radio rate
    size_t radio_cp_len() const;     //!< cp_len at the radio rate
    size_t radio_window_len() const; //!< window_len at the radio rate

    const ofdm_splitter_config& config() const { return _cfg; }

private:
    struct impl;
    ofdm_splitter_config _cfg;
    std::unique_ptr<impl> _impl;
};

} // namespace bonded
