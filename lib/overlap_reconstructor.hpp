//
// Copyright 2026 Gabriel Lando
// SPDX-License-Identifier: GPL-3.0-or-later
//

#pragma once

#include "bonded_api.h"
#include <complex>
#include <cstddef>
#include <memory>
#include <vector>

namespace bonded {

/*!
 * Overlap-aware wideband reconstruction for two bonded radios.
 *
 * Two radios each capture a wider slice of a target channel with an
 * intentional spectral overlap (instead of a hard split):
 *
 *   target: a channel of width target_bw centered at fc
 *   radio1: tuned at fc + radio1_offset_hz, captures per_radio_bw
 *   radio2: tuned at fc + radio2_offset_hz, captures per_radio_bw
 *
 * The shared overlap (width overlap_bw, centered on fc) is used to align the
 * two radios (relative gain, phase and fractional delay) and to crossfade the
 * two complex spectra.  The out-of-band corners are discarded.  The result is a
 * single contiguous complex baseband at target_rate.
 *
 * Pipeline (all complex / phase preserving — suitable for decoding):
 *   1. rational-resample each radio per_radio_rate -> target_rate
 *   2. per FFT frame: translate each radio to the target baseband by an integer
 *      bin shift, estimate/apply the radio2->radio1 alignment over the overlap,
 *      crossfade the overlap, discard corners, IFFT (overlap-save).
 *
 * The class is pure C++ (self-contained FFT, no GNU Radio / FFTW dependency) so
 * it can be unit-tested under ctest and wrapped by a GNU Radio block.
 */

struct overlap_reconstructor_config
{
    double target_rate      = 20e6;  //!< Output sample rate (Hz)
    double per_radio_rate   = 15e6;  //!< Per-radio input sample rate (Hz)
    double target_bw        = 20e6;  //!< Reconstructed channel bandwidth (Hz)
    double overlap_bw       = 5e6;   //!< Shared overlap bandwidth (Hz)
    double radio1_offset_hz = -5e6;  //!< radio1 center - fc (Hz), lower slice
    double radio2_offset_hz = +5e6;  //!< radio2 center - fc (Hz), upper slice
    size_t fft_size         = 4096;  //!< FFT size (power of two, multiple of 4)
    bool estimate_alignment = true;  //!< Estimate gain/phase/delay from overlap
};

//! Alignment estimate of radio2 relative to radio1 (over the overlap band).
struct overlap_alignment
{
    double gain        = 1.0;  //!< |radio1/radio2| amplitude ratio
    double phase_rad   = 0.0;  //!< constant phase offset (rad)
    double delay_samps = 0.0;  //!< fractional delay at target_rate (samples)
    bool valid         = false;
};

class BONDED_API overlap_reconstructor
{
public:
    explicit overlap_reconstructor(const overlap_reconstructor_config& cfg);

    /*!
     * One-shot reconstruction of a finite capture (offline / unit test).
     * rx1, rx2 are equal-duration complex baseband captures at per_radio_rate.
     * Returns the combined baseband at target_rate.  Performs a two-pass
     * alignment estimate over the whole capture.
     */
    std::vector<std::complex<float>> reconstruct(
        const std::vector<std::complex<float>>& rx1,
        const std::vector<std::complex<float>>& rx2);

    /*!
     * Streaming reconstruction (real-time / GNU Radio).
     * Consumes n samples from each radio (at per_radio_rate) and returns the
     * combined samples (at target_rate) ready so far.  Alignment is estimated
     * adaptively from the overlap and stabilises after the first frames.
     * Maintains internal resampler / translation / overlap-save state across
     * calls; call reset() to start a new stream.
     */
    std::vector<std::complex<float>> process(
        const std::complex<float>* rx1, const std::complex<float>* rx2, size_t n);

    //! Clear all streaming state.
    void reset();

    //! Alignment recovered by the last reconstruct()/process() activity.
    const overlap_alignment& alignment() const { return _alignment; }

    //! Override the alignment estimate (e.g. hold a previous estimate).
    void set_alignment(const overlap_alignment& a) { _alignment = a; }

    const overlap_reconstructor_config& config() const { return _cfg; }

private:
    struct impl; //!< All designed filters and streaming DSP state.

    overlap_reconstructor_config _cfg;
    overlap_alignment _alignment;
    std::shared_ptr<impl> _impl;
};

/*!
 * Polyphase rational resampler (one-shot) for complex baseband.
 * Resamples x by the ratio interp/decim with a windowed-sinc anti-alias filter.
 * Shared by the reconstructor and by the test band-splitter.
 */
BONDED_API std::vector<std::complex<float>> rational_resample(
    const std::vector<std::complex<float>>& x, size_t interp, size_t decim);

} // namespace bonded
