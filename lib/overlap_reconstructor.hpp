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
 * Two radios each capture a slice of a target channel, with an intentional
 * spectral overlap around the channel centre fc instead of a hard split:
 *
 *   target: a channel of width target_bw centered at fc, sampled at target_rate
 *   radio1: tuned at fc + radio1_offset_hz (lower slice), at per_radio_rate
 *   radio2: tuned at fc + radio2_offset_hz (upper slice), at per_radio_rate
 *
 * Each radio is resampled to target_rate, translated to the target baseband
 * and band-shaped by complementary complex cos^2/sin^2 crossfade filters whose
 * transition is the overlap; the out-of-band corners are discarded. Radio 2 is
 * aligned to radio 1 (gain as a power ratio, phase from the cross-correlation)
 * over the overlap during a warmup, after which the alignment is frozen (the
 * radios share a 10 MHz reference, so it is static) and a fused path runs one
 * complex polyphase filter per radio. The result is one contiguous,
 * phase-preserving complex baseband suitable for decoding.
 */

struct overlap_reconstructor_config
{
    double target_rate      = 20e6;  //!< Output sample rate (Hz)
    double per_radio_rate   = 15e6;  //!< Per-radio input sample rate (Hz)
    double target_bw        = 20e6;  //!< Reconstructed channel bandwidth (Hz)
    double overlap_bw       = 5e6;   //!< Shared overlap bandwidth (Hz)
    double radio1_offset_hz = -5e6;  //!< radio1 center - fc (Hz), lower slice
    double radio2_offset_hz = +5e6;  //!< radio2 center - fc (Hz), upper slice
    //! Output samples to process before freezing the (static) alignment and
    //! switching to the fused steady-state path.
    size_t warmup_samples = 200000;
};

//! Alignment estimate of radio2 relative to radio1 (over the overlap band).
struct overlap_alignment
{
    double gain      = 1.0;  //!< |radio1/radio2| amplitude ratio
    double phase_rad = 0.0;  //!< constant phase offset (rad)
    bool valid       = false;
};

class BONDED_API overlap_reconstructor
{
public:
    explicit overlap_reconstructor(const overlap_reconstructor_config& cfg);

    /*!
     * One-shot reconstruction of a finite capture (offline / unit test).
     * rx1, rx2 are equal-duration complex baseband captures at per_radio_rate.
     * Returns the combined baseband at target_rate.  Starts from a fresh
     * state, like reset() followed by process().
     */
    std::vector<std::complex<float>> reconstruct(
        const std::vector<std::complex<float>>& rx1,
        const std::vector<std::complex<float>>& rx2);

    /*!
     * Streaming reconstruction (real-time / GNU Radio).
     * Consumes n samples from each radio (at per_radio_rate) and returns the
     * combined samples (at target_rate) ready so far.  Filter and alignment
     * state is kept across calls; call reset() to start a new stream.
     */
    std::vector<std::complex<float>> process(
        const std::complex<float>* rx1, const std::complex<float>* rx2, size_t n);

    /*!
     * Streaming reconstruction without copying the internal output buffer.
     * The returned reference remains valid until the next process/process_view/reset call.
     */
    const std::vector<std::complex<float>>& process_view(
        const std::complex<float>* rx1, const std::complex<float>* rx2, size_t n);

    /*!
     * Like process_view(), but the final combine writes up to `cap` samples
     * straight into `dst` (typically the GNU Radio output buffer). Returns
     * how many samples were written to `dst`; whatever did not fit is left in
     * pending() and must be drained before the next call.  Only the frozen
     * path writes directly; during the warmup everything goes to pending().
     */
    size_t process_into(const std::complex<float>* rx1,
        const std::complex<float>* rx2,
        size_t n,
        std::complex<float>* dst,
        size_t cap);

    //! Output samples from the last process_into() that did not fit in `dst`.
    const std::vector<std::complex<float>>& pending() const;

    //! Clear all streaming state.
    void reset();

    //! Alignment estimated so far (frozen once the warmup is over).
    const overlap_alignment& alignment() const { return _alignment; }

    const overlap_reconstructor_config& config() const { return _cfg; }

private:
    struct impl; //!< All designed filters and streaming DSP state.

    overlap_reconstructor_config _cfg;
    overlap_alignment _alignment;
    std::shared_ptr<impl> _impl;
};

} // namespace bonded
