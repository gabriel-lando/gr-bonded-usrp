# AGENTS.md: gr-bonded-usrp

Guidance for coding agents working on this module. The user-facing
description (blocks, build, tuning) is in `README.md`; read it first.

## 1. What this is

A GNU Radio 3.10 OOT module that bonds two Ettus B210 radios into one wider
channel, for RX and for TX. The radios share an external 10 MHz reference and
PPS: they are frequency-locked and time-aligned, but **not phase-coherent**
(each AD9361 has its own LO, so the relative phase is random after every
tune). The engine is C++ and uses only the public UHD API.

RX and TX use **deliberately different** strategies:

- **RX: intentional overlap plus complex crossfade.** The overlap is used to
  estimate radio B's gain and phase against radio A, and is where the two
  slices are crossfaded.
- **TX: hard split at the channel centre.** Each radio radiates one half,
  and the receiver's per-subcarrier equalizer absorbs each half's gain, phase
  and delay, so no calibration is needed. For OFDM the split is done per
  symbol in the subcarrier domain (`ofdm_splitter`), which keeps the halves
  orthogonal at any relative phase; `band_splitter` is the generic filter
  split.

## 2. Hard rules

- **All DSP stays in C++.** Python/GRC is only flowgraph wiring. An early
  prototype that combined two independent Python `usrp_source` streams desynced
  on the first overflow and produced undecodable output. Real time also needs
  tight SIMD loops.
- **Never patch UHD.** Configuring bonding from inside UHD was tried first and
  abandoned. Everything here goes through `multi_usrp`.
- `ctest` must stay green. Add a test with every behaviour change; the DSP and
  alignment logic are deliberately kept testable without hardware.

## 3. Layout

| file | role |
|---|---|
| `lib/overlap_reconstructor.{hpp,cpp}` | RX reconstruction: warmup path (resample, rotate, align, crossfade) and frozen fused path (one complex polyphase filter per radio, sliced across threads) |
| `lib/ofdm_splitter.{hpp,cpp}` | TX split for CP-OFDM bursts: per-symbol FFT, each radio's subcarriers synthesized at its rate (IFFT, CP, raised-cosine overlap), translation phase from the radio-rate timeline |
| `lib/band_splitter.{hpp,cpp}` | TX split for any signal: NCO + polyphase resampler per radio |
| `lib/bonded_resample.hpp` | `poly_resampler` (real prototype), `poly_resampler_c` (complex prototype, prepare/compute/commit), filter design, dot kernels, `nco` |
| `lib/bonded_receiver.{hpp,cpp}` | device setup and PPS alignment, one recv thread and timestamped ring per device, aligned extraction |
| `lib/bonded_alignment.hpp` | `extract_aligned_samples()`: places chunks on the time grid by timestamp |
| `lib/bonded_transmitter.{hpp,cpp}` | device setup, scheduler thread + one sender thread per device, zero-fill, underrun re-alignment, phase correction |
| `lib/bonded_csv.hpp` | CSV parameter parsing (throws on malformed numbers) |
| `lib/*_impl.{cc,h}` | GR wrappers: forecast, buffering, env knobs only |
| `tests/` | Boost.Test unit tests; `python/bonded_usrp/qa_*.py` flowgraph tests |

Engine code is in `namespace bonded`; header-only helpers are in
`bonded::detail`. The engine headers under `lib/` are not installed. Only the
block headers in `include/gnuradio/bonded_usrp/` are public.

## 4. Build and test

```bash
mkdir -p build && cd build && cmake .. && make -j$(nproc)
ctest --output-on-failure
sudo make install && sudo ldconfig
```

`overlap_reconstructor.cpp` is built with `-O3 -ffast-math -funroll-loops
-march=native`. Bit-exact comparisons between differently vectorized loops
can therefore differ by an ULP; compare with a relative tolerance. Reinstall
after changing any header used by the bindings, or Python crashes with
`free(): invalid pointer`.

## 5. DSP invariants (do not regress)

- **Complex (IQ) reconstruction.** Never stitch magnitude spectra: that
  discards phase and cannot be decoded.
- **Gain = power ratio** `sqrt(ΣA²/ΣB²)` over the overlap, and phase =
  `arg(Σ A·conj(B))`. Never use `|Σ A·conj(B)|/Σ|B|²` for the gain: it
  collapses toward 0 on bursty traffic, which once decoded 0 packets. The
  power ratio assumes similar content in both radios. Each radio's resampler
  attenuates the far side of the overlap, which cancels for symmetric content
  (802.11) but biases the gain for lopsided test signals; tests that measure
  flatness use symmetric overlap content.
- **Alignment is estimated during a warmup, then the fused path takes over**
  (the shared reference makes it static). The fused path must continue the
  warmup's timeline exactly: the fused filters are fed the warmup input
  (prepare/commit without computing) and both paths share `nco1/nco2`. If the
  fused path restarts at the freeze, each radio's translation phase jumps by
  a different amount and the alignment becomes wrong (a static seam notch).
  `test_overlap_frozen_path_continues_warmup_path` guards this.
- **The fused path tracks the phase.** `arg(Σ a·b*)` of the two fused
  outputs gets contributions from the overlap only (their spectra share
  nothing else), so it measures the alignment phase for ~2 % of the cost.
  Only calls whose normalized correlation `ρ·n` exceeds 3 count (noise
  averages ~0.03): a warmup on a noise-only overlap freezes an arbitrary
  phase, and noise must not walk a good one. The gain stays the warmup's.
- **The right RX phase depends on the emitter.** Each RX radio has its own
  antenna, so the relative channel of the pair changes with where the signal
  comes from (live: a weak continuous source 40° away from our frames, which
  is what the warmup froze on). Tracking therefore sums `Σ a·b*` weighted by
  its own magnitude (energy squared) and never forgets: it locks onto the
  strongest coherent signal, the one dense constellations need. A short
  memory (tried: forgetting 0.8 per call) hops between emitters inside frames
  and cost 64-QAM 20 points. A phase error Δ only costs `cos(Δ/2)` at the
  centre of the overlap, so tens of degrees are harmless; near 180° is not.
- **Translation phases are exact.** Every frequency shift uses `nco`
  (phase in double-precision cycles, float increment only within a call). A
  float phasor advanced by a float increment runs at the increment's rounded
  frequency: that was the "RX pair drifts 9–20 °/s after the freeze" (two
  radios shifted in opposite directions, ~6 °/s each at 20 MHz) and it cost
  64-QAM through the seam. The `*_no_phase_drift` tests guard it.
- **Overlap sets the cost.** The resampling ratio depends only on O/T. 25 %
  (ratio 8/5) is the cheapest practical choice; smaller overlaps need longer
  filters and cost more CPU, not less.
- **TX split for OFDM is per symbol, not a filter.** The receiver sees the
  two halves summed with the radios' relative phase φ. A split filter longer
  than the cyclic prefix only sums back to an impulse at φ = 0; at 180° the
  subcarriers next to DC drop to about −8 dB EVM (offline model, confirmed
  live: 64-QAM 22 % at −114°, ~97 % near 0°). `ofdm_splitter` synthesizes
  each radio's own subcarriers as a complete CP-OFDM signal, so the FFT at
  the receiver separates the radios exactly at any φ (≤ −37 dB EVM in
  `test_ofdm_split_any_relative_phase`). Its invariants:
  - the translation phase of a symbol is computed from that symbol's start on
    the radio-rate output (`exp(-j2π·bin·q/NR)`, integer arithmetic), so
    consecutive symbols of a burst stay coherent;
  - symbols overlap by a raised-cosine window inside the cyclic prefix
    (`window_len` ≤ `cp_len`): without it each radio's spectrum ends as
    sharply as the symbol boundary and the radios' own interpolation filters
    (DUC) turn that into ISI;
  - a burst is released as soon as its last whole symbol is consumed (same
    no-holdback rule as `band_splitter`);
  - the burst start comes from the length tag: the 802.11 flowgraph must not
    put padding between the PHY and the splitter (it would shift the frame
    off its tag).
- **`band_splitter` seam width** (generic split). The split filter's
  transition is the only spectrum radiated by both radios, and it fades with
  their relative phase. Keep `taps_per_phase` high (64). Each radio's cutoff
  is |offset|, so the radios can sit at `fc ∓ 14/64·T`, away from their own
  band edge.

## 6. Streaming invariants (each of these was a real bug)

- **`band_splitter` must not hold back output.** `forecast()` returns 0
  required input while leftover output is pending, and `general_work()`
  consumes only what fits in the output. Otherwise 802.11's bursty input
  leaves each frame's tail parked until the next frame, hundreds of ms later.
- **`overlap_reconstructor` consumes only what fits in the output room.**
  Consuming everything available makes it and the PHY downstream run in
  series.
- **`bonded_source` needs a large `set_min_output_buffer`** (262144). With
  GNU Radio's default, scheduler handoffs cap throughput while no thread is
  saturated. *Nothing at 100 % CPU while throughput is pinned* means a
  buffer/handoff problem, not compute.
- **RX rings: holes stay holes.** After an overflow, extraction must place each
  chunk by its own timestamp and mark the block unsuccessful; `bonded_source`
  then drops it on every device. Concatenating across the hole shifts one
  radio for the rest of the session.
- **RX ring bookkeeping is O(1)** (`extract_result::removed`). Rescanning the
  deque under its mutex at tens of MS/s causes a feedback spiral.
- **A full ring must not stall the source.** Once the flowgraph falls 4 s
  behind, every push drops the oldest chunk of that ring. Choosing the start
  time and extracting therefore happen under *all* ring locks (otherwise a
  front moves in between and every fetch starts on a hole), and
  `bonded_source::work()` fetches again right away after a degraded block
  instead of returning 0: a source that produced nothing is retried by the
  scheduler only about twice a second, which kept the rings full and every
  fetch degraded for the rest of the session (live at 60 MHz: decode 0).
- **TX streams continuously** with zero-fill. A B210 needs ~300–600 µs to ramp
  TX power after idle, longer than an 802.11 preamble.
- **TX: one sender thread per device**, all fed the identical block sequence by
  a scheduler, plus a common timed start. A single thread sending to both in
  turn underruns past ~2×20 MS/s. After an underrun a B210 resumes immediately,
  so its timeline slips. Any underrun therefore restarts every device on a new
  common time.
- **TX device setup runs in the `bonded_sink` constructor**, not in `start()`.
  GR does not order `start()` across blocks, so upstream message sources would
  queue seconds of frames that then go out as one burst.
- **TX linger.** After data the scheduler waits up to 2 ms for the rest of a
  frame before zero-filling, so one frame split across several `work()` calls
  is not cut.
- **Zero-fill moves the splitter's timeline against the radios'.** The
  splitters number their output samples, but the radios also send the
  zero-fill in between, so each frame meets the radio LOs at a different
  phase: the relative phase of the two halves changes from frame to frame
  (seen live). Harmless for `ofdm_splitter` (each frame is self-contained);
  for `band_splitter` it means the seam quality varies per frame.

## 7. Performance notes (measured; do not retry the dead ends)

- The fused per-radio filter is ~97 % of reconstruction time; optimize that,
  not the rest. Its hot loop runs on de-interleaved re/im data
  (`dot_cc_split`), about 1.9× faster than interleaved complex.
- Parallelism that works: slicing each call's output range across persistent
  workers with atomic counters and spin → yield → sleep backoff. Condition
  variables cost ~200 µs per wake against ~0.2 ms slices. OpenMP fork-join
  inside the block was measured and does not help live.
- Things measured and found *not* to matter: the output memcpy
  (`BONDED_DIRECT_OUT`), CPU pinning/SMT placement, and the ring mutex hold
  time.
- Two B210s on one USB 3 controller deliver 45–50 MS/s each without overflow.
  At high bandwidths the limit is the host DSP and the flowgraph, not the
  radios.
- Check the CPU governor (`performance`) before any benchmark.

## 8. Known open points

- **60 MHz decode on an 8-core host with stock gr-ieee802-11 is marginal.**
  Reconstruction alone has margin, but the stock PHY plus the two receive
  threads (each near a full core at 37.5 MS/s) and the reconstructor's
  spinning workers leave little headroom. When ambient 802.11 traffic in the
  60 MHz span is heavy, the receive threads fall behind, every block has a
  hole and the source drops them for the rest of the session (most sessions
  are fine; some decode nothing). A faster PHY front-end (a fused
  Schmidl-Cox detector in gr-ieee802-11, outside this module) is what made
  80 MHz decode real time on the same CPU class.
- Each TX radio's LO leakage sits on the subcarrier at its centre (±12 in
  the `ofdm_splitter` plan): about 3 dB worse EVM there than elsewhere.
- `band_splitter` runs at roughly 33 MS/s input on one core. That is enough
  for bursty 802.11 traffic at 60 MHz, but not for continuous full-rate
  traffic.
