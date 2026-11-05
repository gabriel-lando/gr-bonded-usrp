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
- **TX: hard split at the channel centre.** The seam falls on the unused DC
  subcarrier of 802.11. The receiver's per-subcarrier equalizer absorbs each
  half's gain, phase and delay, so no calibration is needed.

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
| `lib/band_splitter.{hpp,cpp}` | TX split: NCO + polyphase resampler per radio |
| `lib/bonded_resample.hpp` | `poly_resampler` (real prototype), `poly_resampler_c` (complex prototype, prepare/compute/commit), filter design, dot kernels |
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
- **Alignment is estimated during a warmup, then frozen** (the shared
  reference makes it static). The frozen path must continue the warmup's
  timeline exactly: the fused filters are fed the warmup input (prepare/commit
  without computing) and both paths share `rot1/rot2`. If the fused path
  restarts at the freeze, each radio's translation phase jumps by a different
  amount and the frozen alignment becomes wrong (a static seam notch).
  `test_overlap_frozen_path_continues_warmup_path` guards this.
- **Overlap sets the cost.** The resampling ratio depends only on O/T. 25 %
  (ratio 8/5) is the cheapest practical choice; smaller overlaps need longer
  filters and cost more CPU, not less.
- **TX seam width.** The split filter's transition is the only spectrum
  radiated by both radios, and it fades with their relative phase. Keep
  `taps_per_phase` high (64). Each radio's cutoff is |offset|, so the radios
  can sit at `fc ∓ 14/64·T`, away from their own band edge.

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

- **Warmup content.** The alignment freezes after `warmup_samples` output
  samples (~10 ms at 20 MHz). If the overlap holds only uncorrelated noise
  during that window, the frozen phase is arbitrary and the subcarriers around
  the RX seam are attenuated for the whole session. Low-order modulations
  survive this (the equalizer absorbs a static channel), denser ones may not.
  Not yet measured live; a fix would freeze only once the overlap shows
  coherent content.
- The RX pair's relative phase drifts after the freeze (about 9 °/s at
  20 MHz and 20 °/s at 40 MHz were measured, against ~0.3 °/s for a TX pair).
  This costs 64-QAM through the RX seam. Tracking it would mean not freezing,
  which conflicts with the real-time budget of the warmup path.
- `band_splitter` runs at roughly 33 MS/s input on one core. That is enough
  for bursty 802.11 traffic at 60 MHz, but not for continuous full-rate
  traffic.
