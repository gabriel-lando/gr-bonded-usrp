# AGENTS.md — gr-bonded-usrp (bonded wideband TX/RX, self-contained module)

> This is the **product** repo: a self-contained GNU Radio OOT module that bonds
> two low-cost USB SDRs (Ettus B210) into one wideband channel. The C++ DSP +
> acquisition engine is compiled **into this module** and links only the
> **public** UHD API — it no longer lives inside the UHD source tree. Test and
> capture scripts live in a **separate sibling repo** (see §8).

## 1. Goal

Bond two B210s via an **intentional frequency overlap** (NOT a hard split):

- **RX (thesis Part 1):** reconstruct one contiguous **20 MHz** wideband baseband
  from two overlapping ~15 MHz slices and decode real **802.11a/g** packets.
- **TX (Part 2):** transmit a wideband signal split across the two radios. Needs
  inter-radio **phase calibration** — two B210s on a shared 10 MHz are frequency-
  locked but **not** phase-coherent (see `[[bonding-tx-coherence]]`).

Each radio tunes to `fc ± T/4` and captures a wider slice sharing a configurable
overlap centred on `fc`; the overlap aligns the radios (gain/phase) and crossfades
them into one channel; out-of-band corners are discarded.

## 2. Hard rule (do not violate)

**All signal-processing/DSP logic lives in C++** (the engine below). Python / GRC
is ONLY a thin flowgraph/test interface. See `[[bonding-dsp-in-cpp]]`.

## 3. Architecture

```
Radio A (B210 @ fc-T/4) ┐ shared 10 MHz + PPS    Radio B (B210 @ fc+T/4) ┘
        │ bonded::bonded_receiver      — synchronized acquisition (lock 10 MHz,
        │                                detect PPS, latch common time=0; <1 µs)
        ▼ bonded::overlap_reconstructor — complex/IQ reconstruction DSP → fc baseband
        ▼ gr-ieee802-11 receiver        — Schmidl-Cox → FFT → equalize → Viterbi → MAC
TX path: bonded::bonded_transmitter (bonded_sink block) ← split wideband → 2 radios
```

GR blocks (this module): `bonded_source`, `overlap_reconstructor`, `bonded_sink`.
The C++ engine and the block wrappers are all under `lib/`.

## 4. Layout (everything in `lib/`)

| File | Role |
|---|---|
| `overlap_reconstructor.{hpp,cpp}` | **THE** reconstruction DSP (time-domain VOLK) |
| `bonded_resample.hpp` | shared polyphase resampler + VOLK dot helpers (header-only) |
| `bonded_receiver.{hpp,cpp}` | synchronized multi-device acquisition (RX) |
| `bonded_transmitter.{hpp,cpp}` | synchronized multi-device transmit (TX) |
| `bonded_alignment.hpp` | burst-alignment / ring-buffer helpers (header-only, `bonded::detail`) |
| `bonded_api.h` | `BONDED_API` export macro (replaced UHD's `UHD_API`) |
| `*_impl.{cc,h}` | GR block wrappers (thin: forecast + buffer plumbing only) |
| `tests/*` | hardware-free Boost.Test ctests (see §6) |

**Namespace:** the engine is `namespace bonded` (was `uhd::usrp::bonded` when it
lived inside UHD). It depends on UHD only through the public API
(`uhd/usrp/multi_usrp.hpp`, `uhd/types/*`, `uhd/utils/thread.hpp`).

## 5. Build / install (fast — no libuhd rebuild)

Requires an **installed UHD** (official is fine) with dev headers + `UHDConfig.cmake`,
GNU Radio 3.10, Boost, and optionally VOLK / FFTW3f / OpenMP (auto-detected,
accelerate the reconstructor; built-in fallbacks otherwise).

```bash
cd gr-bonded-usrp && mkdir build && cd build
cmake .. && make -j16
ctest --output-on-failure                 # 3 hardware-free C++ tests
sudo make install && sudo ldconfig         # installs lib + python + grc
```

Editing the engine = rebuild **only this module** (seconds). No more
"install libuhd first", no root-owned UHD build dir. Rebuild + reinstall whenever
a C++ header changes (class layout) to avoid `free(): invalid pointer` from a
stale Python binding.

## 6. Validation

- `ctest` → `overlap_reconstructor_test` (alignment recovery, ±band-flatness, seam,
  streaming), `bonded_receiver_test` (burst alignment / ring extraction),
  `bonded_transmitter_test` (config + freq-plan geometry).
- `python3 ../samples/test_overlap_block.py` (or sibling repo) — GR flowgraph smoke
  test, no hardware.
- **Headline (still pending on real B210s):** offline record→replay vs the
  single-radio 20 MHz baseline — compare decoded frame counts.

## 7. DSP invariants — DO NOT regress

- **Complex (IQ) reconstruction** — phase preserved (OFDM). Power-spectrum
  stitching is unusable for decode.
- **Gain = power ratio** `sqrt(Σ|A|²/Σ|B|²)`, NOT `|Σ A·conj(B)|/Σ|B|²` (the latter
  collapses to 0 on bursty WiFi → "0 packets decoded"; this was a real bug).
- **Alignment is static** (shared 10 MHz) — estimate during a warmup, then **freeze**.
- **Configurable overlap:** smaller overlap lowers USB bandwidth but RAISES
  reconstruction CPU (sharper transition → longer FIRs). The CPU knob is `T`, not `O`.
- **Real-time:** ~0.3× at 20 MHz on a 16-core CPU; marginal at 5 MHz. The 20 MHz
  scientific claim is validated **offline**. Next RT avenue: a fused single
  complex-polyphase filter per radio.

## 8. Sibling test repo (separate, side by side)

Test/capture/flowgraph scripts and the test-only `band_splitter` block are **not**
in this module. They live in a sibling repo (was `host/samples/`):
`wifi_bonded_rx.py`, `wifi_bonded_tx.py`, `wifi_trx.py`, capture/analyze scripts,
`OFFLINE_20MHZ_DECODE.md`, and `band_splitter/` (a dev block that splits one
capture into two simulated slices — see its README; depends on the installed
module's `rational_resample()`). **Captures (`*.fc32` ~248 MB, `*.pcap`) are
gitignored — never commit them.**

## 9. Hardware

4× B210. **TX** `30DBC3D`; **bonded RX pair** `30B56D6` (A, lower) + `30DBC3C`
(B, upper) — MUST share an external 10 MHz ref + PPS; **single-radio baseline**
`30EDB63`. Default `fc = 2.412 GHz`.

## 10. Dropped on extraction — do not resurrect

- Hard-split pipeline and the power-spectrum stitcher (display-only, discards phase).
- Legacy `bonded=true` device-args path (`bonded_usrp.*` + a `multi_usrp.cpp`
  patch) — it was the only edit to UHD core; the real RX/TX engine never used it.
- `spectrum_stitcher`, `frequency_plan` — unused by any block/engine.

Memory: `[[bonding-dsp-in-cpp]]`, `[[bonding-overlap-design]]`, `[[bonding-tx-coherence]]`.
