# gr-bonded-usrp

GNU Radio out-of-tree module that bonds two Ettus USRP B210 radios into one
wider channel, for receive and for transmit. Each radio covers one half of the
channel. The module keeps both radios on a common time base and merges or
splits the signal in C++.

The two radios must share an external **10 MHz reference and PPS**. That makes
them frequency-locked and lets them start on the same device time. Each B210
still has its own local oscillator, so the two radios are **not
phase-coherent**: their relative phase is arbitrary after every tune. The RX
and TX paths handle this in different ways.

## How it works

### Receive: overlap and complex crossfade

```
           radio A (fc - offset)  -.
 bonded_source                      >-- overlap_reconstructor --> channel at fc
           radio B (fc + offset)  -'
```

The two radios capture slices of the channel that **overlap** around the
centre `fc`. `overlap_reconstructor` resamples each slice to the output rate
and shifts it to its place in the channel. It then joins the two slices with
complementary cos²/sin² crossfade filters, whose transition band is the
overlap, and discards whatever falls outside the channel.

During a short warmup the block estimates radio B's gain and phase relative to
radio A from the overlap. The gain is a power ratio and the phase comes from
the cross-correlation. After the warmup the block switches to a fused path:
one complex polyphase filter per radio, with the output split across threads.
The shared reference keeps the alignment constant, but the fused path keeps
measuring the phase whenever the overlap holds a signal both radios see, so a
warmup that saw only noise is corrected on the first frames. Because each
radio has its own antenna, the exact phase differs slightly from emitter to
emitter; the tracking settles on the strongest signal. The output is
phase-preserving complex baseband that a standard receiver can decode, for
example an 802.11 OFDM receiver.

`bonded_source` streams both radios and aligns their samples by device
timestamp. If one radio drops samples (overflow), the block is dropped on
**both** radios, so they never slip against each other.

Frequency plan without wasted spectrum, for a channel of width `T` with
overlap `O`:

| | value | example, T = 20 MHz, O = 5 MHz |
|---|---|---|
| per-radio rate | `(T + O) / 2` | 12.5 MS/s |
| radio centres | `fc ∓ (T − O) / 4` | fc ∓ 3.75 MHz |
| output rate | `T` | 20 MS/s |

The CPU cost depends on `O` *as a fraction of* `T`, because that sets the
resampling ratio. 25 % (as above) gives the ratio 8/5 and is the cheapest
practical choice. Smaller overlaps need longer filters.

### Transmit: hard split at the centre

```
 wideband baseband --> ofdm_splitter --> radio A (fc - offset) -.
                     (or band_splitter) \--> radio B (fc + offset) --> bonded_sink
```

Each radio transmits one half of the channel: radio A the lower half, radio B
the upper half. The radios are not phase-coherent, so the receiver sees each
half with its own gain, phase and delay. A standard OFDM receiver's
per-subcarrier channel estimate absorbs that, as long as every subcarrier is
radiated by one radio only and the two halves stay orthogonal. No phase
calibration is needed.

`ofdm_splitter` (for OFDM, recommended) splits every OFDM symbol in the
subcarrier domain. It takes the FFT of the symbol, gives each radio the
subcarriers of its half, and synthesizes them directly at the radio's rate
with their own cyclic prefix. Each radio therefore sends a complete CP-OFDM
signal, and the receiver's FFT separates the two exactly, whatever their
relative phase. Consecutive symbols overlap by a short raised-cosine window,
which keeps each radio's spectrum compact. The block needs to know where the
bursts start: a length tag (`packet_len` from gr-ieee802-11 or gr-digital OFDM
transmitters) on the first sample of each burst.

`band_splitter` (generic) cuts the wideband signal at `fc` with a sharp
filter: each half is shifted to DC and resampled through a lowpass of
half-width |offset|. It works for any signal, but its filter is much longer
than an OFDM cyclic prefix. When the two radios arrive with a relative phase
near 180°, the subcarriers next to the seam get inter-symbol interference.
BPSK to 16-QAM do not mind; 64-QAM fails unless the phase is brought near 0°
(`set_phase_offset` on the sink, from a receiver's measurement).

`bonded_sink` streams both radios **continuously** for the whole session. A
B210 needs several hundred microseconds to ramp its TX power after idling,
which is longer than an 802.11 preamble. So the gaps between input bursts are
filled with zeros, and every burst goes out at full power. A scheduler thread
hands the same sample sequence to one sender thread per radio, both starting
at the same device time. If either radio underruns, both restart on a new
common time.

Recommended TX plans for 802.11 with a channel of width `T`:

| splitter | per-radio rate | radio offsets | example, T = 20 MHz |
|---|---|---|---|
| `ofdm_splitter` | `5/8 · T` | `∓ 12/64 · T` (12 subcarriers) | 12.5 MS/s at fc ∓ 3.75 MHz |
| `band_splitter` | `T/2` | `∓ 14/64 · T` | 10 MS/s at fc ∓ 4.375 MHz |

The `ofdm_splitter` plan is the RX plan: each radio covers its 26 subcarriers
with room for the window to spill a few subcarriers past the seam. `T/2` at
`∓ 13/64 · T` also works with it, with less margin at the radios' band
edges. The `band_splitter` plan keeps the seam away from the edge of each
radio's band, where its DUC rolls off.

Drive level: the radios' DACs clip at full scale (|x| = 1) and OFDM has a high
peak-to-average ratio. Keep each radio's RMS around 0.3 or below and use the
TX gain for power; clipping shows up first as 64-QAM losses.

## Blocks

All blocks live in `gnuradio.bonded_usrp`. GRC definitions are installed under
the *[bonded_usrp]* category.

### `bonded_source` (RX)

```python
bonded_usrp.bonded_source(serials_csv, rate, freq, gain,
    clock_source="external", time_source="external", strict_lock=False,
    lock_timeout=5.0, fetch_timeout=1.0, subdev="", freq_plan_csv="",
    stream_args="")
```

| parameter | meaning |
|---|---|
| `serials_csv` | B210 serials, radio A first (`"SERIAL_A,SERIAL_B"`) |
| `rate` | per-radio sample rate (Hz) |
| `freq` | centre frequency for every radio, unless `freq_plan_csv` is set |
| `gain` | RX gain (dB) |
| `clock_source`, `time_source` | `external` for the shared 10 MHz and PPS |
| `strict_lock` | throw instead of warn if the reference does not lock within `lock_timeout` s |
| `fetch_timeout` | seconds `work()` waits for aligned data |
| `freq_plan_csv` | per-radio centre frequencies (Hz), e.g. `"2.40825e9,2.41575e9"` |
| `stream_args` | extra UHD RX stream args |

Outputs: two ports per radio (port `2*d + ch` = radio `d`, channel `ch`). Use
ports 0 and 2 for `overlap_reconstructor`.

### `overlap_reconstructor` (RX)

```python
bonded_usrp.overlap_reconstructor(rate_in, rate_out, target_bw, overlap_bw,
    radio1_offset, radio2_offset)
```

Inputs: radio A, radio B (at `rate_in`). Output: the channel at `rate_out`.
`radio1_offset` and `radio2_offset` are each radio's centre minus `fc` and must
match the tuning. For the plan above: `(12.5e6, 20e6, 20e6, 5e6, -3.75e6,
3.75e6)`.

### `band_splitter` (TX)

```python
bonded_usrp.band_splitter(rate_in, rate_out, radio1_offset, radio2_offset,
    taps_per_phase=64)
```

Input: the wideband signal at `rate_in`. Outputs: radio A, radio B (at
`rate_out`). Requires `2·|offset| ≤ rate_out`.

### `ofdm_splitter` (TX, OFDM)

```python
bonded_usrp.ofdm_splitter(rate_in, rate_out, radio1_offset, radio2_offset,
    fft_len=64, cp_len=16, window_len=8, len_tag_key="packet_len")
```

Input: CP-OFDM bursts at `rate_in`, each starting at a `len_tag_key` length
tag, made of back-to-back symbols of `fft_len + cp_len` samples (a trailing
partial symbol is ignored; gr-ieee802-11 frames end with one extra sample).
Outputs: radio A (the subcarriers below the midpoint between the radio
centres), radio B (those above; a subcarrier exactly on the midpoint is
shared). Input outside bursts becomes silence on the same timeline.

| parameter | meaning |
|---|---|
| `rate_in`, `rate_out` | wideband and per-radio rates; `fft_len` and `cp_len` must be whole numbers of samples at `rate_out` |
| `radio1_offset`, `radio2_offset` | radio centres minus `fc`, whole subcarriers (`rate_in / fft_len`), radio 1 below radio 2 |
| `window_len` | raised-cosine overlap between symbols, in input samples, at most `cp_len`. The receiver's FFT window must start at least this far into the cyclic prefix; 8 of 16 leaves the other half for the channel and timing |

For the plan above: `(20e6, 12.5e6, -3.75e6, 3.75e6)`.

### `bonded_sink` (TX)

```python
bonded_usrp.bonded_sink(serials_csv, rate, freq, gain,
    clock_source="external", time_source="external", strict_lock=False,
    lock_timeout=5.0, subdev="", freq_plan_csv="", stream_args="",
    delay_trim_a=0, delay_trim_b=0, gain_plan_csv="", freq_offset_csv="",
    phase_offset_csv="")
```

| parameter | meaning |
|---|---|
| `serials_csv` | exactly two B210 serials, radio A (lower half) first |
| `rate` | per-radio sample rate (Hz) |
| `gain` | normalized TX gain, 0 to 1 (`gain_plan_csv` sets it per radio) |
| `freq_plan_csv` | per-radio centre frequencies (Hz); empty tunes both to `freq` |
| `delay_trim_a/b` | per-radio delay in samples |
| `freq_offset_csv`, `phase_offset_csv` | per-radio digital frequency (Hz) and phase (rad) correction |

Methods: `set_phase_offset(device, rad)` and `phase_offset(device)` change and
read a radio's phase while running. `frames_sent()` counts the input bursts
handed to the radios.

The device setup (opening the radios, locking, PPS alignment) runs in the
constructor and takes a few seconds. This way no upstream block can queue a
backlog before the radios start streaming. Streaming starts 0.5 s after the
flowgraph starts. On stop, whatever is still queued is transmitted before
the bursts end, so a finite source is sent completely.

## Minimal flowgraphs

```python
from gnuradio import blocks, gr
from gnuradio import bonded_usrp

fc, T, O = 2.412e9, 20e6, 5e6
per_radio, off = (T + O) / 2, (T - O) / 4

# Receive one 20 MHz channel with two radios.
rx = gr.top_block()
src = bonded_usrp.bonded_source("SERIAL_A,SERIAL_B", per_radio, fc, 40.0,
                                freq_plan_csv=f"{fc - off},{fc + off}")
rec = bonded_usrp.overlap_reconstructor(per_radio, T, T, O, -off, off)
rx.connect((src, 0), (rec, 0))
rx.connect((src, 2), (rec, 1))
rx.connect(rec, blocks.file_sink(gr.sizeof_gr_complex, "channel.fc32"))
for port in (1, 3):
    rx.connect((src, port), blocks.null_sink(gr.sizeof_gr_complex))

# Transmit a 20 MHz baseband signal (here: from a file) with two radios.
tx = gr.top_block()
off_tx = 14 / 64 * T
wide = blocks.file_source(gr.sizeof_gr_complex, "baseband_20msps.fc32", False)
spl = bonded_usrp.band_splitter(T, T / 2, -off_tx, off_tx)
snk = bonded_usrp.bonded_sink("SERIAL_C,SERIAL_D", T / 2, fc, 0.7,
                              freq_plan_csv=f"{fc - off_tx},{fc + off_tx}")
tx.connect(wide, spl)
tx.connect((spl, 0), (snk, 0))
tx.connect((spl, 1), (snk, 1))
tx.run()  # returns once the whole file has been transmitted
```

For OFDM bursts carrying length tags (for example straight from a
gr-ieee802-11 transmitter), use the subcarrier split instead:

```python
off_tx = 12 / 64 * T
spl = bonded_usrp.ofdm_splitter(T, 5 / 8 * T, -off_tx, off_tx)
snk = bonded_usrp.bonded_sink("SERIAL_C,SERIAL_D", 5 / 8 * T, fc, 0.7,
                              freq_plan_csv=f"{fc - off_tx},{fc + off_tx}")
```

## Requirements

- Two (per direction) Ettus B210, all sharing one external 10 MHz reference and
  PPS (for example from an OctoClock). Ideally each radio has its own USB 3
  port.
- UHD 4.x with development headers (`UHDConfig.cmake`). Tested with 4.10.
- GNU Radio 3.10 with development headers and pybind11. Tested with 3.10.12.
- CMake ≥ 3.16 and a C++17 compiler.
- GNU Radio's FFT component (`gnuradio-fft`, part of a normal GNU Radio
  install) for `ofdm_splitter`.
- Optional: VOLK (normally installed with GNU Radio) for the NCO rotators and
  some dot products; Boost.Test for the C++ unit tests; numpy for the Python
  tests.

The engine uses only the public UHD API. UHD itself is not modified.

## Build and install

```bash
git clone https://github.com/gabriel-lando/gr-bonded-usrp.git
cd gr-bonded-usrp
mkdir build && cd build
cmake ..                      # -DCMAKE_INSTALL_PREFIX=... to install elsewhere
make -j$(nproc)
ctest --output-on-failure     # hardware-free tests
sudo make install
sudo ldconfig
```

Check the installation:

```bash
python3 -c "from gnuradio import bonded_usrp; print(bonded_usrp.bonded_source)"
```

The reconstructor is compiled with `-march=native`, so build the module on the
machine that will run it. After changing a C++ header, reinstall before running
Python flowgraphs again: a stale binding usually crashes with
`free(): invalid pointer`.

## Tests

`ctest` runs, without hardware:

| test | covers |
|---|---|
| `overlap_reconstructor_test` | gain and phase recovery, band flatness including the seam, both frequency plans, warmup and frozen paths (and that the frozen path continues the warmup timeline), phase tracking after a noise-only warmup and its immunity to noise, no phase drift over a long run, streaming, thread count invariance, direct-output API, config checks |
| `bonded_resample_test` | filter design, dot-product kernels, streaming polyphase resamplers, sliced computation, oscillator phase over long runs |
| `bonded_receiver_test` | timestamp alignment of the receive rings: stale data, holes after an overflow, multi-channel, two-device alignment, bookkeeping |
| `bonded_transmitter_test` | configuration checks, phase-offset API, frequency/phase correction across chunks |
| `band_splitter_test` | flatness, half assignment, seam width, pulled-in plan, streaming vs one-shot, output length, no phase drift, config checks |
| `ofdm_splitter_test` | every subcarrier recovered at any relative phase of the radios (both plans), subcarrier assignment, radio-rate timeline and burst release, chunked vs one-call streaming, config checks |
| `bonded_csv_test` | parsing of the comma-separated block parameters |
| `qa_band_splitter`, `qa_ofdm_splitter`, `qa_overlap_reconstructor` | the GNU Radio blocks in flowgraphs (Python) |

Device setup, streaming and the over-the-air behaviour require radios and are
not covered by `ctest`.

## Runtime tuning

These environment variables are read when a block starts. None needs a
rebuild.

| variable | default | effect |
|---|---|---|
| `BONDED_THREADS` | 4 | reconstruction threads (frozen path) |
| `BONDED_PARALLEL=0` | on | single-threaded reconstruction |
| `BONDED_CHUNK` | 8192 | reconstructor output multiple; its output buffer is 8× this |
| `BONDED_SRC_BUFFER` | 262144 | `bonded_source` output buffer (items); lowering it costs throughput |
| `BONDED_DIRECT_OUT=0` | on | reconstructor copies through an internal buffer instead of writing into the output buffer (for comparisons) |
| `BONDED_TRACK=0` | on | reconstructor keeps the warmup's phase instead of tracking it on the fused path |
| `BONDED_PROFILE=1` | off | reconstructor profile once per second on stderr |
| `BONDED_RX_PROFILE=1` | off | acquisition profile once per second on stderr |

The two profile lines:

```text
[bonded_profile]    recon calls=N parallel_calls=N avg_ms=X in_msps=X out_msps=X stage_pct={...}
[bonded_rx_profile] calls=N wait=X% ringlock=X% | dev0 recv=X served=X dropped=X ring=Xms | dev1 ...
```

To diagnose a throughput problem:

- `recv` below the per-radio rate means the radios or USB are not delivering.
- `recv` at the rate but `served` below it means the flowgraph is behind. A
  growing `ring` is the early warning; drops start once it reaches 4 s.
- A high `wait` with `ring` near zero is healthy.
- For the reconstructor, duty cycle = `calls × avg_ms / 1000`. Near 100 % means
  it is compute-bound. Well below 100 % while `in_msps` is under the per-radio
  rate means it is starved by the source.

The reconstructor also reports its alignment on stderr, once at the end of
the warmup and whenever tracking moves the phase by more than 2° (at most
once per second):

```text
[overlap_reconstructor] alignment frozen: gain G, phase P deg, coherence C
[overlap_reconstructor] alignment phase P1 -> P2 deg (tracked; this call P3)
```

A coherence near 1 means the warmup saw a signal common to both radios; near
0 means it saw only noise, and the first frames through the overlap will
correct the phase.

Set the CPU frequency governor to `performance` for high bandwidths.

## Measured capacity

On a Ryzen 7 5700G (8 cores), with the frequency plans above:

| | bandwidth |
|---|---|
| RX, reconstruction only | 80 MHz real time (64 MHz with margin) |
| RX with a GNU Radio 802.11 receiver attached | 40 MHz decoding at parity with a single-radio receiver |
| TX, bonded, 802.11 frames | 20, 40 and 60 MHz without underruns |
| 64-QAM 3/4 through a bonded TX and a bonded RX (`ofdm_splitter`) | 97–100 % of the frames at 20 and 40 MHz; ~85 % at 60 MHz, where the same TX power is spread over three times the bandwidth |

## Source layout

| path | content |
|---|---|
| `lib/overlap_reconstructor.{hpp,cpp}` | RX reconstruction engine |
| `lib/ofdm_splitter.{hpp,cpp}` | TX split engine for OFDM (subcarrier domain) |
| `lib/band_splitter.{hpp,cpp}` | TX split engine (filter) |
| `lib/bonded_resample.hpp` | polyphase resamplers, filter design, dot kernels |
| `lib/bonded_receiver.{hpp,cpp}` | synchronized multi-radio receive (rings, alignment) |
| `lib/bonded_alignment.hpp` | timestamp alignment of the receive rings |
| `lib/bonded_transmitter.{hpp,cpp}` | synchronized multi-radio transmit (scheduler, per-radio senders) |
| `lib/bonded_csv.hpp` | block parameter parsing |
| `lib/*_impl.{cc,h}` | GNU Radio block wrappers |
| `include/gnuradio/bonded_usrp/` | public block headers |
| `python/bonded_usrp/` | Python bindings and flowgraph tests |
| `grc/` | GRC block definitions |
| `tests/` | C++ unit tests |

## License

GPL-3.0-or-later.
