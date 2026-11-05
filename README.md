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
the cross-correlation. After the warmup the estimate is frozen (the shared
reference keeps it constant) and the block switches to a fused path: one
complex polyphase filter per radio, with the output split across threads. The
output is phase-preserving complex baseband that a standard receiver can
decode, for example an 802.11 OFDM receiver.

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
 wideband baseband --> band_splitter --> radio A (fc - offset) -.
                                    \--> radio B (fc + offset) --> bonded_sink
```

`band_splitter` cuts the wideband signal in two at `fc`. Each half is shifted
to DC and resampled through a lowpass of half-width |offset|, so the −6 dB
points of both filters meet at the centre. For OFDM signals with an unused DC
subcarrier, such as 802.11, the cut falls on that empty subcarrier. A standard
receiver's per-subcarrier channel estimate absorbs each radio's own gain,
phase and delay, so the transmitter needs no phase calibration.

Only the few subcarriers inside the split filter's transition are radiated by
both radios. They can fade, depending on the relative phase. `taps_per_phase`
(default 64) keeps that zone within ±1 subcarrier of 802.11. If the relative
phase is near 180°, the subcarriers next to the seam get some ISI. That only
matters for dense constellations (64-QAM). The phase can be changed while
running (`set_phase_offset`) to bring it near 0°, for example from a receiver's
measurement.

`bonded_sink` streams both radios **continuously** for the whole session. A
B210 needs several hundred microseconds to ramp its TX power after idling,
which is longer than an 802.11 preamble. So the gaps between input bursts are
filled with zeros, and every burst goes out at full power. A scheduler thread
hands the same sample sequence to one sender thread per radio, both starting
at the same device time. If either radio underruns, both restart on a new
common time.

Recommended TX plan for 802.11 with a channel of width `T`: per-radio rate
`T/2`, radio offsets `∓ 14/64 · T`. Each radio then carries the 26 used
subcarriers of its side, and the seam stays away from the edge of each radio's
own band, where its DUC rolls off. `∓ T/4` also works.

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

## Requirements

- Two (per direction) Ettus B210, all sharing one external 10 MHz reference and
  PPS (for example from an OctoClock). Ideally each radio has its own USB 3
  port.
- UHD 4.x with development headers (`UHDConfig.cmake`). Tested with 4.10.
- GNU Radio 3.10 with development headers and pybind11. Tested with 3.10.12.
- CMake ≥ 3.16 and a C++17 compiler.
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
| `overlap_reconstructor_test` | gain and phase recovery, band flatness including the seam, both frequency plans, warmup and frozen paths (and that the frozen path continues the warmup timeline), streaming, thread count invariance, direct-output API, config checks |
| `bonded_resample_test` | filter design, dot-product kernels, streaming polyphase resamplers, sliced computation |
| `bonded_receiver_test` | timestamp alignment of the receive rings: stale data, holes after an overflow, multi-channel, two-device alignment, bookkeeping |
| `bonded_transmitter_test` | configuration checks, phase-offset API, frequency/phase correction across chunks |
| `band_splitter_test` | flatness, half assignment, seam width, pulled-in plan, streaming vs one-shot, output length, config checks |
| `bonded_csv_test` | parsing of the comma-separated block parameters |
| `qa_band_splitter`, `qa_overlap_reconstructor` | the GNU Radio blocks in flowgraphs (Python) |

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

Set the CPU frequency governor to `performance` for high bandwidths.

## Measured capacity

On a Ryzen 7 5700G (8 cores), with the frequency plans above:

| | bandwidth |
|---|---|
| RX, reconstruction only | 80 MHz real time (64 MHz with margin) |
| RX with a GNU Radio 802.11 receiver attached | 40 MHz decoding at parity with a single-radio receiver |
| TX, bonded, 802.11 frames | 20, 40 and 60 MHz without underruns |

## Source layout

| path | content |
|---|---|
| `lib/overlap_reconstructor.{hpp,cpp}` | RX reconstruction engine |
| `lib/band_splitter.{hpp,cpp}` | TX split engine |
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
