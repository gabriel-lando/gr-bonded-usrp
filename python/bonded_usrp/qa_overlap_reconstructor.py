#!/usr/bin/env python3
"""Flowgraph tests of the overlap_reconstructor block (no hardware)."""

import os
import sys
import unittest

import numpy as np
import pmt
from gnuradio import blocks, gr

if os.environ.get("BONDED_USRP_PYTHON_DIR"):
    sys.path.insert(0, os.environ["BONDED_USRP_PYTHON_DIR"])
    import bonded_usrp_python as bonded_usrp
else:
    from gnuradio import bonded_usrp


def rx_time_tag(offset, t):
    secs = int(t)
    return gr.tag_utils.python_to_tag(
        (offset, pmt.intern("rx_time"),
         pmt.make_tuple(pmt.from_uint64(secs), pmt.from_double(t - secs))))


def run_recon(rx1, rx2, rate_in, rate_out, target_bw, overlap, off, tags=()):
    tb = gr.top_block()
    s1 = blocks.vector_source_c(rx1.tolist(), False, 1, list(tags))
    s2 = blocks.vector_source_c(rx2.tolist(), False, 1, list(tags))
    rec = bonded_usrp.overlap_reconstructor(rate_in, rate_out, target_bw, overlap, -off, off)
    snk = blocks.vector_sink_c()
    tb.connect(s1, (rec, 0))
    tb.connect(s2, (rec, 1))
    tb.connect(rec, snk)
    tb.run()
    return np.array(snk.data())


def tone(f, fs, n, phase=0.0):
    return np.exp(1j * (2 * np.pi * f * np.arange(n) / fs + phase)).astype(np.complex64)


class qa_overlap_reconstructor(unittest.TestCase):
    # Edge-aligned plan for a 20 MHz channel and a 5 MHz overlap.
    T, O = 20e6, 5e6
    PR, OFF = 12.5e6, 3.75e6

    def slices(self, f, n, phase2=0.0, pilot=1.1e6, pilot_amp=0.3):
        """What the two radios capture of a tone at fc + f plus a weaker pilot
        at fc + pilot (each at PR; a tone outside a radio's band is removed by
        its anti-alias filter). The pilot sits in the overlap, so the warmup
        always has common content to align on, as a real channel does."""
        def radio(center, phase):
            x = np.zeros(n, dtype=np.complex64)
            for ff, a in ((f, 1.0), (pilot, pilot_amp)):
                if abs(ff - center) < 0.5 * self.PR:
                    x += a * tone(ff - center, self.PR, n, phase)
            return x
        return radio(-self.OFF, 0.0), radio(self.OFF, phase2)

    def amplitude(self, x, f):
        k = np.arange(len(x))
        return abs(np.mean(x * np.exp(-2j * np.pi * f * k / self.T)))

    def test_output_rate(self):
        n = 400000
        rx1, rx2 = self.slices(-6e6, n)
        out = run_recon(rx1, rx2, self.PR, self.T, self.T, self.O, self.OFF)
        # The block works in 8192-sample output chunks; at the end of a
        # finite stream the scheduler may stop it a few chunks early.
        self.assertGreater(len(out), n * self.T / self.PR - 4 * 8192)
        self.assertLessEqual(len(out), n * self.T / self.PR + 8192)

    def test_tone_lands_at_its_frequency(self):
        for f in (-7e6, -1e6, 0.5e6, 6.5e6):
            rx1, rx2 = self.slices(f, 400000)
            out = run_recon(rx1, rx2, self.PR, self.T, self.T, self.O, self.OFF)
            seg = out[300000:-20000]
            spec = np.abs(np.fft.fft(seg))
            freqs = np.fft.fftfreq(len(seg), 1 / self.T)
            self.assertAlmostEqual(freqs[np.argmax(spec)], f, delta=2 * self.T / len(seg))
            self.assertAlmostEqual(self.amplitude(seg, f), 1.0, delta=0.03)

    def test_radio2_phase_is_aligned(self):
        # Overlap content with radio 2 rotated by an arbitrary phase: the two
        # contributions must add coherently (no notch at the seam). The content
        # is symmetric about fc, like 802.11: the gain estimate is a power
        # ratio over the overlap, whose outer parts each radio's resampler
        # attenuates, and that only cancels for symmetric content.
        for f in (0.0, 0.3e6, 1.5e6, 2.0e6):
            pilot = -f if f else 1e9  # 1e9: no pilot (out of every band)
            rx1, rx2 = self.slices(f, 400000, phase2=2.1, pilot=pilot, pilot_amp=1.0)
            out = run_recon(rx1, rx2, self.PR, self.T, self.T, self.O, self.OFF)
            self.assertAlmostEqual(self.amplitude(out[300000:-20000], f), 1.0, delta=0.03)

    def test_rx_time_gap_is_stepped_over(self):
        # bonded_source drops a block both radios have a hole in and tags
        # the next sample with its rx_time. The reconstructor must continue
        # the radios' timeline across the gap: closing it up would move each
        # radio's translation phase differently and notch the overlap.
        os.environ["BONDED_TRACK"] = "0"  # deterministic alignment
        try:
            n, at, gap = 450000, 150000, 131072  # gap after the warmup
            k = np.arange(n)
            # Content across the whole channel, overlap included.
            rx1, rx2 = self.slices(0.4e6, n, phase2=1.3, pilot=-0.4e6, pilot_amp=1.0)
            for f in (-7e6, 6e6):
                a, b = self.slices(f, n, phase2=1.3, pilot=1e9)
                rx1, rx2 = rx1 + a, rx2 + b
            ref = run_recon(rx1, rx2, self.PR, self.T, self.T, self.O, self.OFF)
            keep = np.r_[0:at, at + gap:n]
            tags = [rx_time_tag(0, 5.0), rx_time_tag(at, 5.0 + (at + gap) / self.PR)]
            out = run_recon(rx1[keep], rx2[keep], self.PR, self.T, self.T, self.O, self.OFF, tags)
            closed = run_recon(rx1[keep], rx2[keep], self.PR, self.T, self.T, self.O, self.OFF)
        finally:
            del os.environ["BONDED_TRACK"]
        shift = (gap // 5) * 8  # outputs of the skipped whole resampler periods
        lo, hi = int(at * self.T / self.PR) + 8192, len(out) - 20000

        def rel_err(x, s):
            seg, r = x[lo:hi], ref[lo + s:hi + s]
            return np.sum(np.abs(seg - r) ** 2) / np.sum(np.abs(r) ** 2)

        self.assertLess(rel_err(out, shift), 1e-6)
        self.assertGreater(min(rel_err(closed, s) for s in range(shift - 2, shift + 3)), 1e-2)


if __name__ == "__main__":
    unittest.main()
