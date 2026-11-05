#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Flowgraph tests of the band_splitter block (no hardware)."""

import os
import sys
import unittest

import numpy as np
from gnuradio import blocks, gr

if os.environ.get("BONDED_USRP_PYTHON_DIR"):
    sys.path.insert(0, os.environ["BONDED_USRP_PYTHON_DIR"])
    import bonded_usrp_python as bonded_usrp
else:
    from gnuradio import bonded_usrp


def run_splitter(x, rate_in, rate_out, offset, taps=64):
    tb = gr.top_block()
    src = blocks.vector_source_c(x.tolist(), False)
    spl = bonded_usrp.band_splitter(rate_in, rate_out, -offset, offset, taps)
    s0, s1 = blocks.vector_sink_c(), blocks.vector_sink_c()
    tb.connect(src, spl)
    tb.connect((spl, 0), s0)
    tb.connect((spl, 1), s1)
    tb.run()
    return np.array(s0.data()), np.array(s1.data())


def tone(f, fs, n):
    return np.exp(2j * np.pi * f * np.arange(n) / fs).astype(np.complex64)


class qa_band_splitter(unittest.TestCase):
    def test_every_input_sample_comes_out(self):
        # The block must not hold the tail of its input back waiting for more
        # (it would only leave with the next burst): the whole input of a
        # finite source must come out, at exactly the resampling ratio.
        for rate_in, rate_out, offset in ((20e6, 10e6, 4.375e6), (40e6, 25e6, 8.75e6)):
            for n in (1000, 30000, 262144 + 17):
                x = (np.random.default_rng(n).standard_normal(n)
                     + 1j * np.random.default_rng(n + 1).standard_normal(n)).astype(np.complex64)
                o0, o1 = run_splitter(x, rate_in, rate_out, offset)
                expect = n * rate_out / rate_in
                self.assertEqual(len(o0), len(o1))
                self.assertLessEqual(abs(len(o0) - expect), 2, (rate_in, n, len(o0)))

    def test_halves_go_to_their_radio(self):
        fs, n, off = 20e6, 40000, 4.375e6
        o0, o1 = run_splitter(tone(-6e6, fs, n), fs, 10e6, off)
        lo, hi = np.mean(np.abs(o0[2000:-2000]) ** 2), np.mean(np.abs(o1[2000:-2000]) ** 2)
        self.assertGreater(lo, 0.9)
        self.assertLess(hi, 1e-3)
        o0, o1 = run_splitter(tone(+6e6, fs, n), fs, 10e6, off)
        lo, hi = np.mean(np.abs(o0[2000:-2000]) ** 2), np.mean(np.abs(o1[2000:-2000]) ** 2)
        self.assertGreater(hi, 0.9)
        self.assertLess(lo, 1e-3)

    def test_lower_slice_is_shifted_to_dc(self):
        # A tone at fc - 6 MHz lands at -6 + 4.375 = -1.625 MHz in radio A.
        fs, rate, off = 20e6, 10e6, 4.375e6
        o0, _ = run_splitter(tone(-6e6, fs, 40000), fs, rate, off)
        seg = o0[2000:-2000]
        spec = np.abs(np.fft.fftshift(np.fft.fft(seg)))
        freqs = np.fft.fftshift(np.fft.fftfreq(len(seg), 1 / rate))
        self.assertAlmostEqual(freqs[np.argmax(spec)], -6e6 + off, delta=rate / len(seg) * 2)

    def test_rejects_bad_config(self):
        with self.assertRaises(Exception):
            bonded_usrp.band_splitter(20e6, 10e6, -6e6, 6e6, 64)  # 12 MHz slice in 10 MS/s


if __name__ == "__main__":
    unittest.main()
