#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Flowgraph tests of the ofdm_splitter block (no hardware)."""

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

FFT, CP = 64, 16
T, R, OFF = 20e6, 12.5e6, 3.75e6  # default plan: 5/8 rate, 12 subcarriers


def ofdm_burst(nsym, seed, bins=None):
    """Back-to-back CP-OFDM symbols, random QPSK on `bins` (default +-1..26)."""
    rng = np.random.default_rng(seed)
    ks = [k for k in range(-26, 27) if k != 0] if bins is None else bins
    out = []
    for _ in range(nsym):
        X = np.zeros(FFT, complex)
        X[np.array(ks) % FFT] = (rng.choice([-1, 1], len(ks)) + 1j * rng.choice([-1, 1], len(ks)))
        body = np.fft.ifft(X) * 0.3
        out.append(np.concatenate([body[-CP:], body]))
    return np.concatenate(out).astype(np.complex64)


def run_splitter(x, bursts):
    """bursts: [(offset, length)], tagged as packet_len on the input."""
    tags = [gr.tag_utils.python_to_tag((o, pmt.intern("packet_len"), pmt.from_long(n)))
            for o, n in bursts]
    tb = gr.top_block()
    src = blocks.vector_source_c(x.tolist(), False, 1, tags)
    spl = bonded_usrp.ofdm_splitter(T, R, -OFF, OFF)
    s0, s1 = blocks.vector_sink_c(), blocks.vector_sink_c()
    tb.connect(src, spl)
    tb.connect((spl, 0), s0)
    tb.connect((spl, 1), s1)
    tb.run()
    return np.array(s0.data()), np.array(s1.data())


class qa_ofdm_splitter(unittest.TestCase):
    def test_tagged_bursts_keep_the_timeline(self):
        # Silence keeps the radio rate, and each burst comes out where its tag
        # put it, in full, without waiting for more input.
        b1, b2 = ofdm_burst(6, 1), ofdm_burst(4, 2)
        x = np.zeros(8000, np.complex64)
        x[1600:1600 + len(b1)] = b1
        x[5600:5600 + len(b2)] = b2
        o0, o1 = run_splitter(x, [(1600, len(b1)), (5600, len(b2))])
        self.assertEqual(len(o0), len(o1))
        self.assertEqual(len(o0), round(len(x) * R / T))
        nz = np.nonzero(np.abs(o0) > 0)[0]
        self.assertEqual(nz[0], 1000)  # 1600 * 5/8
        self.assertGreater(nz[-1], 3500 + 4 * 50 - 2)  # 5600 * 5/8 + 4 symbols

    def test_last_burst_is_not_held_back(self):
        # A burst that ends the input is released whole: its last symbol and
        # the window tail after it.
        b = ofdm_burst(5, 3)
        o0, _ = run_splitter(b, [(0, len(b))])
        self.assertEqual(len(o0), 5 * 50 + 5)  # 5 symbols of 50 + 5-sample tail

    def test_halves_go_to_their_radio(self):
        lo = ofdm_burst(4, 4, bins=[-20, -9, -2])
        hi = ofdm_burst(4, 5, bins=[2, 9, 20])
        for x, low in ((lo, True), (hi, False)):
            o0, o1 = run_splitter(x, [(0, len(x))])
            e0, e1 = np.sum(np.abs(o0) ** 2), np.sum(np.abs(o1) ** 2)
            big, small = (e0, e1) if low else (e1, e0)
            self.assertGreater(big, 1e-3)
            self.assertLess(small, 1e-9 * big)

    def test_untagged_input_is_silence(self):
        o0, o1 = run_splitter(ofdm_burst(5, 6), [])
        self.assertTrue(np.all(o0 == 0) and np.all(o1 == 0))

    def test_rejects_bad_config(self):
        with self.assertRaises(Exception):
            bonded_usrp.ofdm_splitter(T, 12e6, -OFF, OFF)  # 64 * 12/20 not whole
        with self.assertRaises(Exception):
            bonded_usrp.ofdm_splitter(T, R, -3.7e6, 3.7e6)  # not whole subcarriers


if __name__ == "__main__":
    unittest.main()
