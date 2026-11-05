//
// Copyright 2026 Gabriel Lando
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Subcarrier-domain TX split of CP-OFDM bursts. Per input symbol:
//   1. FFT of the symbol body (the fft_len samples after its cyclic prefix);
//   2. per radio, its subcarriers moved onto a radio-rate FFT grid centred on
//      that radio, inverse FFT, cyclic prefix and suffix, raised-cosine edges;
//   3. overlap-added at the symbol's place on the radio-rate timeline, with
//      the translation phase of that place so consecutive symbols stay
//      coherent.
//

#include "ofdm_splitter.hpp"
#include "bonded_resample.hpp"

#include <gnuradio/fft/fft.h>

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <string>

namespace bonded {

namespace {
using namespace detail; // cf, kPi

struct bin_map
{
    size_t in;  // input FFT bin
    size_t out; // radio FFT bin
    float w;
};
} // namespace

struct ofdm_splitter::impl
{
    size_t N, CP, S;     // input symbol geometry
    size_t NR, CPR, SR;  // radio symbol geometry
    size_t WR;           // window length (radio samples)
    uint64_t interp, decim;
    long bins[2];        // radio centres, in input subcarriers
    std::vector<bin_map> maps[2];
    std::vector<float> ramp;
    gr::fft::fft_complex_fwd fwd;
    gr::fft::fft_complex_rev rev;

    // Timeline: in_pos input samples consumed; output acc[r] starts at
    // absolute radio sample acc_start, and [acc_start, final_end) is complete
    // (no later symbol can add to it).
    uint64_t in_pos = 0, acc_start = 0, final_end = 0;
    std::vector<cf> acc[2];

    bool in_burst      = false;
    uint64_t burst_left = 0; // input samples of the current burst not consumed
    uint64_t q_next     = 0; // radio sample where the next symbol starts
    std::vector<cf> sym;
    size_t sym_fill = 0;

    impl(const ofdm_splitter_config& c, size_t nr, size_t cpr, size_t wr,
        uint64_t L, uint64_t M, long b1, long b2)
        : N(c.fft_len), CP(c.cp_len), S(c.fft_len + c.cp_len), NR(nr), CPR(cpr),
          SR(nr + cpr), WR(wr), interp(L), decim(M), bins{b1, b2},
          fwd(static_cast<int>(c.fft_len)), rev(static_cast<int>(nr)), sym(S)
    {
        // Split at the midpoint between the radio centres: radio 0 takes the
        // subcarriers below it, radio 1 those above, half each on it.
        const long mid2 = b1 + b2;
        const long half = static_cast<long>(NR / 2);
        for (size_t j = 0; j < N; j++) {
            const long k = (j < (N + 1) / 2) ? static_cast<long>(j)
                                             : static_cast<long>(j) - static_cast<long>(N);
            for (int r = 0; r < 2; r++) {
                const long side = (r == 0) ? mid2 - 2 * k : 2 * k - mid2;
                const float w   = side > 0 ? 1.0f : (side == 0 ? 0.5f : 0.0f);
                const long rb   = k - bins[r];
                if (w == 0.0f || rb <= -half || rb >= half) {
                    continue; // other radio's, or beyond this radio's Nyquist
                }
                const long nrl = static_cast<long>(NR);
                maps[r].push_back(
                    { j, static_cast<size_t>(((rb % nrl) + nrl) % nrl),
                        w / static_cast<float>(N) });
            }
        }
        ramp.resize(WR);
        for (size_t t = 0; t < WR; t++) {
            ramp[t] = static_cast<float>(0.5
                - 0.5 * std::cos(kPi * (static_cast<double>(t) + 0.5)
                                 / static_cast<double>(WR)));
        }
    }

    uint64_t out_index(uint64_t in) const
    {
        return static_cast<uint64_t>(
            (static_cast<unsigned __int128>(in) * interp + decim / 2) / decim);
    }

    void ensure(uint64_t end)
    {
        const size_t need = static_cast<size_t>(end - acc_start);
        if (acc[0].size() < need) {
            acc[0].resize(need, cf(0, 0));
            acc[1].resize(need, cf(0, 0));
        }
    }

    void advance_final(uint64_t end)
    {
        if (end > final_end) {
            ensure(end);
            final_end = end;
        }
    }

    //! Input outside any burst: silence, keeping the timeline.
    void idle(size_t k)
    {
        in_pos += k;
        advance_final(out_index(in_pos));
    }

    void synthesize(uint64_t q)
    {
        ensure(q + SR + WR);
        std::copy(sym.begin() + static_cast<std::ptrdiff_t>(CP), sym.end(), fwd.get_inbuf());
        fwd.execute();
        const cf* X = fwd.get_outbuf();
        for (int r = 0; r < 2; r++) {
            cf* Y = rev.get_inbuf();
            std::fill(Y, Y + NR, cf(0, 0));
            for (const auto& m : maps[r]) {
                Y[m.out] += X[m.in] * m.w;
            }
            rev.execute();
            const cf* body = rev.get_outbuf();
            // Translation phase at radio sample q: exp(-j 2 pi bins[r] q / NR).
            const long nrl     = static_cast<long>(NR);
            const uint64_t bm  = static_cast<uint64_t>(((bins[r] % nrl) + nrl) % nrl);
            const uint64_t cyc = (bm * (q % NR)) % NR;
            const cf ph        = std::polar(1.0f,
                static_cast<float>(-2.0 * kPi * static_cast<double>(cyc) / static_cast<double>(NR)));
            cf* dst = acc[r].data() + (q - acc_start);
            for (size_t m = 0; m < SR + WR; m++) {
                float w = 1.0f;
                if (m < WR) {
                    w = ramp[m];
                } else if (m >= SR) {
                    w = ramp[SR + WR - 1 - m];
                }
                dst[m] += body[(m + NR - CPR) % NR] * (ph * w);
            }
        }
    }

    void end_burst()
    {
        in_burst = false;
        sym_fill = 0;
        advance_final(q_next + WR); // the last symbol's suffix is final now
    }

    //! Ends the burst as soon as no whole symbol is left in it, so its
    //! output is released without waiting for the trailing samples.
    void maybe_end()
    {
        if (in_burst && sym_fill == 0 && burst_left < S) {
            end_burst();
        }
    }

    void start_burst(uint64_t len)
    {
        if (in_burst) {
            end_burst();
        }
        in_burst   = true;
        burst_left = len;
        sym_fill   = 0;
        q_next     = std::max(out_index(in_pos), final_end);
        advance_final(q_next);
        maybe_end();
    }

    //! Consumes up to k samples of the current burst; returns how many.
    size_t feed_burst(const cf* p, size_t k)
    {
        k = static_cast<size_t>(std::min<uint64_t>(k, burst_left));
        size_t used = 0;
        while (in_burst && used < k) {
            const size_t c = std::min(k - used, S - sym_fill);
            std::copy(p + used, p + used + c, sym.begin() + static_cast<std::ptrdiff_t>(sym_fill));
            sym_fill += c;
            used += c;
            burst_left -= c;
            in_pos += c;
            if (sym_fill == S) {
                synthesize(q_next);
                q_next += SR;
                advance_final(q_next);
                sym_fill = 0;
            }
            maybe_end();
        }
        return used;
    }

    void process(const cf* in, size_t n, const std::vector<burst>& bursts,
        std::vector<cf>& out1, std::vector<cf>& out2)
    {
        size_t i = 0, bi = 0;
        while (i < n) {
            if (bi < bursts.size() && bursts[bi].offset <= i) {
                start_burst(bursts[bi].length);
                bi++;
                continue;
            }
            const size_t lim = (bi < bursts.size()) ? std::min(bursts[bi].offset, n) : n;
            if (in_burst) {
                i += feed_burst(in + i, lim - i);
            } else {
                // Outside a burst, or the trailing partial symbol of one.
                idle(lim - i);
                i = lim;
            }
        }
        const size_t ready = static_cast<size_t>(final_end - acc_start);
        for (int r = 0; r < 2; r++) {
            auto& dst = (r == 0) ? out1 : out2;
            dst.insert(dst.end(), acc[r].begin(), acc[r].begin() + static_cast<std::ptrdiff_t>(ready));
            acc[r].erase(acc[r].begin(), acc[r].begin() + static_cast<std::ptrdiff_t>(ready));
        }
        acc_start = final_end;
    }
};

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

namespace {
struct geometry
{
    uint64_t L, M;
    size_t NR, CPR, WR;
    long b1, b2;
};

geometry check(const ofdm_splitter_config& c)
{
    auto fail = [](const std::string& m) {
        throw std::invalid_argument("ofdm_splitter: " + m);
    };
    if (c.input_rate <= 0 || c.per_radio_rate <= 0) {
        fail("sample rates must be > 0");
    }
    if (c.per_radio_rate > c.input_rate) {
        fail("per_radio_rate must be <= input_rate");
    }
    if (c.fft_len < 2 || c.cp_len == 0) {
        fail("fft_len must be >= 2 and cp_len > 0");
    }
    geometry g{};
    const auto ri = static_cast<uint64_t>(std::llround(c.input_rate));
    const auto ro = static_cast<uint64_t>(std::llround(c.per_radio_rate));
    const uint64_t d = std::gcd(ri, ro);
    g.L = ro / d;
    g.M = ri / d;
    if ((c.fft_len * g.L) % g.M || (c.cp_len * g.L) % g.M) {
        fail("fft_len and cp_len must be whole numbers of radio samples");
    }
    g.NR  = static_cast<size_t>(c.fft_len * g.L / g.M);
    g.CPR = static_cast<size_t>(c.cp_len * g.L / g.M);
    if (c.window_len > c.cp_len) {
        fail("window_len must not exceed cp_len");
    }
    g.WR = static_cast<size_t>(c.window_len * g.L / g.M);
    auto to_bins = [&](double f) {
        const double b = f * static_cast<double>(c.fft_len) / c.input_rate;
        if (std::abs(b - std::round(b)) > 1e-6) {
            fail("radio offsets must be whole subcarriers (input_rate/fft_len)");
        }
        return static_cast<long>(std::lround(b));
    };
    g.b1 = to_bins(c.radio1_offset_hz);
    g.b2 = to_bins(c.radio2_offset_hz);
    if (g.b1 >= g.b2) {
        fail("radio1_offset must be below radio2_offset");
    }
    return g;
}
} // namespace

ofdm_splitter::ofdm_splitter(const ofdm_splitter_config& cfg) : _cfg(cfg)
{
    reset();
}

ofdm_splitter::~ofdm_splitter() = default;

void ofdm_splitter::reset()
{
    const geometry g = check(_cfg);
    _impl = std::make_unique<impl>(_cfg, g.NR, g.CPR, g.WR, g.L, g.M, g.b1, g.b2);
}

void ofdm_splitter::process(const std::complex<float>* in,
    size_t n,
    const std::vector<burst>& bursts,
    std::vector<std::complex<float>>& out1,
    std::vector<std::complex<float>>& out2)
{
    _impl->process(in, n, bursts, out1, out2);
}

size_t ofdm_splitter::radio_fft_len() const { return _impl->NR; }
size_t ofdm_splitter::radio_cp_len() const { return _impl->CPR; }
size_t ofdm_splitter::radio_window_len() const { return _impl->WR; }

} // namespace bonded
