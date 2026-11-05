//
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Time-domain reconstruction: each radio is resampled (polyphase), translated
// to the target baseband, and combined through complementary complex crossfade
// FIR filters (which also discard the out-of-band corners).  Radio 2 is aligned
// to radio 1 (gain + phase) using a running estimate over the shared overlap.
// After the warmup the alignment is frozen and each radio runs through one
// fused complex polyphase filter, with the output range of every call split
// across worker threads.
//

#include "overlap_reconstructor.hpp"
#include "bonded_resample.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <thread>

#ifdef BONDED_HAVE_VOLK
#    include <volk/volk.h>
#endif

namespace bonded {

namespace {

using namespace detail; // cf, cd, kPi, hamming, design_lowpass_real, dot_cr,
                        // dot_cc_hot, poly_resampler, poly_resampler_c

bool env_flag_enabled(const char* name)
{
    const char* value = std::getenv(name);
    if (value == nullptr || *value == '\0') {
        return false;
    }
    return std::string(value) != "0";
}

//! Same, but defaulting to ON when the variable is unset.
bool env_flag_enabled_default_on(const char* name)
{
    const char* value = std::getenv(name);
    if (value == nullptr || *value == '\0') {
        return true;
    }
    return std::string(value) != "0";
}

//! Hint the CPU we are in a spin-wait (cheap; avoids a full yield syscall).
inline void cpu_relax()
{
#if defined(__x86_64__) || defined(__i386__)
    __builtin_ia32_pause();
#else
    std::this_thread::yield();
#endif
}

// --- minimal radix-2 FFT (double), used only for one-time filter design ------
void design_fft(std::vector<cd>& a, bool inverse)
{
    const size_t n = a.size();
    for (size_t i = 1, j = 0; i < n; i++) {
        size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) {
            j ^= bit;
        }
        j ^= bit;
        if (i < j) {
            std::swap(a[i], a[j]);
        }
    }
    for (size_t len = 2; len <= n; len <<= 1) {
        const double ang = (inverse ? 2.0 : -2.0) * kPi / static_cast<double>(len);
        const cd wlen(std::cos(ang), std::sin(ang));
        for (size_t i = 0; i < n; i += len) {
            cd w(1.0, 0.0);
            for (size_t k = 0; k < len / 2; k++) {
                const cd u = a[i + k];
                const cd v = a[i + k + len / 2] * w;
                a[i + k]           = u + v;
                a[i + k + len / 2] = u - v;
                w *= wlen;
            }
        }
    }
    if (inverse) {
        for (auto& x : a) {
            x /= static_cast<double>(n);
        }
    }
}

// Complex FIR designed by frequency sampling of a (complex) target response,
// windowed to L taps. Returns taps in convolution order b[0..L-1].
std::vector<cf> design_complex_fir(size_t L, double fs,
    const std::function<cd(double)>& response)
{
    size_t N = 1;
    while (N < 8 * L) {
        N <<= 1;
    }
    std::vector<cd> H(N);
    const double df = fs / static_cast<double>(N);
    for (size_t k = 0; k < N; k++) {
        const double f = (k < N / 2 ? static_cast<double>(k)
                                    : static_cast<double>(k) - static_cast<double>(N))
                         * df;
        H[k] = response(f);
    }
    design_fft(H, true); // -> impulse response (periodic)
    std::vector<cf> b(L);
    const long half = static_cast<long>(L / 2);
    for (size_t i = 0; i < L; i++) {
        const long n = static_cast<long>(i) - half; // centered, may be negative
        const size_t idx = static_cast<size_t>((n % static_cast<long>(N) + static_cast<long>(N))
                                               % static_cast<long>(N));
        const cd v = H[idx] * hamming(i, L);
        b[i] = cf(static_cast<float>(v.real()), static_cast<float>(v.imag()));
    }
    return b;
}

// Complex FIR taps stored reversed (window order) for the warmup crossfade.
struct cfir
{
    std::vector<cf> taps;
    size_t K = 0;

    void init(std::vector<cf> conv_taps)
    {
        K = conv_taps.size();
        std::reverse(conv_taps.begin(), conv_taps.end());
        taps = std::move(conv_taps);
    }
};

// Streaming FIR (decim=1), real taps, complex data (for overlap analysis).
struct rfir
{
    std::vector<float> taps; // window order
    size_t K = 0;
    std::vector<cf> buf;     // persistent [tail | block] working buffer

    void init(std::vector<float> conv_taps)
    {
        K = conv_taps.size();
        std::reverse(conv_taps.begin(), conv_taps.end());
        taps = std::move(conv_taps);
        buf.assign(K - 1, cf(0, 0));
    }

    void filter(const cf* x, size_t n, std::vector<cf>& out)
    {
        buf.resize(K - 1 + n);
        std::copy(x, x + n, buf.begin() + static_cast<std::ptrdiff_t>(K - 1));
        const size_t base = out.size();
        out.resize(base + n);
        for (size_t i = 0; i < n; i++) {
            out[base + i] = dot_cr(buf.data() + i, taps.data(), K);
        }
        std::copy(buf.end() - static_cast<std::ptrdiff_t>(K - 1), buf.end(),
            buf.begin());
        buf.resize(K - 1);
    }
};

} // namespace

struct overlap_reconstructor::impl
{
    using clock = std::chrono::steady_clock;

    struct profile_counters {
        clock::duration total{};
        clock::duration resample{};
        clock::duration rotate{};
        clock::duration align{};
        clock::duration crossfade{};
        size_t calls = 0;
        size_t parallel_calls = 0;
        size_t input_samples = 0;
        size_t output_samples = 0;
    };

    overlap_reconstructor_config cfg;
    size_t interp, decim;
    // Translation phase of each radio at the next output sample; shared by
    // both paths so the frozen path continues the warmup's phase exactly.
    cf inc1, inc2, rot1{1, 0}, rot2{1, 0};

    poly_resampler rs1, rs2;     // warmup path (real proto + separate rotate/crossfade)
    poly_resampler_c rs1c, rs2c; // frozen path (fused complex polyphase per radio)
    cfir hlow, hhigh;
    rfir ov1, ov2; // overlap-band lowpass for alignment

    // Running alignment statistics over the overlap.
    double Saa = 0.0, Sbb = 0.0;
    cd Sxy{0, 0};
    cf alpha{1, 0};
    overlap_alignment align;
    // Alignment is static (shared clock); estimate during a warmup window then
    // freeze it and skip the overlap-analysis filters (a real-time saving).
    bool frozen = false;
    size_t processed = 0;
    size_t kWarmup = 200000;

    // Reusable scratch buffers (retain capacity across calls).
    std::vector<cf> r1_, r2_, o1_, o2_, out_;
    // Warmup crossfade history buffers ([tail | rotated block]) and the K-1
    // sample tails carried across calls.
    std::vector<cf> ch1_, ch2_, tail1_, tail2_;
    // Fused-path resampler outputs (already band-shaped + translated per radio).
    std::vector<cf> a_, b_;
    bool profile_enabled = false;
    profile_counters profile_stats;
    clock::time_point last_profile_log{};

    // Optional caller-supplied output buffer for the frozen path's combine
    // (see process_into); null means "write everything to out_".
    cf* direct_dst_        = nullptr;
    size_t direct_cap_     = 0;
    size_t direct_written_ = 0;

    // --- Parallelism (frozen/fused path only) -------------------------------
    // The fused per-radio filters dominate the cost (~98% of process()) and
    // every output sample is independent of the others, so each call's output
    // range is cut into `n_threads_` slices; each slice runs both radios' FIR,
    // their frequency rotation (closed-form start phase) and the combine for
    // its range. The caller runs slice 0, persistent workers the rest.
    // BONDED_THREADS sets the count (default 4); BONDED_PARALLEL=0 forces 1.
    //
    // Handoff is atomic sequence/done counters with a bounded spin, not a
    // condition variable: a slice is only ~0.1-0.2 ms of work, while a condvar
    // wake costs ~200 us. The waiters back off spin -> yield -> short sleep, so
    // an idle flowgraph does not pin cores.
    std::vector<std::thread> workers_;
    std::atomic<std::uint32_t> job_seq_{0};  // bumped by the caller to submit
    std::atomic<std::uint32_t> done_{0};     // slices finished by the workers
    std::atomic<bool> worker_stop_{false};
    size_t n_threads_ = 4;
    bool parallel_enabled = true;

    // The current call's batch, published before job_seq_ (release).
    struct batch {
        size_t p1 = 0, p2 = 0, mf = 0, direct = 0, slices = 1;
        cf rotA{1, 0}, rotB{1, 0}, alpha{1, 0};
    } job_;

    //! Spin -> yield -> sleep backoff. `spins` is the caller's iteration count.
    static void backoff(int spins)
    {
        constexpr int kSpin  = 20000; // ~0.2 ms of pauses: covers one slice
        constexpr int kYield = 25000;
        if (spins < kSpin) {
            cpu_relax();
        } else if (spins < kYield) {
            std::this_thread::yield();
        } else {
            std::this_thread::sleep_for(std::chrono::microseconds(100));
        }
    }

    //! Slice k of the current batch: outputs [i0, i1) of both radios, rotated,
    //! then combined into the caller's buffer / out_.
    void run_slice(size_t k)
    {
        const batch& j = job_;
        const size_t total = std::max(j.p1, j.p2);
        const size_t i0 = total * k / j.slices;
        const size_t i1 = total * (k + 1) / j.slices;
        cf* a = a_.data();
        cf* b = b_.data();
        const size_t e1 = std::min(i1, j.p1), e2 = std::min(i1, j.p2);
        if (i0 < e1) {
            rs1c.compute(a, i0, e1);
            rotate_tile(a + i0, e1 - i0, j.rotA * inc_pow(inc1, i0), inc1);
        }
        if (i0 < e2) {
            rs2c.compute(b, i0, e2);
            rotate_tile(b + i0, e2 - i0, j.rotB * inc_pow(inc2, i0), inc2);
        }
        const size_t c1 = std::min(i1, j.mf);
        const size_t split = std::min(std::max(i0, j.direct), c1);
        for (size_t i = i0; i < split; i++) {
            direct_dst_[i] = a[i] + j.alpha * b[i];
        }
        cf* o = out_.data() - j.direct;
        for (size_t i = split; i < c1; i++) {
            o[i] = a[i] + j.alpha * b[i];
        }
    }

    void worker_loop(size_t k)
    {
        std::uint32_t seen = 0;
        while (true) {
            int spins = 0;
            std::uint32_t s;
            while ((s = job_seq_.load(std::memory_order_acquire)) == seen) {
                if (worker_stop_.load(std::memory_order_relaxed)) {
                    return;
                }
                backoff(spins++);
            }
            if (worker_stop_.load(std::memory_order_relaxed)) {
                return;
            }
            seen = s;
            run_slice(k);
            done_.fetch_add(1, std::memory_order_release);
        }
    }

    ~impl()
    {
        worker_stop_.store(true, std::memory_order_relaxed);
        job_seq_.fetch_add(1, std::memory_order_release); // wake the waiters
        for (auto& w : workers_) {
            w.join();
        }
    }

    explicit impl(const overlap_reconstructor_config& c) : cfg(c)
    {
        parallel_enabled = env_flag_enabled_default_on("BONDED_PARALLEL");
        if (const char* env = std::getenv("BONDED_THREADS")) {
            const long v = std::strtol(env, nullptr, 10);
            if (v >= 1 && v <= 16) {
                n_threads_ = static_cast<size_t>(v);
            }
        }
        if (!parallel_enabled) {
            n_threads_ = 1;
        }
        if (c.warmup_samples > 0) {
            kWarmup = c.warmup_samples;
        }
        profile_enabled = env_flag_enabled("BONDED_PROFILE");
        if (profile_enabled) {
            last_profile_log = clock::now();
        }

        const size_t g =
            std::gcd(static_cast<size_t>(std::llround(c.target_rate)),
                static_cast<size_t>(std::llround(c.per_radio_rate)));
        interp = static_cast<size_t>(std::llround(c.target_rate)) / g;
        decim  = static_cast<size_t>(std::llround(c.per_radio_rate)) / g;
        rs1.init(interp, decim);
        rs2.init(interp, decim);

        const double fs        = c.target_rate;
        const double half_ov   = 0.5 * c.overlap_bw;
        const double half_band = 0.5 * c.target_bw;
        // Complementary complex crossfade filters (target baseband). The
        // transition band equals the overlap, so the FIR length scales with
        // target/overlap: a wide overlap -> short, cheap filters; a narrow
        // overlap -> a sharper, longer (costlier) filter. Bounded for safety.
        const double frac = c.overlap_bw / c.target_rate; // normalized transition
        auto odd = [](long v) { return static_cast<size_t>(v | 1L); };
        const size_t L = odd(std::clamp<long>(
            std::lround(8.0 / std::max(frac, 1e-6)), 17, 257));
        const std::vector<cf> blow = design_complex_fir(L, fs, [=](double f) -> cd {
            if (std::abs(f) > half_band) return cd(0, 0);
            if (f <= -half_ov) return cd(1, 0);
            if (f >= half_ov)  return cd(0, 0);
            const double x = (f + half_ov) / (2.0 * half_ov);
            const double c = std::cos(0.5 * kPi * x);
            return cd(c * c, 0);
        });
        const std::vector<cf> bhigh = design_complex_fir(L, fs, [=](double f) -> cd {
            if (std::abs(f) > half_band) return cd(0, 0);
            if (f <= -half_ov) return cd(0, 0);
            if (f >= half_ov)  return cd(1, 0);
            const double x = (f + half_ov) / (2.0 * half_ov);
            const double s = std::sin(0.5 * kPi * x);
            return cd(s * s, 0);
        });
        hlow.init(blow);
        hhigh.init(bhigh);
        // Overlap-band lowpass (|f| < half_ov) for the alignment estimate.
        const size_t Lo = odd(std::clamp<long>(
            std::lround(6.0 / std::max(frac, 1e-6)), 17, 195));
        ov1.init(design_lowpass_real(Lo, half_ov / fs, 1.0));
        ov2.init(design_lowpass_real(Lo, half_ov / fs, 1.0));

        // Persistent FIR history for the unified crossfade loop (K-1 zeros).
        tail1_.assign(hlow.K - 1, cf(0, 0));
        tail2_.assign(hhigh.K - 1, cf(0, 0));

        const double w1 = 2.0 * kPi * c.radio1_offset_hz / fs;
        const double w2 = 2.0 * kPi * c.radio2_offset_hz / fs;
        inc1 = cf(static_cast<float>(std::cos(w1)), static_cast<float>(std::sin(w1)));
        inc2 = cf(static_cast<float>(std::cos(w2)), static_cast<float>(std::sin(w2)));

        // --- Fused per-radio prototype (frozen/steady-state path) -------------
        // rotate-then-FIR identity: (h * (r·e^{jwn}))[n] = e^{jwn}(h'*r)[n] with
        // h'[k]=h[k]e^{-jwk}. So fuse the resampler anti-alias P with the
        // (modulated) crossfade h' into ONE complex polyphase filter per radio;
        // the residual output rotation e^{jwn} is applied cheaply afterwards.
        // By the noble identity, the fused prototype on the upsampled grid is
        //   g = conv(P, upsample(h', decim)).
        const size_t maxLM = std::max(interp, decim);
        const std::vector<float> P = design_lowpass_real(
            12 * maxLM + 1, 0.5 / static_cast<double>(maxLM),
            static_cast<double>(interp)); // == poly_resampler's internal proto
        auto build_fused = [&](const std::vector<cf>& b, double w) {
            // modulated crossfade taps h'[k] = b[k] e^{-j w k}
            std::vector<cf> hp(b.size());
            for (size_t k = 0; k < b.size(); k++) {
                const double a = -w * static_cast<double>(k);
                hp[k] = b[k] * cf(static_cast<float>(std::cos(a)),
                                  static_cast<float>(std::sin(a)));
            }
            // upsample by decim (insert decim-1 zeros between taps)
            std::vector<cf> hup((b.size() - 1) * decim + 1, cf(0, 0));
            for (size_t k = 0; k < b.size(); k++) {
                hup[k * decim] = hp[k];
            }
            // g = conv(P, hup)
            std::vector<cf> g(P.size() + hup.size() - 1, cf(0, 0));
            for (size_t u = 0; u < P.size(); u++) {
                const float pu = P[u];
                for (size_t v = 0; v < hup.size(); v++) {
                    g[u + v] += pu * hup[v];
                }
            }
            return g;
        };
        rs1c.init(interp, decim, build_fused(blow, w1));
        rs2c.init(interp, decim, build_fused(bhigh, w2));
    }

    void maybe_log_profile()
    {
        if (!profile_enabled) {
            return;
        }

        const auto now = clock::now();
        if (now - last_profile_log < std::chrono::seconds(1) || profile_stats.calls == 0) {
            return;
        }

        const auto total_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(
            profile_stats.total);
        if (total_ns.count() <= 0) {
            last_profile_log = now;
            return;
        }

        const auto to_ms = [](clock::duration duration) {
            return std::chrono::duration<double, std::milli>(duration).count();
        };
        const double total_ms = to_ms(profile_stats.total);
        const double avg_ms = total_ms / static_cast<double>(profile_stats.calls);
        const double in_msps = static_cast<double>(profile_stats.input_samples) / 1e6;
        const double out_msps = static_cast<double>(profile_stats.output_samples) / 1e6;
        const double window_s = std::chrono::duration<double>(now - last_profile_log).count();
        const auto pct = [&](clock::duration duration) {
            return 100.0 * to_ms(duration) / std::max(total_ms, 1e-9);
        };

        std::cerr << "[bonded_profile] recon calls=" << profile_stats.calls
                  << " parallel_calls=" << profile_stats.parallel_calls
                  << " avg_ms=" << avg_ms
                  << " in_msps=" << (in_msps / window_s)
                  << " out_msps=" << (out_msps / window_s)
                  << " stage_pct={resample:" << pct(profile_stats.resample)
                  << ", rotate:" << pct(profile_stats.rotate)
                  << ", align:" << pct(profile_stats.align)
                  << ", crossfade:" << pct(profile_stats.crossfade)
                  << "}" << std::endl;

        profile_stats = {};
        last_profile_log = now;
    }

    // x[i] *= start * inc^i (frequency translation). The phase at any index is
    // closed form, so each slice gets its own start phase (see inc_pow) and
    // the slices can rotate concurrently.
    static void rotate_tile(cf* x, size_t n, cf start, cf inc)
    {
#ifdef BONDED_HAVE_VOLK
        cf p = start;
        volk_32fc_s32fc_x2_rotator2_32fc(reinterpret_cast<lv_32fc_t*>(x),
            reinterpret_cast<const lv_32fc_t*>(x),
            reinterpret_cast<const lv_32fc_t*>(&inc),
            reinterpret_cast<lv_32fc_t*>(&p), static_cast<unsigned>(n));
#else
        cf p = start;
        for (size_t i = 0; i < n; i++) {
            x[i] *= p;
            p *= inc;
            const float mg = std::abs(p);
            if (mg > 0) p /= mg;
        }
#endif
    }

    //! inc^k for unit-modulus inc, in closed form: (cos(k*w), sin(k*w)).
    static cf inc_pow(cf inc, size_t k)
    {
        const double w = std::atan2(static_cast<double>(inc.imag()),
            static_cast<double>(inc.real()));
        return cf(static_cast<float>(std::cos(w * static_cast<double>(k))),
            static_cast<float>(std::sin(w * static_cast<double>(k))));
    }

    //! phase *= inc^N, renormalized.
    static void advance_phase(cf& phase, cf inc, size_t N)
    {
        phase = phase * inc_pow(inc, N);
        const float mg = std::abs(phase);
        if (mg > 0) phase /= mg;
    }

    void rotate(std::vector<cf>& x, cf& phase, cf inc)
    {
        const size_t N = x.size();
        if (N == 0) {
            return;
        }
        rotate_tile(x.data(), N, phase, inc);
        advance_phase(phase, inc, N);
    }

    const std::vector<cf>& process(const cf* rx1, const cf* rx2, size_t n)
    {
        const auto total_start = profile_enabled ? clock::now() : clock::time_point{};
        const bool use_parallel = n >= 2048;
        direct_written_ = 0; // only the frozen path can write directly

        if (frozen) {
            // Steady-state fast path: ONE fused complex polyphase filter per
            // radio (resample + translate + crossfade band-shape in a single
            // pass), then alignment scale + sum, cut into slices of the output
            // range that run concurrently (see run_slice).
            const auto resample_start_f =
                profile_enabled ? clock::now() : clock::time_point{};
            job_.p1 = rs1c.prepare(rx1, n);
            job_.p2 = rs2c.prepare(rx2, n);
            job_.mf = std::min(job_.p1, job_.p2);
            a_.resize(job_.p1);
            b_.resize(job_.p2);
            // Combine straight into the caller's buffer when one was supplied
            // (process_into), so the GR thread does not copy every output
            // sample a second time. Anything past `cap` spills to out_.
            job_.direct =
                (direct_dst_ != nullptr) ? std::min(job_.mf, direct_cap_) : 0;
            direct_written_ = job_.direct;
            out_.resize(job_.mf - job_.direct);
            job_.rotA  = rot1;
            job_.rotB  = rot2;
            job_.alpha = alpha;
            // Only hand off when the call is big enough to amortize the
            // handoff; tiny calls stay serial.
            job_.slices = use_parallel ? n_threads_ : 1;
            if (job_.slices > 1) {
                while (workers_.size() < job_.slices - 1) {
                    const size_t k = workers_.size() + 1;
                    workers_.emplace_back([this, k] { worker_loop(k); });
                }
                done_.store(0, std::memory_order_relaxed);
                job_seq_.fetch_add(1, std::memory_order_release);
                run_slice(0);
                int spins = 0;
                while (done_.load(std::memory_order_acquire) != job_.slices - 1) {
                    backoff(spins++);
                }
            } else {
                run_slice(0);
            }
            rs1c.commit();
            rs2c.commit();
            advance_phase(rot1, inc1, job_.p1);
            advance_phase(rot2, inc2, job_.p2);
            if (profile_enabled) {
                profile_stats.resample += clock::now() - resample_start_f;
                profile_stats.total += clock::now() - total_start;
                profile_stats.calls += 1;
                // Count calls actually dispatched to the workers, not merely
                // ones large enough to be eligible.
                profile_stats.parallel_calls += job_.slices > 1 ? 1 : 0;
                profile_stats.input_samples += n;
                profile_stats.output_samples += direct_written_ + out_.size();
                maybe_log_profile();
            }
            return out_;
        }

        const auto resample_start = profile_enabled ? clock::now() : clock::time_point{};
        r1_.clear();
        r2_.clear();
        rs1.feed(rx1, n, r1_);
        rs2.feed(rx2, n, r2_);
        // Keep the fused filters' input history and output grid in step with
        // the warmup resamplers (no outputs computed), so that after the
        // freeze they continue the same timeline instead of restarting: a
        // restart would change each radio's translation phase differently
        // and leave the frozen alignment wrong.
        rs1c.prepare(rx1, n);
        rs1c.commit();
        rs2c.prepare(rx2, n);
        rs2c.commit();

        const size_t m = std::min(r1_.size(), r2_.size());
        r1_.resize(m);
        r2_.resize(m);
        if (profile_enabled) {
            profile_stats.resample += clock::now() - resample_start;
        }

        // Frequency-translate each radio to the target baseband.
        const auto rotate_start = profile_enabled ? clock::now() : clock::time_point{};
        rotate(r1_, rot1, inc1);
        rotate(r2_, rot2, inc2);
        if (profile_enabled) {
            profile_stats.rotate += clock::now() - rotate_start;
        }

        // Warmup: estimate gain+phase over the shared overlap, then freeze
        // (the clock is shared, so the alignment is static).
        const auto align_start = profile_enabled ? clock::now() : clock::time_point{};
        o1_.clear();
        o2_.clear();
        ov1.filter(r1_.data(), m, o1_);
        ov2.filter(r2_.data(), m, o2_);
        for (size_t i = 0; i < m; i++) {
            Saa += std::norm(o1_[i]);
            Sbb += std::norm(o2_[i]);
            Sxy += static_cast<cd>(o1_[i]) * std::conj(static_cast<cd>(o2_[i]));
        }
        if (Sbb > 0.0) {
            const double gain = std::sqrt(Saa / Sbb);
            const double ph   = std::arg(Sxy);
            alpha = cf(static_cast<float>(gain * std::cos(ph)),
                       static_cast<float>(gain * std::sin(ph)));
            align.gain      = gain;
            align.phase_rad = ph;
            align.valid     = true;
        }
        processed += m;
        if (processed >= kWarmup && Sbb > 0.0) {
            frozen = true;
        }
        if (profile_enabled) {
            profile_stats.align += clock::now() - align_start;
        }

        // Apply the alignment to radio 2 and crossfade the two slices into the
        // contiguous target channel: each output is a pair of dot products over
        // [tail | rotated block]. The alignment scale is folded into the radio-2
        // history copy (FIR linearity), so there is no separate scaling pass.
        const auto crossfade_start = profile_enabled ? clock::now() : clock::time_point{};
        const size_t K = hlow.K; // == hhigh.K (same design length)
        ch1_.resize(K - 1 + m);
        ch2_.resize(K - 1 + m);
        std::copy(tail1_.begin(), tail1_.end(), ch1_.begin());
        std::copy(tail2_.begin(), tail2_.end(), ch2_.begin());
        std::copy(r1_.begin(), r1_.begin() + static_cast<std::ptrdiff_t>(m),
            ch1_.begin() + static_cast<std::ptrdiff_t>(K - 1));
        cf* c2 = ch2_.data() + (K - 1);
        for (size_t i = 0; i < m; i++) {
            c2[i] = r2_[i] * alpha;
        }
        out_.resize(m);
        const cf* h1  = ch1_.data();
        const cf* h2  = ch2_.data();
        const cf* tlo = hlow.taps.data();
        const cf* thi = hhigh.taps.data();
        cf* o         = out_.data();
        for (size_t i = 0; i < m; i++) {
            o[i] = dot_cc_hot(h1 + i, tlo, K) + dot_cc_hot(h2 + i, thi, K);
        }
        tail1_.assign(ch1_.end() - static_cast<std::ptrdiff_t>(K - 1), ch1_.end());
        tail2_.assign(ch2_.end() - static_cast<std::ptrdiff_t>(K - 1), ch2_.end());
        if (profile_enabled) {
            profile_stats.crossfade += clock::now() - crossfade_start;
        }

        if (profile_enabled) {
            profile_stats.total += clock::now() - total_start;
            profile_stats.calls += 1;
            profile_stats.input_samples += n;
            profile_stats.output_samples += out_.size();
            maybe_log_profile();
        }
        return out_;
    }
};

overlap_reconstructor::overlap_reconstructor(const overlap_reconstructor_config& cfg)
    : _cfg(cfg)
{
    if (cfg.target_rate <= 0 || cfg.per_radio_rate <= 0) {
        throw std::invalid_argument("overlap_reconstructor: sample rates must be > 0");
    }
    if (cfg.overlap_bw <= 0 || cfg.target_bw <= 0) {
        throw std::invalid_argument(
            "overlap_reconstructor: target_bw and overlap_bw must be > 0");
    }
}

std::vector<std::complex<float>> overlap_reconstructor::reconstruct(
    const std::vector<std::complex<float>>& rx1,
    const std::vector<std::complex<float>>& rx2)
{
    reset();
    const size_t n = std::min(rx1.size(), rx2.size());
    auto out = _impl->process(rx1.data(), rx2.data(), n);
    _alignment = _impl->align;
    return out;
}

std::vector<std::complex<float>> overlap_reconstructor::process(
    const std::complex<float>* rx1, const std::complex<float>* rx2, size_t n)
{
    if (!_impl) {
        _impl = std::make_shared<impl>(_cfg);
    }
    auto out = _impl->process(rx1, rx2, n);
    _alignment = _impl->align;
    return out;
}

size_t overlap_reconstructor::process_into(const std::complex<float>* rx1,
    const std::complex<float>* rx2,
    size_t n,
    std::complex<float>* dst,
    size_t cap)
{
    if (!_impl) {
        _impl = std::make_shared<impl>(_cfg);
    }
    _impl->direct_dst_ = dst;
    _impl->direct_cap_ = cap;
    _impl->process(rx1, rx2, n);
    _impl->direct_dst_ = nullptr;
    _impl->direct_cap_ = 0;
    _alignment = _impl->align;
    return _impl->direct_written_;
}

const std::vector<std::complex<float>>& overlap_reconstructor::pending() const
{
    static const std::vector<std::complex<float>> kEmpty;
    return _impl ? _impl->out_ : kEmpty;
}

const std::vector<std::complex<float>>& overlap_reconstructor::process_view(
    const std::complex<float>* rx1,
    const std::complex<float>* rx2,
    size_t n)
{
    if (!_impl) {
        _impl = std::make_shared<impl>(_cfg);
    }
    const auto& out = _impl->process(rx1, rx2, n);
    _alignment = _impl->align;
    return out;
}

void overlap_reconstructor::reset()
{
    _impl = std::make_shared<impl>(_cfg);
}

} // namespace bonded
