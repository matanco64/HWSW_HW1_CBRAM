// ============================================================
// Stage 3 — DBM, SoA + time-skewing + VECTORIZED stencil (single core)
// Changes from stage 2 (dbm_stage2.cpp):
//
//   Stage 2 made the kernel cache-resident (time-skewing) but turned
//   it compute-bound: it added ~40% instructions and ran 1.16× slower
//   than plain SoA. The remaining question is purely a *compute* one —
//   can we issue that arithmetic faster?
//
//   The hot 4-point Jacobi stencil never auto-vectorized in stage 1/2:
//   accessed through std::vector::operator[] with no aliasing guarantee,
//   GCC bailed ("complicated access pattern") and emitted SCALAR code,
//   leaving ~8× of AVX2 int32 throughput on the table. Stage 3 fixes
//   exactly that, with three changes confined to advance_tile():
//
//     1. __restrict raw pointers into the scratch buffers so GCC
//        auto-vectorizes the stencil to 32-byte (8×int32) AVX2 vectors.
//     2. The per-sweep metal re-pin is FUSED into the stencil store
//        (branchless select), removing a whole extra pass over the tile.
//     3. The tile copy-in/out loops are __restrict'd so they vectorize.
//
//   This pays off most *here*, not in stage 1: stage 2's tiles are
//   L2-resident, so the stencil is throughput-bound and SIMD converts
//   almost directly to wall-clock; a vectorized stage 1 would still
//   stall on DRAM latency. Cache-residency (time-skewing) and compute
//   throughput (SIMD) are complementary — this is what lets time-skewing
//   finally beat plain SoA on a single core.
//
//   All arithmetic is unchanged fixed-point int32 add + arithmetic
//   shift; the SIMD path is exact and order-independent, so the result
//   is byte-for-byte identical to stage 0/1/2 →
//   V_final_stage3.bin == V_final_stage0.bin.
// ============================================================

#include "physics_dbm.h"
#include "io.h"
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <sys/stat.h>
#include <vector>

// Blocking parameters. (TILE + 2*T_BLOCK)^2 int32 * 2 buffers should fit the
// P-core's 2 MB private L2: (256+20)^2 * 4 * 2 ≈ 0.6 MB. Blocking for L1 (the
// first attempt: TILE=64, T_BLOCK=4) loses badly — the tiny tiles make the
// halo + copy-in/out overhead ~2x the useful stencil work, and the data was
// already L2-resident anyway. Both are tunable knobs for the report's sweep.
static const int TILE    = 256;
static const int T_BLOCK = 10;

// Fixed conditions on the grid edges: Dirichlet anode/cathode rows, Neumann
// (insulating) side walls. O(N). Used to refresh a chunk's edges between blocks.
static void apply_edges(std::vector<int32_t>& V, int N) {
    for (int c = 0; c < N; ++c) V[c] = V_APPLIED;                   // anode (row 0)
    for (int c = 0; c < N; ++c) V[(N - 1) * N + c] = 0;             // cathode (row N-1)
    for (int r = 0; r < N; ++r) {                                   // insulating side walls
        V[r * N]             = V[r * N + 1];
        V[r * N + (N - 1)]   = V[r * N + (N - 2)];
    }
}

// Uniform draw in [0, n) from the deterministic PRNG, full 128-bit range.
static __int128 rand_below(Rng& rng, __int128 n) {
    unsigned __int128 r = ((unsigned __int128)rng.next() << 64) | rng.next();
    return (__int128)(r % (unsigned __int128)n);
}

// ── Temporal-blocking core ───────────────────────────────────────────────────

// Advance one TILE-owned region [ro0,ro1)×[co0,co1) by T Jacobi sweeps, reading
// the chunk-start field from `src` and writing the time+T owned interior into
// `dst`. A T-cell halo (clamped to the grid) is loaded into L1-resident scratch;
// the valid region shrinks one cell per side each sweep so no stale data leaks in.
static void advance_tile(const std::vector<int32_t>& src, std::vector<int32_t>& dst,
                         const std::vector<uint8_t>& metal, int N, int T,
                         int ro0, int ro1, int co0, int co1) {
    static std::vector<int32_t> a, b;   // L2-resident scratch, reused across tiles/calls
    static std::vector<uint8_t> m;      // tile-local metal copy (fixed for all T sweeps)

    const int R0 = std::max(0, ro0 - T), R1 = std::min(N, ro1 + T);
    const int C0 = std::max(0, co0 - T), C1 = std::min(N, co1 + T);
    const int LH = R1 - R0, LW = C1 - C0;
    a.resize((size_t)LH * LW);
    b.resize((size_t)LH * LW);
    m.resize((size_t)LH * LW);

    // Load the chunk-start state into scratch `a`. Buffer `b` only needs the
    // Dirichlet anode/cathode rows: every other cell a sub-sweep reads was
    // rewritten by the previous sub-sweep (the valid region shrinks one cell
    // per side per sweep, and the Neumann walls are refreshed every sub-sweep),
    // but rows 0 and N-1 are never recomputed and must be valid in BOTH buffers.
    // metal[] is fixed for the whole growth step, so one local copy per chunk
    // replaces T strided reads of the global array (a full-grid DRAM stream
    // per sub-sweep otherwise — the dominant remaining traffic at large N).
    {
        const int32_t* __restrict srcp = src.data();
        const uint8_t* __restrict metp = metal.data();
        int32_t* __restrict ap = a.data();
        uint8_t* __restrict mp = m.data();
        for (int r = R0; r < R1; ++r) {
            size_t go = (size_t)r * N, lo = (size_t)(r - R0) * LW;
            for (int c = C0; c < C1; ++c) {
                ap[lo + (c - C0)] = srcp[go + c];
                mp[lo + (c - C0)] = metp[go + c];
            }
        }
    }
    if (R0 == 0)
        for (int c = C0; c < C1; ++c)
            b[c - C0] = src[c];
    if (R1 == N)
        for (int c = C0; c < C1; ++c)
            b[(size_t)(N - 1 - R0) * LW + (c - C0)] = src[(size_t)(N - 1) * N + c];

    std::vector<int32_t>* cur = &a;
    std::vector<int32_t>* nxt = &b;
    for (int s = 1; s <= T; ++s) {
        const int rr0 = std::max(1, ro0 - (T - s)), rr1 = std::min(N - 1, ro1 + (T - s));
        const int cc0 = std::max(1, co0 - (T - s)), cc1 = std::min(N - 1, co1 + (T - s));
        // Raw __restrict pointers: tell GCC the buffers don't alias so it
        // auto-vectorizes the stencil to 32-byte AVX2 (8×int32). Reached
        // through std::vector::operator[] it bailed ("complicated access
        // pattern") and emitted scalar code — stage 1/2 ran the stencil scalar.
        const int32_t* __restrict cv = cur->data();
        int32_t*       __restrict nv = nxt->data();
        const uint8_t* __restrict mv = m.data();

        // 4-point Jacobi stencil with the metal re-pin FUSED into the store
        // (branchless select). The stencil and re-pin both cover cols
        // [cc0,cc1) (>=1 .. <=N-2) and the re-pin ran last in stage 2, so
        // folding it here is bit-identical — and removes a whole pass.
        for (int r = rr0; r < rr1; ++r) {
            size_t base = (size_t)(r - R0) * LW - C0;
            // omp simd forces the AVX2 vectorization GCC's -O2 cost model
            // otherwise vetoes ("not profitable"). The buffers are __restrict
            // and the stencil reads `cv` / writes `nv` (distinct), so there is
            // no loop-carried dependence — the assertion is sound.
            #pragma omp simd
            for (int c = cc0; c < cc1; ++c) {
                size_t li = base + c;
                int32_t avg = (cv[li - LW] + cv[li + LW] + cv[li - 1] + cv[li + 1]) >> 2;
                nv[li] = mv[li] ? 0 : avg;
            }
        }
        // Neumann side walls, only if this tile's halo reaches a global wall.
        // These touch cols 0 / N-1 only — disjoint from the stencil/re-pin
        // range above — so they stay exactly as in stage 2.
        if (C0 == 0)
            for (int r = rr0; r < rr1; ++r) {
                size_t row = (size_t)(r - R0) * LW - C0;
                nv[row + 0] = nv[row + 1];                 // left wall = its neighbour
            }
        if (C1 == N)
            for (int r = rr0; r < rr1; ++r) {
                size_t row = (size_t)(r - R0) * LW - C0;
                nv[row + (N - 1)] = nv[row + (N - 2)];     // right wall = its neighbour
            }
        std::swap(cur, nxt);
    }

    // Write the owned interior (time+T) back to the global destination buffer.
    const int32_t* __restrict resp = cur->data();
    int32_t* __restrict dstp = dst.data();
    for (int r = ro0; r < ro1; ++r) {
        size_t go = (size_t)r * N, lo = (size_t)(r - R0) * LW;
        for (int c = co0; c < co1; ++c)
            dstp[go + c] = resp[lo + (c - C0)];
    }
}

// Advance the whole grid by T sweeps, tiling the interior so each tile + halo
// stays hot in L1 across all T sweeps.
static void block_advance(const std::vector<int32_t>& src, std::vector<int32_t>& dst,
                          const std::vector<uint8_t>& metal, int N, int T) {
    for (int ro0 = 1; ro0 < N - 1; ro0 += TILE) {
        int ro1 = std::min(ro0 + TILE, N - 1);
        for (int co0 = 1; co0 < N - 1; co0 += TILE) {
            int co1 = std::min(co0 + TILE, N - 1);
            advance_tile(src, dst, metal, N, T, ro0, ro1, co0, co1);
        }
    }
}

// Phase 1: warm-started Jacobi solve, time-blocked. JACOBI_ITERS sweeps done in
// chunks of T_BLOCK; edges refreshed on the destination between chunks so the
// next chunk loads a consistent full state (and so the final V matches stage 0).
static void solve_potential(std::vector<int32_t>& V, std::vector<int32_t>& V_next,
                            const std::vector<uint8_t>& metal, int N) {
    std::vector<int32_t>* src = &V;
    std::vector<int32_t>* dst = &V_next;
    int done = 0;
    while (done < JACOBI_ITERS) {
        int T = std::min(T_BLOCK, JACOBI_ITERS - done);
        block_advance(*src, *dst, metal, N, T);
        apply_edges(*dst, N);
        std::swap(src, dst);
        done += T;
    }
    if (src != &V)               // result landed in V_next; move it into V
        V.swap(V_next);
}

// ── Progress bar ─────────────────────────────────────────────────────────────

static void print_bar(int step, int tip_row, int N, double elapsed) {
    int span   = N - 3;
    int done   = N - 2 - tip_row;
    if (done < 0) done = 0;
    double frac = span > 0 ? (double)done / span : 0.0;
    if (frac > 1.0) frac = 1.0;
    int filled  = (int)(frac * 30);
    double eta  = (frac > 0.005 && elapsed > 0) ? elapsed * (1.0 - frac) / frac : 0.0;
    fprintf(stderr, "\r  [");
    for (int i = 0; i < 30; ++i) fputc(i < filled ? '#' : ' ', stderr);
    fprintf(stderr, "] row %-4d  step %-7d  %.0f step/s  ETA %.0fs   ",
            tip_row, step, step / (elapsed + 1e-9), eta);
    fflush(stderr);
}

// ── Per-step phases (unchanged from stage 1) ────────────────────────────────

static void collect_candidates(const std::vector<uint8_t>& metal, int N,
                               std::vector<int>& cand) {
    cand.clear();
    for (int r = 1; r < N - 1; ++r) {
        for (int c = 1; c < N - 1; ++c) {
            int idx = r * N + c;
            if (metal[idx]) continue;
            if (metal[idx - N] || metal[idx + N] ||
                metal[idx - 1] || metal[idx + 1])
                cand.push_back(idx);
        }
    }
}

static int pick_candidate(const std::vector<int32_t>& V,
                          const std::vector<int>& cand, Rng& rng) {
    __int128 total = 0;
    for (int idx : cand) {
        int64_t v = V[idx] > 0 ? V[idx] : 0;
        total += (__int128)v * v * v;
    }
    int chosen = cand.back();
    if (total > 0) {
        __int128 thr = rand_below(rng, total), acc = 0;
        for (int idx : cand) {
            int64_t v = V[idx] > 0 ? V[idx] : 0;
            acc += (__int128)v * v * v;
            if (thr < acc) { chosen = idx; break; }
        }
    }
    return chosen;
}

static void snapshot_fields(const std::vector<int32_t>& V,
                            const std::vector<uint8_t>& metal, int N,
                            std::vector<int32_t>& V_buf,
                            std::vector<int32_t>& sigma_buf) {
    size_t nn = (size_t)N * N;
    for (size_t i = 0; i < nn; ++i) {
        V_buf[i]     = V[i];
        sigma_buf[i] = metal[i] ? SIGMA_MAX : SIGMA_LOW;
    }
}

static void write_frame(int step, int N,
                        const std::vector<int32_t>& V_buf,
                        const std::vector<int32_t>& sigma_buf) {
    char path[256];
    snprintf(path, sizeof(path), "frames_stage3/frame_%06d_V.bin", step);
    dump_binary(path, V_buf.data(), N);
    snprintf(path, sizeof(path), "frames_stage3/frame_%06d_S.bin", step);
    dump_binary(path, sigma_buf.data(), N);
}

int main(int argc, char* argv[]) {
    int N = DEFAULT_N;
    bool verbose = false;
    bool dump_frames = false;   // opt-in: -f enables per-frame dumps (for the video)
    int max_steps = 0;          // -s K caps growth steps; 0 = run until bridged.
    int frame_int = FRAME_INTERVAL;  // overridable with -F K
    for (int i = 1; i < argc; ++i) {
        if (argv[i][0] == '-') {
            if (argv[i][1] == 'v') verbose = true;
            else if (argv[i][1] == 'f') dump_frames = true;
            else if (argv[i][1] == 'n') dump_frames = false;   // legacy no-op (off is the default)
            else if (argv[i][1] == 's' && i + 1 < argc) max_steps = atoi(argv[++i]);
            else if (argv[i][1] == 'F' && i + 1 < argc) frame_int = atoi(argv[++i]);
        }
        else N = atoi(argv[i]);
    }
    size_t nn = (size_t)N * N;

    std::vector<int32_t> V(nn), V_next(nn);
    std::vector<uint8_t> metal(nn, 0);
    Rng rng(PRNG_SEED);

    // Init: linear potential ramp V_APPLIED (anode) → 0 (cathode); a single
    // metallic seed at the cathode surface.
    for (int r = 0; r < N; ++r) {
        int32_t v = (int32_t)((int64_t)V_APPLIED * (N - 1 - r) / (N - 1));
        for (int c = 0; c < N; ++c) V[r * N + c] = v;
    }
    metal[(N - 2) * N + (N / 2)] = 1;

    mkdir("frames_stage3", 0755);
    FILE* ts_fp = nullptr;
    if (dump_frames) {
        ts_fp = fopen("frames_stage3/frame_timestamps.csv", "w");
        if (ts_fp) fprintf(ts_fp, "step,elapsed_ms\n");
    }
    std::vector<int32_t> sigma_buf(nn), V_buf(nn);
    std::vector<int>     cand;

    bool bridged = false;
    int  step    = 0;
    int  limit   = (max_steps > 0 && max_steps < (int)nn) ? max_steps : (int)nn;
    auto t_start = std::chrono::steady_clock::now();
    for (step = 0; step < limit; ++step) {
        solve_potential(V, V_next, metal, N);      // phase 1: time-blocked Jacobi

        collect_candidates(metal, N, cand);
        if (cand.empty()) break;

        int chosen = pick_candidate(V, cand, rng);
        metal[chosen] = 1;

        if (verbose && step % frame_int == 0) {
            double elapsed = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - t_start).count();
            print_bar(step, chosen / N, N, elapsed);
        }

        if (dump_frames && step % frame_int == 0) {
            snapshot_fields(V, metal, N, V_buf, sigma_buf);
            write_frame(step, N, V_buf, sigma_buf);
            if (ts_fp) {
                double ms = std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - t_start).count() * 1000.0;
                fprintf(ts_fp, "%d,%.1f\n", step, ms);
                fflush(ts_fp);
            }
        }

        if (chosen / N <= 1) { bridged = true; break; }
    }

    if (verbose)
        fprintf(stderr, "\n  %s at step %d\n",
                bridged ? "Bridged" : "Stopped", step);

    snapshot_fields(V, metal, N, V_buf, sigma_buf);
    if (dump_frames) {
        write_frame(step, N, V_buf, sigma_buf);       // closing frame for the video
        if (ts_fp) {
            double ms = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - t_start).count() * 1000.0;
            fprintf(ts_fp, "%d,%.1f\n", step, ms);
            fclose(ts_fp);
        }
    }
    dump_binary("V_final_stage3.bin", V_buf.data(), N);
    dump_binary("sigma_final_stage3.bin", sigma_buf.data(), N);

    printf("Stage 3 (SoA + time-skewing + SIMD) done. N=%d, steps=%d, bridged=%d, Jacobi_iters=%d\n",
           N, step, (int)bridged, JACOBI_ITERS);
    return 0;
}
