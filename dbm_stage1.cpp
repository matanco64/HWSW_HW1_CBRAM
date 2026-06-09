// ============================================================
// Stage 1 — DBM, Struct-of-Arrays (SoA)
// Changes from stage 0 (dbm_stage0.cpp):
//
//   1. Layout: struct Cell{int32 V,metal}  →  two flat arrays,
//      int32_t V[] and uint8_t metal[]. The Jacobi hot loop now
//      streams a contiguous int32 V[] — 16 values per 64-byte
//      cache line (~100% utilization) instead of 8 (~50%), because
//      the interleaved `metal` field no longer rides along in cache.
//
//   2. Only V is double-buffered (V / V_next). `metal` is FIXED for
//      all JACOBI_ITERS sweeps of a growth step, so it is a single
//      shared array that is never copied. This deletes the per-sweep
//      full-grid `grid_next = grid` copy — the ~20% memmove the stage-0
//      flamegraph exposed — outright.
//
//   3. metal[] is uint8_t, so pin_cluster scans 1 byte/cell instead
//      of dragging an 8-byte Cell through cache to read one flag.
//
// Numerics, PRNG sequence and candidate order are byte-for-byte
// identical to stage 0 → V_final_stage1.bin == V_final_stage0.bin.
// ============================================================

#include "physics_dbm.h"
#include "io.h"
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <sys/stat.h>
#include <vector>

// Fixed conditions on the grid edges: Dirichlet anode/cathode rows, Neumann
// (insulating) side walls. O(N) — touches only the border, negligible cost.
static void apply_edges(std::vector<int32_t>& V, int N) {
    for (int c = 0; c < N; ++c) V[c] = V_APPLIED;                   // anode (row 0)
    for (int c = 0; c < N; ++c) V[(N - 1) * N + c] = 0;             // cathode (row N-1)
    for (int r = 0; r < N; ++r) {                                   // insulating side walls
        V[r * N]             = V[r * N + 1];
        V[r * N + (N - 1)]   = V[r * N + (N - 2)];
    }
}

// Internal Dirichlet condition: pin every metallic (filament) cell to the
// cathode potential V=0. Still O(N²) every sweep, but now a tight scan over a
// separate uint8_t metal[] — 1 byte/cell, no struct drag.
static void pin_cluster(std::vector<int32_t>& V, const std::vector<uint8_t>& metal) {
    for (size_t i = 0; i < metal.size(); ++i)
        if (metal[i]) V[i] = 0;
}

// One Jacobi sweep of Laplace: V_new = (V_N + V_E + V_S + V_W) / 4.
// Uniform medium, so a plain 4-point average; integer-exact via >> 2.
// Reads and writes contiguous int32 arrays — the SoA cache win lives here.
static void jacobi_sweep(const std::vector<int32_t>& V,
                         std::vector<int32_t>& Vn, int N) {
    for (int r = 1; r < N - 1; ++r) {
        for (int c = 1; c < N - 1; ++c) {
            int idx = r * N + c;
            Vn[idx] = (V[idx - N] + V[idx + N] +
                       V[idx - 1] + V[idx + 1]) >> 2;
        }
    }
}

// Uniform draw in [0, n) from the deterministic PRNG, full 128-bit range.
static __int128 rand_below(Rng& rng, __int128 n) {
    unsigned __int128 r = ((unsigned __int128)rng.next() << 64) | rng.next();
    return (__int128)(r % (unsigned __int128)n);
}

// ── Per-step phases ──────────────────────────────────────────────────────────

// Phase 1: warm-started Jacobi solve of the Laplace potential for the current
// cluster. No buffer copy: every cell of V_next is written each sweep — interior
// by jacobi_sweep, edges by apply_edges, metal cells by pin_cluster — and metal[]
// is read-only and shared, so there is nothing to copy.
static void solve_potential(std::vector<int32_t>& V, std::vector<int32_t>& V_next,
                            const std::vector<uint8_t>& metal, int N) {
    for (int k = 0; k < JACOBI_ITERS; ++k) {
        jacobi_sweep(V, V_next, N);
        apply_edges(V_next, N);
        pin_cluster(V_next, metal);
        std::swap(V, V_next);
    }
}

// Phase 2: row-major scan for empty interior cells with ≥1 metal 4-neighbor.
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

// Phase 3: DBM rule — choose a candidate with probability ∝ V^ETA (=V^3).
// Full-precision __int128 weight: a Q16.16 cube would underflow to 0 for the
// small V near the cathode and stall growth.
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

// Phase 4a: copy the current V and σ fields into flat output buffers.
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

// Phase 4b: write one animation frame (V and σ) to disk. Pure I/O — skipped
// entirely under -n so the compute profile isn't polluted by it.
static void write_frame(int step, int N,
                        const std::vector<int32_t>& V_buf,
                        const std::vector<int32_t>& sigma_buf) {
    char path[256];
    snprintf(path, sizeof(path), "frames_stage1/frame_%06d_V.bin", step);
    dump_binary(path, V_buf.data(), N);
    snprintf(path, sizeof(path), "frames_stage1/frame_%06d_S.bin", step);
    dump_binary(path, sigma_buf.data(), N);
}

int main(int argc, char* argv[]) {
    int N = DEFAULT_N;
    bool verbose = false;
    bool dump_frames = true;   // -n disables per-frame dumps (compute-only profiling)
    for (int i = 1; i < argc; ++i) {
        if (argv[i][0] == '-') {
            if (argv[i][1] == 'v') verbose = true;
            else if (argv[i][1] == 'n') dump_frames = false;
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

    mkdir("frames_stage1", 0755);
    std::vector<int32_t> sigma_buf(nn), V_buf(nn);
    std::vector<int>     cand;              // candidate cell indices (row-major order)

    bool bridged = false;
    int  step    = 0;
    for (step = 0; step < (int)nn; ++step) {
        solve_potential(V, V_next, metal, N);      // phase 1: Jacobi (hot loop)

        collect_candidates(metal, N, cand);        // phase 2: O(N²) frontier scan
        if (cand.empty()) break;

        int chosen = pick_candidate(V, cand, rng);      // phase 3: V³ weighted pick
        metal[chosen] = 1;

        if (verbose && step % FRAME_INTERVAL == 0)
            fprintf(stderr, "\r  step %d  cells %d  tip_row %d   ",
                    step, step + 1, chosen / N);

        if (dump_frames && step % FRAME_INTERVAL == 0) {   // phase 4: output
            snapshot_fields(V, metal, N, V_buf, sigma_buf);
            write_frame(step, N, V_buf, sigma_buf);
        }

        if (chosen / N <= 1) { bridged = true; break; }   // reached the anode
    }

    if (verbose)
        fprintf(stderr, "\n  %s at step %d\n",
                bridged ? "Bridged" : "Stopped", step);

    // Final state → buffers, used for both the closing frame and the dumps.
    snapshot_fields(V, metal, N, V_buf, sigma_buf);
    if (dump_frames)                                  // closing frame for the video
        write_frame(step, N, V_buf, sigma_buf);
    dump_binary("V_final_stage1.bin", V_buf.data(), N);
    dump_binary("sigma_final_stage1.bin", sigma_buf.data(), N);

    printf("Stage 1 (SoA) done. N=%d, steps=%d, bridged=%d, Jacobi_iters=%d\n",
           N, step, (int)bridged, JACOBI_ITERS);
    return 0;
}
