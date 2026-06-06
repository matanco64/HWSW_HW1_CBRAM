// ============================================================
// Stage 0 — DBM, AoS Naive
// Changes from previous: (none — this is the baseline)
//
// Dielectric Breakdown Model. A metallic filament grows one cell
// at a time; each growth step warm-starts JACOBI_ITERS Jacobi
// sweeps of Laplace, then adds the empty cluster-adjacent cell
// chosen with probability ∝ V^ETA.
//
// Memory layout: struct Cell { int32_t V, metal; } — 8 bytes.
// A 64-byte cache line holds 8 cells. The Jacobi hot loop reads
// only V from each neighbor, so half of every line fetched is the
// unused `metal` field. Cache-line utilization ~50%. This is the
// bottleneck stage 1 (SoA) targets.
// ============================================================

#include "physics_dbm.h"
#include "io.h"
#include <cstdio>
#include <cstdlib>
#include <sys/stat.h>
#include <vector>

struct Cell {
    int32_t V;
    int32_t metal;   // 1 = part of the filament, 0 = electrolyte
};

// Re-impose the boundary conditions after a Jacobi sweep.
static void apply_boundary(std::vector<Cell>& grid, int N) {
    for (int c = 0; c < N; ++c) grid[c].V = V_APPLIED;              // anode (row 0)
    for (int c = 0; c < N; ++c) grid[(N - 1) * N + c].V = 0;        // cathode (row N-1)
    for (int r = 0; r < N; ++r) {                                   // insulating side walls
        grid[r * N].V         = grid[r * N + 1].V;
        grid[r * N + (N - 1)].V = grid[r * N + (N - 2)].V;
    }
    for (size_t i = 0; i < grid.size(); ++i)                        // cluster ≡ cathode (V=0)
        if (grid[i].metal) grid[i].V = 0;
}

// One Jacobi sweep of Laplace: V_new = (V_N + V_E + V_S + V_W) / 4.
// Uniform medium, so a plain 4-point average; integer-exact via >> 2.
static void jacobi_sweep(const std::vector<Cell>& grid,
                         std::vector<Cell>& next, int N) {
    for (int r = 1; r < N - 1; ++r) {
        for (int c = 1; c < N - 1; ++c) {
            int idx = r * N + c;
            next[idx].V = (grid[idx - N].V + grid[idx + N].V +
                           grid[idx - 1].V + grid[idx + 1].V) >> 2;
        }
    }
}

// Uniform draw in [0, n) from the deterministic PRNG, full 128-bit range.
static __int128 rand_below(Rng& rng, __int128 n) {
    unsigned __int128 r = ((unsigned __int128)rng.next() << 64) | rng.next();
    return (__int128)(r % (unsigned __int128)n);
}

int main(int argc, char* argv[]) {
    int N = DEFAULT_N;
    bool verbose = false;
    for (int i = 1; i < argc; ++i) {
        if (argv[i][0] == '-') { if (argv[i][1] == 'v') verbose = true; }
        else N = atoi(argv[i]);
    }
    size_t nn = (size_t)N * N;

    std::vector<Cell> grid(nn), grid_next(nn);
    Rng rng(PRNG_SEED);

    // Init: linear potential ramp V_APPLIED (anode) → 0 (cathode); a single
    // metallic seed at the cathode surface.
    for (int r = 0; r < N; ++r) {
        int32_t v = (int32_t)((int64_t)V_APPLIED * (N - 1 - r) / (N - 1));
        for (int c = 0; c < N; ++c) { grid[r * N + c].V = v; grid[r * N + c].metal = 0; }
    }
    grid[(N - 2) * N + (N / 2)].metal = 1;

    mkdir("frames_stage0", 0755);
    std::vector<int32_t> sigma_buf(nn), V_buf(nn);
    std::vector<int>     cand;              // candidate cell indices (row-major order)
    char path[256];

    bool bridged = false;
    int  step    = 0;
    for (step = 0; step < (int)nn; ++step) {
        // ── Hot loop: warm-started Jacobi solve of the potential ──
        for (int k = 0; k < JACOBI_ITERS; ++k) {
            grid_next = grid;
            jacobi_sweep(grid, grid_next, N);
            apply_boundary(grid_next, N);
            std::swap(grid, grid_next);
        }

        // ── Candidate sites: empty interior cells 4-adjacent to the cluster ──
        cand.clear();
        for (int r = 1; r < N - 1; ++r) {
            for (int c = 1; c < N - 1; ++c) {
                int idx = r * N + c;
                if (grid[idx].metal) continue;
                if (grid[idx - N].metal || grid[idx + N].metal ||
                    grid[idx - 1].metal || grid[idx + 1].metal)
                    cand.push_back(idx);
            }
        }
        if (cand.empty()) break;

        // ── DBM rule: pick one candidate with probability ∝ V^ETA (=V^3) ──
        // Full-precision integer weight: a Q16.16 cube would underflow to 0 for
        // the small V near the cathode and stall growth.
        __int128 total = 0;
        for (int idx : cand) {
            int64_t v = grid[idx].V > 0 ? grid[idx].V : 0;
            total += (__int128)v * v * v;
        }
        int chosen = cand.back();
        if (total > 0) {
            __int128 thr = rand_below(rng, total), acc = 0;
            for (int idx : cand) {
                int64_t v = grid[idx].V > 0 ? grid[idx].V : 0;
                acc += (__int128)v * v * v;
                if (thr < acc) { chosen = idx; break; }
            }
        }
        grid[chosen].metal = 1;

        if (verbose && step % FRAME_INTERVAL == 0)
            fprintf(stderr, "\r  step %d  cells %d  tip_row %d   ",
                    step, step + 1, chosen / N);

        if (step % FRAME_INTERVAL == 0) {
            for (size_t i = 0; i < nn; ++i)
                sigma_buf[i] = grid[i].metal ? SIGMA_MAX : SIGMA_LOW;
            snprintf(path, sizeof(path), "frames_stage0/frame_%04d.ppm", step);
            write_ppm(path, sigma_buf.data(), N, SIGMA_MAX);
        }

        if (chosen / N <= 1) { bridged = true; break; }   // reached the anode
    }

    if (verbose)
        fprintf(stderr, "\n  %s at step %d\n",
                bridged ? "Bridged" : "Stopped", step);

    for (size_t i = 0; i < nn; ++i) {
        V_buf[i]     = grid[i].V;
        sigma_buf[i] = grid[i].metal ? SIGMA_MAX : SIGMA_LOW;
    }
    dump_binary("V_final_stage0.bin", V_buf.data(), N);
    dump_binary("sigma_final_stage0.bin", sigma_buf.data(), N);

    printf("Stage 0 (DBM) done. N=%d, steps=%d, bridged=%d, Jacobi_iters=%d\n",
           N, step, (int)bridged, JACOBI_ITERS);
    return 0;
}
