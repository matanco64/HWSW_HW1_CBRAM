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

// Fixed conditions on the grid edges: Dirichlet anode/cathode rows, Neumann
// (insulating) side walls. O(N) — touches only the border, negligible cost.
static void apply_edges(std::vector<Cell>& grid, int N) {
    for (int c = 0; c < N; ++c) grid[c].V = V_APPLIED;              // anode (row 0)
    for (int c = 0; c < N; ++c) grid[(N - 1) * N + c].V = 0;        // cathode (row N-1)
    for (int r = 0; r < N; ++r) {                                   // insulating side walls
        grid[r * N].V         = grid[r * N + 1].V;
        grid[r * N + (N - 1)].V = grid[r * N + (N - 2)].V;
    }
}

// Internal Dirichlet condition: pin every metallic (filament) cell to the
// cathode potential V=0. O(N²) — a full-grid scan re-run on EVERY Jacobi
// iteration, and in this AoS layout each check pulls a whole 8-byte Cell to
// read one flag. This is what dominates apply_boundary in the profile.
static void pin_cluster(std::vector<Cell>& grid) {
    for (size_t i = 0; i < grid.size(); ++i)
        if (grid[i].metal) grid[i].V = 0;
}

// Re-impose the boundary conditions after a Jacobi sweep.
static void apply_boundary(std::vector<Cell>& grid, int N) {
    apply_edges(grid, N);
    pin_cluster(grid);
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

// ── Per-step phases (split out so perf can charge time to each) ──────────────

// Seed the write buffer from the read buffer before a sweep. In this AoS layout
// the WHOLE Cell (V *and* metal) is copied, even though Jacobi only writes V —
// the metal copy is pure AoS tax that SoA eliminates (separate metal[] array,
// never copied during the solve). Shows up as a fat memmove in the profile.
static void copy_grid(const std::vector<Cell>& src, std::vector<Cell>& dst) {
    dst = src;
}

// Phase 1: warm-started Jacobi solve of the Laplace potential for the current
// cluster. JACOBI_ITERS double-buffered sweeps; never reset between steps.
static void solve_potential(std::vector<Cell>& grid,
                            std::vector<Cell>& grid_next, int N) {
    for (int k = 0; k < JACOBI_ITERS; ++k) {
        copy_grid(grid, grid_next);
        jacobi_sweep(grid, grid_next, N);
        apply_boundary(grid_next, N);
        std::swap(grid, grid_next);
    }
}

// Phase 2: row-major scan for empty interior cells with ≥1 metal 4-neighbor.
// O(N²) every step — a hidden cost that grows with the grid.
static void collect_candidates(const std::vector<Cell>& grid, int N,
                               std::vector<int>& cand) {
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
}

// Phase 3: DBM rule — choose a candidate with probability ∝ V^ETA (=V^3).
// Full-precision __int128 weight: a Q16.16 cube would underflow to 0 for the
// small V near the cathode and stall growth.
static int pick_candidate(const std::vector<Cell>& grid,
                          const std::vector<int>& cand, Rng& rng) {
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
    return chosen;
}

// Phase 4a: copy the current V and σ fields into flat output buffers.
static void snapshot_fields(const std::vector<Cell>& grid, int N,
                            std::vector<int32_t>& V_buf,
                            std::vector<int32_t>& sigma_buf) {
    size_t nn = (size_t)N * N;
    for (size_t i = 0; i < nn; ++i) {
        V_buf[i]     = grid[i].V;
        sigma_buf[i] = grid[i].metal ? SIGMA_MAX : SIGMA_LOW;
    }
}

// Phase 4b: write one animation frame (V and σ) to disk. Pure I/O — skipped
// entirely under -n so the compute profile isn't polluted by it.
static void write_frame(int step, int N,
                        const std::vector<int32_t>& V_buf,
                        const std::vector<int32_t>& sigma_buf) {
    char path[256];
    snprintf(path, sizeof(path), "frames_stage0/frame_%06d_V.bin", step);
    dump_binary(path, V_buf.data(), N);
    snprintf(path, sizeof(path), "frames_stage0/frame_%06d_S.bin", step);
    dump_binary(path, sigma_buf.data(), N);
}

int main(int argc, char* argv[]) {
    int N = DEFAULT_N;
    bool verbose = false;
    bool dump_frames = false;   // opt-in: -f enables per-frame dumps (for the video)
    int max_steps = 0;          // -s K caps growth steps; 0 = run until bridged.
                                // Lets large-N profiling runs do a fixed, identical
                                // amount of work per stage instead of a full bridge.
    for (int i = 1; i < argc; ++i) {
        if (argv[i][0] == '-') {
            if (argv[i][1] == 'v') verbose = true;
            else if (argv[i][1] == 'f') dump_frames = true;
            else if (argv[i][1] == 'n') dump_frames = false;   // legacy no-op (off is the default)
            else if (argv[i][1] == 's' && i + 1 < argc) max_steps = atoi(argv[++i]);
        }
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

    bool bridged = false;
    int  step    = 0;
    int  limit   = (max_steps > 0 && max_steps < (int)nn) ? max_steps : (int)nn;
    for (step = 0; step < limit; ++step) {
        solve_potential(grid, grid_next, N);       // phase 1: Jacobi (hot loop)

        collect_candidates(grid, N, cand);         // phase 2: O(N²) frontier scan
        if (cand.empty()) break;

        int chosen = pick_candidate(grid, cand, rng);   // phase 3: V³ weighted pick
        grid[chosen].metal = 1;

        if (verbose && step % FRAME_INTERVAL == 0)
            fprintf(stderr, "\r  step %d  cells %d  tip_row %d   ",
                    step, step + 1, chosen / N);

        if (dump_frames && step % FRAME_INTERVAL == 0) {   // phase 4: output
            snapshot_fields(grid, N, V_buf, sigma_buf);
            write_frame(step, N, V_buf, sigma_buf);
        }

        if (chosen / N <= 1) { bridged = true; break; }   // reached the anode
    }

    if (verbose)
        fprintf(stderr, "\n  %s at step %d\n",
                bridged ? "Bridged" : "Stopped", step);

    // Final state → buffers, used for both the closing frame and the dumps.
    snapshot_fields(grid, N, V_buf, sigma_buf);
    if (dump_frames)                                  // closing frame for the video
        write_frame(step, N, V_buf, sigma_buf);
    dump_binary("V_final_stage0.bin", V_buf.data(), N);
    dump_binary("sigma_final_stage0.bin", sigma_buf.data(), N);

    printf("Stage 0 (DBM) done. N=%d, steps=%d, bridged=%d, Jacobi_iters=%d\n",
           N, step, (int)bridged, JACOBI_ITERS);
    return 0;
}
