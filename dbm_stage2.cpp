// ============================================================
// Stage 2 — DBM, SoA + cache blocking / time-skewing
// Changes from stage 1 (dbm_stage1.cpp):
//
//   The Jacobi solve does JACOBI_ITERS (=30) sweeps over the whole
//   grid per growth step. In stage 1 each sweep streams both V
//   buffers (N*N int32 each) through the cache from scratch, so the
//   working set is re-read from L2/DRAM 30 times — bandwidth bound
//   once N*N*4 exceeds the cache.
//
//   Stage 2 applies *temporal blocking* (a.k.a. time skewing) via
//   OVERLAPPED TILING: the 30 sweeps are processed in chunks of
//   T_BLOCK. For each chunk, the grid is cut into TILE×TILE tiles;
//   each tile is loaded together with a T-cell halo into a small
//   scratch pair that fits the P-core's private L2, marched forward
//   T sweeps (the valid region shrinks one cell per side each sweep),
//   then the owned interior is written back. A tile's whole T-step
//   history stays hot in L2 — the big arrays are streamed from
//   L3/DRAM ~30/T_BLOCK times instead of 30, at the cost of
//   recomputing the halos. This only pays off once the working set
//   (2 * N² * 4 B) exceeds the 24 MB L3, i.e. N ≳ 1800; below that
//   stage 1 never touches DRAM and the halo work is pure overhead.
//
//   Tiles are independent (all read the chunk-start state, write
//   disjoint outputs), which also makes stage 3 (OpenMP) trivial.
//
//   The per-sweep boundary conditions are reproduced inside each
//   tile in global coordinates: Dirichlet anode/cathode rows are
//   constant (loaded once), Neumann side walls are refreshed each
//   sub-sweep when the tile touches a wall, and metal cells are
//   re-pinned to 0 each sub-sweep. Result: byte-for-byte identical
//   to stage 0/1 → V_final_stage2.bin == V_final_stage0.bin.
// ============================================================

#include "physics_dbm.h"
#include "io.h"
#include <algorithm>
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
    for (int r = R0; r < R1; ++r)
        for (int c = C0; c < C1; ++c) {
            size_t li = (size_t)(r - R0) * LW + (c - C0);
            a[li] = src[(size_t)r * N + c];
            m[li] = metal[(size_t)r * N + c];
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
        std::vector<int32_t>& cv = *cur;
        std::vector<int32_t>& nv = *nxt;

        // 4-point Jacobi stencil over the (shrinking) interior region.
        for (int r = rr0; r < rr1; ++r) {
            size_t base = (size_t)(r - R0) * LW - C0;
            for (int c = cc0; c < cc1; ++c) {
                size_t li = base + c;
                nv[li] = (cv[li - LW] + cv[li + LW] + cv[li - 1] + cv[li + 1]) >> 2;
            }
        }
        // Neumann side walls, only if this tile's halo reaches a global wall.
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
        // Internal Dirichlet: re-pin metal cells to V=0 (overrides the stencil).
        // Branchless select so the compiler vectorizes it like stage 1's
        // pin_cluster, instead of a per-cell compare-and-branch.
        for (int r = rr0; r < rr1; ++r) {
            size_t base = (size_t)(r - R0) * LW - C0;
            for (int c = cc0; c < cc1; ++c)
                nv[base + c] = m[base + c] ? 0 : nv[base + c];
        }
        std::swap(cur, nxt);
    }

    // Write the owned interior (time+T) back to the global destination buffer.
    std::vector<int32_t>& res = *cur;
    for (int r = ro0; r < ro1; ++r)
        for (int c = co0; c < co1; ++c)
            dst[(size_t)r * N + c] = res[(size_t)(r - R0) * LW + (c - C0)];
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
    snprintf(path, sizeof(path), "frames_stage2/frame_%06d_V.bin", step);
    dump_binary(path, V_buf.data(), N);
    snprintf(path, sizeof(path), "frames_stage2/frame_%06d_S.bin", step);
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

    mkdir("frames_stage2", 0755);
    std::vector<int32_t> sigma_buf(nn), V_buf(nn);
    std::vector<int>     cand;

    bool bridged = false;
    int  step    = 0;
    int  limit   = (max_steps > 0 && max_steps < (int)nn) ? max_steps : (int)nn;
    for (step = 0; step < limit; ++step) {
        solve_potential(V, V_next, metal, N);      // phase 1: time-blocked Jacobi

        collect_candidates(metal, N, cand);
        if (cand.empty()) break;

        int chosen = pick_candidate(V, cand, rng);
        metal[chosen] = 1;

        if (verbose && step % FRAME_INTERVAL == 0)
            fprintf(stderr, "\r  step %d  cells %d  tip_row %d   ",
                    step, step + 1, chosen / N);

        if (dump_frames && step % FRAME_INTERVAL == 0) {
            snapshot_fields(V, metal, N, V_buf, sigma_buf);
            write_frame(step, N, V_buf, sigma_buf);
        }

        if (chosen / N <= 1) { bridged = true; break; }
    }

    if (verbose)
        fprintf(stderr, "\n  %s at step %d\n",
                bridged ? "Bridged" : "Stopped", step);

    snapshot_fields(V, metal, N, V_buf, sigma_buf);
    if (dump_frames)
        write_frame(step, N, V_buf, sigma_buf);
    dump_binary("V_final_stage2.bin", V_buf.data(), N);
    dump_binary("sigma_final_stage2.bin", sigma_buf.data(), N);

    printf("Stage 2 (SoA + time-skewing) done. N=%d, steps=%d, bridged=%d, Jacobi_iters=%d\n",
           N, step, (int)bridged, JACOBI_ITERS);
    return 0;
}
