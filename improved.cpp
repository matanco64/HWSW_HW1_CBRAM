// ============================================================
// improved.cpp — CBRAM physics-improvement series
// Baseline (Step 0): identical physics to naive.cpp.
// Each fix step is applied cumulatively in-place.
// Current step: 7 — Joule heating feedback (all 7 fixes active)
// Output: frames_improved/, sigma_final_improved.bin
// ============================================================

#include "physics_improved.h"
#include "io.h"
#include <cstdio>
#include <sys/stat.h>
#include <vector>

struct Cell {
    int32_t V, sigma, ion, temp;
};

static void apply_boundary(std::vector<Cell>& grid, int N) {
    for (int i = 0; i < N; ++i) grid[i].V = V_APPLIED;
    for (int i = 0; i < N; ++i) grid[(N - 1) * N + i].V = 0;
    for (int r = 0; r < N; ++r) {
        grid[r * N].V         = grid[r * N + 1].V;
        grid[r * N + (N-1)].V = grid[r * N + (N-2)].V;
    }
}

static void jacobi_sweep(const std::vector<Cell>& grid,
                         std::vector<Cell>& next, int N)
{
    for (int r = 1; r < N - 1; ++r) {
        for (int c = 1; c < N - 1; ++c) {
            int idx = r * N + c;

            int32_t s_c = grid[idx].sigma;
            int32_t s_e = (s_c + grid[idx + 1].sigma) >> 1;
            int32_t s_w = (s_c + grid[idx - 1].sigma) >> 1;
            int32_t s_n = (s_c + grid[idx + N].sigma) >> 1;
            int32_t s_s = (s_c + grid[idx - N].sigma) >> 1;

            int32_t denom = s_e + s_w + s_n + s_s;
            if (denom == 0) { next[idx].V = 0; continue; }

            int64_t num = (int64_t)s_e * grid[idx + 1].V
                        + (int64_t)s_w * grid[idx - 1].V
                        + (int64_t)s_n * grid[idx + N].V
                        + (int64_t)s_s * grid[idx - N].V;

            next[idx].V = (int32_t)(num / denom);
        }
    }
}

static void drift_diffusion(std::vector<Cell>& grid, int N) {
    for (int r = N - 2; r >= 1; --r) {
        for (int c = 1; c < N - 1; ++c) {
            int src = r * N + c;
            int dst = (r + 1) * N + c;
            int32_t E_down = grid[src].V - grid[dst].V;
            if (E_down <= 0) continue;
            // Step 5: sinh field-enhanced hopping; Step 7: Arrhenius temperature-scaled mobility
            int32_t arg   = q_mul(SINH_COEFF, E_down);
            int32_t dT    = grid[src].temp - TEMP_INIT;
            int32_t mob_t = q_mul(ION_MOBILITY, q_exp_approx(q_mul(EA_OVER_K, dT)));
            int32_t flux  = q_mul(q_mul(mob_t, q_sinh(arg)), grid[src].ion);
            int32_t avail = grid[src].ion - ION_LOW;
            if (flux > avail) flux = avail;
            if (flux <= 0) continue;
            grid[src].ion -= flux;
            if (r + 1 < N - 1)
                grid[dst].ion = q_clamp(grid[dst].ion + flux, ION_LOW, ION_HIGH * 4);
        }
    }
    // Step 1: lateral diffusion — J_x = -D_LAT · ∂²C/∂x²
    for (int r = 1; r < N - 1; ++r) {
        for (int c = 1; c < N - 1; ++c) {
            int idx = r * N + c;
            int32_t lap = grid[idx - 1].ion - 2 * grid[idx].ion + grid[idx + 1].ion;
            grid[idx].ion = q_clamp(grid[idx].ion + q_mul(D_LAT, lap), ION_LOW, ION_HIGH * 4);
        }
    }
    // Step 6: Butler-Volmer anode injection — inject ∝ exp(α·F·η/RT)
    // No sigma factor: exchange current i₀ is a material constant, not local conductivity.
    // (Step 3 moved seeds to row N-2 so row-1 sigma stays at SIGMA_LOW throughout.)
    for (int c = 1; c < N - 1; ++c) {
        int32_t eta = grid[N + c].V - V_EQ;
        if (eta <= 0) continue;
        int32_t bv  = q_exp_approx(q_mul(BV_ALPHA, eta));
        int32_t inj = q_mul(ION_INJECT, bv);
        grid[N + c].ion = q_clamp(grid[N + c].ion + inj, ION_LOW, ION_HIGH);
    }
}

// Step 4: tip-only deposition — sigma grows only adjacent to existing metal (σ > SIGMA_HIGH).
// Prevents broad electrolyte deposition; channels all growth to the advancing filament tip.
// Uses SIGMA_GROWTH_IMP (10× baseline) so the tip traverses N=200 rows within 3000 steps.
static void update_sigma(std::vector<Cell>& grid, int N, int t) {
    uint32_t t_mix = (uint32_t)t * 1664525u + 1013904223u;
    for (int r = 1; r < N - 1; ++r) {
        for (int c = 1; c < N - 1; ++c) {
            int idx = r * N + c;
            bool adj_metal =
                grid[idx - 1].sigma > SIGMA_HIGH ||
                grid[idx + 1].sigma > SIGMA_HIGH ||
                grid[idx - N].sigma > SIGMA_HIGH ||
                grid[idx + N].sigma > SIGMA_HIGH;
            if (!adj_metal) continue;
            uint32_t rng = (uint32_t)idx ^ t_mix;
            rng ^= rng << 13; rng ^= rng >> 17; rng ^= rng << 5;
            int32_t prob_q16 = q_mul(SIGMA_GROWTH_IMP, grid[idx].ion);
            // guard against overflow when prob ≥ 1.0 in Q16.16
            uint32_t threshold = (prob_q16 >= ONE) ? 0xFFFFFFFFu : ((uint32_t)prob_q16 << 16);
            if (rng < threshold)
                grid[idx].sigma = q_clamp(grid[idx].sigma + SIGMA_GROWTH_IMP, SIGMA_LOW, SIGMA_MAX);
        }
    }
}

// Step 7: update temperature field — T = T_amb + σ·|∇V|² · R_TH
static void update_temperature(std::vector<Cell>& grid, int N) {
    for (int r = 1; r < N - 1; ++r) {
        for (int c = 1; c < N - 1; ++c) {
            int idx = r * N + c;
            int32_t dV_r = grid[idx + N].V - grid[idx - N].V;
            int32_t dV_c = grid[idx + 1].V - grid[idx - 1].V;
            int32_t E2   = q_mul(dV_r, dV_r) + q_mul(dV_c, dV_c);
            int32_t P    = q_mul(grid[idx].sigma, E2);
            int32_t dT   = q_mul(P, R_TH);
            grid[idx].temp = q_clamp(TEMP_INIT + dT, TEMP_INIT, TEMP_MAX);
        }
    }
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

    {
        std::vector<int32_t> V(nn), sigma(nn), ion(nn), temp(nn);
        init_fields(N, V.data(), sigma.data(), ion.data(), temp.data());
        for (size_t i = 0; i < nn; ++i)
            grid[i] = {V[i], sigma[i], ion[i], temp[i]};
    }

    // Step 3: move seeds from anode row 1 to cathode row N-2.
    // Clear the seeds init_fields placed at row 1, then replant at row N-2
    // using the same RNG sequence so the same columns are chosen.
    {
        for (int c = 1; c < N - 1; ++c) {
            grid[N + c].sigma = SIGMA_LOW;
            grid[N + c].ion   = ION_LOW;
        }
        std::mt19937 seed_rng(PRNG_SEED);
        std::uniform_int_distribution<int> col_dist(1, N - 2);
        std::uniform_real_distribution<float> str_dist(0.2f, 1.0f);
        for (int s = 0; s < NUM_SEEDS; ++s) {
            int   c = col_dist(seed_rng);
            (void)str_dist(seed_rng);  // consume same RNG calls as init_fields
            int idx = (N - 2) * N + c;
            // Step 4: seeds fully metallic so tip-only condition fires immediately
            grid[idx].sigma = SIGMA_MAX;
            grid[idx].ion   = ION_HIGH;
        }
        // Step 4: pre-fill electrolyte with ions — physically, dissolved metal is present before SET
        for (int r = 1; r < N - 1; ++r)
            for (int c = 1; c < N - 1; ++c)
                if (grid[r * N + c].sigma < SIGMA_HIGH)
                    grid[r * N + c].ion = ION_HIGH;
    }

    mkdir("frames_improved", 0755);

    std::vector<int32_t> sigma_buf(nn);

    char path[256];
    for (int t = 0; t < TOTAL_TIMESTEPS; ++t) {
        if (verbose) {
            fprintf(stderr, "\r  timestep %d/%d", t + 1, TOTAL_TIMESTEPS);
            fflush(stderr);
        }
        for (int k = 0; k < JACOBI_ITERS; ++k) {
            grid_next = grid;
            jacobi_sweep(grid, grid_next, N);
            apply_boundary(grid_next, N);
            std::swap(grid, grid_next);
        }

        update_temperature(grid, N);  // Step 7: Joule heating
        drift_diffusion(grid, N);
        update_sigma(grid, N, t);

        if (t % FRAME_INTERVAL == 0) {
            for (size_t i = 0; i < nn; ++i) sigma_buf[i] = grid[i].sigma;
            snprintf(path, sizeof(path), "frames_improved/frame_%04d.ppm", t);
            write_ppm(path, sigma_buf.data(), N, SIGMA_MAX);
        }

        // Step 4: SET complete when filament reaches anode (row 1), not cathode
        bool bridged = false;
        for (int c = 1; c < N - 1; ++c)
            if (grid[N + c].sigma >= SIGMA_MAX) { bridged = true; break; }
        if (bridged) {
            if (verbose) fprintf(stderr, "\n  Filament bridged at t=%d — SET complete\n", t + 1);
            break;
        }
    }

    if (verbose) fprintf(stderr, "\n");

    for (size_t i = 0; i < nn; ++i) sigma_buf[i] = grid[i].sigma;
    dump_binary("sigma_final_improved.bin", sigma_buf.data(), N);

    printf("Improved done. Step=7 (Joule heating + BV injection). N=%d, timesteps=%d, Jacobi_iters=%d\n",
           N, TOTAL_TIMESTEPS, JACOBI_ITERS);
    return 0;
}
