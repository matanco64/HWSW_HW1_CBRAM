// ============================================================
// Stage 0 — AoS Naive
// Changes from previous: (none — this is the baseline)
//
// Memory layout: struct Cell { int32_t V, sigma, ion, temp; }
// One Cell per grid point, 16 bytes. A 64-byte cache line holds
// 4 cells. The Poisson hot loop reads only V and sigma from each
// neighbor, so only 8 of the 16 bytes fetched per line are used.
// Cache-line utilization: ~50%. This is the bottleneck the plan
// targets with the SoA transformation in stage 1.
// ============================================================

#include "physics.h"
#include "io.h"
#include <cstdio>
#include <sys/stat.h>
#include <vector>

struct Cell {
    int32_t V, sigma, ion, temp;
};

static void apply_boundary(std::vector<Cell>& grid, int N) {
    // Dirichlet: top = V_APPLIED, bottom = 0
    for (int i = 0; i < N; ++i) grid[i].V = V_APPLIED;
    for (int i = 0; i < N; ++i) grid[(N - 1) * N + i].V = 0;
    // Neumann: insulating left/right walls (dV/dx = 0)
    for (int r = 0; r < N; ++r) {
        grid[r * N].V         = grid[r * N + 1].V;
        grid[r * N + (N-1)].V = grid[r * N + (N-2)].V;
    }
}

// Jacobi iteration for Poisson equation with variable conductivity.
// V_new[i,j] = (s_e*V[i,j+1] + s_w*V[i,j-1] + s_n*V[i+1,j] + s_s*V[i-1,j])
//              / (s_e + s_w + s_n + s_s)
// All arithmetic in Q16.16; accumulator in Q32.32 (int64_t).
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

// Drift-diffusion: positive metal ions drift downward (anode row 0 → cathode row N-1).
// Upwind explicit Euler swept bottom-to-top so destination row is not revisited
// this step — avoids double-transporting the same ions in a single sweep.
static void drift_diffusion(std::vector<Cell>& grid, int N) {
    for (int r = N - 2; r >= 1; --r) {
        for (int c = 1; c < N - 1; ++c) {
            int src = r * N + c;
            int dst = (r + 1) * N + c;
            // J-driven drift: flux = µ * σ * E * ion — localizes to high-conductivity paths
            int32_t E_down = grid[src].V - grid[dst].V;
            if (E_down <= 0) continue;
            int32_t flux = q_mul(q_mul(q_mul(ION_MOBILITY, grid[src].sigma), E_down), grid[src].ion);
            int32_t avail = grid[src].ion - ION_LOW;
            if (flux > avail) flux = avail;
            if (flux <= 0) continue;
            grid[src].ion -= flux;
            if (r + 1 < N - 1)
                grid[dst].ion = q_clamp(grid[dst].ion + flux, ION_LOW, ION_HIGH * 4);
        }
    }
    // Active anode (row 1): electrodissolution proportional to local conductivity
    // — supplies ions where current density is highest, localizing the filament
    for (int c = 1; c < N - 1; ++c) {
        int32_t inject = q_mul(ION_INJECT, grid[N + c].sigma);
        grid[N + c].ion = q_clamp(grid[N + c].ion + inject, ION_LOW, ION_HIGH);
    }
}

// Conductivity update: sigma grows with ion concentration.
static void update_sigma(std::vector<Cell>& grid, int N) {
    for (int r = 1; r < N - 1; ++r) {
        for (int c = 1; c < N - 1; ++c) {
            int idx = r * N + c;
            int32_t growth = q_mul(SIGMA_GROWTH, grid[idx].ion);
            grid[idx].sigma = q_clamp(grid[idx].sigma + growth, SIGMA_LOW, SIGMA_MAX);
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

    // Init via shared routine then scatter into AoS
    {
        std::vector<int32_t> V(nn), sigma(nn), ion(nn), temp(nn);
        init_fields(N, V.data(), sigma.data(), ion.data(), temp.data());
        for (size_t i = 0; i < nn; ++i)
            grid[i] = {V[i], sigma[i], ion[i], temp[i]};
    }

    mkdir("frames_stage0", 0755);

    // Pre-allocate scratch buffers for frame/dump extraction
    std::vector<int32_t> sigma_buf(nn), V_buf(nn);

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

        drift_diffusion(grid, N);
        update_sigma(grid, N);

        if (t % FRAME_INTERVAL == 0) {
            for (size_t i = 0; i < nn; ++i) sigma_buf[i] = grid[i].sigma;
            snprintf(path, sizeof(path), "frames_stage0/frame_%04d.ppm", t);
            write_ppm(path, sigma_buf.data(), N, SIGMA_MAX);
        }

        // Stop when any filament bridges to the cathode (SET complete)
        bool bridged = false;
        for (int c = 1; c < N - 1; ++c)
            if (grid[(N - 2) * N + c].sigma >= SIGMA_MAX) { bridged = true; break; }
        if (bridged) {
            if (verbose) fprintf(stderr, "\n  Filament bridged at t=%d — SET complete\n", t + 1);
            break;
        }
    }

    if (verbose) fprintf(stderr, "\n");
    for (size_t i = 0; i < nn; ++i) V_buf[i] = grid[i].V;
    dump_binary("V_final_stage0.bin", V_buf.data(), N);

    printf("Stage 0 done. N=%d, timesteps=%d, Jacobi_iters=%d\n",
           N, TOTAL_TIMESTEPS, JACOBI_ITERS);
    return 0;
}
