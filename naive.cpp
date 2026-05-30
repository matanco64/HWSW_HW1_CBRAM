// ============================================================
// Stage 0 — AoS Naive
// Changes from previous: (none — this is the baseline)
//
// Memory layout: struct Cell { int32_t V, sigma, ion, temp; }
// One Cell per grid point, 16 bytes. A 64-byte cache line holds
// 4 cells. The Poisson hot loop reads only V and sigma from each
// neighbor, so only 8 of the 16 floats fetched per line are used.
// Cache-line utilization: ~50%. This is the bottleneck the plan
// targets with the SoA transformation in stage 1.
// ============================================================

#include "physics.h"
#include "io.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sys/stat.h>

struct Cell {
    int32_t V, sigma, ion, temp;
};

static void apply_boundary(Cell* grid, int N) {
    // Top row: V = V_APPLIED (Dirichlet)
    for (int i = 0; i < N; ++i)
        grid[i].V = V_APPLIED;
    // Bottom row: V = 0 (Dirichlet)
    for (int i = 0; i < N; ++i)
        grid[(N - 1) * N + i].V = 0;
}

// Jacobi iteration for Poisson equation with variable conductivity.
// V_new[i,j] = (s_e*V[i,j+1] + s_w*V[i,j-1] + s_n*V[i+1,j] + s_s*V[i-1,j])
//              / (s_e + s_w + s_n + s_s)
// where s_e = arithmetic mean of sigma at i,j and i,j+1, etc.
// All arithmetic in Q16.16; accumulator in Q32.32 (int64_t).
static void jacobi_sweep(const Cell* __restrict__ grid,
                               Cell* __restrict__ next, int N)
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

// Drift-diffusion: ions move toward lower potential (down the field).
// Simple upwind explicit Euler: ion flows in direction of E = -dV/dy.
static void drift_diffusion(Cell* grid, int N) {
    for (int r = 1; r < N - 1; ++r) {
        for (int c = 1; c < N - 1; ++c) {
            int idx = r * N + c;
            int32_t dV = grid[idx - N].V - grid[idx + N].V; // E points down
            int32_t flux = q_mul(ION_DRIFT, dV > 0 ? dV : -dV);
            if (dV > 0)
                grid[idx].ion = q_clamp(grid[idx].ion + flux, 0, ION_HIGH * 2);
        }
    }
}

// Conductivity update: sigma grows logistically with ion concentration.
static void update_sigma(Cell* grid, int N) {
    for (int r = 1; r < N - 1; ++r) {
        for (int c = 1; c < N - 1; ++c) {
            int idx = r * N + c;
            int32_t growth = q_mul(SIGMA_GROWTH, grid[idx].ion);
            grid[idx].sigma = q_clamp(grid[idx].sigma + growth, SIGMA_LOW, SIGMA_MAX);
        }
    }
}

int main(int argc, char* argv[]) {
    int N = (argc > 1) ? atoi(argv[1]) : DEFAULT_N;

    Cell* grid      = (Cell*)malloc((size_t)N * N * sizeof(Cell));
    Cell* grid_next = (Cell*)malloc((size_t)N * N * sizeof(Cell));
    if (!grid || !grid_next) { fprintf(stderr, "malloc failed\n"); return 1; }

    // Init via shared routine — extract field pointers from AoS
    // We init separate arrays then scatter into the struct.
    int32_t* tmp_V     = (int32_t*)malloc((size_t)N * N * sizeof(int32_t));
    int32_t* tmp_sigma = (int32_t*)malloc((size_t)N * N * sizeof(int32_t));
    int32_t* tmp_ion   = (int32_t*)malloc((size_t)N * N * sizeof(int32_t));
    int32_t* tmp_temp  = (int32_t*)malloc((size_t)N * N * sizeof(int32_t));
    init_fields(N, tmp_V, tmp_sigma, tmp_ion, tmp_temp);
    for (int i = 0; i < N * N; ++i) {
        grid[i].V     = tmp_V[i];
        grid[i].sigma = tmp_sigma[i];
        grid[i].ion   = tmp_ion[i];
        grid[i].temp  = tmp_temp[i];
    }
    free(tmp_V); free(tmp_sigma); free(tmp_ion); free(tmp_temp);

    mkdir("frames_stage0", 0755);

    char path[256];
    for (int t = 0; t < TOTAL_TIMESTEPS; ++t) {
        // Jacobi solve: K iterations
        for (int k = 0; k < JACOBI_ITERS; ++k) {
            memcpy(grid_next, grid, (size_t)N * N * sizeof(Cell));
            jacobi_sweep(grid, grid_next, N);
            apply_boundary(grid_next, N);
            Cell* tmp = grid; grid = grid_next; grid_next = tmp;
        }

        drift_diffusion(grid, N);
        update_sigma(grid, N);

        if (t % FRAME_INTERVAL == 0) {
            // Extract sigma for PPM writer
            int32_t* sigma_buf = (int32_t*)malloc((size_t)N * N * sizeof(int32_t));
            for (int i = 0; i < N * N; ++i) sigma_buf[i] = grid[i].sigma;
            snprintf(path, sizeof(path), "frames_stage0/frame_%04d.ppm", t);
            write_ppm(path, sigma_buf, N, SIGMA_MAX);
            free(sigma_buf);
        }
    }

    // Final state dump for correctness verification
    int32_t* V_buf = (int32_t*)malloc((size_t)N * N * sizeof(int32_t));
    for (int i = 0; i < N * N; ++i) V_buf[i] = grid[i].V;
    dump_binary("V_final_stage0.bin", V_buf, N);
    free(V_buf);

    free(grid);
    free(grid_next);
    printf("Stage 0 done. N=%d, timesteps=%d, Jacobi_iters=%d\n",
           N, TOTAL_TIMESTEPS, JACOBI_ITERS);
    return 0;
}
