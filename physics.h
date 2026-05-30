#pragma once
#include <cstdint>
#include <random>

// Q16.16 fixed-point: integer part in upper 16 bits, fractional in lower 16
static constexpr int SHIFT = 16;
static constexpr int32_t ONE = 1 << SHIFT;  // 1.0 in Q16.16

// Grid parameters (can be overridden via CLI)
static constexpr int DEFAULT_N         = 4096;
static constexpr int JACOBI_ITERS      = 50;   // K: Jacobi iterations per timestep
static constexpr int TOTAL_TIMESTEPS   = 200;  // outer simulation steps
static constexpr int FRAME_INTERVAL    = 2;    // dump PPM every N timesteps
static constexpr uint32_t PRNG_SEED    = 42;

// Physical parameters in Q16.16
static constexpr int32_t V_APPLIED     = 1 * ONE;       // top electrode voltage (1.0 V)
static constexpr int32_t SIGMA_LOW     = ONE / 100;     // background conductivity (0.01)
static constexpr int32_t SIGMA_HIGH    = 10 * ONE;      // seed defect conductivity (10.0)
static constexpr int32_t ION_LOW       = ONE / 1000;    // background ion concentration
static constexpr int32_t ION_HIGH      = ONE;           // seed ion concentration (1.0)
static constexpr int32_t ION_DRIFT     = ONE / 200;     // drift-diffusion step size
static constexpr int32_t SIGMA_GROWTH  = ONE / 50;      // conductivity growth rate
static constexpr int32_t SIGMA_MAX     = 20 * ONE;      // conductivity ceiling
static constexpr int32_t TEMP_INIT     = 1 * ONE;       // uniform initial temperature

static constexpr int NUM_SEEDS = 8; // number of defect seed sites

// Clamp a Q16.16 value to [lo, hi]
inline int32_t q_clamp(int32_t x, int32_t lo, int32_t hi) {
    return x < lo ? lo : (x > hi ? hi : x);
}

// Multiply two Q16.16 values, result in Q16.16
inline int32_t q_mul(int32_t a, int32_t b) {
    return (int32_t)(((int64_t)a * b) >> SHIFT);
}

// Shared init: fills V[], sigma[], ion[], temp[] from a fixed seed.
// Works for both AoS (pass field pointers from struct) and SoA flat arrays.
inline void init_fields(int N,
                        int32_t* V, int32_t* sigma, int32_t* ion, int32_t* temp)
{
    // Zero / uniform init
    for (int i = 0; i < N * N; ++i) {
        V[i]     = 0;
        sigma[i] = SIGMA_LOW;
        ion[i]   = ION_LOW;
        temp[i]  = TEMP_INIT;
    }

    // Top electrode: V = V_APPLIED (Dirichlet)
    for (int i = 0; i < N; ++i)
        V[i] = V_APPLIED;

    // Place defect seeds with fixed PRNG
    std::mt19937 rng(PRNG_SEED);
    std::uniform_int_distribution<int> row_dist(1, N - 2); // avoid electrodes
    std::uniform_int_distribution<int> col_dist(1, N - 2);

    for (int s = 0; s < NUM_SEEDS; ++s) {
        int r = row_dist(rng);
        int c = col_dist(rng);
        int idx = r * N + c;
        sigma[idx] = SIGMA_HIGH;
        ion[idx]   = ION_HIGH;
    }
}
