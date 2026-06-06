#pragma once
#include <algorithm>
#include <cstdint>
#include <random>

// Q16.16 fixed-point: integer part in upper 16 bits, fractional in lower 16
static constexpr int SHIFT = 16;
static constexpr int32_t ONE = 1 << SHIFT;  // 1.0 in Q16.16

// Grid parameters (can be overridden via CLI)
static constexpr int DEFAULT_N         = 4096;
static constexpr int JACOBI_ITERS      = 50;
static constexpr int TOTAL_TIMESTEPS   = 3000;
static constexpr int FRAME_INTERVAL    = 2;
static constexpr uint32_t PRNG_SEED    = 42;

// Physical parameters in Q16.16 (unchanged from physics.h)
static constexpr int32_t V_APPLIED     = 1 * ONE;
static constexpr int32_t SIGMA_LOW     = ONE / 100;
static constexpr int32_t SIGMA_HIGH    = 10 * ONE;
static constexpr int32_t ION_LOW       = ONE / 1000;
static constexpr int32_t ION_HIGH      = ONE;
static constexpr int32_t ION_MOBILITY  = 4 * ONE;
static constexpr int32_t ION_INJECT    = ONE / 10;
static constexpr int32_t SIGMA_GROWTH  = ONE / 10;
static constexpr int32_t SIGMA_MAX     = 20 * ONE;
static constexpr int32_t TEMP_INIT     = 1 * ONE;

static constexpr int NUM_SEEDS = 3;

// ── Constants added by fix steps ─────────────────────────────────────────────
// Step 1 — lateral diffusion
static constexpr int32_t D_LAT         = ONE / 20;    // lateral diffusion coeff (0.05)

// Step 4 — tip-only deposition: faster per-deposit growth so filament advances in ≤3000 steps
static constexpr int32_t SIGMA_GROWTH_IMP = ONE;          // 10× SIGMA_GROWTH; ~10 deposits to metallic

// Step 5 — sinh field-enhanced hopping
static constexpr int32_t SINH_COEFF    = 2 * ONE;     // q·a / 2kT scale factor

// Step 6 — Butler-Volmer injection
static constexpr int32_t V_EQ          = 0;           // equilibrium potential (0 V)
static constexpr int32_t BV_ALPHA      = 2 * ONE;     // α·F/RT (scaled)

// Step 7 — Joule heating
static constexpr int32_t R_TH          = 5 * ONE;     // thermal resistance (normalized)
static constexpr int32_t TEMP_MAX      = 10 * ONE;    // temperature ceiling (normalized)
static constexpr int32_t EA_OVER_K     = ONE / 4;     // Eₐ/kT scale factor (~0.25)

// ─────────────────────────────────────────────────────────────────────────────

inline int32_t q_clamp(int32_t x, int32_t lo, int32_t hi) {
    return x < lo ? lo : (x > hi ? hi : x);
}

inline int32_t q_mul(int32_t a, int32_t b) {
    return (int32_t)(((int64_t)a * b) >> SHIFT);
}

// Step 6: exp approximation in Q16.16 (Taylor series up to x^5; good for |x| ≤ 3)
inline int32_t q_exp_approx(int32_t x) {
    int32_t result = ONE + x;
    int32_t term = x;
    for (int i = 2; i <= 5; ++i) {
        term = (int32_t)((int64_t)term * x / ((int64_t)ONE * i));
        result += term;
    }
    return result < ONE ? ONE : result;
}

// Step 5: sinh approximation in Q16.16
// |x| ≤ 1.5: cubic  sinh(x) ≈ x + x³/6  (error < 0.4%)
// |x| > 1.5: Taylor e^x up to x^6, then sinh ≈ e^x / 2
inline int32_t q_sinh(int32_t x) {
    if (x < 0) return -q_sinh(-x);
    if (x > 3 * ONE / 2) {
        int32_t ex = ONE;
        int32_t term = x;
        for (int i = 1; i <= 6; ++i) {
            ex += term;
            term = (int32_t)((int64_t)term * x / ((int64_t)ONE * (i + 1)));
        }
        return ex >> 1;
    }
    int32_t x3 = (int32_t)(((int64_t)x * x >> SHIFT) * (int64_t)x >> SHIFT);
    return x + (int32_t)((int64_t)x3 / 6);
}

// Shared init: seeds are placed at ANODE side (row 1) by default.
// improved.cpp overrides the seed row via the cathode-nucleation fix (Step 3).
inline void init_fields(int N,
                        int32_t* V, int32_t* sigma, int32_t* ion, int32_t* temp)
{
    for (int i = 0; i < N * N; ++i) {
        V[i]     = 0;
        sigma[i] = SIGMA_LOW;
        ion[i]   = ION_LOW;
        temp[i]  = TEMP_INIT;
    }
    for (int i = 0; i < N; ++i)
        V[i] = V_APPLIED;

    std::mt19937 rng(PRNG_SEED);
    std::uniform_int_distribution<int> col_dist(1, N - 2);
    std::uniform_real_distribution<float> strength_dist(0.2f, 1.0f);

    for (int s = 0; s < NUM_SEEDS; ++s) {
        int c   = col_dist(rng);
        float k = strength_dist(rng);
        int idx = 1 * N + c;
        sigma[idx] = (int32_t)(SIGMA_HIGH * k);
        ion[idx]   = (int32_t)(ION_HIGH   * k);
    }
}
