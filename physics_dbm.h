#pragma once
#include <cstdint>

// ============================================================
// physics_dbm.h — shared constants + helpers for the DBM stages.
//
// Dielectric Breakdown Model: a metallic filament grows one cell
// at a time. Each growth step solves Laplace for the potential V
// (the cluster pinned to the cathode potential, the anode held at
// V_APPLIED), then adds one empty cell adjacent to the cluster,
// chosen with probability ∝ V^ETA.
//
// Everything is Q16.16 fixed-point and integer arithmetic so that
// every optimization stage produces a BIT-IDENTICAL V_final — the
// correctness gate the HW optimization story rests on.
// ============================================================

// Q16.16 fixed-point: integer part in upper 16 bits, fraction in lower 16.
static constexpr int     SHIFT = 16;
static constexpr int32_t ONE   = 1 << SHIFT;        // 1.0

// Grid / solver parameters (N overridable via CLI).
static constexpr int      DEFAULT_N    = 200;
static constexpr int      JACOBI_ITERS = 30;        // warm-started sweeps per growth step
static constexpr int      ETA          = 3;         // DBM exponent (integer → exact cube)
static constexpr uint64_t PRNG_SEED    = 42;

// Physical constants in Q16.16.
static constexpr int32_t V_APPLIED = 1 * ONE;       // anode potential (1.0 V)
static constexpr int32_t SIGMA_LOW = ONE / 100;     // electrolyte conductivity (0.01)
static constexpr int32_t SIGMA_MAX = 20 * ONE;      // metallic conductivity (20.0)

// Output cadence.
static constexpr int FRAME_INTERVAL = 2;           // dump a PPM every N growth steps

// Multiply two Q16.16 values, result in Q16.16.
inline int32_t q_mul(int32_t a, int32_t b) {
    return (int32_t)(((int64_t)a * b) >> SHIFT);
}

// Clamp a Q16.16 value to [lo, hi].
inline int32_t q_clamp(int32_t x, int32_t lo, int32_t hi) {
    return x < lo ? lo : (x > hi ? hi : x);
}

// Deterministic 64-bit PRNG (xorshift64*). A fixed sequence of draws plus the
// bit-identical V field gives an identical filament across every stage.
struct Rng {
    uint64_t s;
    explicit Rng(uint64_t seed) : s(seed ? seed : 0x9E3779B97F4A7C15ULL) {}
    uint64_t next() {
        s ^= s >> 12; s ^= s << 25; s ^= s >> 27;
        return s * 0x2545F4914F6CDD1DULL;
    }
    // Uniform integer in [0, n) for n > 0.
    uint64_t below(uint64_t n) { return next() % n; }
};
