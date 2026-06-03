# CBRAM Naive Simulation — Deep Analysis & Fix Plan

> Reference: `cbram_memristor_guide.md`, `naive.cpp`, `physics.h`, `simulation_loop.md`

---

## 1. What the Simulation Does

The grid is N×N. Anode at **row 0** (V = 1.0), cathode at **row N−1** (V = 0). Three defect
seeds with random column positions and strengths are placed at **row 1** (anode-adjacent).

Each outer timestep:
1. K = 50 Jacobi sweeps of `∇·(σ∇V) = 0`
2. Ion drift downward: `flux = µ · σ · E · C`
3. σ grows where ions accumulate: `σ += SIGMA_GROWTH · ion`
4. Bridge check: if any σ[N−2, c] ≥ SIGMA_MAX → SET complete, stop

---

## 2. What is Correct

| Aspect | Why correct |
|---|---|
| Poisson eq. `∇·(σ∇V) = 0`, Jacobi discretization | Correct current-continuity equation with variable conductivity; weighted-average formula is the exact FD discretization |
| Ion drift direction (anode → cathode) | Positive ions drift along E field (high V → low V); `E_down = V[r,c] − V[r+1,c]`, flux skipped when E_down ≤ 0 |
| J-driven flux `σ · E · C` | Creates the positive feedback loop: high-σ paths carry more current → more ion drift → more deposition → σ grows; 1000× selectivity over background at init |
| Upwind sweep (bottom-to-top) | Each ion packet transported once per timestep; top-to-bottom would double-transport |
| Anode injection ∝ σ | Physically: electrode dissolves proportional to current density J = σ·E; ties ion supply to developing filament |
| Bridge detection at row N−2 | Correct: cell just above cathode reaches SIGMA_MAX → HRS→LRS transition |
| Q16.16 integer arithmetic | Guarantees bit-identical results across optimization stages |

---

## 3. What is Wrong

### 3.1 Seeds at the ANODE, not the CATHODE — *moderate*

**Guide (§3.1 Step 3):** "At the inert electrode (cathode), metal ions are reduced and
electrodeposited." Real CBRAM nucleates at the cathode because that is where electrons are
available for M^n+ + ne⁻ → M. The filament grows from cathode **toward** anode.

The simulation inverts this: seeds at row 1 (near anode), growth toward cathode. This changes
which tip is the "growing tip" and prevents symmetry from breaking at the cathode surface where
real stochastic nucleation begins.

---

### 3.2 No stochastic nucleation — **CRITICAL (primary cause of no lightning)**

**Guide (§3.2):** "Nucleation is stochastic at multiple sites simultaneously."
**Guide (§6.3):** "C2C variability originates from probabilistic Poisson-distributed random events
in ion hopping."

The simulation has **zero runtime randomness**. PRNG seed 42 fixes seeds deterministically;
every update is deterministic. The state at timestep T is a pure function of initial conditions.

Lightning/dendritic morphology is **Diffusion-Limited Aggregation (DLA)**. DLA produces fractal,
branching patterns only when:
- Growth probability is field-enhanced at tips (partially present), **AND**
- There is **randomness in each deposition event**

Without noise, DLA degenerates: every ion follows the steepest-field path → smooth, straight
column directly below each seed. No branching is possible, because branching requires an ion to
deviate from the highest-probability path, which requires thermal noise.

---

### 3.3 No lateral diffusion — **CRITICAL (second cause of no lightning)**

**Guide (§9.1) Nernst-Planck:**
```
J = -D · ∂C/∂x  −  (D·z·q / k_BT) · C · ∂φ/∂x
    ↑ diffusion          ↑ drift
```

The simulation implements **only** the drift term, and even that only in the **downward direction**.
No `J_x = -D · ∂C/∂x` lateral spreading exists anywhere.

Without lateral diffusion, all ions released from a seed at column `c` travel in a perfectly
straight vertical column. Result: three vertical stripes, not a tree. Lightning shape requires
ions to spread horizontally — diffusion is the only mechanism for this, as the lateral electric
field in a 1D stack geometry is zero at the center.

---

### 3.4 Linear flux instead of sinh field-enhanced hopping — *moderate*

**Guide (§9.1):** The field-enhanced hopping rate from transition state theory is:
```
ν = ν₀ · exp(−Eₐ/kT) · sinh(q·a·E / 2kT)
```
At high fields (q·a·E >> 2kT), sinh grows as exp(E) — **exponential in field**. This creates
the sharp V_SET threshold and explosive tip runaway as the gap to cathode approaches zero
(E → ∞ → ν → ∞).

The simulation uses `flux ∝ E` (linear). The growth rate difference between the filament tip
(high E) and the sides (low E) is only 2–3×. With sinh it would be orders of magnitude.

---

### 3.5 Linear anode injection instead of Butler-Volmer — *moderate*

**Guide (§9.2):** Anodic dissolution follows:
```
i = i₀ · [exp(α·F·η / RT) − exp(−(1−α)·F·η / RT)]
```
Exponential in overpotential η. Creates a threshold voltage below which essentially no ions
are injected.

The simulation uses `inject = ION_INJECT · σ` — linear, no voltage threshold. Ions are
injected at constant rate proportional to σ regardless of V_applied, so there is no V_SET.

---

### 3.6 σ grows everywhere along the ion path, not just at the tip — *moderate*

Real CBRAM: metal reduction (M^n+ + ne⁻ → M) occurs only at the **metallic surface** (the
growing tip). Ions in transit through the electrolyte remain ionic — they don't deposit until
they contact already-deposited metal.

Simulation: `sigma[r,c] += SIGMA_GROWTH · ion[r,c]` — σ grows everywhere ions drift through.
This creates a broad conductivity cone from the seed all the way to the tip, rather than a
sharp metallic wire with a clear tip.

---

### 3.7 Joule heating unused — *minor for SET, critical for RESET*

**Guide (§9.3):** `ΔT = I²R · R_th` at filament neck. Elevated T increases ion mobility via
`D(T) = D₀ · exp(−Eₐ/kT)`. This concentrates growth at the narrowest neck (highest current
density), making the filament needle-like rather than tapered.

The `temp` field is allocated in `struct Cell` but is **never updated**. Ion mobility is
constant everywhere.

---

### 3.8 Jacobi convergence severely insufficient for large N — *severe for N ≥ 512*

Jacobi spectral radius ρ = cos(π/N) ≈ 1 − π²/(2N²). Full convergence on N = 4096 requires
~7 billion iterations. K = 50 propagates information at most 50 cells per sweep — a tiny
fraction of the grid. For N = 1024, the electric field at the filament tip is barely
distinguishable from the background until the filament is within 50 cells of the cathode.

This doesn't affect the optimization stages comparison (all use K = 50), but physically
the tip-field enhancement mechanism barely functions for large N.

---

## 4. Why There Is No Lightning Shape — Summary

| Root cause | Physical consequence |
|---|---|
| **No stochastic deposition** | DLA degenerates → no branching at all |
| **No lateral diffusion** | Ions travel straight down → three vertical stripes |
| Linear flux (∝ E, not sinh) | Tip advantage too small → filament widens not elongates |
| σ grows along entire path | Broad diffuse cone instead of sharp wire |
| Seeds at anode not cathode | Symmetry not broken at cathode; wrong nucleation point |
| No Butler-Volmer | No V_SET threshold; injection constant regardless of voltage |
| No Joule heating | No thermal focusing at neck |

The two dominant causes are **(1)** and **(2)**. Without fixing these, all other improvements
are cosmetic.

---

## 5. Fix Plan

`improved.cpp` starts as a copy of `naive.cpp`. Each step below is applied **cumulatively**.
After each step: build → run (`N=200 -v`) → `/grade-filament improved` → inspect frames.

### Step 0 — Baseline (naive.cpp)
- Add `dump_binary("sigma_final_stage0.bin", ...)` to `naive.cpp` so it can be graded.
- Expected grade: **F (~10/100)** — three vertical stripes, no bridge.

---

### Step 1 — Lateral Diffusion
**File:** `improved.cpp`  **Header:** `physics_improved.h`

Add explicit-Euler lateral diffusion pass in `drift_diffusion`, after the downward sweep:

```cpp
// Lateral diffusion: J_x = -D_LAT · ∂²C/∂x²  (central differences)
// Must use a copy of ion[] to avoid order-dependence; process on grid_temp
for (int r = 1; r < N - 1; ++r) {
    for (int c = 1; c < N - 1; ++c) {
        int idx = r * N + c;
        int32_t lap = grid[idx - 1].ion - 2 * grid[idx].ion + grid[idx + 1].ion;
        grid[idx].ion = q_clamp(grid[idx].ion + q_mul(D_LAT, lap), ION_LOW, ION_HIGH * 4);
    }
}
```

New constant in `physics_improved.h`:
```cpp
static constexpr int32_t D_LAT = ONE / 20;   // lateral diffusion coefficient (0.05)
```

> Note: The diffusion pass must read from a snapshot of the `ion` field (copy before
> the loop or use two alternating buffers) to avoid propagation artifacts. Because the
> diffusion stencil is symmetric and uses central differences, single-pass on the same
> buffer introduces a mild left-to-right bias but is acceptable for this study.

**Expected grade:** ~20–30 — ions spread laterally; blobs begin to overlap; still no
branching because no randomness.

---

### Step 2 — Stochastic Deposition
**File:** `improved.cpp`

Replace deterministic `update_sigma` with probabilistic per-cell deposition. Use a fast
inline xorshift32 seeded per-timestep (no `<random>` overhead):

```cpp
// Global LCG state, seeded once; updated each timestep
static uint32_t lcg_state = 12345;

inline uint32_t lcg_next() {
    lcg_state ^= lcg_state << 13;
    lcg_state ^= lcg_state >> 17;
    lcg_state ^= lcg_state << 5;
    return lcg_state;
}

static void update_sigma(std::vector<Cell>& grid, int N) {
    for (int r = 1; r < N - 1; ++r) {
        for (int c = 1; c < N - 1; ++c) {
            int idx = r * N + c;
            if (grid[idx].ion <= ION_LOW) continue;
            // Deposit with probability ∝ SIGMA_GROWTH · ion
            // threshold = SIGMA_GROWTH * ion in [0, ONE*ONE/65536] range
            uint32_t threshold = (uint32_t)((uint64_t)(uint32_t)SIGMA_GROWTH
                                          * (uint32_t)grid[idx].ion >> SHIFT);
            if ((lcg_next() & 0xFFFF) < (threshold & 0xFFFF))
                grid[idx].sigma = q_clamp(grid[idx].sigma + SIGMA_GROWTH,
                                          SIGMA_LOW, SIGMA_MAX);
        }
    }
}
```

**Expected grade:** ~35–45 — first visible branching; stochastic fingers begin to deviate
from the vertical baseline.

---

### Step 3 — Cathode Nucleation
**File:** `improved.cpp`

Move seeds from row 1 to row **N−2** (cathode-adjacent). Keep anode injection at row 1
unchanged. Add a nucleation-boost: when `ion[N−2, c] > ION_HIGH/2`, the deposition
probability is 4× higher (models preferential nucleation at the inert cathode surface):

In `init_fields` local copy: change `int idx = 1 * N + c` → `int idx = (N-2) * N + c`.

In `update_sigma`, add after the existing loop:
```cpp
// Cathode nucleation boost (row N-2): higher deposition probability where ions accumulate
for (int c = 1; c < N - 1; ++c) {
    int idx = (N - 2) * N + c;
    if (grid[idx].ion <= ION_HIGH / 2) continue;
    uint32_t threshold = (uint32_t)((uint64_t)(uint32_t)(SIGMA_GROWTH * 4)
                                  * (uint32_t)grid[idx].ion >> SHIFT);
    if ((lcg_next() & 0xFFFF) < (threshold & 0xFFFF))
        grid[idx].sigma = q_clamp(grid[idx].sigma + SIGMA_GROWTH * 2,
                                  SIGMA_LOW, SIGMA_MAX);
}
```

**Expected grade:** ~45–55 — filament nucleates from the bottom and grows upward; initial
asymmetry is at the cathode as in real CBRAM.

---

### Step 4 — Tip-Only Deposition
**File:** `improved.cpp`

Real deposition only occurs at the metallic surface. Gate `update_sigma` so σ only grows
at cells adjacent to at least one already-metallic neighbor (σ > SIGMA_HIGH):

```cpp
static void update_sigma(std::vector<Cell>& grid, int N) {
    for (int r = 1; r < N - 1; ++r) {
        for (int c = 1; c < N - 1; ++c) {
            int idx = r * N + c;
            if (grid[idx].ion <= ION_LOW) continue;
            // Only deposit adjacent to existing metal
            bool near_metal = grid[idx - 1].sigma > SIGMA_HIGH ||
                               grid[idx + 1].sigma > SIGMA_HIGH ||
                               grid[idx - N].sigma > SIGMA_HIGH ||
                               grid[idx + N].sigma > SIGMA_HIGH;
            if (!near_metal) continue;
            // Stochastic deposition (from Step 2)
            uint32_t threshold = (uint32_t)((uint64_t)(uint32_t)SIGMA_GROWTH
                                          * (uint32_t)grid[idx].ion >> SHIFT);
            if ((lcg_next() & 0xFFFF) < (threshold & 0xFFFF))
                grid[idx].sigma = q_clamp(grid[idx].sigma + SIGMA_GROWTH,
                                          SIGMA_LOW, SIGMA_MAX);
        }
    }
}
```

**Expected grade:** ~55–65 — filament becomes noticeably narrower; tip advances rapidly
while sides are suppressed. Morphology sharpens significantly.

---

### Step 5 — Sinh Field-Enhanced Hopping
**File:** `improved.cpp`

Replace linear `E_down` in drift flux with `sinh(q·a·E / 2kT)`:

```cpp
// sinh approximation in Q16.16: sinh(x) ≈ x + x³/6 for |x| < 1.5
//                                        ≈ exp(x)/2   for x > 1.5
inline int32_t q_sinh(int32_t x) {
    if (x > 3 * ONE / 2) {                    // x > 1.5 → sinh ≈ exp(x)/2
        // exp via: e^x = e^floor(x) * e^frac(x); use LUT or series
        // Simple: clamp to prevent overflow, then approximate
        int32_t ex = ONE;
        int32_t term = x;
        for (int i = 1; i <= 6; ++i) {        // Taylor series e^x
            ex += term;
            term = (int32_t)((int64_t)term * x / (ONE * (i + 1)));
        }
        return ex >> 1;
    }
    int32_t x3 = (int32_t)(((int64_t)x * x >> SHIFT) * (int64_t)x >> SHIFT);
    return x + (int32_t)((int64_t)x3 / 6);   // x + x³/6
}

// In drift_diffusion, replace:
//   int32_t flux = q_mul(q_mul(q_mul(ION_MOBILITY, grid[src].sigma), E_down), grid[src].ion);
// with:
    int32_t arg   = q_mul(SINH_COEFF, E_down);  // q·a·E / 2kT; SINH_COEFF ≈ 2·ONE
    int32_t sinhE = q_sinh(arg);
    int32_t flux  = q_mul(q_mul(q_mul(ION_MOBILITY, grid[src].sigma), sinhE), grid[src].ion);
```

New constant: `static constexpr int32_t SINH_COEFF = 2 * ONE;` in `physics_improved.h`.

**Expected grade:** ~65–72 — sharp threshold behavior; filament snaps through the last
gap rather than advancing gradually.

---

### Step 6 — Butler-Volmer Injection
**File:** `improved.cpp`

Replace linear anode injection with exponential overpotential-dependent injection:

```cpp
// BV anodic branch: inject ∝ exp(α·F·η / RT)
// η = V[row1, c] - V_EQ  (overpotential)
// Approximate exp in Q16.16 via short Taylor series (sufficient for η ∈ [0, 1V])
inline int32_t q_exp_approx(int32_t x) {
    int32_t result = ONE + x;
    int32_t term = x;
    for (int i = 2; i <= 5; ++i) {
        term = (int32_t)((int64_t)term * x / (ONE * i));
        result += term;
    }
    return result < ONE ? ONE : result;  // exp(x) >= 1 for x >= 0
}

// In anode injection loop:
for (int c = 1; c < N - 1; ++c) {
    int32_t eta = grid[N + c].V - V_EQ;
    if (eta <= 0) continue;
    int32_t bv  = q_exp_approx(q_mul(BV_ALPHA, eta));
    int32_t inj = q_mul(q_mul(ION_INJECT, grid[N + c].sigma), bv);
    grid[N + c].ion = q_clamp(grid[N + c].ion + inj, ION_LOW, ION_HIGH);
}
```

New constants in `physics_improved.h`:
```cpp
static constexpr int32_t V_EQ     = 0;            // equilibrium potential (0 V)
static constexpr int32_t BV_ALPHA = 2 * ONE;       // α·F/RT (scaled; ~20 in SI units)
```

**Expected grade:** ~72–80 — ion injection now starts abruptly above threshold V,
mimicking the real SET onset.

---

### Step 7 — Joule Heating Feedback
**File:** `improved.cpp`

Activate the `temp` field. After the Jacobi sweeps, compute local Joule power and update T.
In `drift_diffusion`, scale ion mobility by Arrhenius factor `exp(Eₐ·(T−T₀)/(kT·kT₀))`:

```cpp
// After jacobi_sweep block, before drift_diffusion:
static void update_temperature(std::vector<Cell>& grid, int N) {
    for (int r = 1; r < N - 1; ++r) {
        for (int c = 1; c < N - 1; ++c) {
            int idx = r * N + c;
            // Local power density ∝ σ · |∇V|²
            int32_t dV_r = grid[idx + N].V - grid[idx - N].V;  // central diff
            int32_t dV_c = grid[idx + 1].V - grid[idx - 1].V;
            int32_t E2   = q_mul(dV_r, dV_r) + q_mul(dV_c, dV_c);
            int32_t P    = q_mul(grid[idx].sigma, E2);          // σ·|∇V|²
            // Thermal balance: T = T_amb + P · R_TH
            int32_t dT   = q_mul(P, R_TH);
            grid[idx].temp = q_clamp(TEMP_INIT + dT, TEMP_INIT, TEMP_MAX);
        }
    }
}

// In drift_diffusion, replace ION_MOBILITY with temperature-scaled version:
int32_t dT      = grid[src].temp - TEMP_INIT;
int32_t mob_t   = q_mul(ION_MOBILITY, q_exp_approx(q_mul(EA_OVER_K, dT)));
int32_t flux    = q_mul(q_mul(q_mul(mob_t, grid[src].sigma), sinhE), grid[src].ion);
```

New constants in `physics_improved.h`:
```cpp
static constexpr int32_t R_TH      = 5 * ONE;       // thermal resistance (normalized)
static constexpr int32_t TEMP_MAX  = 10 * ONE;       // temperature ceiling (normalized)
static constexpr int32_t EA_OVER_K = ONE / 4;        // Eₐ/kT scale factor (≈ 0.25)
```

**Expected grade:** ~80–88 — hottest neck grows fastest; filament develops a sharp waist
as seen in TEM images; overall morphology approaches lightning-bolt shape.

---

## 6. Grading Criteria

See `grade_filament.py` for implementation. Summary:

| Metric | Max | What it measures |
|---|---|---|
| Bridging | 30 | Connected bright-σ path from row 1 to row N−2 |
| Narrowness | 25 | Fraction of interior cells that are "bright" (σ > SIGMA_MAX/2) |
| Aspect ratio | 20 | Height/width of bright bounding box |
| Tortuosity | 15 | Bridge path length / straight-line distance |
| Branch count | 10 | Number of junction points in bright skeleton |

**Grade:** A ≥ 90, B ≥ 75, C ≥ 55, D ≥ 35, F < 35

---

## 7. References

- §3.1–3.4: ECM mechanism, filament morphology, SET/RESET — `cbram_memristor_guide.md`
- §6.3: C2C variability, stochastic nucleation — `cbram_memristor_guide.md`
- §9.1: Nernst-Planck, sinh hopping rate — `cbram_memristor_guide.md`
- §9.2: Butler-Volmer electrode kinetics — `cbram_memristor_guide.md`
- §9.3: Joule heating, electrothermal coupling — `cbram_memristor_guide.md`
