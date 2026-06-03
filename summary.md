# CBRAM Filament Physics — Session Summary

**Date:** 2026-06-04  
**Grader:** `grade_filament.py` · **Grid:** N=200 · **Steps:** 3000 timesteps, 50 Jacobi iters

---

## What We Did

### All 7 physics fixes implemented in `improved.cpp` + `physics_improved.h`

| Step | Change | File(s) | Score |
|---|---|---|---|
| 0 | Baseline clone of naive.cpp | improved.cpp | 0/100 F |
| 1 | Lateral diffusion: `J_x = -D_LAT · ∂²C/∂x²` (D_LAT = 0.05) | improved.cpp | 5/100 F |
| 2 | Stochastic deposition: xorshift32 per-cell probability | improved.cpp | 20/100 F |
| 3 | Cathode nucleation: seeds moved to row N-2, sigma decoupled from drift | improved.cpp | 25/100 F |
| 4 | Tip-only deposition: σ grows only when adjacent to σ > SIGMA_HIGH; ION_HIGH pre-fill | improved.cpp | 30/100 F |
| 5 | Sinh field-enhanced hopping: flux ∝ sinh(SINH_COEFF · E) | improved.cpp | 30/100 F |
| 6 | Butler-Volmer injection: anode inject = ION_INJECT · exp(BV_ALPHA · η) | improved.cpp | 30/100 F |
| 7 | Joule heating: temp = TEMP_INIT + σ·|∇V|²·R_TH; mob_t ∝ exp(EA/kT · ΔT) | improved.cpp | 30/100 F |

### Infrastructure built
- `grade_filament.py` — 5-metric 100-pt grader (bridging, narrowness, aspect ratio, tortuosity, branches)
- `cbram_analysis.md` — deep root-cause analysis identifying all 7 physics defects
- `cbram_progress.md` — per-step grade table
- `prompts.md` — full AI prompt log as required by the assignment

---

## Where We Stand

**Current grade: 30/100 (F)**

```
Bridging   : 30/30  ← bridge exists (length 198, nearly straight)
Narrowness :  0/25  ← 90.9% bright — entire grid is metallic blob
Aspect     :  0/20  ← ratio 1.0 — blob fills full width
Tortuosity :  0/15  ← ratio 1.01 — path is basically a straight line
Branches   :  0/10  ← 35 486 junctions (too wide: needs < 10% bright)
```

The bridge forms quickly and robustly. The problem is that it is a **blob**, not a wire.

### Root cause of the blob

Step 4 introduced two changes that interact badly:
1. **`SIGMA_GROWTH_IMP = ONE` (1.0)** — when `prob_q16 = SIGMA_GROWTH_IMP * ion` and `ion = ION_HIGH = 1.0`, the deposit probability = 100%. Every cell adjacent to metal deposits on every timestep.
2. **ION_HIGH pre-fill** — all electrolyte cells are initialized to `ION_HIGH = 1.0` ions. This removes the diffusion-limited supply that is the physical engine of DLA branching.

Together these cause the metal front to advance at full speed in ALL directions simultaneously — there is no selectivity between tip and sides. The result is a blob rather than a dendritic tree.

Steps 5, 6, and 7 are correctly implemented but cannot rescue a process that has 100% deposition probability everywhere. Joule heating can only focus growth at a narrow neck; when 91% of cells are metallic there is no neck for E to concentrate across.

---

## Future Work — How to Get Lightning

The guide (`cbram_memristor_guide.md` §3.2, §6.3, §9.1) says the CF is fractal/dendritic, arising from **Diffusion-Limited Aggregation (DLA)**. Real DLA requires:

> "Nucleation is stochastic at multiple sites simultaneously … C2C variability originates from probabilistic Poisson-distributed random events in ion hopping." — §6.3

DLA only produces branching when **both** conditions hold:
1. Growth probability is genuinely < 1 (thermal noise matters)
2. Ion supply is diffusion-limited (tip competes with sides for ions)

### Fix A — Reduce deposition probability (critical)

In `physics_improved.h`, change:
```cpp
static constexpr int32_t SIGMA_GROWTH_IMP = ONE;       // 100% probability — TOO AGGRESSIVE
```
to something probabilistic:
```cpp
static constexpr int32_t SIGMA_GROWTH_IMP = ONE / 5;   // 20% base probability
```

With `ion = ION_HIGH = 1.0` and `SIGMA_GROWTH_IMP = ONE/5`:
- `prob_q16 = ONE/5 * ONE = ONE/5 → threshold = (ONE/5 << 16) = 20% chance`
- Now stochastic: the tip (higher E, more ions) deposits ~20% of steps; sides deposit less

### Fix B — Reduce or eliminate ION_HIGH pre-fill (critical)

In `improved.cpp` (the pre-fill loop at step 4 init):
```cpp
// Change: all electrolyte to ION_HIGH
grid[r * N + c].ion = ION_HIGH;          // current (too much)
// To: a fraction, e.g. ION_HIGH / 4
grid[r * N + c].ion = ION_HIGH / 4;      // diffusion-limited regime
```

With lower initial ion concentration, the filament tip (where E is highest) depletes local ions faster than diffusion can replenish them from the sides → the tip's local ion deficit creates lateral asymmetry → branching emerges.

### Fix C — Verify Joule heating magnitude (tune R_TH)

Currently `R_TH = 5` but in the blob state `ΔT ≈ 0.009` (nearly zero) because |∇V|² is tiny across the uniform metallic blob. Once Fixes A+B create a narrow filament, E concentrates at the neck and Joule heating will activate naturally. No code change needed — just verify it fires correctly after A+B.

### Fix D — Reduce BV_ALPHA or ION_INJECT if over-injecting

With fix A+B, the anode row may re-flood if BV injection `ION_INJECT * exp(BV_ALPHA * 1.0) ≈ 0.74` is too large. May need to tune:
```cpp
static constexpr int32_t ION_INJECT = ONE / 20;   // reduce from ONE/10 to ONE/20
```

### Expected grade trajectory after A+B

| Metric | Expected |
|---|---|
| Bridging | 30/30 (bridge still forms, but slower — may need more timesteps) |
| Narrowness | 20–25/25 (< 5% bright if DLA is working) |
| Aspect | 15–20/20 (tall thin column with side branches) |
| Tortuosity | 10–15/15 (dendritic path meanders) |
| Branches | 6–10/10 (3+ junction points in narrow bridged structure) |
| **Total** | **~75–90/100 (B or A)** |

---

## Ramp-Up Prompt for Next Session

Copy this verbatim to start the next Claude Code session:

---

> We are working on a CBRAM filament simulation in `/csl/ece882-017/cbram/HWSW_HW1_CBRAM/`. Read `summary.md` to understand the full state. The short version:
>
> All 7 physics fixes are implemented in `improved.cpp` + `physics_improved.h`. The current grade is **30/100 F** — bridge works but the entire grid is a metallic blob (90.9% bright). The root cause is that `SIGMA_GROWTH_IMP = ONE` + ION_HIGH pre-fill gives 100% deposition probability everywhere adjacent to metal, so the filament expands as a blob instead of a dendritic tree.
>
> The next task is: **make the lightning shape emerge** by tuning two parameters:
> 1. Reduce `SIGMA_GROWTH_IMP` from `ONE` to `ONE/5` in `physics_improved.h` (make deposition genuinely probabilistic)
> 2. Reduce the ION_HIGH pre-fill to `ION_HIGH/4` in the step-4 init block in `improved.cpp` (restore diffusion-limited ion supply)
>
> After each change: rebuild with `g++ -O2 -march=native -std=c++17 -I. improved.cpp io.cpp -o build/improved`, run with `./build/improved 200 -v`, grade with `python3 grade_filament.py sigma_final_improved.bin`, and record in `cbram_progress.md`. Target: **B (75/100)** or better. If the bridge breaks (too slow), increase `TOTAL_TIMESTEPS` in `physics_improved.h` before reducing growth rate.

---

*End of summary. Files ready to commit: `improved.cpp`, `physics_improved.h`, `grade_filament.py`, `cbram_analysis.md`, `cbram_progress.md`, `prompts.md`, `summary.md`.*
