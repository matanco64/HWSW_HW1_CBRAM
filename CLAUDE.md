# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build

```bash
# Quick build (all available stages + improved)
bash build.sh

# CMake build (alternative)
cmake -B build -S . && cmake --build build

# Build a single stage manually
g++ -O2 -march=native -std=c++17 -Wall -Wno-unused-result -I. naive.cpp io.cpp -o build/cbram_stage0
g++ -O2 -march=native -std=c++17 -Wall -Wno-unused-result -I. improved.cpp io.cpp -o build/cbram_improved
```

## Run

```bash
# Run all stages, verify correctness, and profile (N=1024 takes ~10-30s for stage 0)
bash run.sh 200        # N=200 for quick iteration
bash run.sh 1024       # N=1024 for profiling

# Run a single binary (must cd to build/ so output files land there)
cd build && ./cbram_stage0 200 -v
cd build && ./cbram_improved 200 -v
```

Outputs per run: `sigma_final_<stage>.bin` and `V_final_<stage>.bin` (N×N int32 row-major), plus `frames_<stage>/frame_NNNN.ppm`.

## Grade

```bash
# Grade the improved simulation (reads sigma_final_improved.bin in CWD)
python3 grade_filament.py sigma_final_improved.bin

# Grade stage 0 baseline
python3 grade_filament.py sigma_final_stage0.bin
```

The grader infers N from file size. The `/grade-filament` skill wraps this and runs the binary first.

## Architecture

### Fixed-point arithmetic
All physics uses **Q16.16 fixed-point** (`int32_t`): `ONE = 1 << 16`. Multiply two Q16.16 values with `q_mul(a, b)` which shifts right 16. Accumulate products in `int64_t` before truncating. Every constant in the headers is expressed in Q16.16.

### Grid layout
N×N grid, row 0 = anode (V = V_APPLIED), row N−1 = cathode (V = 0). Interior rows 1..N−2 are active. Four fields per cell: `V` (electric potential), `sigma` (conductivity), `ion` (ion concentration), `temp` (temperature, unused in stage 0).

**Stage 0 (`naive.cpp`):** Array-of-Structs — `struct Cell { int32_t V, sigma, ion, temp; }`. Cache-line utilization ~50% in the Jacobi hot loop (bottleneck for optimization stages).

**Stage 1–3 (not yet written):** SoA, SoA + time skewing, SoA + time skewing + OpenMP. Must produce bit-identical `V_final.bin` to stage 0.

### Simulation loop (per timestep)
1. `K=50` Jacobi sweeps: solves `∇·(σ∇V) = 0` by weighted-neighbor average in Q16.16
2. `drift_diffusion`: ions drift downward (upwind, bottom-to-top sweep) with flux `µ·σ·E·C`; anode injects ions at row 1 proportional to σ
3. `update_sigma`: conductivity grows where ions accumulate
4. Bridge check: stop when any `sigma[N−2, c] ≥ SIGMA_MAX`

### Two header families
- `physics.h` — shared by stage 0 and all optimization stages; must not change (bit-identity constraint)
- `physics_improved.h` — used only by `improved.cpp`; adds constants for steps 1–7 and Q16.16 helpers `q_sinh`, `q_exp_approx`

### Shared I/O (`io.h` / `io.cpp`)
`write_ppm` renders σ field with viridis colormap; `dump_binary` writes raw `int32_t` rows.

## Physics improvement series (`improved.cpp`)

Goal: evolve `improved.cpp` through 7 cumulative steps to produce a lightning-bolt filament (target ≥ 90/100). Start from a clone of `naive.cpp`; apply each step in-place. After each step: build → `cd build && ./cbram_improved 200 -v` → grade → record in `cbram_progress.md`.

| Step | Fix | Key change |
|------|-----|-----------|
| 0 | Baseline | Clone of naive.cpp |
| 1 | Lateral diffusion | Add `J_x = -D_LAT · ∂²C/∂x²` in `drift_diffusion` using `D_LAT = ONE/20` |
| 2 | Stochastic deposition | Replace deterministic σ growth with xorshift32 probabilistic deposit |
| 3 | Cathode nucleation | Move seeds from row 1 → row N−2; add 4× nucleation boost at cathode |
| 4 | Tip-only deposition | Gate σ growth on adjacency to metallic cell (σ > SIGMA_HIGH) |
| 5 | Sinh hopping | Replace linear `E_down` flux with `q_sinh(SINH_COEFF · E)` |
| 6 | Butler-Volmer injection | Anode injection ∝ `q_exp_approx(BV_ALPHA · η)` instead of linear |
| 7 | Joule heating | Activate `temp` field; Arrhenius mobility scaling via `EA_OVER_K` |

Full code snippets for every step are in `cbram_analysis.md`. All new constants are pre-declared (commented by step) in `physics_improved.h`.

## Grading rubric

| Metric | Points | Target |
|--------|--------|--------|
| Bridging (BFS path row 1 → row N−2) | 30 | bridge present |
| Narrowness (fraction of bright cells) | 25 | < 2% bright |
| Aspect ratio (height/width bbox) | 20 | > 5 |
| Tortuosity (path len / straight-line) | 15 | > 1.5 |
| Branch count (junction cells) | 10 | ≥ 3 junctions |

Bright threshold: σ > SIGMA_MAX/2 = 10·ONE. Grade: A ≥ 90, B ≥ 75, C ≥ 55, D ≥ 35, F < 35.
