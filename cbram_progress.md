# CBRAM Filament Improvement — Progress Tracker

Graded with `grade_filament.py` at N=200. Higher is better; target ≥ B (75/100).

## Scoring rubric (100 pts)

| Metric | Max | What a lightning filament looks like |
|---|---|---|
| Bridging | 30 | Connected bright path (σ > SIGMA_MAX/2) from row 1 to row N-2 |
| Narrowness | 25 | < 2% of interior cells bright (thin wire, not a blob) |
| Aspect ratio | 20 | Bounding box height/width > 5 (tall and thin) |
| Tortuosity | 15 | Bridge path length / straight-line distance > 1.5 (zigzag) |
| Branches | 10 | ≥ 3 junction points in narrow (< 10%) bridged filament |

Grade thresholds: **A** ≥ 90 · **B** ≥ 75 · **C** ≥ 55 · **D** ≥ 35 · **F** < 35

---

## Results

| Step | Name | Key Change | Bridging | Narrow | Aspect | Tortuous | Branches | Total | Grade |
|---|---|---|---|---|---|---|---|---|---|
| 0 | stage0 | baseline (no improvements) | 0/30 | 0/25 (34.7%) | 0/20 (1.0) | 0/15 | 0/10 | **0/100** | F |
| 1 | improved | lateral diffusion (D_LAT=ONE/20) | 0/30 | 5/25 (15.8%) | 0/20 (0.2) | 0/15 | 0/10 | **5/100** | F |
| 2 | improved | stochastic deposition (xorshift32) | 0/30 | 20/25 (2.3%) | 0/20 (0.0) | 0/15 | 0/10 | **20/100** | F |
| 3 | improved | cathode seeds (row N-2) + sigma-decoupled drift | 0/30 | 25/25 (1.8%) | 0/20 (0.0) | 0/15 | 0/10 | **25/100** | F |
| 4 | improved | tip-only deposition (gate on metal-adjacent cells) | 30/30 | 0/25 (91.2%) | 0/20 (1.0) | 0/15 | 0/10 | **30/100** | F |
| 5 | improved | sinh field-enhanced hopping (SINH_COEFF=2) | 30/30 | 0/25 (92.9%) | 0/20 (1.0) | 0/15 | 0/10 | **30/100** | F |
| 6+7 | improved | BV injection (exp(αη)) + Joule heating (Arrhenius mob) | 30/30 | 0/25 (90.9%) | 0/20 (1.0) | 0/15 | 0/10 | **30/100** | F |

---

## What step 0 reveals

The naive simulation produces a **broad saturated base**, not three stripes as expected.
After 3000 timesteps the ions have accumulated uniformly at the bottom rows, flooding 34.7% of
all interior cells above the bright threshold. No bridge is detected because the seeds (row 1)
are not connected to the saturated base via a continuously bright path — the intermediate rows
are below threshold. Both primary root causes are confirmed:

- **No bridge** → the filament never closes the gap between anode and cathode.
- **34.7% bright / aspect 1.0** → wide, diffuse ion accumulation, not a narrow wire.

Fix 1 (lateral diffusion) will help ions spread more realistically but won't fix the bridge.
Fix 2 (stochastic deposition) is the first step expected to produce a bridge.

---

## What step 1 reveals

Lateral diffusion (D_LAT = ONE/20 = 0.05) reduced bright-cell fraction from 34.7% → 15.8%,
earning 5 narrowness points. However, the aspect ratio collapsed from 1.0 → 0.2 (wide horizontal
band), because the diffusion is strong enough to homogenize ion concentration across the full width
of the anode rows before they can drift downward. The three seed columns spread into a uniform
stripe instead of narrow fingers. No bridge, no vertical structure.

**Root cause of aspect regression:** D_LAT competes with ION_MOBILITY·σ·E drift. At the anode
rows where σ is highest, the lateral flux overwhelms the downward drift at early timesteps, locking
ions into a horizontal band. Fix 2 (stochastic deposition) and Fix 4 (tip-only growth) should
counteract this by localizing sigma growth to narrow tips rather than the whole band.
