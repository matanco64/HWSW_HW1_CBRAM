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

---

## Session 4 — Model pivot: Dielectric Breakdown Model (DBM)

The continuum drift-diffusion line above topped out at 30/100 (a saturated blob, never narrow).
We abandoned it for the **Dielectric Breakdown Model**: grow a metallic cluster one cell at a
time; each step solve Laplace (cluster pinned to cathode V=0, anode at V_APPLIED) and add one
cluster-adjacent empty cell with probability ∝ V^η. One knob, η, dials morphology
(0→bush, 1→branched, large→needle).

### Python ground truth (`python_imp/cbram_dbm.py`)

| η | Bridging | Narrow | Aspect | Tortuous | Branches | Total | Grade |
|---|---|---|---|---|---|---|---|
| 1 | 30/30 | 12/25 (8.1%) | 0/20 (1.9) | 10/15 | 10/10 | 62 | C |
| 2 | 30/30 | — | — | — | — | 85 | B |
| **3** | **30/30** | **25/25 (1.2%)** | **20/20 (6.8)** | **10/15 (1.33)** | **10/10** | **95** | **A** |
| 4 | — | — | — | — | — | 95 | A |
| 6 | — | — | — | — | — | 90 | A |

η=3 locked in as the sweet spot.

### C++ naive baseline (`dbm_stage0.cpp`, Q16.16, N=200)

| Step | Name | Key Change | Bridging | Narrow | Aspect | Tortuous | Branches | Total | Grade |
|---|---|---|---|---|---|---|---|---|---|
| 0 | dbm_stage0 | naive AoS DBM, Q16.16, warm-started Jacobi | 30/30 | 25/25 (1.5%) | 15/20 (3.6) | 10/15 (1.34) | 10/10 (143) | **90/100** | A |

Bridged at step 588 in ~0.76 s. Aspect is lower than Python's 6.8 (fixed-point V + xorshift PRNG
differ from float/PCG64), but it is a clear A-grade vertical branched bolt. **`V_final` is
byte-identical across runs** — the determinism the optimization ladder requires.

Key port detail: a Q16.16 `V^3` underflows to 0 for the small V near the cathode and stalls
growth, so the DBM weight is computed at full `__int128` precision (still pure deterministic
integer math → bit-identity holds).

### Optimization ladder (future, all bit-identical to stage 0 via `cmp V_final`)

| Stage | File | Change | Cache story |
|---|---|---|---|
| 1 | `dbm_stage1.cpp` | AoS → SoA (separate `V[]`, `metal[]`) | Jacobi streams V at ~100% cache-line use (vs ~50% AoS) |
| 2 | `dbm_stage2.cpp` | Cache blocking + time-skewing | March a tile through several of the 30 sweeps while hot in L1/L2 |
| 3 | `dbm_stage3.cpp` | OpenMP | Double-buffer Jacobi is race-free → parallel + bit-identical |

(Red-Black Gauss-Seidel was considered and dropped — it reads updated values mid-sweep, so it
cannot be bit-identical to the Jacobi baseline.)

### Visualization

`dbm_stage0.cpp` dumps per-frame V+σ binaries to `frames_stage0/frame_NNNNNN_{V,S}.bin`
(every `FRAME_INTERVAL` growth steps, plus a final frame on the bridged state).
`python_imp/render_dbm_video.py` reads the pairs and encodes a 3-panel MP4 — σ (filament) |
φ (Jacobi solve) | |∇φ| (drives growth):

```bash
cd build && ./cbram_stage0 200                     # writes frames_stage0/*.bin
cd ../python_imp && python3 render_dbm_video.py 200 ../build/frames_stage0 filament_stage0.mp4
```

Output `python_imp/filament_stage0.mp4` (1920×720, 295 frames, ~10 s @ 30 fps).
(Bump `FRAME_INTERVAL` in `physics_dbm.h` for large-N perf runs — frame I/O is ~160 KB/field/frame.)

### Determinism (verified)

Independent runs at N=200 are bit-for-bit reproducible — `V_final`, `sigma_final`, and all 590
frame files are byte-identical across runs (`diff -rq` clean, matching md5s). Sources are all
deterministic: fixed-seed `xorshift64`, row-major candidate order, pure `__int128` integer pick,
integer Jacobi, single-threaded, no float in the compute path. This is the reproducibility every
optimization stage must preserve; the stage0-vs-stageN bit-identity gate (`run.sh` `cmp`) builds
on it. The C++ is intentionally *not* bit-identical to the Python float ground truth.
