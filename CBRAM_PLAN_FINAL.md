# CBRAM_PLAN_FINAL.md — HW1: CBRAM Device Simulation, Cache & Parallelism Optimization

> Build spec for Claude Code. Two-person team (Matan + Yuval).
> Merges CBRAM_PLAN.md (V1) + CBRAM_PLAN_v2.md (V2).
> Key decisions:
> - **Separate source file per optimization stage** (not a single evolving file).
> - **Integer Q16.16 arithmetic** for guaranteed bit-identical `cmp` across all stages (from V2).
> - **Three layered optimizations** (SoA → time skewing → OpenMP), each in its own file (from V2).
> - **Diffs made explicit** via a standard header block in each file and a diff table in this plan.

---

## 0. Summary

A 2D CBRAM device simulator in C++ using Q16.16 fixed-point integer arithmetic. Each grid cell
carries four physical fields (V, σ, ion, T). The dominant kernel is an iterative Jacobi Poisson
solve, with drift-diffusion of metal ions producing a conductive filament between the electrodes.

We deliver **four source files**, each a complete standalone simulator, differing only in how they
use hardware resources. Each stage gates on `cmp`-byte-identical output against the naive baseline
before we move on. The report's central argument is that each optimization peels away one
bottleneck to expose the next — bandwidth → reuse distance → parallelism.

---

## 1. Hardware mechanisms (one per file)

| File | HW mechanism | Why it helps |
|---|---|---|
| `naive.cpp` | none — AoS baseline | reference |
| `opt1_soa.cpp` | cache-line utilization | hot kernels load only the fields they actually read |
| `opt2_blocked.cpp` | reuse distance / arithmetic intensity | tile stays in L2 across multiple Jacobi sweeps |
| `opt3_omp.cpp` | thread-level parallelism + per-core cache | parallel scaling, only viable after escaping bandwidth wall |

**The report's narrative:**
> Each layer attacks the bottleneck exposed by the previous one. SoA reduces wasted bandwidth →
> exposes that reuse distance is still forcing DRAM round-trips → time skewing reduces those round-
> trips → exposes that we're now compute-bound and can scale threads → OpenMP. The cache work does
> not just speed up serial execution — it **unlocks** parallel scaling. An OpenMP version of
> `naive.cpp` would hit the bandwidth wall almost immediately.

---

## 2. Why integer arithmetic (Q16.16)

- **Bit-identical `cmp` is trivially achievable.** Floats break under auto-vectorization (different
  FMA fusion, reordered reductions). `int32_t` arithmetic with the same iteration order produces
  exactly the same bits on every stage — no tolerance gymnastics needed.
- **Same memory footprint.** `int32_t` = 4 bytes = `float`. AoS struct stays 16 B; all cache-line
  math is unchanged.
- **Sufficient range.** Q16.16 in `int32_t` covers ±32,768 with ~1.5×10⁻⁵ resolution. Use
  `int64_t` for the accumulator inside the Poisson update (Q32.32 product), then shift back.
  Real σ ~10, V ~1 → products ~10¹⁰, well inside int64.
- **Caveat for the report:** integer Jacobi converges to a slightly different steady-state than
  float Jacobi (per-step rounding differs). That is fine — `cmp` is always integer-vs-integer,
  never int-vs-float.

**Poisson update (Q16.16, SHIFT = 16):**
```cpp
int32_t s_c = sigma[idx];
int32_t s_e = (s_c + sigma[idx+1]) >> 1;
int32_t s_w = (s_c + sigma[idx-1]) >> 1;
int32_t s_n = (s_c + sigma[idx+N]) >> 1;
int32_t s_s = (s_c + sigma[idx-N]) >> 1;
int32_t denom = s_e + s_w + s_n + s_s;           // Q16.16

int64_t num = (int64_t)s_e * V[idx+1] + (int64_t)s_w * V[idx-1]
            + (int64_t)s_n * V[idx+N] + (int64_t)s_s * V[idx-N]; // Q32.32

V_next[idx] = (int32_t)(num / denom);            // Q32.32 / Q16.16 = Q16.16
```

---

## 3. Deliverables (assignment map)

1. `naive.cpp` — AoS baseline (stage 0).
2. `opt1_soa.cpp` — SoA layout (stage 1).
3. `opt2_blocked.cpp` — SoA + time-skewed Jacobi (stage 2).
4. `opt3_omp.cpp` — SoA + time-skewing + OpenMP (stage 3).
5. `physics.h` — shared constants, Q16.16 helpers, init routine (identical across all files).
6. `io.cpp` / `io.h` — PPM frame dumper + binary state dumper (never changes across stages).
7. `run.sh` — compile all four, verify correctness, profile all, produce ablation table.
8. `viz/make_video.sh` — ffmpeg side-by-side filament race (naive vs opt3).
9. `report_draft.md` — ≤3-page draft (Matan finalizes to PDF).
10. `prompts.md` — log of all AI prompts (log continuously, not at the end).
11. (Separately) names + IDs PDF.

---

## 4. Repository structure

```
cbram-hw1/
├── CBRAM_PLAN_FINAL.md
├── physics.h            (shared constants, Q16.16 macros, init routine)
├── io.h / io.cpp        (PPM dump, binary state dump — unchanged across stages)
├── naive.cpp            (stage 0 — AoS)
├── opt1_soa.cpp         (stage 1 — SoA)
├── opt2_blocked.cpp     (stage 2 — SoA + time skewing)
├── opt3_omp.cpp         (stage 3 — SoA + time skewing + OpenMP)
├── run.sh
├── viz/
│   ├── make_video.sh
│   └── colormap.h
├── results/             (perf output dumps — gitignore this)
├── report_draft.md
└── prompts.md
```

---

## 5. Making diffs explicit between files

Each source file opens with a standard header block that explicitly lists what changed from the
previous stage and why. This makes the diff story immediately visible to a grader without needing
to run `diff`.

```cpp
// ============================================================
// Stage N — <name>
// Changes from stage N-1:
//   - <change 1 and the HW reason>
//   - <change 2 and the HW reason>
// Unchanged: arithmetic, iteration order, physics, init seed.
// ============================================================
```

Additionally, `run.sh` prints a `diff --unified` between consecutive pairs at the end of its
output, and the report includes the diff summary table from §1.

To inspect any particular change:
```bash
diff -u naive.cpp opt1_soa.cpp        # AoS → SoA
diff -u opt1_soa.cpp opt2_blocked.cpp # SoA → time skewing
diff -u opt2_blocked.cpp opt3_omp.cpp # skewing → OpenMP
```

---

## 6. Shared physics design (identical across all four files)

**Domain.** 2D grid, N×N cells. Top row = active electrode (Dirichlet V = V_applied; ion source).
Bottom row = inert electrode (V = 0). Left/right = Neumann (zero-flux). Interior cells only are
updated; boundary values are clamped each iteration.

**Fields per cell (4 × int32_t = 16 B in Q16.16):**
- `V` — electric potential
- `sigma` — local conductivity (grows with ion/atom density)
- `ion` — metal ion concentration in electrolyte
- `temp` — local temperature (uniform init; Joule heating optional)

**Initial conditions.** Deterministic. Uniform low σ everywhere except a few seed defect sites
placed by a fixed-seed PRNG (`std::mt19937` seed = 42). **All four files call the same `init()`
from `physics.h`** — this is mandatory for `cmp` to work.

**Simulation loop (per outer timestep Δt):**
1. **Poisson solve for V** — Jacobi iteration, K iterations per timestep (K=50 default). *Hotspot.*
2. **Drift-diffusion for ions** — explicit Euler update of `ion` using V and σ.
3. **Conductivity update** — `sigma` = smooth (logistic) function of local ion concentration.
4. *(Optional for v1)* **Heat equation for T** — Joule source σ|∇V|². Skip unless time allows.
5. **Periodic frame dump** — write σ as PPM every M timesteps.

**Termination.** Fixed number of outer timesteps (tuned to produce a clear bridging filament).
Optionally detect bridge via a BFS/DFS connectivity check on high-σ cells.

**Output for correctness.** At end: dump the full `V` field as binary (`V_final_stageN.bin`).
`naive.cpp` produces the reference; every other stage is `cmp`'d byte-for-byte against it.

---

## 7. Per-file specs

### `naive.cpp` — Stage 0: AoS baseline

```
Stage 0 — AoS Naive
Changes from previous: (none — this is the baseline)
```

- `struct Cell { int32_t V, sigma, ion, temp; };` — 16 B, 4 cells per 64-byte cache line.
- `Cell* grid; Cell* grid_next;` — double-buffered Jacobi.
- Poisson hot loop reads `grid[idx].V`, `grid[idx].sigma` from neighbors, but each cache-line
  fetch carries all 4 fields. **Cache-line utilization ≈ 50%** — `ion` and `temp` ride along
  unused in the Poisson kernel.
- Frames → `frames_stage0/`. Final dump → `V_final_stage0.bin` (reference for all `cmp` checks).

---

### `opt1_soa.cpp` — Stage 1: Struct of Arrays

```
Stage 1 — SoA (Struct of Arrays)
Changes from stage 0:
  - Four flat arrays replace the Cell struct: int32_t *V, *sigma, *ion, *temp; int32_t *V_next.
  - Allocated with posix_memalign(..., 64, ...) for 64-byte alignment.
  HW reason: Poisson kernel fetches 16 contiguous V values per cache line (100% utilization)
  instead of 4 cells * 4 fields (50% utilization). ~2x reduction in DRAM bandwidth for the
  hot kernel. SoA contiguous-int32 layout also auto-vectorizes more aggressively.
Unchanged: arithmetic, iteration order, K, N, init seed.
```

- **No algorithmic changes.** Same Jacobi formula, same drift-diffusion, same σ update.
- `cmp V_final_stage1.bin V_final_stage0.bin` must pass byte-for-byte.
- Frames → `frames_stage1/`.

---

### `opt2_blocked.cpp` — Stage 2: Time Skewing (Temporal Blocking)

```
Stage 2 — Time Skewing (Temporal Blocking)
Changes from stage 1:
  - The K-iteration Jacobi inner loop is restructured into trapezoidal spatial-temporal tiles.
  - Tile parameters: spatial block B (chosen so B² × 8 B × 2 buffers fits in ~½ L2),
    temporal depth M (typically 4–8; must satisfy B - 2M > 0).
  HW reason: Without blocking, each cell's V value is evicted from L2 between Jacobi iterations
  (reuse distance > L2 capacity). With a trapezoid of depth M, the same data block is reused
  M times while still in L2, reducing DRAM round-trips by ~M× for the Poisson kernel.
Unchanged: arithmetic per cell, total K updates per cell, init seed.
```

One trapezoid over spatial block `(jj, ii)` of width `B` and temporal depth `M`:
- For iteration `m = 0..M-1`, update interior region shrinking by 1 cell per side per iteration.
- Each cell receives exactly K updates total, same arithmetic as stage 1.
- Boundary tiles (grid edges) truncate the trapezoid; corners handled similarly.
- **This is the fiddly part** — write it carefully, verify `cmp` before profiling.
- `cmp V_final_stage2.bin V_final_stage0.bin` must pass.
- Frames → `frames_stage2/`.

---

### `opt3_omp.cpp` — Stage 3: OpenMP Parallelization

```
Stage 3 — OpenMP (Thread-Level Parallelism)
Changes from stage 2:
  - #pragma omp parallel for collapse(2) over the (jj, ii) tile loop.
  - Build flag: add -fopenmp.
  HW reason: After time skewing, the Poisson kernel is compute-bound (not bandwidth-bound).
  Multiple cores each have their own L1/L2, and the tiles are independent — no read/write
  conflicts under Jacobi double-buffering (read V[], write V_next[]).
  Without the cache work of stages 1–2, OpenMP would immediately hit the shared DRAM bandwidth
  wall and show near-zero scaling; the cache optimizations unlock it.
Unchanged: arithmetic, tile logic, init seed.
```

- Inner iteration loop within a tile stays serial (no intra-tile parallelism needed).
- `OMP_PROC_BIND=close OMP_PLACES=cores` — bind threads to physical cores, not SMT.
- `cmp V_final_stage3.bin V_final_stage0.bin` must pass (Jacobi parallelization is deterministic
  given the same tile decomposition).
- Scaling sweep: T = 1, 2, 4, up to the machine's physical core count (confirm with `lscpu`).

---

## 8. Correctness verification

**Primary gate — byte-identical at every stage:**
```bash
cmp V_final_stage1.bin V_final_stage0.bin && echo "Stage 1: BIT-IDENTICAL OK"
cmp V_final_stage2.bin V_final_stage0.bin && echo "Stage 2: BIT-IDENTICAL OK"
cmp V_final_stage3.bin V_final_stage0.bin && echo "Stage 3: BIT-IDENTICAL OK"
```
Because all arithmetic is integer with identical per-cell iteration order, this is genuinely
achievable — not a tolerance approximation.

**If `cmp` fails, do not proceed.** Diagnose:
- Stage 1: stale read of `V[]` instead of `V_next[]`, or wrong stride.
- Stage 2: trapezoid bounds wrong; a cell receiving ≠ K updates. Check edge tiles.
- Stage 3: race condition (should not happen with Jacobi double-buffering; if it does, the
  tile decomposition was wrong).

**Secondary check:** final σ PPM frames should be visually identical across all stages. Include a
4-up comparison frame in the report.

---

## 9. Build & environment

```bash
CXXFLAGS="-O2 -march=native -std=c++17 -Wall -Wno-unused-result"

g++ $CXXFLAGS naive.cpp        io.cpp -o cbram_stage0
g++ $CXXFLAGS opt1_soa.cpp     io.cpp -o cbram_stage1
g++ $CXXFLAGS opt2_blocked.cpp io.cpp -o cbram_stage2
g++ $CXXFLAGS -fopenmp opt3_omp.cpp io.cpp -o cbram_stage3
```

**Hard rules — violating any of these invalidates the comparison:**
- Same `-O2 -march=native` for stages 0–2; stage 3 adds `-fopenmp` only.
- Same N, K, timestep count, init seed across all four binaries.
- Time only the simulation kernel; exclude allocation and I/O.

**Environment setup:**
```bash
lscpu > results/env.txt                          # log cache sizes + core count
taskset -c 2 ./cbram_stage0 4096               # pin to one core for serial stages
# for stage 3:
OMP_PROC_BIND=close OMP_PLACES=cores ./cbram_stage3 4096
```

Tune N (default 4096, 256 MB AoS working set >> LLC) and timestep count so `cbram_stage0`
takes 10–30 seconds. This puts the working set firmly in the DRAM-streaming regime where the
bandwidth differential shows up cleanly.

---

## 10. Profiling — ablation study

Same `perf` invocation on every stage:

```bash
for STAGE in 0 1 2 3; do
  perf stat -r 5 \
    -e cycles,instructions,\
cache-references,cache-misses,\
L1-dcache-loads,L1-dcache-load-misses,\
LLC-loads,LLC-load-misses,\
dTLB-loads,dTLB-load-misses \
    ./cbram_stage$STAGE 4096 2> results/stage${STAGE}.perf
done
```

Confirm the hotspot is the Poisson loop:
```bash
perf record -g ./cbram_stage0 4096 && perf report --stdio | head -40
```

**Thread scaling sweep on stage 3:**
```bash
for T in 1 2 4 8; do   # cap at physical core count
  OMP_NUM_THREADS=$T OMP_PROC_BIND=close OMP_PLACES=cores \
    perf stat -r 5 -e cycles,instructions,LLC-load-misses \
    ./cbram_stage3 4096 2> results/stage3_t${T}.perf
done
```

**Ablation table for the report:**

| Metric | Stage 0 (AoS) | Stage 1 (SoA) | Stage 2 (+Blocking) | Stage 3 (+OMP) |
|---|---|---|---|---|
| Wall-clock (s) | | | | |
| Speedup vs stage 0 | 1.0× | | | |
| IPC | | | | |
| LLC-load-miss rate | | | | |
| Est. DRAM bytes (LLC-misses × 64 B) | | | | |

Plus a stage-3 thread-scaling plot (speedup vs threads, with ideal-linear reference).

**Honest caveats to include:**
- Stage 1's IPC jump is partly from auto-vectorization on contiguous int32 (SoA is more
  SIMD-friendly). Attribute both effects explicitly — bandwidth reduction *and* vectorization.
  Optionally run with `-fno-tree-vectorize` on stages 0 and 1 to isolate the pure bandwidth
  component; include this as a secondary experiment if time allows.
- Stage 3 scaling will be sub-linear — that is expected and interesting. The report should show
  it scales *better* than an OpenMP naive would (which would hit the bandwidth wall immediately).
- `LLC-load-misses` includes HW prefetcher behavior; treat it as a proxy, corroborated by
  wall-clock and IPC.

---

## 11. Visualization

**Per-stage frame dumps:**
- Each binary writes σ PPM frames to `frames_stage{0,1,2,3}/` every M timesteps.
- Colormap via viridis (single-header `viz/colormap.h`). Aim for ~100 frames per run.
- Use simulation-timestep index for filenames (not wall-clock) so frames align across stages.

**Side-by-side video (`viz/make_video.sh`):**
- Stage 0 vs stage 3 side-by-side via ffmpeg.
- Overlay a wall-clock timestamp on each panel — the viewer sees stage 3's clock advancing
  faster while both filaments grow identically.
- Output `filament_race.mp4`.

**Static figures for the report:**
1. **Memory layout diagram** — one 64 B cache line under AoS (4 cells × 4 fields, V+σ
   highlighted = 50%) vs SoA (16 V values = 100%). Makes the stage-1 mechanism instantly visible.
2. **Trapezoid time-skewing diagram** — one 2D tile showing 4 iterations, each shrinking by 1
   cell per side. Makes the stage-2 mechanism visible.
3. **Ablation bar chart** — wall-clock for all 4 stages (log scale).
4. **Thread scaling plot** — stage 3 speedup vs threads with ideal-linear reference.

---

## 12. Report structure (≤3 pages)

1. **Approach** (½ page) — Why CBRAM (research-relevant; multi-field physics naturally motivates
   the layout question). The layered-bottlenecks insight: each optimization exposes the next one.
2. **Naive (stage 0)** (½ page) — What it computes. AoS is the natural physical grouping but
   wastes 50% of bandwidth in the Poisson kernel. Profile numbers. Cache-line utilization framing.
3. **Optimizations** (1 page) — Three subsections (SoA / time skewing / OpenMP). For each:
   the change made → HW/SW insight → profiling numbers → was the result expected or surprising?
   Include the layout diagram and trapezoid diagram here.
4. **Ablation comparison** (½ page) — The headline table + bar chart + scaling plot. Argue why
   each layer's gain matches its mechanism (cache-line size, reuse distance ÷ M, cores × per-core L2).
5. **What we tried that didn't help / next steps** (½ page) — e.g., `__builtin_prefetch` (HW
   prefetcher already handles streaming), manual SIMD intrinsics (auto-vectorizer gets most of it
   on SoA int32), huge pages (small effect on streaming), naive + OMP (bandwidth wall kills scaling).

Front-load the correct mechanism: *cache-line utilization*, not "every access is a miss."

---

## 13. Work sequence

**Gate rule: `cmp` must pass before moving to the next stage. No exceptions.**

1. **Scaffold** — `physics.h`, `io.cpp`/`io.h` (frame + state dump), `run.sh` skeleton,
   `viz/` directory.
2. **Stage 0 (`naive.cpp`)** — implement full physics + Poisson + drift-diffusion + σ update.
   Tune N and timestep count for 10–30 s runtime. Save reference: `cp V_final_stage0.bin ref.bin`.
3. **Stage 1 (`opt1_soa.cpp`)** — port to SoA. `cmp` against `ref.bin` must pass.
4. **Stage 2 (`opt2_blocked.cpp`)** — add time skewing. `cmp` must pass.
   *This is the hardest stage — budget 2–3× the effort of others.*
5. **Stage 3 (`opt3_omp.cpp`)** — add OpenMP. `cmp` must pass.
6. **Visualization** — PPM dumper (already wired in stage 0), `make_video.sh`.
7. **`run.sh`** — full ablation across all 4 stages + thread scaling sweep.
8. **Report draft** — fill in ablation table with real numbers; generate 4 figures.
9. **`prompts.md`** — must be maintained *throughout*, not filled retroactively.

**Two-person split suggestion:**
- One owns physics + correctness (stages 0–1, init, drift-diffusion, σ update, `cmp`).
- The other owns systems (stages 2–3, `perf`, `run.sh`, viz pipeline, scaling experiment).
- Both review the trapezoid math in stage 2 together — it's the bug-prone part.

---

## 14. Open question before stage 3

**How many physical cores does the test machine have?** Confirm with `lscpu` before starting
stage 3, and plan the scaling sweep accordingly (T = 1, 2, 4 on a 4-core laptop; up to 8–16 on a
lab machine). Use physical cores only — do not count SMT/hyper-threading threads.

---

## 15. Pitfalls

- ❌ Skipping `cmp` between stages. Integer arithmetic makes byte-equality genuinely achievable
  — use it as the hard gate.
- ❌ Different compiler flags across stages (except `-fopenmp` for stage 3). Confounded comparison.
- ❌ Algorithmic changes alongside the layout/scheduling change (e.g., switching to red-black
  Gauss-Seidel in stage 2). The story is "same algorithm, better HW use." Keep it that way.
- ❌ Profiling before `cmp` passes. Speed of wrong code is meaningless.
- ❌ Reporting stage-3 scaling with SMT threads (they share L1/L2 — the scaling story is muddied).
- ❌ Claiming pure-bandwidth speedup for stage 1 when vectorization is partly responsible.
- ❌ Padding the AoS struct artificially to inflate the bandwidth difference. Use honest,
  physically-motivated fields only.
- ❌ Trapezoid off-by-one at grid boundaries causing edge cells to receive ≠ K updates (silent
  correctness break — which is exactly why `cmp` is the gate, not a visual check).
- ❌ Timing allocation and I/O inside the measured kernel window. Time only the simulation loop.
- ❌ Reusing a class example program. AoS→SoA as a technique is allowed; applying it to a CBRAM
  simulator is novel and that's the correct framing.
