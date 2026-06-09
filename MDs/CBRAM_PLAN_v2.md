# CBRAM_PLAN_v2.md — HW1: CBRAM Cache + Threading Optimization

> Build spec for Claude Code. Two-person team (Matan + Yuval).
> **Supersedes CBRAM_PLAN.md.** Key changes in v2:
> - Integer (Q16.16 fixed-point) arithmetic → trivially bit-identical across all stages.
> - Scope expanded to **C** (int + SoA + time skewing + OpenMP).
> - Project structured as a **linear git history with one commit per optimization**, so each
>   stage can be checked out and profiled independently → ablation table in the report.

---

## 0. Summary

A 2D CBRAM device simulator in C++ using Q16.16 fixed-point integer arithmetic. Each grid cell
carries V, σ, ion concentration, and temperature; the dominant kernel is an iterative Jacobi
Poisson solve, with simple drift-diffusion producing a conductive filament between the
electrodes.

We build the project as a **linear sequence of commits**, each adding exactly one optimization.
Every commit produces a working binary that simulates the same physics and writes identical
output. The grader (or us) can `git checkout` any stage, build, run, and `perf stat` it.

---

## 1. Hardware mechanisms in play (one per commit)

Each stage targets one HW concept; the report's strength is showing they **stack** correctly.

| Stage | HW mechanism | Why it helps |
|---|---|---|
| 0 — naive | none | baseline |
| 1 — SoA | cache-line utilization | hot kernels load only the fields they read |
| 2 — time skewing | reuse distance / arithmetic intensity | tile stays in L2 across multiple Jacobi sweeps |
| 3 — OpenMP | thread-level parallelism + per-core cache | parallel scaling, only possible once off the bandwidth wall |

**The narrative for the report:**
> Each layer attacks a different bottleneck exposed by the previous one. SoA reduces wasted
> bandwidth → exposes that we're still bandwidth-bound by reuse distance → time skewing reduces
> DRAM traffic → exposes that we're now compute-bound and can scale to more cores → OpenMP.
> **The cache optimizations do not just make the serial code faster — they unlock parallel
> scaling.** Without them, threading hits the DRAM bandwidth wall almost immediately.

---

## 2. Why integer arithmetic (Q16.16)

- **Bit-identical comparison is trivial** — `cmp` works across all four stages with no
  tolerance gymnastics. Floats break this under auto-vectorization (different FMA fusion,
  reordered reductions); integers don't.
- **Same memory footprint** — `int32_t` is 4 bytes, same as `float`; the AoS struct stays 16 B
  and all cache-line math is unchanged.
- **Range is sufficient.** Q16.16 in `int32_t` represents ±32,768 with ~1.5×10⁻⁵ resolution.
  Use `int64_t` for the product accumulator inside the Poisson update (numerator is Q32.32),
  then shift back to Q16.16 at store. Numerical headroom is comfortable for this problem
  (real σ ~10, V ~1 → products ~10¹⁰, well inside int64).
- **Caveat to mention in report:** integer Jacobi converges to a slightly different
  steady-state than float Jacobi (different per-step rounding). That is fine — `cmp` is
  always between *integer* runs, never int-vs-float.

Poisson update template (Q16.16 throughout; SHIFT = 16):

```cpp
int32_t s_c = sigma[idx];
int32_t s_e = (s_c + sigma[idx+1]) >> 1;
int32_t s_w = (s_c + sigma[idx-1]) >> 1;
int32_t s_n = (s_c + sigma[idx+N]) >> 1;
int32_t s_s = (s_c + sigma[idx-N]) >> 1;
int32_t denom = s_e + s_w + s_n + s_s;        // Q16.16

int64_t num = (int64_t)s_e * V[idx+1] + (int64_t)s_w * V[idx-1]
            + (int64_t)s_n * V[idx+N] + (int64_t)s_s * V[idx-N];   // Q32.32

V_next[idx] = (int32_t)(num / denom);          // Q32.32 / Q16.16 = Q16.16
```

---

## 3. Deliverables

1. A git repository with one branch (`main`) and at least 4 tagged commits (see §5).
2. `naive.cpp` and `optimized.cpp` as separate top-level files for the assignment submission
   (copied from the relevant commits — see §6).
3. `run.sh` — build, verify correctness, and profile all 4 stages; produce the ablation table.
4. `viz/` — PPM dumper + `make_video.sh` (ffmpeg) for the side-by-side filament-race video.
5. `report_draft.md` — ≤3-page draft (Matan finalizes to PDF).
6. `prompts.md` — log of AI prompts used.
7. (Provided separately) names + IDs PDF.

---

## 4. Repository structure

Flat layout, single source file evolving across commits — easier to diff between stages than
splitting into many files.

```
cbram-hw1/
├── CBRAM_PLAN_v2.md       (this file)
├── src/
│   ├── cbram.cpp          (the simulator; evolves across commits)
│   ├── physics.h          (shared constants, init routine, Q16.16 helpers)
│   └── io.cpp             (PPM dump, binary state dump — never changes)
├── naive.cpp              (copy of cbram.cpp at tag stage-0, for submission)
├── optimized.cpp          (copy of cbram.cpp at tag stage-3, for submission)
├── run.sh                 (build + verify + profile all stages)
├── viz/
│   ├── make_video.sh
│   └── colormap.h
├── results/               (perf output dumps, gitignored)
├── report_draft.md
└── prompts.md
```

---

## 5. Git commit sequence (THE PLAN)

Each stage = one commit on `main`, tagged `stage-N`. **After each commit, the binary must
build, run, and produce a final-state dump that `cmp`'s bit-identical against `stage-0`.**
If `cmp` fails, do not advance to the next stage — fix or revert.

### Commit / tag `stage-0-naive` — AoS baseline

- `struct Cell { int32_t V, sigma, ion, temp; };` (16 B; 4 cells per 64 B line)
- `Cell* grid; Cell* grid_next;` — double-buffered Jacobi
- Outer simulation loop: K Jacobi iterations + drift-diffusion + σ update per timestep
- Deterministic init (fixed `std::mt19937` seed), Dirichlet top/bottom electrodes, Neumann sides
- PPM frame dumper writes σ every M timesteps to `frames_stage0/`
- Final-state binary dump: `V_final.bin` (just the V field, in row-major Q16.16)
- This is the reference output for `cmp` at every later stage.

### Commit / tag `stage-1-soa` — Struct of Arrays

- Replace `Cell* grid` with four aligned arrays: `int32_t *V, *sigma, *ion, *temp;` plus
  `int32_t* V_next`
- Allocate with `posix_memalign(..., 64, ...)` for cache-line alignment
- **No algorithmic changes.** Same iteration order, same arithmetic, same K per timestep.
- `cmp V_final.bin stage0_V_final.bin` must pass byte-for-byte (this is why we went integer)
- Frames go to `frames_stage1/`; should be visually identical to `frames_stage0/`

### Commit / tag `stage-2-blocked` — Time skewing (temporal blocking)

- Restructure the K-iteration Jacobi loop into trapezoidal spatial-temporal tiles
- Tile parameters: spatial block `B` (chosen so `B² × 8 B × 2 buffers` fits in ~½ L2),
  temporal depth `M` (typically 4–8); choose so `B - 2M > 0`
- One trapezoid: for iteration `m = 0..M-1`, update interior region shrinking by 1 cell per
  side per iteration. Each cell still receives exactly K updates total, each using the
  iteration-(m-1) values of its neighbors — same arithmetic per cell as stage 1.
- Boundary tiles (touching grid edges) handled by truncating the trapezoid; corner tiles
  similarly. This is the fiddly part — write it carefully and test before profiling.
- `cmp` must still pass. If it fails, the trapezoid bookkeeping is wrong; fix it.

### Commit / tag `stage-3-omp` — OpenMP parallelization

- `#pragma omp parallel for collapse(2)` over the outer two loops of the trapezoid scheduler
  (the (jj, ii) tile loop). Within a trapezoid there is no parallelism — keep the inner
  iteration loops serial.
- Build with `-fopenmp`. Default to `OMP_NUM_THREADS=` all physical cores; report scaling at
  T = 1, 2, 4, 8 (capped at the test machine's physical core count, see §8).
- Jacobi has no read/write conflicts (read `V`, write `V_next`), so no atomics or barriers
  needed beyond the implicit barrier at the end of the parallel-for.
- `cmp` must still pass — Jacobi parallelization is deterministic given the same trapezoid
  decomposition.

### Optional `stage-4-extras` (only if time allows, do not block submission on this)

- Huge pages via `madvise(MADV_HUGEPAGE)` for the V/σ arrays.
- First-touch NUMA placement (only matters if test machine is multi-socket; check `lscpu`).
- Mention as future work if not implemented.

---

## 6. Submission artifacts vs commits

The assignment asks for `naive.cpp` and `optimized.cpp` as single files. Resolve by:

```bash
git checkout stage-0-naive -- src/cbram.cpp && cp src/cbram.cpp naive.cpp
git checkout stage-3-omp   -- src/cbram.cpp && cp src/cbram.cpp optimized.cpp
git checkout main
```

`run.sh` automates this so the two top-level files always reflect the right commits.

---

## 7. Correctness verification

**Primary check at every stage:** `cmp V_final.bin stage0_V_final.bin` must pass byte-for-byte.
Because everything is integer arithmetic with identical iteration order per cell, this is
genuinely achievable — not a tolerance approximation. If it fails:
- Stage 1: likely a stale read of `V[]` instead of `V_next[]`, or wrong stride.
- Stage 2: trapezoid bounds wrong, or a cell receiving ≠ K updates.
- Stage 3: race condition (a thread writing where another reads — should not happen with
  Jacobi double-buffering, so if it does the trapezoid decomposition was wrong).

**Secondary check:** final filament image (PPM of σ) should be visually identical across
all stages; include a 4-up frame in the report.

---

## 8. Build & environment

Identical flags across all four stages (stage 3 also needs `-fopenmp`):

```bash
CXXFLAGS="-O2 -march=native -std=c++17 -Wall -Wno-unused-result"
# stage 0–2:
g++ $CXXFLAGS src/cbram.cpp src/io.cpp -o cbram_stageN
# stage 3:
g++ $CXXFLAGS -fopenmp src/cbram.cpp src/io.cpp -o cbram_stage3
```

- **Pin and stabilize:** `taskset` for serial stages; for OpenMP stage, set
  `OMP_PROC_BIND=close` and `OMP_PLACES=cores`.
- **Auto-detect:** `lscpu` → cache sizes + physical core count → log into `results/env.txt`
  so the report cites the real hardware.
- **Same N (grid size) and K (Jacobi iters per timestep) across all stages.** Tune them once
  so stage-0 takes 10–30 s.

**Hard rules** (do not violate, or comparisons are invalid):
- Same compiler flags across stages 0–2 (stage 3 adds `-fopenmp` only).
- Same N, K, timestep count, init seed.
- Time only the simulation kernel; exclude allocation/IO.

---

## 9. Profiling — the ablation study

Run every stage with the same `perf` invocation. This produces the headline table for the
report — a true ablation showing what each layer buys.

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

Stage 3 also gets a **scaling sweep**:

```bash
for T in 1 2 4 8; do
  OMP_NUM_THREADS=$T OMP_PROC_BIND=close OMP_PLACES=cores \
    perf stat -r 5 -e cycles,instructions,LLC-load-misses \
    ./cbram_stage3 4096 2> results/stage3_t${T}.perf
done
```

Cap T at the machine's physical core count (not logical/SMT — SMT threads share L1/L2 and
muddy the scaling story).

**Metrics for the ablation table:**

| Metric | Stage 0 | Stage 1 | Stage 2 | Stage 3 |
|---|---|---|---|---|
| Wall-clock (s) | | | | |
| Speedup vs stage 0 | 1.0× | | | |
| IPC | | | | |
| LLC-load-miss rate | | | | |
| Est. DRAM bytes (LLC-load-misses × 64) | | | | |

Plus a stage-3 scaling plot: threads vs speedup, with the ideal-linear line for reference.

**Honest framing for the report:**
- Stage 1's IPC will jump partly from auto-vectorization on SoA (compiler vectorizes
  contiguous int32 loops more easily). Attribute that effect explicitly; do not claim "pure
  bandwidth." (Mention can run with `-fno-tree-vectorize` to confirm the bandwidth component
  in isolation — optional secondary experiment.)
- Stage 2 should show a sharp drop in LLC-load-misses and estimated DRAM bytes.
- Stage 3 scaling will be sub-linear; that's fine and expected. The interesting fact for the
  report is that scaling is *better* than it would be from stage-0 (which would hit the
  bandwidth wall almost immediately) — i.e., the cache work enabled the thread work. Worth
  proving by also running an OMP version of stage 0 as a comparison data point if time allows.

---

## 10. Visualization

**Per-stage frame dumps** for the filament race:
- Each stage's PPM frames go to `frames_stage{0,1,2,3}/`
- Color-map σ values with viridis (single-header in `viz/colormap.h`)
- Aim for ~100 frames per run, mapped to simulation time (not wall-clock)

**Side-by-side video** (`viz/make_video.sh`):
- Stage 0 vs stage 3, side-by-side via ffmpeg
- Wall-clock timestamp overlay on each panel — viewer sees stage 3's clock running ahead
  while both filaments grow identically
- Output `filament_race.mp4`

**Static figures for the report:**
1. **Memory layout diagram** — one 64 B cache line under AoS (4 cells × 4 fields with V+σ
   highlighted = 50% utilization) vs SoA (16 V's = 100%). Makes stage-1 mechanism instantly
   visible.
2. **Trapezoid time-skewing diagram** — one 2D trapezoid showing 4 iterations of a tile,
   each shrinking by 1 cell per side. Makes stage-2 mechanism visible.
3. **Headline ablation chart** — bar chart of wall-clock for the 4 stages on a log scale.
4. **Thread scaling plot** — stage 3 speedup vs threads, with ideal-linear reference.

---

## 11. Report structure (≤3 pages)

1. **Approach** (½ page) — CBRAM motivation; the layered-bottlenecks insight: each
   optimization exposes the next one.
2. **Naive (stage 0)** (½ page) — what it computes; AoS is the natural physical grouping;
   bandwidth bound with 50% line utilization. Profile numbers.
3. **Optimizations** (1 page) — three stacked subsections (SoA, time skewing, OpenMP) with
   the layout / trapezoid figures. Each subsection: change → HW insight → numbers → was the
   result expected?
4. **Ablation comparison** (½ page) — the headline table + bar chart + scaling plot. Argue
   why each layer's gain matches its mechanism (line size × N, reuse distance ÷ M, threads
   times per-core L2).
5. **What we tried that didn't help** (½ page) — e.g., `__builtin_prefetch` (HW prefetcher
   already optimal), manual SIMD intrinsics (auto-vectorizer already gets most of it on SoA),
   huge pages (small effect on streaming workloads).

---

## 12. Work sequence

For each stage in order: implement → `cmp` against stage-0 dump must pass → `perf stat` →
commit and tag → only then move to next stage. **Do not skip ahead** — without a stage's
correctness gate passing, the next stage's numbers are not trustworthy.

1. Scaffold repo, `physics.h`, `io.cpp` (frame + state dumper), `run.sh` skeleton.
2. Implement stage 0 (AoS naive). Tune N (default 4096) and timestep count for 10–30 s runtime.
   Tag `stage-0-naive`. Save reference output: `cp V_final.bin stage0_V_final.bin`.
3. Implement stage 1 (SoA). `cmp` must pass. Tag `stage-1-soa`.
4. Implement stage 2 (time skewing). `cmp` must pass. Tag `stage-2-blocked`. *This is the
   hardest stage — budget 2–3× the time of the others.*
5. Implement stage 3 (OpenMP). `cmp` must pass. Tag `stage-3-omp`.
6. Implement PPM dumper + ffmpeg side-by-side video script.
7. `run.sh`: full ablation across all 4 stages + thread scaling sweep on stage 3.
8. Draft `report_draft.md` with the ablation table + 4 figures.
9. Log every AI prompt into `prompts.md` continuously, not at the end.

**Two-person split suggestion:**
- One owns physics + correctness (stages 0–1, init, drift-diffusion, σ update, `cmp`)
- The other owns systems (stages 2–3, `perf`, `run.sh`, viz pipeline, scaling experiment)
- Both review the trapezoid math in stage 2 together — it's the bug-prone part.

---

## 13. Open question to resolve before stage 3

**How many physical cores does the test machine have?** OpenMP scaling at 8 cores tells a
much richer story than at 2 cores. If it's a laptop with 4 physical cores, the scaling plot
will be shorter (T = 1, 2, 4) and that's fine. If it's a Technion lab machine with 8–16
cores, plan the scaling sweep accordingly. **Confirm before starting stage 3.**

---

## 14. Pitfalls (the v1 list, refreshed)

- ❌ Skipping the per-stage `cmp` check. Integer arithmetic makes byte-equality genuinely
  achievable — use it as the gate.
- ❌ Different flags across stages. Same `-O2 -march=native` for 0–2; add `-fopenmp` only at 3.
- ❌ Algorithmic changes alongside the layout/scheduling change (e.g., switching to red-black
  in stage 2 because "it parallelizes better"). The report's story is "same algorithm, better
  HW use." Keep it that way.
- ❌ Profiling before `cmp` passes. Untrusted speed of wrong code.
- ❌ Reporting stage-3 scaling with SMT threads enabled. Use physical cores only.
- ❌ Claiming pure-bandwidth speedup for stage 1 when vectorization is partly responsible.
- ❌ Using padding to inflate the AoS struct artificially. Honest physical fields only.
- ❌ Trapezoid off-by-one at grid boundaries causing edge cells to receive ≠ K updates
  (silent correctness break unless `cmp` is checked — which is why `cmp` is the gate).
