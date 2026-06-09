# CBRAM_PLAN.md — HW1: CBRAM Device Simulation, Cache Optimization

> Build spec for Claude Code. Two-person team (Matan + Yuval). Goal: a clean, honest
> AoS vs SoA cache story on a CBRAM device simulator, with bit-identical (or
> tolerance-equivalent) output and large, explainable `perf` deltas.

---

## 0. Summary

A 2D CBRAM device simulator in C++. Each grid cell carries multiple physical fields
(potential V, conductivity σ, ion concentration, temperature, …). The dominant kernel
is an iterative Poisson solve for V on the device cross-section, with simple
drift-diffusion of ions producing a conductive filament. Both versions compute the
**same physics**, differing only in memory layout:

- **Naive:** Array of Structs (AoS) — one `Cell` struct per grid point.
- **Optimized:** Struct of Arrays (SoA) — one flat array per physical field.

---

## 1. The honest hardware mechanism (this is the report's argument)

The Poisson hot loop reads only **two fields per neighbor**: V and σ. Under AoS, each
cache line fetched from DRAM contains 4 cells × 4 fields each = 16 floats, of which
only 8 (the V's and σ's) are used by the kernel. **Cache-line utilization is 50%** —
half the DRAM bandwidth is spent transporting `ion` and `temp` values the kernel
never reads.

Under SoA, a cache line of `V[]` holds 16 consecutive V values, all useful in a row
sweep. Same for `σ[]`. Cache-line utilization is ~100%. The Poisson kernel's DRAM
bandwidth drops by ~2× — and the *aggregate* win across all kernels is bigger because
drift-diffusion uses V+ion, heat-eq uses σ+T, etc., each kernel paying only for its
hot fields.

A secondary effect: SoA is far more SIMD-friendly. With `-O2 -march=native`, the
compiler will auto-vectorize the SoA loops more aggressively than AoS, giving an
additional IPC gain. **The report must attribute both effects honestly** — the layout
change buys (a) better cache-line utilization → fewer DRAM bytes, and (b) better
vectorization → more useful work per cycle. Both are HW-level wins from the same
software change.

For a 4096×4096 grid (256 MB AoS), the working set is far larger than any LLC, so
each sweep streams from DRAM and the bandwidth differential shows up cleanly in
`LLC-load-misses` and wall-clock.

**Originality:** CBRAM simulation maps directly to Matan's research on emerging
non-volatile memory and is essentially absent from architecture coursework. The
*technique* (AoS→SoA) is well-known and explicitly allowed by the assignment; only
reusing class *example programs* is forbidden.

---

## 2. Deliverables (map to assignment)

1. `naive.cpp` — AoS implementation.
2. `optimized.cpp` — SoA implementation, identical I/O.
3. `run.sh` — compile (identical flags), run, verify correctness, profile both.
4. `viz/` — PPM dump utility + `make_video.sh` (ffmpeg) for side-by-side animation.
5. `report_draft.md` — ≤3-page report draft (Matan finalizes to PDF).
6. `prompts.md` — log of AI prompts used (assignment requires this).
7. (Provided separately) names + IDs PDF.

---

## 3. Decision points (resolve before building)

- **Grid size N.** Default 4096 (256 MB AoS, well > LLC). If runtime is too long,
  drop to 2048 (64 MB AoS) — still > most L3s. Tune so naive runtime is 5–15 s.
- **Number of fields per cell.** Default 4 (V, σ, ion, T) → 16B struct, 2× expected
  speedup. If we want a bigger ratio later, add fields (J_x, J_y, defect density,
  filament order parameter) — keep them physically motivated, not padding.
- **Vectorization policy.** Default: same flags for both (`-O2 -march=native`),
  honestly report combined bandwidth+IPC win. Alternative: also run with
  `-fno-tree-vectorize` on both to isolate the pure bandwidth effect — *do this as a
  secondary experiment in the report* if time allows.

---

## 4. Shared design (physics & data, identical across versions)

**Domain.** 2D grid, N×N cells. Top row = active electrode (Dirichlet V = V_applied,
ion source). Bottom row = inert electrode (V = 0). Left/right = Neumann (zero-flux).

**Fields per cell (4 floats, 16 B total):**
- `V` — electric potential
- `sigma` — local conductivity (scales with local ion/atom density)
- `ion` — metal ion concentration in electrolyte
- `temp` — local temperature (initialized uniform; Joule heating optional)

**Initial conditions.** Deterministic. Uniform low σ everywhere except a few seed
defect sites placed by a *fixed-seed PRNG* (e.g., `std::mt19937` with seed 42).
**Both versions must use identical initial state** — write a shared init routine
that fills either layout from the same seed.

**Simulation loop (per outer timestep Δt):**
1. **Poisson solve for V** — Jacobi iteration to convergence (or fixed K iterations
   per timestep, e.g., K=50). *This is the hotspot.*
2. **Drift-diffusion for ions** — explicit Euler update of `ion` using V and σ.
3. **Conductivity update** — `sigma` becomes a smooth function of local ion
   concentration (e.g., logistic / Arrhenius).
4. *(Optional)* **Heat equation for T** — another stencil with Joule source σ|∇V|².
   Skip for v1 to keep code small.
5. **Periodic frame dump** — every M timesteps, write σ as PPM to `frames_naive/`
   or `frames_opt/`.

**Outer loop runs until the filament bridges the electrodes** — detect by checking
if a connected high-σ path reaches the bottom electrode, OR just run a fixed number
of timesteps tuned to produce a clear filament.

**Output for correctness check.** At end of simulation, dump the full `V` field and
`sigma` field as binary files (`V_final_naive.bin`, etc.).

---

## 5. `naive.cpp` (AoS) spec

- `struct Cell { float V, sigma, ion, temp; };`
- Two buffers: `Cell* grid; Cell* grid_next;` swap pointers each Jacobi iteration.
- Poisson hot loop: nested `for (j) for (i)` over interior cells, computing
  `grid_next[idx].V` from neighbors' V and σ via the arithmetic-mean interface
  conductivities (see code sketch in chat). After updating V, carry the other three
  fields forward unchanged.
- All updates pass through the `Cell` struct, so every kernel pays for full struct
  loads even when it only needs 1–2 fields.

## 6. `optimized.cpp` (SoA) spec

- Four separate aligned arrays: `float *V, *sigma, *ion, *temp;` plus a `V_next`
  double-buffer. Allocate with `posix_memalign` to 64 B for line alignment.
- Identical arithmetic per cell, identical iteration order, identical neighbor
  indexing — *only the data layout changes*. Each kernel reads only the field
  arrays it actually needs.
- Do **not** introduce any algorithmic change (no SOR, no red-black Gauss-Seidel,
  no Chebyshev acceleration). The report's story is "same algorithm, better layout"
  and must stay that way.

---

## 7. Correctness verification

**Primary check — bit-identical:**
```bash
cmp V_final_naive.bin V_final_opt.bin && echo "BIT-IDENTICAL OK"
```
This should hold because per-cell arithmetic order is unchanged. If it fails, the
likely culprit is the compiler emitting different FMA / reordered instructions for
the two layouts.

**Fallback — tolerance check:** if bit-identity fails despite identical arithmetic
(can happen with aggressive auto-vectorization), use a max-relative-error check:
```
max |V_naive - V_opt| / max |V_naive|  < 1e-5
```
If this trips, investigate — don't lower the tolerance silently.

**Sanity check:** the final filament images (PPM dumps of σ) should be visually
identical between the two versions — include a side-by-side frame in the report.

---

## 8. Build & environment (fair, reproducible)

```bash
# Identical flags — this is non-negotiable
CXXFLAGS="-O2 -march=native -std=c++17 -Wall"
g++ $CXXFLAGS naive.cpp     -o naive
g++ $CXXFLAGS optimized.cpp -o optimized
```

- Pin to one core: `taskset -c 2 ./naive ...`
- If permitted on the machine: `sudo cpupower frequency-set -g performance` to
  stop turbo from drifting between runs.
- Auto-detect and log cache sizes (`lscpu`, `sysconf(_SC_LEVEL{1,2,3}_CACHE_SIZE)`)
  so the report cites the real cache hierarchy.
- Warm up once, then `perf stat -r 5` to get mean ± stdev.

**Hard rule:** do not make naive `-O0` and optimized `-O3`. The measured win must
come from the layout, not the flags.

---

## 9. Profiling methodology

Drop any unsupported event (`perf list` to check). Run both binaries identically:

```bash
perf stat -r 5 \
  -e cycles,instructions,\
cache-references,cache-misses,\
L1-dcache-loads,L1-dcache-load-misses,\
LLC-loads,LLC-load-misses,\
dTLB-loads,dTLB-load-misses \
  ./naive 4096

perf stat -r 5 ... ./optimized 4096

# confirm hotspot is the Poisson loop
perf record -g ./naive 4096 && perf report --stdio | head -40
```

**Metrics to report:**
- Wall-clock per outer timestep (or total).
- IPC = instructions / cycles (expect SoA ↑ from vectorization).
- LLC-load-miss rate (expect SoA ↓ ~2×).
- Estimated DRAM bytes ≈ `LLC-load-misses × 64 B` (expect SoA ↓ ~2×).
- Optional: compute arithmetic intensity (FLOPs / DRAM byte) for both and locate
  each on a roofline plot.

**Honest caveats to include in the report:**
- `LLC-load-misses` includes prefetcher behavior; treat as a proxy, corroborated
  by wall-clock and IPC.
- If running in a VM/container, document any counters that are unavailable and why.

---

## 10. Visualization

**Per-version frame dumps.**
Each version periodically writes the σ field as a PPM image to `frames_{naive,opt}/`,
mapped to a perceptually-uniform colormap (viridis or just blue→white→red). Aim for
~100 frames per run.

**Side-by-side speedup video** (`make_video.sh`):
- Use ffmpeg to compose `frames_naive/` and `frames_opt/` side-by-side.
- Overlay a wall-clock timer on each panel (use the timestamp recorded when each
  frame was written) so the viewer sees the optimized panel's clock advancing faster
  while both filaments grow identically.
- Output `filament_race.mp4` for the submission (or a GIF if size matters).

**Two static figures for the report:**
1. **Memory layout diagram** (SVG/PNG, drawn once) showing one 64 B cache line under
   AoS (4 cells × 4 fields, with V+σ highlighted = 50% of the line) vs SoA (16 V's
   in one line = 100% used). Makes the mechanism immediately visible to a grader.
2. **Headline bar chart** of wall-clock, IPC, LLC misses, and DRAM bytes — naive vs
   optimized. Or a roofline plot with both versions as points.

Skip cache-miss-as-overlay-on-grid — needs kernel perf integration or a software
cache model, and the cost isn't worth it for this assignment.

---

## 11. Report draft structure (≤3 pages)

1. **Approach.** Why a CBRAM simulator (research-relevant + multi-field physics
   naturally motivates the layout question). The cache-line-utilization insight.
2. **Naive version.** What it computes, why AoS is the natural choice (groups
   physically-related quantities), why it's slow (50% bandwidth wasted on
   non-hot fields). Profiling numbers.
3. **Optimization.** AoS→SoA, the HW insight (line size, working set vs LLC,
   vectorizability). Profiling numbers. Note both effects: bandwidth ↓ and IPC ↑.
4. **Comparison.** Side-by-side table + the memory-layout figure + the filament
   side-by-side frame proving identical physics.
5. **Things tried that didn't help / next steps.** E.g., why naive padding to 64 B
   doesn't fix it, why we didn't add temporal blocking (would have muddied the
   single-cause story).

Front-load the corrected mechanism: cache-line utilization, not "every access misses."

---

## 12. Work sequence (suggested for the agent)

1. Scaffold repo with both source files, shared header for physics constants, a
   shared init routine, `run.sh`, `report_draft.md`, `prompts.md`.
2. Implement physics + Poisson + drift-diffusion in `naive.cpp`. Verify it runs,
   dumps frames, produces a recognizable filament.
3. Port to SoA in `optimized.cpp`. Run `cmp` on final state dumps — must pass.
4. Tune N and timestep count so naive runtime is 5–15 s and the filament is clear.
5. Implement PPM writer + `make_video.sh`. Generate the side-by-side video.
6. Write `run.sh` (build → run → cmp → perf stat both → perf record).
7. Collect numbers; draft `report_draft.md`; produce the memory-layout figure.
8. Optional secondary experiment: `-fno-tree-vectorize` runs to isolate pure
   bandwidth effect; report both numbers.
9. Log every prompt into `prompts.md`.

---

## 13. Pitfalls to avoid

- ❌ Different compiler flags for the two versions (confounded comparison).
- ❌ Algorithmic changes alongside the layout change (SOR, red-black, etc.) — muddies
  the story. Same algorithm, same iteration order, only layout differs.
- ❌ Claiming "every access is a cache miss" — wrong. The story is cache-line
  utilization and bandwidth efficiency.
- ❌ Padding the AoS struct artificially to make AoS look worse. Use honest,
  physically-motivated fields only.
- ❌ Timing allocation/IO instead of just the simulation kernel.
- ❌ Trusting `perf` numbers before `cmp` (or tolerance check) passes.
- ❌ Reusing an example program shown in class. AoS→SoA as a *technique* is fine;
  applying it to a CBRAM sim is novel.
