// ============================================================
// HW1 - CBRAM Device Simulation: Profiling & Optimization
// Build:  typst compile report/cbram_report.typ      (VSCode/tinymist: just open it)
//   `results/` is a symlink to ../results so flame-graph SVGs resolve inside the
//   Typst sandbox with no --root flag.
//
// All measured numbers live in the `data` block below - to refresh from a new
// profiling run, edit that ONE block and recompile. Current numbers are the
// course-server run (single-core Xeon E5-2630 v3, N=6144 / 80 steps).
// ============================================================

#set document(title: "CBRAM Simulation - Profiling & Optimization")
#set page(paper: "a4", margin: (x: 1.7cm, y: 1.5cm), numbering: "1 / 1")
#set text(font: "New Computer Modern", size: 9.7pt)
#set par(justify: true, leading: 0.55em)
#show heading.where(level: 1): set text(size: 11.5pt)
#show heading: set block(above: 0.85em, below: 0.45em)
#set heading(numbering: "1.")
#show raw: set text(font: "DejaVu Sans Mono", size: 8.6pt)

// ── DATA (edit here to refresh from a new run) ──────────────────────────────
// machine: where the numbers were taken. N = 6144, 80 growth steps.
// l1m = L1-dcache miss %, cmiss = cache (last-level) misses in billions,
// dtlb = dTLB-load miss %, instr/cyc = retired instructions / cycles (billions).
#let machine = "course server · 1-core Intel Xeon E5-2630 v3 @ 2.40 GHz · L1d 32 KiB, L2 4 MiB, L3 16 MiB · N = 6144, 80 growth steps (working set ≈ 302 MB, ≈ 19× the L3)"
#let s0 = (name: "Stage 0 - AoS (baseline)",     time: 375.82, ipc: 1.746, l1m: 19.00, cmiss: 3.012, dtlb: 0.18, instr: 1508, cyc: 864)
#let s1 = (name: "Stage 1 - SoA",                time: 245.88, ipc: 2.359, l1m: 7.30,  cmiss: 0.462, dtlb: 0.05, instr: 1343, cyc: 569)
#let s2 = (name: "Stage 2 - SoA + time-skewing", time: 284.92, ipc: 2.845, l1m: 4.25,  cmiss: 0.162, dtlb: 0.01, instr: 1881, cyc: 661)
#let s3 = (name: "Stage 3 - SoA + skewing + SIMD", time: 83.53, ipc: 2.110, l1m: 15.24, cmiss: 0.153, dtlb: 0.06, instr: 408, cyc: 193)

#let sp(a, b) = calc.round(a / b, digits: 2)        // ratio helper
#let pct(x)   = calc.round(x * 100, digits: 0)

// ── Title ───────────────────────────────────────────────────────────────────
#align(center)[
  #text(15pt, weight: "bold")[Profiling and Optimizing a CBRAM Filament-Growth Simulation]
  #v(0.2em)
  #text(9.7pt)[
    Matan Cohen (ID: #sys.inputs.at("matan-id", default: "<ID>")) #h(1.2em) · #h(1.2em) Yuval Kogan (ID: #sys.inputs.at("yuval-id", default: "<ID>")) \
    HW/SW Co-Design - HW1 · #datetime.today().display("[month repr:long] [year]") \
    Code & data: #link("https://github.com/matanco64/HWSW_HW1_CBRAM")[`github.com/matanco64/HWSW_HW1_CBRAM`]
  ]
]
#v(0.15em)
#line(length: 100%, stroke: 0.4pt)

// ── Hero figure (a still frame extracted from comparison_0vs1.mp4) ───────────
// A real-time race frame: the instant SoA has bridged but AoS has not.
// To re-render this PNG after regenerating the video, run:  report/regen_assets.sh
#v(0.2em)
#figure(
  kind: image, supplement: [Figure],
  image("assets/filament_0vs1.png", width: 14.5cm),
  caption: [A single real-time instant from a side-by-side race ($eta=3$,
    $N=1024$). *Top:* Stage 0 (AoS); *bottom:* Stage 1 (SoA) - each row shows the
    filament $sigma$, the potential $phi$, and the field $|nabla phi|$. Because the
    stages are *bit-identical* they trace the same growth sequence, so this is a
    fair race: at this same wall-clock moment SoA has already *bridged* (step 3209)
    while AoS is only halfway (step 1800). The on-screen gap *is* the speed-up.],
)
#v(0.2em)

= Problem and approach

We simulate conductive-bridge RAM (CBRAM) forming, a metallic filament grows across the dielectric until it bridges
the two electrodes. CBRAM is one physical realization of a *memristor* - a
two-terminal device whose resistance encodes its history of current flow. We model it with the *Dielectric Breakdown Model* (DBM) on an
$N times N$ grid (1024 for visualization and 6144 for simulation). Each *growth step*:
#set enum(numbering: "i)")
+ warm-starts a Jacobi solve of the Laplace potential $nabla^2 V = 0$ with the existing filament pinned to the cathode ($V=0$) and the anode held at $V_"app"$, then 
+ adds one filament-adjacent empty cell, chosen at random with probability $prop V^eta$ ($eta = 3$). The field concentrates at the filament tip, so growth is self-reinforcing and a narrow branched filament emerges - the morphology CBRAM exhibits (Fig. 1).


The simulation is overwhelmingly dominated by step (i): every growth step runs
`JACOBI_ITERS = 30` double-buffered stencil sweeps over the whole grid, and there
are $O(N^2)$ growth steps. This makes the *Jacobi sweep a low-arithmetic-intensity,
memory-bound stencil*: it moves far more bytes than it performs arithmetic, so what
limits it is *how much data it touches and how close that data sits* - not raw
compute. That is exactly the kernel where data-layout and cache-locality
optimizations are decisive, which is why we chose it.

*Methodology.*
- *One source file per stage* (`dbm_stage0..2.cpp`), identical compiler flags (`-O2 -march=native -std=c++17`).
- *Correctness gate:* all arithmetic is fixed-point `int32` (Q16.16) with a deterministic xorshift64 PRNG; every stage must produce a *bit-identical* `V_final` (verified with `cmp`) - optimizations cannot "cheat" by changing the numerics.
- *Grid size $N = 6144$:* working set $approx 19 times$ the 16 MB L3, so the full memory hierarchy is genuinely exercised.
- *Identical workload:* each run capped at 80 growth steps; profiled with `perf stat -r 3` plus `perf record` flame graphs (Fig. 2).

#block(fill: luma(245), inset: 6pt, radius: 3pt, width: 100%)[
  *Measurement note.* Numbers are from the #machine. runs are highly repeatable ($lt.eq 1.3%$
  spread over 3 runs). The virtualized PMU exposes cycles, instructions, and
  cache/L1/dTLB counters; LLC and top-down events are unavailable, so DRAM traffic
  is read from the last-level `cache-misses` event.
]

= Baseline: Array-of-Structs (Stage 0)

The natural first implementation stores the grid as an array of structs,
`struct Cell { int32_t V, metal; }` - one object per cell holding both the
potential and the filament flag. It is the obvious, readable layout and we
deliberately left it unoptimized as the reference.

*Why it is slow.* The Jacobi sweep reads only the four neighbours' `V`, but in
AoS each `Cell` is 8 bytes, so a 64-byte cache line holds 8 cells and *half of
every line fetched is the unused `metal` field* - cache-line utilization $approx
50%$. With the grid $approx 19 times$ the L3 that wasted half *doubles the DRAM
traffic* the whole hierarchy must move. Two further taxes are visible in the flame
graph (Fig. 2a): `apply_boundary` - the $O(N^2)$ `pin_cluster` scan re-pinning
every metal cell to $V=0$ on each of the 30 sweeps - is *27%* of the run, and the
double-buffer `grid_next = grid` copy (`copy_grid`) duplicates `metal` every sweep
even though only `V` changes - a fat memmove at *20%*.

*Baseline profile.* IPC just *#s0.ipc*, an L1-dcache miss rate of *#s0.l1m%*,
*#s0.cmiss billion* last-level cache misses streamed from DRAM, and a runtime of
*#s0.time s*. The machine spends its time waiting on memory, not computing - this
is the bottleneck Stage 1 attacks.

= Optimization: Struct-of-Arrays (Stage 1)

*The change.* Split the AoS grid into two flat arrays: a contiguous
`int32_t V[]` and a separate `uint8_t metal[]`. The hardware insight is direct:
the Jacobi sweep now streams `V[]` densely - *16 useful values per 64-byte line
(~100% utilization)* instead of 8 - so it moves half as many bytes for the same
work. Two corollaries fall out for free: (1) `metal[]` is *fixed* for all 30
sweeps of a step, so it is never double-buffered - the per-sweep `metal` memmove
*disappears entirely* (visible in Fig. 2b: the `copy_grid` block is simply gone);
(2) `metal[]` is `uint8_t`, so `pin_cluster` scans 1 byte/cell, an 8#sym.times
denser scan.

*Result.* L1-miss falls *#s0.l1m% #sym.arrow #s1.l1m%* (#pct(1 - s1.l1m/s0.l1m)% fewer),
last-level misses drop *#s0.cmiss B #sym.arrow #s1.cmiss B* - #sp(s0.cmiss, s1.cmiss)#sym.times
less DRAM traffic - and with the stalls gone IPC rises *#s0.ipc #sym.arrow #s1.ipc
(+#pct(s1.ipc/s0.ipc - 1)%)*. Deleting the copy even cuts instructions
(#s0.instr B #sym.arrow #s1.instr B). Runtime improves *#s0.time s #sym.arrow
#s1.time s - a #sp(s0.time, s1.time)#sym.times speed-up*, bit-identical to the
baseline. This is the clean win: a pure data-layout change that halves the bytes
moved per sweep without touching the math - the single most effective
optimization in the project.

= Stage 2 - temporal blocking: a memory win that lost on compute

This is the most instructive part of the project, so we document it in full.

*Hypothesis.* The 30 sweeps per step re-stream the entire $N^2$ grid 30 times;
with the grid far larger than the L3 that is $30 times$ the minimum DRAM traffic.
*Temporal blocking* (time-skewing) should fix it: load a `TILE`#sym.times`TILE`
tile plus a halo into a cache-resident scratch buffer, march it `T_BLOCK` sweeps
locally (the valid region shrinks one cell per side per sweep), write it back -
cutting global passes from $30$ to $30 \/ T_"block"$. On paper, a large win.

*What actually happened.* Time-skewing did *exactly* what it promised on memory:
Stage 2 has the *best memory behaviour of any stage* - L1-miss #s2.l1m%, just
#s2.cmiss B last-level misses (#sp(s1.cmiss, s2.cmiss)#sym.times fewer than SoA,
#sp(s0.cmiss, s2.cmiss)#sym.times fewer than the baseline), and the *highest* IPC
(#s2.ipc). The restructured tile-marcher (`block_advance`, Fig. 2c) is genuinely
cache-resident. And yet it is *#sp(s2.time, s1.time)#sym.times slower than plain
SoA* (#s2.time s vs #s1.time s). The counters say why exactly:

#align(center)[#table(
  columns: (auto, 1fr, 1fr, 1fr, 1fr, 1fr),
  align: (left, center, center, center, center, center),
  inset: (x: 6pt, y: 3pt), stroke: 0.4pt + luma(170),
  table.header([*Stage*], [*Time (s)*], [*vs base*], [*IPC*], [*L1-miss %*], [*Cache-miss (B)*]),
  [#s0.name], [#s0.time], [1.00#sym.times], [#s0.ipc], [#s0.l1m], [#s0.cmiss],
  [#s1.name], [#s1.time], [#sp(s0.time, s1.time)#sym.times], [#s1.ipc], [#s1.l1m], [#s1.cmiss],
  [#s2.name], [#s2.time], [#sp(s0.time, s2.time)#sym.times], [#s2.ipc], [#s2.l1m], [#s2.cmiss],
)]
#v(-0.2em)
#align(center)[#text(8.3pt, style: "italic")[
  Instructions retired: Stage 1 #s1.instr B #sym.arrow Stage 2 #s2.instr B
  (+#pct(s2.instr/s1.instr - 1)%) - the halo recomputation that makes the most
  memory-efficient stage the slower one.
]]

*Root cause.* Each tile recomputes a halo of redundant border cells every chunk;
for our `TILE`/`T_BLOCK` that overhead *adds #pct(s2.instr/s1.instr - 1)% instructions*
(#s1.instr B #sym.arrow #s2.instr B). The higher IPC absorbs some of it - cycles
rise only +#pct(s2.cyc/s1.cyc - 1)% - but not enough, leaving Stage 2
#pct(s2.time/s1.time - 1)% slower than SoA. Time-skewing converted a memory-bound
kernel into a lighter compute-bound one, and the compute it *added* outweighed the
traffic it *saved* - traffic that, single-core, was no longer the limiter anyway.

*The lesson.* An optimization that *succeeds at its stated target* can still be a
net loss. Time-skewing is unambiguously the most cache-friendly stage; it is also
slower, because a net win needs (a) the resource it saves - here, memory traffic -
to be what is actually *limiting* the kernel, *and* (b) the price paid elsewhere -
the extra instructions - to be *smaller* than the saving. Neither held for us:
Stage 1 had already pulled the working set close enough that cutting traffic
further bought little, while the halo recomputation was paid in full. The lesson
also points the way out: if the loss is *added compute*, make that compute cheap
enough to pay back by vectorizing the now cache-resident stencil. *We did exactly
this* (Stage 3, Appendix A): SIMD over the L2-resident tile collapses instructions
*#sp(s2.instr, s3.instr)#sym.times* and runtime to *#s3.time s* - *#sp(s2.time, s3.time)#sym.times
faster than Stage 2* and *#sp(s1.time, s3.time)#sym.times faster than SoA*,
bit-identical. Time-skewing wins once its compute is vectorized away.

= Comparison and conclusions

SoA (Stage 1) is the unambiguous win - *#sp(s0.time, s1.time)#sym.times* faster
than the baseline from a layout change alone, by raising cache-line utilization
from $approx$50% to $approx$100% and deleting a per-sweep copy, cutting DRAM
misses #sp(s0.cmiss, s1.cmiss)#sym.times and lifting IPC #pct(s1.ipc/s0.ipc - 1)%,
all while staying bit-identical. Temporal blocking (Stage 2) is the more subtle and
more valuable lesson: it is the *most cache-efficient* stage yet
*#sp(s2.time, s1.time)#sym.times slower than SoA*, because it traded a memory
bottleneck SoA had already neutralized for a compute bottleneck it then lost - +#pct(s2.instr/s1.instr - 1)%
instructions for #sp(s1.cmiss, s2.cmiss)#sym.times less traffic that was no longer
the limiter. The takeaway is the core of HW/SW co-design: profile to find the
*binding* constraint, and spend complexity only where the resource saved is the
one actually limiting you, for less than it costs elsewhere. Time-skewing's
cache-residency was real headroom, and *Stage 3 cashes it in*: vectorizing the
now-cache-resident stencil makes the added compute pay for itself, turning the
slowest "optimization" into the fastest stage overall - *#s3.time s,
#sp(s0.time, s3.time)#sym.times over the baseline* and #sp(s1.time, s3.time)#sym.times
over SoA, still bit-identical (full analysis in Appendix A). The remaining headroom,
making memory traffic bind again across many cores on one bus, is the natural next step.

// ── Flame-graph figure ──────────────────────────────────────────────────────
#v(0.4em)
#figure(
  kind: image, supplement: [Figure],
  stack(dir: ttb, spacing: 5pt,
    align(left)[#text(8.4pt)[*(a) Stage 0 - AoS:* `solve_potential` 98%, split into the `copy_grid` memmove (20%) and `apply_boundary`/`pin_cluster` scan (27%).]],
    image("results/stage0_flamegraph.svg", width: 12.8cm),
    align(left)[#text(8.4pt)[*(b) Stage 1 - SoA:* the `copy_grid` block is gone; `solve_potential` is now almost pure stencil.]],
    image("results/stage1_flamegraph.svg", width: 12.8cm),
    align(left)[#text(8.4pt)[*(c) Stage 2 - time-skewing:* restructured into the cache-resident tile marcher `block_advance`/`advance_tile`.]],
    image("results/stage2_flamegraph.svg", width: 12.8cm),
  ),
  caption: [Self-time flame graphs (`perf record -F 999`, $N=6144$). The layout
    change (a→b) erases the memmove; time-skewing (c) reorganizes the whole solve
    around tiles.],
)

= Reproducing these results
Every number and figure above regenerates from the bundled sources (file map in
`README.txt`). From the submission root:
- `./build.sh` - compiles `cbram_stage0/1/2` with *identical* flags (`-O2 -march=native -std=c++17`) into `build/`.
- `./run.sh 6144 80` - the exact pipeline used here: build, log the CPU to `env.txt`, run all three stages at $N=6144$ / 80 steps, *gate* each `V_final` bit-identical to Stage 0 with `cmp`, then `perf stat -r 3` (#sym.arrow `results/stage{0,1,2}.perf`) and `perf record` flame graphs (#sym.arrow `results/stage*_flamegraph.svg`, Fig. 2).
- One stage + correctness check: `./build/cbram_stage1 6144 -s 80`, then `cmp build/V_final_stage1.bin build/V_final_stage0.bin` (silent = bit-identical).

Figure 1 is a frame of `comparison_0vs1.mp4`; the numbers above are from the #machine.
