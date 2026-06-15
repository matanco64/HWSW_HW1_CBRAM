// ============================================================
// HW1 — CBRAM Device Simulation: Profiling & Optimization
// Build:  typst compile report/cbram_report.typ      (VSCode/tinymist: just open it)
//   `results/` is a symlink to ../results so flame-graph SVGs resolve inside the
//   Typst sandbox with no --root flag.
//
// All measured numbers live in the `data` block below — to refresh from a new
// profiling run, edit that ONE block and recompile. Current numbers are the
// course-server run (single-core Xeon E5-2630 v3, N=6144 / 80 steps).
// ============================================================

#set document(title: "CBRAM Simulation — Profiling & Optimization")
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
#let s0 = (name: "Stage 0 — AoS (baseline)",     time: 358.43, ipc: 1.835, l1m: 18.91, cmiss: 2.844, dtlb: 0.17, instr: 1508, cyc: 822)
#let s1 = (name: "Stage 1 — SoA",                time: 243.58, ipc: 2.436, l1m: 7.25,  cmiss: 0.436, dtlb: 0.04, instr: 1343, cyc: 551)
#let s2 = (name: "Stage 2 — SoA + time-skewing", time: 281.54, ipc: 2.881, l1m: 4.26,  cmiss: 0.087, dtlb: 0.01, instr: 1878, cyc: 652)

#let sp(a, b) = calc.round(a / b, digits: 2)        // ratio helper
#let pct(x)   = calc.round(x * 100, digits: 0)

// ── Title ───────────────────────────────────────────────────────────────────
#align(center)[
  #text(15pt, weight: "bold")[Profiling and Optimizing a CBRAM Filament-Growth Simulation]
  #v(0.2em)
  #text(9.7pt)[
    Matan Cohen (ID: <ID>) #h(1.2em) · #h(1.2em) Yuval Kogan (ID: <ID>) \
    HW/SW Co-Design — HW1 · #datetime.today().display("[month repr:long] [year]")
  ]
]
#v(0.15em)
#line(length: 100%, stroke: 0.4pt)

= Problem and approach

We simulate conductive-bridge RAM (CBRAM) forming, the physical event that
stores a bit: a metallic filament grows across the dielectric until it bridges
the two electrodes. We model it with the *Dielectric Breakdown Model* (DBM) on an
$N times N$ grid. Each *growth step* (i) warm-starts a Jacobi solve of the Laplace
potential $nabla^2 V = 0$ with the existing filament pinned to the cathode
($V=0$) and the anode held at $V_"app"$, then (ii) adds one filament-adjacent
empty cell, chosen at random with probability $prop V^eta$ ($eta = 3$). The
field concentrates at the filament tip, so growth is self-reinforcing and a
narrow branched filament emerges — the morphology CBRAM exhibits.

The simulation is overwhelmingly dominated by step (i): every growth step runs
`JACOBI_ITERS = 30` double-buffered stencil sweeps over the whole grid, and there
are $O(N^2)$ growth steps. This makes the *Jacobi sweep a textbook
memory-bandwidth-bound stencil* — exactly the kernel where data-layout and
cache-locality optimizations are decisive, which is why we chose it.

*Methodology.* We keep one source file per stage (`dbm_stage0..2.cpp`), identical
compiler flags (`-O2 -march=native -std=c++17`), and a strict correctness gate:
all arithmetic is fixed-point `int32` (Q16.16) with a deterministic xorshift64
PRNG, so every stage must produce a *bit-identical* `V_final` (verified with
`cmp`). Optimizations therefore cannot "cheat" by changing the numerics. We
profile at $N = 6144$ — the working set is $approx 19 times$ the 16 MB L3, so the
full memory hierarchy is genuinely exercised — capping each run at 80 growth
steps so every stage does identical work, with `perf stat -r 3` plus `perf
record` flame graphs (Fig. 1).

#block(fill: luma(245), inset: 6pt, radius: 3pt, width: 100%)[
  *Measurement note.* Numbers are from the #machine. Being a single core with no
  hybrid-PMU split or thermal throttling, runs are highly repeatable ($lt.eq 1.3%$
  spread over 3 runs). We cross-checked every trend on a laptop (Core Ultra 7
  155H, 24 MB L3); the relationships between stages were identical.
]

= Baseline: Array-of-Structs (Stage 0)

The natural first implementation stores the grid as an array of structs,
`struct Cell { int32_t V, metal; }` — one object per cell holding both the
potential and the filament flag. It is the obvious, readable layout and we
deliberately left it unoptimized as the reference.

*Why it is slow.* The Jacobi sweep reads only the four neighbours' `V`, but in
AoS each `Cell` is 8 bytes, so a 64-byte cache line holds 8 cells and *half of
every line fetched is the unused `metal` field* — cache-line utilization $approx
50%$. With the grid $approx 19 times$ the L3 that wasted half *doubles the DRAM
traffic* the whole hierarchy must move. Two further taxes are visible in the flame
graph (Fig. 1a): `apply_boundary` — the $O(N^2)$ `pin_cluster` scan re-pinning
every metal cell to $V=0$ on each of the 30 sweeps — is *27%* of the run, and the
double-buffer `grid_next = grid` copy (`copy_grid`) duplicates `metal` every sweep
even though only `V` changes — a fat memmove at *20%*.

*Baseline profile.* IPC just *#s0.ipc*, an L1-dcache miss rate of *#s0.l1m%*,
*#s0.cmiss billion* last-level cache misses streamed from DRAM, and a runtime of
*#s0.time s*. The machine spends its time waiting on memory, not computing — this
is the bottleneck Stage 1 attacks.

= Optimization: Struct-of-Arrays (Stage 1)

*The change.* Split the AoS grid into two flat arrays: a contiguous
`int32_t V[]` and a separate `uint8_t metal[]`. The hardware insight is direct:
the Jacobi sweep now streams `V[]` densely — *16 useful values per 64-byte line
(~100% utilization)* instead of 8 — so it moves half as many bytes for the same
work. Two corollaries fall out for free: (1) `metal[]` is *fixed* for all 30
sweeps of a step, so it is never double-buffered — the per-sweep `metal` memmove
*disappears entirely* (visible in Fig. 1b: the `copy_grid` block is simply gone);
(2) `metal[]` is `uint8_t`, so `pin_cluster` scans 1 byte/cell, an 8#sym.times
denser scan.

*Result.* L1-miss falls *#s0.l1m% #sym.arrow #s1.l1m%* (#pct(1 - s1.l1m/s0.l1m)% fewer),
last-level misses drop *#s0.cmiss B #sym.arrow #s1.cmiss B* — #sp(s0.cmiss, s1.cmiss)#sym.times
less DRAM traffic — and with the stalls gone IPC rises *#s0.ipc #sym.arrow #s1.ipc
(+#pct(s1.ipc/s0.ipc - 1)%)*. Deleting the copy even cuts instructions
(#s0.instr B #sym.arrow #s1.instr B). Runtime improves *#s0.time s #sym.arrow
#s1.time s — a #sp(s0.time, s1.time)#sym.times speed-up*, bit-identical to the
baseline. This is the clean win: a pure data-layout change that lifts the
bandwidth ceiling without touching the math — the single most effective
optimization in the project.

= Stage 2 — temporal blocking: a memory win that lost on compute

This is the most instructive part of the project, so we document it in full.

*Hypothesis.* The 30 sweeps per step re-stream the entire $N^2$ grid 30 times;
with the grid far larger than the L3 that is $30 times$ the minimum DRAM traffic.
*Temporal blocking* (time-skewing) should fix it: load a `TILE`#sym.times`TILE`
tile plus a halo into a cache-resident scratch buffer, march it `T_BLOCK` sweeps
locally (the valid region shrinks one cell per side per sweep), write it back —
cutting global passes from $30$ to $30 \/ T_"block"$. On paper, a large win.

*What actually happened.* Time-skewing did *exactly* what it promised on memory:
Stage 2 has the *best memory behaviour of any stage* — L1-miss #s2.l1m%, just
#s2.cmiss B last-level misses (#sp(s1.cmiss, s2.cmiss)#sym.times fewer than SoA,
#sp(s0.cmiss, s2.cmiss)#sym.times fewer than the baseline), and the *highest* IPC
(#s2.ipc). The restructured tile-marcher (`block_advance`, Fig. 1c) is genuinely
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
  (+#pct(s2.instr/s1.instr - 1)%) — the halo recomputation that makes the most
  memory-efficient stage the slower one.
]]

*Root cause.* Each tile recomputes a halo of redundant border cells every chunk;
for our `TILE`/`T_BLOCK` that overhead *adds #pct(s2.instr/s1.instr - 1)% instructions*
(#s1.instr B #sym.arrow #s2.instr B). The higher IPC absorbs some of it — cycles
rise only +#pct(s2.cyc/s1.cyc - 1)% — but not enough, leaving Stage 2
#pct(s2.time/s1.time - 1)% slower than SoA. Time-skewing converted a memory-bound
kernel into a lighter compute-bound one, and the compute it *added* slightly
outweighed the bandwidth it *saved*.

*The lesson.* An optimization that *succeeds at its stated target* can still be a
net loss. Time-skewing is unambiguously the most cache-friendly stage; it is also
slower, because a net win needs (a) the resource it saves — DRAM bandwidth — to be
the *binding* constraint, *and* (b) the price paid elsewhere — instructions here —
to be *smaller* than the saving. Both must hold; for us only (a) did, and only
just. Where bandwidth is genuinely the binding constraint — many cores sharing and
saturating one memory bus — the same code should tip to a win; that multicore
regime is the natural next step, left as future work.

= Comparison and conclusions

SoA (Stage 1) is the unambiguous win — *#sp(s0.time, s1.time)#sym.times* faster
than the baseline from a layout change alone, by raising cache-line utilization
from $approx$50% to $approx$100% and deleting a per-sweep copy, cutting DRAM
misses #sp(s0.cmiss, s1.cmiss)#sym.times and lifting IPC #pct(s1.ipc/s0.ipc - 1)%,
all while staying bit-identical. Temporal blocking (Stage 2) is the more subtle and
more valuable lesson: it is the *most cache-efficient* stage yet
*#sp(s2.time, s1.time)#sym.times slower than SoA*, because it traded a bandwidth
bottleneck it had already won for a compute bottleneck it then lost — +#pct(s2.instr/s1.instr - 1)%
instructions for #sp(s1.cmiss, s2.cmiss)#sym.times less traffic that was no longer
the limiter. The takeaway is the core of HW/SW co-design: profile to find the
*binding* constraint, and spend complexity only where the resource saved is the
one actually limiting you, for less than it costs elsewhere. A multicore run,
where shared-bus bandwidth becomes binding, is where we expect time-skewing's
memory advantage to finally turn into wall-clock.

// ── Flame-graph figure ──────────────────────────────────────────────────────
#v(0.4em)
#figure(
  kind: image, supplement: [Figure],
  stack(dir: ttb, spacing: 5pt,
    align(left)[#text(8.4pt)[*(a) Stage 0 — AoS:* `solve_potential` 98%, split into the `copy_grid` memmove (20%) and `apply_boundary`/`pin_cluster` scan (27%).]],
    image("results/stage0_flamegraph.svg", width: 12.8cm),
    align(left)[#text(8.4pt)[*(b) Stage 1 — SoA:* the `copy_grid` block is gone; `solve_potential` is now almost pure stencil.]],
    image("results/stage1_flamegraph.svg", width: 12.8cm),
    align(left)[#text(8.4pt)[*(c) Stage 2 — time-skewing:* restructured into the cache-resident tile marcher `block_advance`/`advance_tile`.]],
    image("results/stage2_flamegraph.svg", width: 12.8cm),
  ),
  caption: [Self-time flame graphs (`perf record -F 999`, $N=6144$). The layout
    change (a→b) erases the memmove; time-skewing (c) reorganizes the whole solve
    around tiles.],
)
