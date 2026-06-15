// ============================================================
// HW1 — CBRAM Device Simulation: Profiling & Optimization
// Build:  typst compile report/cbram_report.typ
//
// All measured numbers live in the `data` block below so that when the
// clean university-server run lands we update ONE place, recompile, done.
// Current numbers are LAPTOP placeholders — see the `machine` field.
// ============================================================

#set document(title: "CBRAM Simulation — Profiling & Optimization")
#set page(paper: "a4", margin: (x: 1.7cm, y: 1.6cm), numbering: "1 / 1")
#set text(font: "New Computer Modern", size: 9.7pt)
#set par(justify: true, leading: 0.56em)
#show heading.where(level: 1): set text(size: 11.5pt)
#show heading.where(level: 2): set text(size: 10.2pt)
#show heading: set block(above: 0.9em, below: 0.5em)
#set heading(numbering: "1.")
#show raw: set text(font: "DejaVu Sans Mono", size: 8.6pt)

// ── DATA (edit here when the university numbers arrive) ──────────────────────
#let machine = "PLACEHOLDER — laptop, Intel Core Ultra 7 155H, N=500 full bridge"
#let s0 = (name: "Stage 0 — AoS (baseline)",     time: 33.939, ipc: 2.194, l1m: 15.735, cm: 0.013, instr: none,   cyc: none)
#let s1 = (name: "Stage 1 — SoA",                time: 27.125, ipc: 2.528, l1m: 3.739,  cm: 0.082, instr: 157.83, cyc: 62.44)
#let s2 = (name: "Stage 2 — SoA + time-skewing", time: 40.719, ipc: 2.349, l1m: 5.662,  cm: 0.034, instr: 222.41, cyc: 94.69)
#let conc = (n: 3072, inst: 8, s1: 15.5, s2: 9.35)   // 8 concurrent instances

#let sp(a, b) = calc.round(a / b, digits: 2)   // speedup helper
#let fmt(x) = if x == none { "—" } else { str(x) }

// ── Title ───────────────────────────────────────────────────────────────────
#align(center)[
  #text(15pt, weight: "bold")[Profiling and Optimizing a CBRAM Filament-Growth Simulation]
  #v(0.2em)
  #text(9.7pt)[
    Matan Cohen (ID: #underline[#h(2.2em)]) #h(1.2em) · #h(1.2em) Yuval #underline[#h(3em)] (ID: #underline[#h(2.2em)]) \
    HW/SW Co-Design — HW1 · #datetime.today().display("[month repr:long] [year]")
  ]
]
#v(0.2em)
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
profile with `perf stat -r 3` (cycles, instructions, cache and L1/dTLB misses)
plus `perf record` flame graphs to locate hot phases.

#block(fill: luma(245), inset: 6pt, radius: 3pt, width: 100%)[
  *Note — placeholder data.* All figures below are from #machine. The grid is
  small enough that its working set fits in cache; the final report will use the
  university server's large-$N$ run, which exercises the full memory hierarchy.
  The _relationships_ between stages are already representative.
]

= Baseline: (Array-of-Structs) AoS (Stage 0)

The natural first implementation stores the grid as an array of structs,
`struct Cell { int32_t V, metal; }` — one object per cell holding both the
potential and the filament flag. It is the obvious, readable layout and we
deliberately left it unoptimized as the reference.

*Why it is slow.* The Jacobi sweep reads only the four neighbours' `V`, but in
AoS each `Cell` is 8 bytes, so a 64-byte cache line holds 8 cells and *half of
every line fetched is the unused `metal` field* — cache-line utilization $approx
50%$. Two further taxes show up in the flame graph: (a) `pin_cluster`, the
$O(N^2)$ scan that re-pins every metal cell to $V=0$ on *each* of the 30 sweeps,
drags a whole 8-byte `Cell` through cache just to test one flag; and (b) the
double-buffer `grid_next = grid` copy duplicates `metal` every sweep even though
only `V` changes — a fat memmove (#sym.tilde 20% of the step in the profile).

*Baseline profile.* IPC #s0.ipc, L1-dcache miss rate *#s0.l1m%*, runtime
*#s0.time s*. The high L1 miss rate is the signature of the wasted half-lines —
this is the bottleneck Stage 1 attacks.

= Optimization: Struct-of-Arrays (Stage 1)

*The change.* Split the AoS grid into two flat arrays: a contiguous
`int32_t V[]` and a separate `uint8_t metal[]`. The hardware insight is direct:
the Jacobi sweep now streams `V[]` densely — *16 useful values per 64-byte line
(~100% utilization)* instead of 8 — so it fetches half as many lines for the same
work. Two corollaries fall out for free: (1) `metal[]` is *fixed* for all 30
sweeps of a step, so it is never double-buffered — the per-sweep `metal` memmove
*disappears entirely*; (2) `metal[]` is `uint8_t`, so `pin_cluster` scans 1
byte/cell, an 8#sym.times denser scan.

*Result.* L1-dcache miss rate drops *#s0.l1m% #sym.arrow #s1.l1m%* (#calc.round((1 - s1.l1m/s0.l1m)*100, digits: 0)% fewer
misses), IPC rises #s0.ipc #sym.arrow #s1.ipc, and runtime improves
*#s0.time s #sym.arrow #s1.time s — a #sp(s0.time, s1.time)#sym.times speed-up*,
with the output still bit-identical to the baseline. This is the expected,
clean win: a pure data-layout change that removes wasted memory traffic without
touching the math. On the large-$N$ server run, where the grid no longer fits in
cache, we expect this gap to *widen*, since bandwidth is then the true ceiling.

= An optimization that did *not* pay: temporal blocking (Stage 2)

This is the most instructive part of the project, so we document it in full.

*Hypothesis.* The 30 sweeps per step re-stream the entire $N^2$ grid 30 times.
If the grid exceeds the cache, that is $30 times$ the minimum DRAM traffic.
*Temporal blocking* (a.k.a. time-skewing) should fix it: load a `TILE`#sym.times`TILE`
tile plus a halo into a cache-resident scratch buffer, march it `T_BLOCK` sweeps
locally (the valid region shrinks one cell per side per sweep), write it back —
cutting global passes from $30$ to $30 \/ T_"block"$. On paper, a large win.

*What actually happened.* Our first version (L1-sized tiles, `TILE=64`,
`T_BLOCK=4`) ran *2.4#sym.times slower* than Stage 1. After resizing the tiles
for the 2 MB L2 and trimming the halo/copy overhead, Stage 2 still does not beat
Stage 1 single-threaded — it runs at *#s2.time s vs #s1.time s
(#sp(s1.time, s2.time)#sym.times*, a regression). The `perf` counters explain
why precisely:

#align(center)[#table(
  columns: (auto, 1fr, 1fr, 1fr, 1fr),
  align: (left, center, center, center, center),
  inset: (x: 7pt, y: 3pt), stroke: 0.4pt + luma(170),
  table.header([*Stage*], [*Time (s)*], [*Speed-up*], [*IPC*], [*L1-miss %*]),
  [#s0.name], [#s0.time], [1.00#sym.times (ref)], [#s0.ipc], [#s0.l1m],
  [#s1.name], [#s1.time], [#sp(s0.time, s1.time)#sym.times], [#s1.ipc], [#s1.l1m],
  [#s2.name], [#s2.time], [#sp(s0.time, s2.time)#sym.times], [#s2.ipc], [#s2.l1m],
)]
#v(-0.2em)
#align(center)[#text(8.3pt, style: "italic")[
  Instructions retired: Stage 1 #s1.instr B #sym.arrow Stage 2 #s2.instr B
  (#calc.round((s2.instr/s1.instr - 1)*100, digits: 0)% more), at lower IPC.
]]

*Root-cause.* At this $N$ the working set is $2 N^2 times 4"B" approx 2$ MB,
which *already fits in one P-core's 2 MB L2*. Stage 1 was therefore *never
memory-bandwidth-bound* (cache-miss rate #s1.cm%, IPC #s1.ipc) — there was no
DRAM traffic for temporal blocking to save. All blocking did was add halo
recomputation and addressing overhead: *+#calc.round((s2.instr/s1.instr - 1)*100, digits: 0)% instructions* at a
*lower* IPC (#s1.ipc #sym.arrow #s2.ipc). We optimized a bottleneck that did not
exist. Temporal blocking only pays once $2 N^2 dot 4"B" > $ L3 (24 MB), i.e.
$N gt.tilde 1800$ — below that it is pure overhead.

*Where it _does_ pay.* The win appears under *bandwidth contention*. Running
#conc.inst concurrent instances at $N=$#conc.n — so the shared memory subsystem
is saturated, as it would be when many cores compete for the same bus — Stage 2
finishes in *#conc.s2 s vs Stage 1's #conc.s1 s, a #sp(conc.s1, conc.s2)#sym.times
win*. The lesson is the core HW/SW-co-design point: *a locality optimization's
value is not intrinsic — it is conditional on the operating point* (working-set
vs cache size, and shared- vs private-bandwidth). Stage 2 is not a dead end; it
is simply the wrong tool for a single core whose working set already fits in L2.

= Comparison and conclusions

SoA (Stage 1) is the unambiguous win — #sp(s0.time, s1.time)#sym.times faster
from a layout change alone, by raising cache-line utilization from $approx$50% to
$approx$100% and deleting a per-sweep copy, all while staying bit-identical to the
baseline. Temporal blocking (Stage 2) taught us the opposite, and more valuable,
lesson: the *same* code can be a regression or a #sp(conc.s1, conc.s2)#sym.times
win depending purely on whether the working set spills the cache and whether
bandwidth is contended. An optimization is only worth its complexity when the
bottleneck it targets actually dominates — and on a single core with an
L2-resident grid, the bandwidth bottleneck it removes simply was not there. The
large-$N$ rerun of all stages on the university server, where the grid spills the
24 MB L3, will replace the placeholder numbers above and is where we expect the
SoA gap to widen and the time-skewing trade-off to shift in its favour.
