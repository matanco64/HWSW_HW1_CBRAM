// ============================================================
// HW1 - CBRAM: APPENDIX A - Stage 3 (vectorizing the cache-resident stencil)
// Build:  typst compile report/cbram_stage3_appendix.typ
// Companion to cbram_report.typ - kept separate to respect the 3-page limit on
// the main report. Numbers are from the SAME course-server run as the main report
// (single-core Xeon E5-2630 v3, N = 6144, 80 growth steps).
// ============================================================

#set document(title: "CBRAM Simulation - Appendix A: Stage 3 (SIMD)")
#set page(paper: "a4", margin: (x: 1.7cm, y: 1.5cm), numbering: "A-1 / 1")
#set text(font: "New Computer Modern", size: 9.7pt)
#set par(justify: true, leading: 0.55em)
#show heading.where(level: 1): set text(size: 11.5pt)
#show heading: set block(above: 0.85em, below: 0.45em)
#set heading(numbering: "A.1")
#show raw: set text(font: "DejaVu Sans Mono", size: 8.6pt)

// ── DATA (same run as the main report) ──────────────────────────────────────
#let s0 = (name: "Stage 0 - AoS",                time: 375.82, ipc: 1.746, l1m: 19.00, cmiss: 3.012, instr: 1508, cyc: 864)
#let s1 = (name: "Stage 1 - SoA",                time: 245.88, ipc: 2.359, l1m: 7.30,  cmiss: 0.462, instr: 1343, cyc: 569)
#let s2 = (name: "Stage 2 - SoA + skewing",      time: 284.92, ipc: 2.845, l1m: 4.25,  cmiss: 0.162, instr: 1881, cyc: 661)
#let s3 = (name: "Stage 3 - SoA + skewing + SIMD", time: 83.53, ipc: 2.110, l1m: 15.24, cmiss: 0.153, instr: 408, cyc: 193)
#let sp(a, b) = calc.round(a / b, digits: 2)
#let pct(x)   = calc.round(x * 100, digits: 0)

// ── Title ───────────────────────────────────────────────────────────────────
#align(center)[
  #text(15pt, weight: "bold")[Appendix A - Stage 3: vectorizing the cache-resident stencil]
  #v(0.2em)
  #text(9.7pt)[
    Matan Cohen (ID: #sys.inputs.at("matan-id", default: "<ID>")) #h(1.2em) · #h(1.2em) Yuval Kogan (ID: #sys.inputs.at("yuval-id", default: "<ID>")) \
    HW/SW Co-Design - HW1 · supplement to the main report (kept separate to respect the 3-page limit) \
    Code & data: #link("https://github.com/matanco64/HWSW_HW1_CBRAM")[`github.com/matanco64/HWSW_HW1_CBRAM`]
  ]
]
#v(0.15em)
#line(length: 100%, stroke: 0.4pt)
#v(0.3em)

The main report ends on an open question: Stage 2 (time-skewing) made the Jacobi
sweep *cache-resident* - the best memory behaviour of the three layouts - but lost
#pct(s2.time/s1.time - 1)% to SoA because the halo recomputation *added*
#pct(s2.instr/s1.instr - 1)% instructions, and on a single core the traffic it saved
was no longer the limiter. That diagnosis names its own cure: if the loss is *added
compute*, make the compute cheap. Stage 3 does exactly that, and it is the fastest
stage in the project - *#s3.time s*, *#sp(s0.time, s3.time)#sym.times* over the
baseline, *#sp(s1.time, s3.time)#sym.times* over SoA, and *#sp(s2.time, s3.time)#sym.times*
over the Stage 2 it builds on - while staying bit-identical to every other stage.

= The lever: the stencil was never vectorized

Inspecting GCC's vectorizer report (`-fopt-info-vec-missed`) on Stages 1 and 2
revealed that *the 4-point Jacobi stencil - the hot loop - ran scalar in both*. GCC
15.2 (`-O2 -march=native`) bailed with "complicated access pattern" / "not
profitable": the grid is reached through `std::vector::operator[]`, so the compiler
cannot prove the read buffer `cur` and the write buffer `nxt` (or the `metal` array)
do not alias, and its cost model then vetoes the transform. The Xeon E5-2630 v3
(Haswell) has AVX2 - *8 `int32` lanes per 256-bit vector* - so the stencil was
leaving #sym.tilde 8#sym.times of integer throughput on the table. Crucially this
matters *most* for Stage 2: its tile is L2-resident, so it is *throughput*-bound, and
vectorized arithmetic converts almost directly to wall-clock. The same SIMD applied
to Stage 1 would still stall on DRAM latency and help far less. *Vectorization is the
missing complement to time-skewing* - cache-residency and compute throughput exploit
each other.

= The change: three edits, all inside `advance_tile()`

`dbm_stage3.cpp` is `dbm_stage2.cpp` with three localized edits; everything else
(`block_advance`, `solve_potential`, candidate selection, PRNG, IO, all constants) is
copied verbatim, so the *only* variable is the stencil codegen.

#set enum(numbering: "1.")
+ *Vectorize the stencil.* Hoist `__restrict` raw pointers and annotate the inner run
  with `#raw("#pragma omp simd")` (compiled with `-fopenmp-simd`, no OpenMP runtime).
  This removes the alias ambiguity *and* overrides the cost-model veto. GCC then
  reports `loop vectorized using 32 byte vectors`. The arithmetic shift `>> 2` maps to
  `vpsrad` (`_mm256_srai_epi32`), exact on these non-negative averages.
+ *Fuse the Dirichlet re-pin into the store.* Stage 2 re-pinned metal cells to $V=0$
  in a *separate pass* over the tile every sweep. Folding it into the stencil store as
  a branchless select writes each cell once and deletes an entire per-sweep pass:
+ *`__restrict` the tile copy-in/out* so they emit `memcpy`-style vector moves.

#v(0.2em)
#block(fill: luma(245), inset: 6pt, radius: 3pt, width: 100%)[
```cpp
const int32_t* __restrict cv = cur->data();   // hoisted, non-aliasing
int32_t*       __restrict nv = nxt->data();
const uint8_t* __restrict mv = m.data();
#pragma omp simd
for (int c = cc0; c < cc1; ++c) {
    size_t li = base + c;
    int32_t avg = (cv[li - LW] + cv[li + LW] + cv[li - 1] + cv[li + 1]) >> 2;
    nv[li] = mv[li] ? 0 : avg;                 // re-pin fused into the store
}
```
]

= Result: the instruction count collapses

#align(center)[#table(
  columns: (auto, 1fr, 1fr, 1fr, 1fr, 1fr),
  align: (left, center, center, center, center, center),
  inset: (x: 6pt, y: 3pt), stroke: 0.4pt + luma(170),
  table.header([*Stage*], [*Time (s)*], [*vs base*], [*Instr (B)*], [*IPC*], [*L1-miss %*]),
  [#s0.name], [#s0.time], [1.00#sym.times], [#s0.instr], [#s0.ipc], [#s0.l1m],
  [#s1.name], [#s1.time], [#sp(s0.time, s1.time)#sym.times], [#s1.instr], [#s1.ipc], [#s1.l1m],
  [#s2.name], [#s2.time], [#sp(s0.time, s2.time)#sym.times], [#s2.instr], [#s2.ipc], [#s2.l1m],
  [*#s3.name*], [*#s3.time*], [*#sp(s0.time, s3.time)#sym.times*], [*#s3.instr*], [#s3.ipc], [#s3.l1m],
)]
#v(0.2em)

The mechanism is visible directly in the counters: instructions retired fall
*#s2.instr B #sym.arrow #s3.instr B* (a *#sp(s2.instr, s3.instr)#sym.times* cut) - the
roughly 8#sym.times{} narrower stencil arithmetic plus the deleted re-pin pass - and cycles
drop *#s2.cyc B #sym.arrow #s3.cyc B* (#sp(s2.cyc, s3.cyc)#sym.times). IPC *decreases*
(#s2.ipc #sym.arrow #s3.ipc), which is exactly right: each vector instruction now does
the work of eight scalar ones, so there are far fewer, heavier instructions in flight
- fewer instructions at slightly lower IPC still finishes in *#sp(s2.cyc, s3.cyc)#sym.times*
the cycles. L1-miss rises (#s2.l1m% #sym.arrow #s3.l1m%) because the vector loads pull
more bytes per instruction, but the tile is L2-resident so those misses are cheap and
the throughput gain dominates. This is the payoff the main report predicted:
time-skewing's cache-residency was latent headroom, and vectorizing its now-resident
stencil *cashes it in* - the slowest "optimization" becomes the fastest stage.

= Correctness preserved

Stage 3 clears the same gate as every other stage: `V_final_stage3.bin` is
*bit-identical* to Stage 0's (`cmp` silent, exit 0; `run.sh` enforces this before
trusting any perf number). This holds because the fixed-point `int32` add + arithmetic
shift is exact and order-independent - `_mm256_srai_epi32` matches scalar `>> 2`
lane-for-lane - and the fused re-pin writes the identical final value. Vectorizing the
math changed *how* the bytes are produced, not *which* bytes, so the SIMD kernel is
provably equivalent, not merely close.

#figure(
  kind: image, supplement: [Figure],
  image("results/stage3_flamegraph.svg", width: 12.8cm),
  caption: [Stage 3 self-time flame graph (`perf record -F 999`, $N=6144$). The
    `advance_tile` stencil that dominated Stage 2 is now a thin, vectorized band - the
    runtime is no longer spent in scalar stencil arithmetic.],
)

#v(0.3em)
*Next step.* Stage 3 is single-core. The remaining headroom is the *other* path the
main report names: parallelizing tiles across cores so the shared memory bus becomes
the binding resource again - at which point Stage 2's traffic savings, dormant on one
core, would finally pay off too. That multicore stage is the natural follow-up.
