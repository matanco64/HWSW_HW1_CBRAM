// ============================================================
// HW1 - CBRAM: BONUS / supplementary note (not part of the 3-page main report)
// Build:  typst compile report/cbram_bonus.typ
// Companion to cbram_report.typ - kept separate to respect the 3-page limit.
// ============================================================

#set document(title: "CBRAM Simulation - Bonus: the bit-identical gate")
#set page(paper: "a4", margin: (x: 1.7cm, y: 1.5cm), numbering: "1 / 1")
#set text(font: "New Computer Modern", size: 9.7pt)
#set par(justify: true, leading: 0.55em)
#show heading.where(level: 1): set text(size: 11.5pt)
#show heading: set block(above: 0.85em, below: 0.45em)
#show raw: set text(font: "DejaVu Sans Mono", size: 8.6pt)

// ── Title ───────────────────────────────────────────────────────────────────
#align(center)[
  #text(15pt, weight: "bold")[Bonus note: what the bit-identical gate actually proves]
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

Throughout the project our correctness gate was a single `cmp` on the raw
`V_final` dump (the final $N times N$ fixed-point `int32` potential grid): every
stage (AoS, SoA, time-skewed) must produce *bit-identical* output, so no
optimization can "cheat" by changing the numerics.
Passing that gate across three *structurally different* kernels is a stronger result
than a typical regression test, on four independent levels - this note unpacks why.

= Algorithmic - why tiling is even allowed

Each Jacobi sweep is double-buffered: every cell reads sweep $k$ from buffer $A$ and
writes sweep $k+1$ to buffer $B$, so cells within a sweep are *mutually
independent*. Spatial tiling (SoA's blocked sweep, Stage 2's tile marcher) only
changes the *order* in which those independent writes happen, never their inputs -
so it is safe by construction. The contrast is Gauss-Seidel, which reads
partially-updated neighbours within the same sweep: there the result *depends* on
traversal order, and tiling would silently change the answer. Choosing Jacobi is
precisely what makes the whole family of reorderings admissible.

= Geometric - the gate verifies the skew

Time-skewing reads neighbours at *earlier* time-levels, so a cell is correct only if
every neighbour it touches is itself valid at the right level. That is exactly why
the tile carries a halo of `T_BLOCK` cells per side and the valid region shrinks one
cell per side per local sweep. An off-by-one in the halo width - or applying the
in-tile boundary re-pinning at the wrong level - would silently feed an edge cell a
stale time-level, producing a *plausible but wrong* filament that no eyeball
inspection of the field plots would catch. Bit-identity is therefore a direct,
machine-checked proof that the skewing geometry (halo width, shrink schedule,
boundary interleaving) is exactly right - the gate doubles as the correctness proof
of the transformation, not just a smoke test.

= Numerical - why it is achievable at all

All arithmetic is fixed-point Q16.16 `int32` with an arithmetic-shift divide, so the
stencil has *zero* rounding freedom. This is what lets `cmp` be a *zero-tolerance*
gate. Under `-O2 -march=native` the three kernels vectorize and FMA-contract their
differently-shaped loops differently; in floating point those compiler choices would
perturb the low mantissa bits - the math would be "the same" yet the bytes would
differ, forcing a tolerance-based comparison that can mask real bugs (a genuine
geometry error hiding under the rounding noise). Integer math makes "same
computation" mean "same bits," independent of how each kernel was compiled,
vectorized, or fused. Picking fixed point was therefore not incidental: it is what
turns the correctness gate from approximate to exact.

= PRNG lock-step - same trajectory, not just same endpoint

Growth-site selection (probability $prop V^eta$, $eta = 3$) draws from a
deterministic xorshift64 stream. Because $V$ is bit-identical at *every* step and the
stream is seeded and consumed in the same order, each stage picks the *same* cell at
every growth step - so the filaments evolve in lock-step, not merely to the same
final picture. Only the Jacobi solve is restructured between stages; the
site-selection logic and its RNG consumption are untouched, and a passing `cmp`
confirms the streams never desynchronized. This is also what makes the side-by-side
race in the main report's Fig. 1 a *fair* one: identical growth sequence, with only
the wall-clock differing - so the on-screen gap is purely the speed-up.
