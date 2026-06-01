# CBRAM Simulation Loop — Technical Reference

## Overview

The simulation models the **SET operation** (filament formation) of a CBRAM device as a 2D grid evolving over discrete timesteps. Each grid cell holds four state variables:

| Field | Symbol | Meaning |
|---|---|---|
| `V` | φ | Electric potential |
| `sigma` | σ | Local conductivity (proxy for metal deposition) |
| `ion` | C | Metal ion concentration |
| `temp` | T | Temperature (allocated, not yet updated) |

All values are stored as **Q16.16 fixed-point integers**: the top 16 bits are the integer part, the bottom 16 bits are the fractional part. `1.0` is represented as `65536` (= `ONE`). This guarantees bit-identical results across all optimization stages (stages 0–3), allowing a simple `cmp` to verify correctness.

---

## Grid Layout

```
col:  0    1    2   ...  N-2   N-1
     ┌────┬────┬────┬────┬────┐
row 0│ V_A│ V_A│ V_A│ V_A│ V_A│  ← Anode (Dirichlet: V = V_APPLIED = 1.0)
     ├────┼────┼────┼────┼────┤
row 1│ ins│seed│    │seed│ ins│  ← Seeds placed here (row 1, random columns)
     ├────┼────┼────┼────┼────┤
  ...│    │    │    │    │    │  ← Interior: all physics computed here
     ├────┼────┼────┼────┼────┤
N-2  │    │    │    │    │    │  ← Bridge detection row
     ├────┼────┼────┼────┼────┤
N-1  │  0 │  0 │  0 │  0 │  0 │  ← Cathode (Dirichlet: V = 0)
     └────┴────┴────┴────┴────┘
      ins                  ins    ← Left/right: Neumann (insulating walls)
```

`ins` = insulating boundary (Neumann: dV/dx = 0, implemented by copying adjacent interior V).

---

## Per-Timestep Loop

Each outer timestep `t` runs the following steps **in order**:

### Step 1 — Jacobi Sweep (repeated K = 50 times)

**What it solves:** The steady-state Poisson equation with variable conductivity:

```
∇·(σ ∇V) = 0
```

This says: current is conserved everywhere. Given the conductivity field σ (which changes slowly due to filament growth), find the potential field V that satisfies this.

**How it works (Jacobi iteration):**

Each interior cell's new potential is the conductivity-weighted average of its four neighbors:

```
V_new[r,c] = (s_E·V[r,c+1] + s_W·V[r,c-1] + s_N·V[r+1,c] + s_S·V[r-1,c])
             / (s_E + s_W + s_N + s_S)
```

where the interface conductances are arithmetic means of adjacent cells:

```
s_E = (σ[r,c] + σ[r,c+1]) / 2
s_W = (σ[r,c] + σ[r,c-1]) / 2
s_N = (σ[r,c] + σ[r+1,c]) / 2
s_S = (σ[r,c] + σ[r-1,c]) / 2
```

**Why repeated 50 times:** A single Jacobi sweep propagates information only one cell. To get a reasonable approximation of the true solution, we repeat K=50 times per timestep. This is not full convergence — it's a fixed-budget approximation that is good enough for the slowly-evolving filament physics while keeping runtime predictable.

**Why this is the hot loop:** For an N×N grid, each sweep is O(N²) work, repeated K times per timestep, for T total timesteps → O(K·T·N²) total. For N=1024, K=50, T=3000: ~157 billion operations. This is the dominant cost and the target of all four optimization stages.

**Implementation detail:** Uses two buffers (`grid` and `grid_next`). Each sweep copies the full grid, updates V in `grid_next` from `grid`, then swaps. This ensures all reads come from the same generation (true Jacobi, not Gauss-Seidel).

---

### Step 2 — Drift-Diffusion (ion transport)

**What it models:** Metal ions (M⁺) drifting from the anode toward the cathode under the electric field. This is the simplified drift term of the Nernst-Planck equation:

```
J_ion = µ · σ · E · C
```

where:
- `µ` = `ION_MOBILITY` — ionic mobility constant  
- `σ` = local conductivity (high-conductivity paths carry more current → more ion drift)
- `E` = local downward electric field = `V[r,c] - V[r+1,c]`
- `C` = `ion` concentration at the source cell

**Why `σ · E` (J-driven):** Using current density J = σ·E rather than just E localizes ion transport to already-conductive paths. Background cells (σ ≈ 0.01) move 1000× fewer ions than seed cells (σ ≈ 10) even at the same field. This creates the positive feedback loop that forms a narrow filament:

```
σ high → more ion flux → more deposition → σ grows → more ion flux → ...
```

**Sweep direction:** Cells are processed bottom-to-top (`r = N-2` down to `r = 1`). This ensures each ion packet is transported exactly once per timestep — if we swept top-to-bottom, an ion moved to `r+1` could be immediately moved again in the same sweep.

**Anode injection:** After the sweep, row 1 (anode-adjacent) is replenished. The injection rate is proportional to local σ — physically, the metal electrode dissolves faster where current density is highest:

```
ion[1,c] += ION_INJECT · σ[1,c]
```

This ties the ion supply to the developing filament, preventing background flooding.

---

### Step 3 — Conductivity Update (sigma growth)

**What it models:** Metal ion deposition increases local conductivity. As ions arrive at a cell, they reduce and form metallic deposits, raising σ:

```
σ[r,c] += SIGMA_GROWTH · ion[r,c]
```

σ is clamped to `[SIGMA_LOW, SIGMA_MAX]` = [0.01, 20.0].

**Why this creates a filament:** Cells along the ion drift path accumulate deposits. The field concentrates at the filament tip (where σ drops sharply from high to low), which drives more ions there → tip extends → filament grows toward the cathode.

---

### Step 4 — Bridge Detection (termination)

After each timestep, the simulation checks whether any cell in row `N-2` (one row above the cathode) has reached `SIGMA_MAX`:

```
for c in 1..N-2:
    if σ[N-2, c] >= SIGMA_MAX → filament bridged → stop
```

This models the compliance current triggering when the filament completes the SET operation (HRS → LRS transition). Without this, the simulation would continue spreading ions laterally after bridging, which is unphysical.

---

## Initialization

Before the loop:

1. All cells: `V=0`, `σ=SIGMA_LOW`, `ion=ION_LOW`, `temp=TEMP_INIT`
2. Row 0 (anode): `V = V_APPLIED = 1.0`
3. `NUM_SEEDS` defect sites placed at row 1 with random columns (PRNG seed 42). Each seed gets a random strength `k ∈ [0.2, 1.0]`:
   - `σ = SIGMA_HIGH × k` (up to 10.0)
   - `ion = ION_HIGH × k` (up to 1.0)

   The random strength causes filaments to compete — stronger seeds grow faster and can suppress weaker neighbors by capturing their share of the electric field.

---

## What the Visualization Shows

Each PPM frame visualizes the `σ` (conductivity) field mapped through the viridis colormap:
- **Purple/dark** = background conductivity (SIGMA_LOW = 0.01) — no filament
- **Yellow/bright** = maximum conductivity (SIGMA_MAX = 20.0) — fully formed filament

The animation shows the filament growing downward from the seeds at the anode (top) toward the cathode (bottom) over time.

---

## Simplifications vs. Full CBRAM Physics

| Full ECM model | This simulation | Impact |
|---|---|---|
| Nernst-Planck: drift + diffusion | Drift only | Diffusion is small at high fields (our regime) — minor error |
| Butler-Volmer electrode kinetics (exponential threshold) | Linear injection ∝ σ | No sharp switching threshold; smoother filament formation |
| Filament nucleates at cathode, grows toward anode | Seeds at anode, filament grows toward cathode | Opposite growth direction — valid alternative model |
| Joule heating / electrothermal coupling | Not implemented (`temp` unused) | RESET behavior not modeled; SET shape slightly off |
| RESET operation | Not implemented | Only SET (forming) is simulated |
| 3D geometry, fractal filament | 2D cross-section, smooth continuum | No dendritic branching; statistical variability absent |
| Exponential field-enhanced hopping: sinh(qaE/2kT) | Linear in E | Underestimates switching nonlinearity at high fields |

---

## Key Parameters

| Constant | Value | Role |
|---|---|---|
| `N` | 200 (CLI) | Grid size (N×N cells) |
| `JACOBI_ITERS` | 50 | Poisson solve budget per timestep |
| `TOTAL_TIMESTEPS` | 3000 | Max simulation steps (early-exit if bridged) |
| `V_APPLIED` | 1.0 | Anode voltage |
| `SIGMA_LOW` | 0.01 | Background conductivity |
| `SIGMA_HIGH` | 10.0 | Peak seed conductivity |
| `SIGMA_MAX` | 20.0 | Conductivity ceiling + bridge threshold |
| `ION_MOBILITY` | 4.0 | Drift mobility µ in flux = µ·σ·E·ion |
| `ION_INJECT` | 0.1 | Anode dissolution rate per unit σ |
| `SIGMA_GROWTH` | 0.1 | Conductivity growth per unit ion per step |
| `NUM_SEEDS` | 3 | Number of defect nucleation sites |
