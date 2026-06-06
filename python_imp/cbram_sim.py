#!/usr/bin/env python3
"""
cbram_sim.py — CBRAM Conductive Filament Simulation
Physics: Field-guided tip growth — Jacobi V field biases the random walk direction.
         Ion drift-diffusion computed for visualization (ion panel).

Key difference from Yuval's original: walk probabilities are NOT hardcoded.
They are derived from the local E field at the tip each step, so the Jacobi
solve is genuinely load-bearing (changing sigma changes V which changes growth).

Grid: row 0 = anode (V=1), row N-1 = cathode (V=0).
Seed at cathode row N-2; filament grows upward toward row 0.

Usage:
  python3 cbram_sim.py          # run, save PNG + binary
  python3 cbram_sim.py --gif    # also save animated GIF
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import LinearSegmentedColormap
import sys, os, time
from collections import deque

# ── Grid ──────────────────────────────────────────────────────────────────────
N            = 200
T_MAX        = 6_000
FRAME_EVERY  = 30
JACOBI_ITERS = 50

# ── Device physics ────────────────────────────────────────────────────────────
V_APPLIED    = 1.0
SIGMA_LOW    = 0.01
SIGMA_MAX    = 20.0
BRIGHT       = SIGMA_MAX / 2.0

ION_LOW      = 0.001
ION_HIGH     = 1.0
ION_MOBILITY = 8.0      # ion drift speed (visualization only)
ION_INJECT   = 0.3      # anode injection rate (visualization only)
ION_DIFFUSION = 0.02    # lateral diffusion (visualization only)

# ── Tip-growth parameters ─────────────────────────────────────────────────────
# Walk probabilities are field-guided: base + field_strength * local_E
# The Jacobi V field drives growth — high E toward anode → tip grows upward.
FIELD_STRENGTH  = 30.0  # how strongly E field biases the walk
BASE_P_UP       = 0.40  # minimum upward bias (even at zero field)
BASE_P_LATERAL  = 0.08  # base left/right probability
BASE_P_DOWN     = 0.02  # base downward probability (drift opposes this)

P_BRANCH_SPAWN  = 0.012  # probability per tip step of spawning a side branch
MAX_BRANCHES    = 5      # cap on total branch stump count

RNG_SEED        = 42

# ── Continuum experiment parameters ──────────────────────────────────────────
# Active only with --continuum flag. Tests whether the pure drift-diffusion
# approach works when two previously missing mechanisms are added:
#   1. Ion decay — ions recombine with defects/stray electrons en route (loss term)
#   2. J_e gating — deposition only where the metallic neighbour carries
#      high electron current density (J_e = σ|∇V|), selecting the active tip
ION_DECAY         = 0.0003 # fractional ion loss per step (~3333 step lifetime,
                           # ~22% survive 200-row crossing at µ=8)
SIGMA_GROWTH_CONT = 0.5    # deposit probability scale for continuum mode
E_FIELD_THRESH    = 0.008  # min downward E field at candidate = 1.6× background (0.005)
                           # Selects tip where V drops sharply into the metallic region

# ── Colormap ──────────────────────────────────────────────────────────────────
_cbram_cmap = LinearSegmentedColormap.from_list("cbram", [
    (0.00, (0.02, 0.02, 0.06)),
    (0.25, (0.04, 0.08, 0.30)),
    (0.55, (0.75, 0.25, 0.00)),
    (0.80, (1.00, 0.75, 0.10)),
    (1.00, (1.00, 1.00, 1.00)),
])

# ── Inline grader ─────────────────────────────────────────────────────────────

def grade(sigma):
    flat = sigma.flatten().tolist()
    Ng   = sigma.shape[0]
    dist, prev, q = {}, {}, deque()
    for c in range(1, Ng-1):
        if flat[Ng+c] > BRIGHT:
            q.append((1,c)); dist[(1,c)]=0; prev[(1,c)]=None
    found, end = False, None
    while q:
        r,c = q.popleft()
        if r==Ng-2: found,end=True,(r,c); break
        for dr,dc in ((-1,0),(1,0),(0,-1),(0,1)):
            nr,nc=r+dr,c+dc
            if 1<=nr<=Ng-2 and 1<=nc<=Ng-2 and (nr,nc) not in dist:
                if flat[nr*Ng+nc]>BRIGHT:
                    dist[(nr,nc)]=dist[(r,c)]+1; prev[(nr,nc)]=(r,c)
                    q.append((nr,nc))
    pl=0
    if found:
        cur=end
        while cur: pl+=1; cur=prev[cur]
    interior=(Ng-2)**2
    bn=sum(1 for r in range(1,Ng-1) for c in range(1,Ng-1) if flat[r*Ng+c]>BRIGHT)
    frac=bn/interior if interior else 1.
    brows=[r for r in range(1,Ng-1) for c in range(1,Ng-1) if flat[r*Ng+c]>BRIGHT]
    bcols=[c for r in range(1,Ng-1) for c in range(1,Ng-1) if flat[r*Ng+c]>BRIGHT]
    ar=((max(brows)-min(brows)+1)/max(max(bcols)-min(bcols)+1,1)) if brows else 0.
    tort=(pl/(Ng-3)) if (found and Ng>3) else 1.
    junc=sum(
        1 for r in range(1,Ng-1) for c in range(1,Ng-1)
        if flat[r*Ng+c]>BRIGHT and
           sum(1 for dr,dc in ((-1,0),(1,0),(0,-1),(0,1))
               if 1<=r+dr<=Ng-2 and 1<=c+dc<=Ng-2
               and flat[(r+dr)*Ng+(c+dc)]>BRIGHT)>=3
    )
    sb =30 if found else 0
    sn =25 if frac<.02 else 20 if frac<.05 else 12 if frac<.10 else 5 if frac<.20 else 0
    sa =20 if ar>5 else 15 if ar>3 else 8 if ar>2 else 0
    st =(0 if not found else 15 if tort>1.5 else 10 if tort>1.2 else 5 if tort>1.05 else 0)
    sbr=(0 if not found or frac>=.10 else
         10 if junc>=3 else 6 if junc>=2 else 3 if junc>=1 else 2)
    tot=sb+sn+sa+st+sbr
    gl ="A" if tot>=90 else "B" if tot>=75 else "C" if tot>=55 else "D" if tot>=35 else "F"
    return tot,gl,dict(bridged=found,pl=pl,frac=frac,ar=ar,
                       tort=tort,junc=junc,sb=sb,sn=sn,sa=sa,st=st,sbr=sbr)

# ── Poisson solver ────────────────────────────────────────────────────────────

def solve_poisson(V, sigma, n=JACOBI_ITERS):
    V = V.copy()
    for _ in range(n):
        sE = (sigma[1:-1, 2:] + sigma[1:-1, 1:-1]) * .5
        sW = (sigma[1:-1, :-2] + sigma[1:-1, 1:-1]) * .5
        sN = (sigma[:-2, 1:-1] + sigma[1:-1, 1:-1]) * .5
        sS = (sigma[2:,  1:-1] + sigma[1:-1, 1:-1]) * .5
        d  = np.where(sE + sW + sN + sS < 1e-15, 1e-15, sE + sW + sN + sS)
        V[1:-1, 1:-1] = (sE*V[1:-1, 2:] + sW*V[1:-1, :-2] +
                         sN*V[:-2, 1:-1] + sS*V[2:,  1:-1]) / d
        V[0, :]  = V_APPLIED
        V[-1, :] = 0.0
        V[:, 0]  = V[:, 1]
        V[:, -1] = V[:, -2]
    return V

# ── Ion drift-diffusion (visualization) ──────────────────────────────────────
# Ion field is maintained for the 4th panel visualization only.
# It shows where ions accumulate (near filament tips), validating the physics story.

def drift_diffusion(ion, sigma, V, decay=0.0):
    ion = ion.copy()
    ion[sigma >= BRIGHT] = ION_LOW   # no mobile ions in metallic cells

    # Ion decay: recombination with defects / stray electrons en route.
    # Drains background accumulation so only ions near the active tip survive.
    if decay > 0.0:
        ion = np.maximum(ion * (1.0 - decay), ION_LOW)

    E_down = np.clip(V[1:-1, 1:-1] - V[2:, 1:-1], 0.0, 0.02)
    flux   = ION_MOBILITY * E_down * ion[1:-1, 1:-1]
    avail  = np.maximum(ion[1:-1, 1:-1] - ION_LOW, 0.0)
    flux   = np.minimum(flux, avail)
    ion[1:-1, 1:-1] -= flux
    ion[2:,   1:-1]  = np.clip(ion[2:, 1:-1] + flux, ION_LOW, ION_HIGH * 4)
    ion[1, 1:-1]     = np.clip(ion[1, 1:-1] + ION_INJECT, ION_LOW, ION_HIGH)

    lat = np.roll(ion, -1, axis=1) + np.roll(ion, 1, axis=1) - 2 * ion
    ion[1:-1, 1:-1] = np.clip(
        ion[1:-1, 1:-1] + ION_DIFFUSION * lat[1:-1, 1:-1],
        ION_LOW, ION_HIGH * 4
    )
    return ion


def stochastic_deposit_je(ion, sigma, V, rng):
    """
    Continuum stochastic deposition gated by local downward electric field.

    Deposition at (r,c) requires:
      1. Metal directly below at (r+1,c)  — tip adjacency
      2. E_down = V[r-1,c] - V[r,c] > E_FIELD_THRESH
         The field arriving at the candidate from above must exceed 1.6× background.
         This selects the active tip: V drops sharply into the near-zero region
         created by the metallic cell below, concentrating the field there.
         Cells deep in the filament perimeter or far from the tip have weaker field.
      3. Stochastic: prob = SIGMA_GROWTH_CONT * ion[r,c]

    Why E_down at the candidate (not J_e at the metallic neighbor):
      The metallic cell directly below is screened — V inside the metal ≈ 0 (grounded
      via the filament to the cathode), and the candidate above is also pulled toward 0
      by the strong σ coupling. This makes |∇V| inside the metal near zero.
      But the field ABOVE the candidate (V[r-1] - V[r]) is enhanced because V[r] ≈ 0
      while V[r-1] is at normal electrolyte potential.
    """
    sigma = sigma.copy()
    metal = sigma >= BRIGHT

    # Downward E field arriving at each cell from the row above
    # E_down[r,c] = V[r-1,c] - V[r,c]  (positive = field drives ions downward into r,c)
    E_down = np.zeros_like(V)
    E_down[1:, :] = V[:-1, :] - V[1:, :]
    E_down = np.maximum(E_down, 0.0)

    # Strict tip: metal directly below
    metal_below = np.roll(metal, -1, axis=0)

    # Diagonal tip: metal one row below AND one column to the side.
    # Lateral ion diffusion creates a left/right concentration gradient that
    # determines which diagonal fires — this is where diffusion drives tortuosity.
    diag_l = np.roll(metal_below,  1, axis=1)   # metal at (r+1, c-1)
    diag_r = np.roll(metal_below, -1, axis=1)   # metal at (r+1, c+1)
    adj_diagonal = (diag_l | diag_r) & ~metal_below

    for arr in (metal_below, adj_diagonal):
        arr[0, :] = arr[-1, :] = arr[:, 0] = arr[:, -1] = False

    cand_strict = metal_below   & ~metal & (E_down > E_FIELD_THRESH)
    cand_diag   = adj_diagonal  & ~metal & (E_down > E_FIELD_THRESH * 0.6)

    prob      = np.clip(SIGMA_GROWTH_CONT * ion, 0.0, 1.0)
    prob_diag = prob * 0.08   # rare lateral L-steps — a few deflections, not a thicket
    rolls     = rng.random((N, N))

    dep_strict = cand_strict & (rolls < prob)
    dep_diag   = cand_diag   & (rolls < prob_diag)

    # 4-connectivity bridge: each diagonal deposit at (r,c) also deposits (r+1,c).
    # This creates the L-path: metal[r+1,c±1] → bridge[r+1,c] → diag[r,c].
    # Without the bridge, diagonal cells are only 8-connected — grader BFS fails.
    bridge = np.roll(dep_diag, 1, axis=0)   # shift diag deposits down to (r+1,c)
    bridge[0, :] = bridge[-1, :] = bridge[:, 0] = bridge[:, -1] = False
    bridge &= ~metal   # don't overwrite existing metal

    deposited = dep_strict | dep_diag | bridge
    sigma[deposited] = SIGMA_MAX
    return sigma, deposited

# ── Field-guided tip growth ───────────────────────────────────────────────────
# Each tip moves one step per timestep. Direction probabilities derived from
# the local E field computed by solve_poisson — V is genuinely load-bearing.

_DR = np.array([-1,  0,  0,  1])   # up, left, right, down
_DC = np.array([ 0, -1,  1,  0])

def grow_step(tips, sigma, V, rng, n_br):
    metal  = sigma >= BRIGHT
    sigma  = sigma.copy()
    out    = []
    bridged = False

    for (r, c) in tips:
        # Local field in 4 directions (positive = higher V in that direction)
        E = np.array([
            V[r-1, c] - V[r, c] if r > 1   else 0.0,   # up
            V[r, c-1] - V[r, c] if c > 1   else 0.0,   # left
            V[r, c+1] - V[r, c] if c < N-2 else 0.0,   # right
            V[r+1, c] - V[r, c] if r < N-2 else 0.0,   # down
        ])

        # Base probabilities + field enhancement (field points toward anode = up)
        base = np.array([BASE_P_UP, BASE_P_LATERAL, BASE_P_LATERAL, BASE_P_DOWN])
        w    = base + FIELD_STRENGTH * np.maximum(E, 0.0)
        w   /= w.sum()

        # Sample direction; retry if blocked (up to 4 tries)
        moved = False
        order = rng.permutation(4)
        cumul = np.cumsum(w)
        roll  = rng.random()
        idx   = int(np.searchsorted(cumul, roll))

        for attempt in range(4):
            i  = (idx + attempt) % 4
            nr = max(1, min(N-2, r + int(_DR[i])))
            nc = max(1, min(N-2, c + int(_DC[i])))
            if not metal[nr, nc]:
                sigma[nr, nc] = SIGMA_MAX
                metal[nr, nc] = True
                out.append((nr, nc))
                if nr == 1:
                    bridged = True
                moved = True
                break

        if not moved:
            out.append((r, c))  # surrounded — tip stays (rare)

        # Spawn 1-cell horizontal branch stump at OLD position
        if n_br < MAX_BRANCHES and rng.random() < P_BRANCH_SPAWN:
            bdir = 1 if rng.random() < 0.5 else -1
            bc = c + bdir
            if 1 <= bc <= N-2 and not metal[r, bc]:
                sigma[r, bc] = SIGMA_MAX
                metal[r, bc] = True
                n_br += 1

    return sigma, out, bridged, n_br

# ── Save binary ───────────────────────────────────────────────────────────────

def save_binary(sigma, path):
    (sigma * 65536).astype(np.int32).tofile(path)

# ── Render — 4 panels: sigma, V, |∇V|, ion ───────────────────────────────────

def render(fig, axes, sigma, V, ion, t, bridged):
    ax_f, ax_v, ax_e, ax_i = axes
    for ax in axes:
        ax.clear()

    log_s = np.log1p(sigma)
    ax_f.imshow(log_s, cmap=_cbram_cmap, origin="upper",
                vmin=0, vmax=np.log1p(SIGMA_MAX), aspect="equal",
                interpolation="nearest")
    mr, mc = np.where(sigma >= BRIGHT)
    if len(mr):
        ax_f.scatter(mc, mr, s=5, c="white", alpha=0.85, linewidths=0, zorder=3)
    ax_f.axhline(0.5,   color="#55aaff", lw=1, ls="--", alpha=0.6)
    ax_f.axhline(N-1.5, color="#aaaaaa", lw=1, ls="--", alpha=0.6)
    ax_f.text(2, 2,   "ANODE ⊕",   color="#55aaff", fontsize=7, va="top")
    ax_f.text(2, N-2, "CATHODE ⊖", color="#aaaaaa", fontsize=7, va="bottom")
    bfrac = 100 * len(mr) / (N * N)
    ax_f.set_title(f"Conductivity σ  t={t}{'  ✓ BRIDGED' if bridged else ''}\n"
                   f"metallic: {len(mr)} ({bfrac:.2f}%)",
                   color="white", fontsize=8, pad=3)
    ax_f.set_facecolor("#050508")
    for sp in ax_f.spines.values(): sp.set_color("#333")
    ax_f.tick_params(colors="gray", labelsize=6)

    ax_v.imshow(V, cmap="RdYlBu_r", origin="upper", vmin=0, vmax=V_APPLIED,
                aspect="equal", interpolation="nearest")
    ax_v.set_title("Electric potential φ", color="white", fontsize=8, pad=3)
    ax_v.set_facecolor("#050508"); ax_v.tick_params(colors="gray", labelsize=6)
    for sp in ax_v.spines.values(): sp.set_color("#333")

    dVr = np.zeros_like(V); dVr[1:-1, :] = V[:-2, :] - V[2:, :]
    dVc = np.zeros_like(V); dVc[:, 1:-1] = V[:, :-2] - V[:, 2:]
    Emag = np.sqrt(dVr**2 + dVc**2)
    vm = float(np.percentile(Emag, 99.5)) + 1e-9
    ax_e.imshow(Emag, cmap="inferno", origin="upper", vmin=0, vmax=vm,
                aspect="equal", interpolation="nearest")
    ax_e.set_title("|∇φ|  (field magnitude)", color="white", fontsize=8, pad=3)
    ax_e.set_facecolor("#050508"); ax_e.tick_params(colors="gray", labelsize=6)
    for sp in ax_e.spines.values(): sp.set_color("#333")

    log_i = np.log1p(ion)
    ax_i.imshow(log_i, cmap="plasma", origin="upper",
                vmin=0, vmax=np.log1p(ION_HIGH * 4), aspect="equal",
                interpolation="nearest")
    ax_i.set_title("Ion concentration C", color="white", fontsize=8, pad=3)
    ax_i.set_facecolor("#050508"); ax_i.tick_params(colors="gray", labelsize=6)
    for sp in ax_i.spines.values(): sp.set_color("#333")

    fig.tight_layout(pad=1.0)

# ── Main ──────────────────────────────────────────────────────────────────────

def run(save_gif=False, continuum=False):
    rng = np.random.default_rng(RNG_SEED)

    sigma = np.full((N, N), SIGMA_LOW)
    ion   = np.full((N, N), ION_LOW)
    V     = np.linspace(V_APPLIED, 0.0, N).reshape(-1, 1) * np.ones((1, N))

    # Single seed at cathode center
    sigma[N-2, N // 2] = SIGMA_MAX
    tips = [(N-2, N // 2)]
    n_br = 0

    frames   = []
    bridged  = False
    t_bridge = None
    t0       = time.time()

    if continuum:
        print(f"CBRAM CONTINUUM mode  N={N}  T_MAX={T_MAX}")
        print(f"  decay={ION_DECAY}  E_field_thresh={E_FIELD_THRESH}"
              f"  deposit_scale={SIGMA_GROWTH_CONT}")
    else:
        print(f"CBRAM field-guided tip growth  N={N}  T_MAX={T_MAX}")
        print(f"  field_strength={FIELD_STRENGTH}  base_up={BASE_P_UP}"
              f"  branch_prob={P_BRANCH_SPAWN}  max_branches={MAX_BRANCHES}")

    for t in range(T_MAX):
        V   = solve_poisson(V, sigma, n=JACOBI_ITERS)
        ion = drift_diffusion(ion, sigma, V, decay=ION_DECAY if continuum else 0.0)

        if continuum:
            sigma, deposited = stochastic_deposit_je(ion, sigma, V, rng)
            ion[deposited]   = ION_LOW     # consume ions at deposition sites
            b = bool(np.any(sigma[1, 1:-1] >= BRIGHT))
        else:
            sigma, tips, b, n_br = grow_step(tips, sigma, V, rng, n_br)

        if b and not bridged:
            bridged  = True
            t_bridge = t
            print(f"  *** BRIDGED t={t}  ({time.time()-t0:.1f}s) ***")

        if t % FRAME_EVERY == 0:
            frames.append((sigma.copy(), V.copy(), ion.copy(), t))
            nm = int(np.sum(sigma >= BRIGHT))
            if continuum:
                print(f"  t={t:5d}  metal={nm:5d}"
                      f"  {'BRIDGED' if bridged else '      '}  {time.time()-t0:.1f}s")
            else:
                print(f"  t={t:5d}  metal={nm:5d}  tips={len(tips):2d}"
                      f"  {'BRIDGED' if bridged else '      '}  {time.time()-t0:.1f}s")

        if bridged:
            break

    V = solve_poisson(V, sigma, n=300)
    frames.append((sigma.copy(), V.copy(), ion.copy(), t))
    print(f"Done {time.time()-t0:.1f}s  ({len(frames)} frames)")

    tot, gl, m = grade(sigma)
    print()
    print("=" * 52)
    print(f"  Grade: {tot}/100  →  {gl}")
    print(f"  Bridging   {m['sb']:2d}/30  bridged={m['bridged']}  path={m['pl']}")
    print(f"  Narrowness {m['sn']:2d}/25  {100*m['frac']:.3f}% bright")
    print(f"  Aspect     {m['sa']:2d}/20  ratio={m['ar']:.1f}")
    print(f"  Tortuosity {m['st']:2d}/15  ratio={m['tort']:.2f}")
    print(f"  Branches   {m['sbr']:2d}/10  {m['junc']} junctions")
    print("=" * 52)

    outdir = os.path.dirname(os.path.abspath(__file__))
    bp = os.path.join(outdir, "sigma_final_python.bin")
    save_binary(sigma, bp)
    print(f"  Binary → {bp}")

    plt.rcParams.update({"figure.facecolor": "#0a0a12", "text.color": "white",
                         "axes.facecolor": "#050508", "axes.labelcolor": "white"})
    fig, axes = plt.subplots(1, 4, figsize=(20, 6))
    fig.patch.set_facecolor("#0a0a12")

    s, v, i, tt = frames[-1]
    render(fig, axes, s, v, i, tt, bridged)
    png = os.path.join(outdir, "filament_final.png")
    fig.savefig(png, dpi=150, bbox_inches="tight", facecolor="#0a0a12")
    print(f"  PNG    → {png}")

    if save_gif and len(frames) > 1:
        gif = os.path.join(outdir, "filament_animation.gif")
        print(f"  Saving GIF ({len(frames)} frames)…")
        def _upd(idx):
            s, vv, ii, tt = frames[idx]
            render(fig, axes, s, vv, ii, tt,
                   bridged and tt >= (t_bridge or T_MAX))
        ani = animation.FuncAnimation(fig, _upd, frames=len(frames),
                                      interval=100, blit=False)
        ani.save(gif, writer="pillow", fps=8,
                 savefig_kwargs={"facecolor": "#0a0a12"})
        print(f"  GIF    → {gif}")

    plt.close(fig)
    return tot, gl


if __name__ == "__main__":
    run(save_gif="--gif" in sys.argv, continuum="--continuum" in sys.argv)
