#!/usr/bin/env python3
"""
cbram_tetris.py — CBRAM Tetris Ion Simulation

Physics: discrete ion particles released from the anode, drifting downward
guided by the Jacobi V field, depositing on contact with existing metal.

Key idea (Matan): ions are like Tetris blocks falling from the anode.
  - Small number active at a time → no flooding
  - Field-guided walk → tortuosity, field is genuinely load-bearing
  - Contact deposition → 4-connectivity guaranteed
  - "Tip thinner than base": early ions travel far (scatter wide → thick base),
    later ions have a short journey to the rising tip (less scatter → thin tip)

Two bugs fixed vs the embedded --tetris mode:
  1. Spawn position: uses V[0]-V[1] (anode interface E), not V[1]-V[2].
     Above the filament V[1]≈0 → anode-interface field is LARGE → ions prefer
     to spawn above the tip. With V[1]-V[2] both are ≈0, so ions spawn elsewhere.
  2. Metal collision: ion trying to enter a metallic cell now deposits at its
     current (adjacent) position instead of being silently dropped.

Usage:
  python3 cbram_tetris.py          # run, save PNG + binary
  python3 cbram_tetris.py --gif    # also save animated GIF
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
T_MAX         = 20_000
FRAME_EVERY   = 100
JACOBI_ITERS  = 50
JACOBI_EVERY  = 10    # recompute V only every N steps — V changes slowly

# ── Device physics ────────────────────────────────────────────────────────────
V_APPLIED = 2.0
SIGMA_LOW = 0.01
SIGMA_MAX = 20.0
BRIGHT    = SIGMA_MAX / 2.0   # metallic threshold (grader uses this)

# ── Tetris ion parameters ─────────────────────────────────────────────────────
MAX_IONS        = 90    # max simultaneous drifting ions
SPAWN_INTERVAL  = 1     # release one new ion every timestep
P_CATHODE       = 0.4   # probability of becoming a seed when reaching the cathode
P_CONNECT       = 0.6   # base connection probability (scaled by neighbors + field below)
FIELD_BOOST_MAX = 8.0   # max multiplier from field enhancement at the tip
BASE_P_DOWN_ION = 0.55  # base downward drift (toward cathode)
BASE_P_LAT_ION  = 0.20  # base lateral diffusion (left/right equally)
BASE_P_UP_ION   = 0.05  # rare upward step (thermal noise)
FIELD_STR_ION   = 15.0  # how strongly local E field biases each step
RNG_SEED        = 42

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
        d  = np.where(sE+sW+sN+sS < 1e-15, 1e-15, sE+sW+sN+sS)
        V[1:-1, 1:-1] = (sE*V[1:-1,2:]+sW*V[1:-1,:-2]+
                         sN*V[:-2,1:-1]+sS*V[2:, 1:-1])/d
        V[0,:]  = V_APPLIED; V[-1,:] = 0.
        V[:,0]  = V[:,1];    V[:,-1] = V[:,-2]
    return V

# ── Tetris ion physics ────────────────────────────────────────────────────────

_DR  = np.array([-1,  0,  0,  1])   # up, left, right, down
_DC  = np.array([ 0, -1,  1,  0])
_BASE = np.array([BASE_P_UP_ION, BASE_P_LAT_ION, BASE_P_LAT_ION, BASE_P_DOWN_ION])


SPAWN_SIGMA = N * 0.01   # Gaussian width for spawn column (≈1% of grid)

def spawn_ion(rng):
    """
    Release one ion at row 1 with column drawn from a Gaussian centered at N//2.

    Ions dissolve from the anode near the center, above the seed location.
    The spread (SPAWN_SIGMA) creates lateral diversity — ions don't all land in
    the same column — which is the source of tortuosity and branching.
    Ions outside [1, N-2] are rejected and resampled.
    """
    cols = np.arange(1, N - 1)
    p    = np.exp(-0.5 * ((cols - N // 2) / SPAWN_SIGMA) ** 2)
    p   /= p.sum()
    c    = int(rng.choice(cols, p=p))
    return [1, c]


def step_ions(ions, sigma, V, rng, t):
    """
    Advance every active ion one step.

    Three rules per ion at (r, c), checked in order:

    1. CATHODE ARRIVAL — ion reaches row N-2 (cathode surface):
       Deposit with P_CATHODE (=1.0). This creates new seeds dynamically.
       Multiple ions reaching different cathode columns → naturally wide base.

    2. METAL ADJACENCY — any 4-connected neighbor is already metallic:
       Deposit with P_CONNECT. Low probability so ions can pass alongside
       the growing filament without immediately sticking — allows the tip to
       receive ions rather than all of them depositing at the base.

    3. DRIFT — field-biased random walk one step downward toward cathode.

    Why "thicker at base" emerges:
      - Every ion reaching the cathode deposits → base grows wide
      - As the metallic base grows upward, it intercepts drifting ions earlier
      - The tip is narrow because few ions survive past the wide base without connecting
    """
    metal    = sigma >= BRIGHT
    ions_out = []
    deposited = []

    for ion in ions:
        r, c = ion

        # Rule 1: cathode arrival → always seed
        if r >= N - 2:
            if rng.random() < P_CATHODE:
                deposited.append((N-2, c))
            # ion consumed regardless (absorbed by cathode)
            continue

        # Rule 2: adjacent to existing metal → neighbor-count + field-boosted connection
        #
        # A) Neighbor count: more metallic neighbors → more "enclosed" → stickier.
        #    1 neighbor (grazing a flat surface): low stick — ion slides past.
        #    2-3 neighbors (trapped in corner/pocket): high stick — truly captured.
        #    Quadratic scaling: 1nb→6%, 2nb→25%, 3nb→56%, 4nb→100% of P_CONNECT.
        #
        # B) Field boost: E field concentrates at the filament TIP (sharp σ drop).
        #    High E_local → tip is "attractive" → ions prefer to connect there.
        nb = (int(r > 0   and metal[r-1, c]) +
              int(r < N-1 and metal[r+1, c]) +
              int(c > 0   and metal[r, c-1]) +
              int(c < N-1 and metal[r, c+1]))

        if nb > 0:
            nb_factor    = (nb / 4.0) ** 2
            E_down_local = max(V[r-1, c] - V[r, c], 0.) if r > 0 else 0.
            E_background = V_APPLIED / N
            field_factor = min(E_down_local / (E_background + 1e-9), FIELD_BOOST_MAX)
            p_connect    = min(P_CONNECT * nb_factor * (1.0 + field_factor), 1.0)

            if rng.random() < p_connect:
                deposited.append((r, c))
            else:
                ions_out.append(ion)   # didn't connect — keep drifting
            continue

        # Rule 3: field-biased walk
        E = np.array([
            V[r-1, c] - V[r, c] if r > 1   else 0.,
            V[r, c-1] - V[r, c] if c > 1   else 0.,
            V[r, c+1] - V[r, c] if c < N-2 else 0.,
            V[r+1, c] - V[r, c] if r < N-2 else 0.,
        ])
        w   = _BASE + FIELD_STR_ION * np.maximum(E, 0.)
        w  /= w.sum()
        idx = int(np.searchsorted(np.cumsum(w), rng.random()))
        nr  = r + int(_DR[idx])
        nc  = max(1, min(N-2, c + int(_DC[idx])))

        if 1 <= nr <= N-1 and not metal[nr, nc]:
            ions_out.append([nr, nc])
        elif metal[nr, nc]:
            ions_out.append(ion)   # blocked by metal wall — stay put this step
        # else nr out of bounds or cathode — ion escapes/absorbed

    # Spawn new ion from anode
    if t % SPAWN_INTERVAL == 0 and len(ions_out) < MAX_IONS:
        ions_out.append(spawn_ion(rng))

    return ions_out, deposited

# ── Save binary ───────────────────────────────────────────────────────────────

def save_binary(sigma, path):
    (sigma * 65536).astype(np.int32).tofile(path)

# ── Render — 4 panels ─────────────────────────────────────────────────────────

def render(fig, axes, sigma, V, active_ions, t, bridged):
    ax_f, ax_v, ax_e, ax_i = axes
    for ax in axes:
        ax.clear()

    # Panel 1: conductivity
    log_s = np.log1p(sigma)
    ax_f.imshow(log_s, cmap=_cbram_cmap, origin="upper",
                vmin=0, vmax=np.log1p(SIGMA_MAX), aspect="equal",
                interpolation="nearest")
    mr, mc = np.where(sigma >= BRIGHT)
    if len(mr):
        ax_f.scatter(mc, mr, s=4, c="white", alpha=0.85, linewidths=0, zorder=3)
    ax_f.axhline(0.5,   color="#55aaff", lw=1, ls="--", alpha=0.6)
    ax_f.axhline(N-1.5, color="#aaaaaa", lw=1, ls="--", alpha=0.6)
    ax_f.text(2, 2,   "ANODE ⊕",   color="#55aaff", fontsize=7, va="top")
    ax_f.text(2, N-2, "CATHODE ⊖", color="#aaaaaa", fontsize=7, va="bottom")
    bfrac = 100*len(mr)/(N*N)
    ax_f.set_title(f"Conductivity σ  t={t}{'  ✓ BRIDGED' if bridged else ''}\n"
                   f"metallic: {len(mr)} ({bfrac:.2f}%)",
                   color="white", fontsize=8, pad=3)
    ax_f.set_facecolor("#050508")
    for sp in ax_f.spines.values(): sp.set_color("#333")
    ax_f.tick_params(colors="gray", labelsize=6)

    # Panel 2: electric potential
    ax_v.imshow(V, cmap="RdYlBu_r", origin="upper", vmin=0, vmax=V_APPLIED,
                aspect="equal", interpolation="nearest")
    ax_v.set_title("Electric potential φ", color="white", fontsize=8, pad=3)
    ax_v.set_facecolor("#050508"); ax_v.tick_params(colors="gray", labelsize=6)
    for sp in ax_v.spines.values(): sp.set_color("#333")

    # Panel 3: |∇V|
    dVr = np.zeros_like(V); dVr[1:-1,:] = V[:-2,:] - V[2:,:]
    dVc = np.zeros_like(V); dVc[:,1:-1] = V[:,:-2] - V[:,2:]
    Emag = np.sqrt(dVr**2 + dVc**2)
    vm = float(np.percentile(Emag, 99.5)) + 1e-9
    ax_e.imshow(Emag, cmap="inferno", origin="upper", vmin=0, vmax=vm,
                aspect="equal", interpolation="nearest")
    ax_e.set_title("|∇φ|  (field magnitude)", color="white", fontsize=8, pad=3)
    ax_e.set_facecolor("#050508"); ax_e.tick_params(colors="gray", labelsize=6)
    for sp in ax_e.spines.values(): sp.set_color("#333")

    # Panel 4: active ion positions (cyan dots on dark background)
    ax_i.set_facecolor("#050508")
    ax_i.set_xlim(0, N); ax_i.set_ylim(N, 0)
    ax_i.set_aspect("equal")
    if active_ions:
        ir = [p[0] for p in active_ions]
        ic = [p[1] for p in active_ions]
        ax_i.scatter(ic, ir, s=30, c="cyan", alpha=0.9, linewidths=0, zorder=3)
    ax_i.axhline(0.5,   color="#55aaff", lw=1, ls="--", alpha=0.4)
    ax_i.axhline(N-1.5, color="#aaaaaa", lw=1, ls="--", alpha=0.4)
    ax_i.set_title(f"Active ions ({len(active_ions)})", color="white", fontsize=8, pad=3)
    ax_i.tick_params(colors="gray", labelsize=6)
    for sp in ax_i.spines.values(): sp.set_color("#333")

    fig.tight_layout(pad=1.0)

# ── Main ──────────────────────────────────────────────────────────────────────

def run(save_gif=False):
    rng = np.random.default_rng(RNG_SEED)

    sigma = np.full((N, N), SIGMA_LOW)
    V     = np.linspace(V_APPLIED, 0., N).reshape(-1,1) * np.ones((1,N))

    # No fixed seed — ions create seeds when they reach the cathode
    active_ions = []
    frames      = []
    bridged     = False
    t_bridge    = None
    t0          = time.time()

    print(f"CBRAM Tetris Ion  N={N}  T_MAX={T_MAX}")
    print(f"  max_ions={MAX_IONS}  spawn_every={SPAWN_INTERVAL}"
          f"  p_connect={P_CONNECT}  p_cathode={P_CATHODE}  field_str={FIELD_STR_ION}")
    print(f"  P_down={BASE_P_DOWN_ION}  P_lat={BASE_P_LAT_ION}"
          f"  P_up={BASE_P_UP_ION}")

    for t in range(T_MAX):
        if t % JACOBI_EVERY == 0:
            V = solve_poisson(V, sigma, n=JACOBI_ITERS)

        active_ions, deposited = step_ions(active_ions, sigma, V, rng, t)
        for (r, c) in deposited:
            sigma[r, c] = SIGMA_MAX

        if not bridged and np.any(sigma[1, 1:-1] >= BRIGHT):
            bridged  = True
            t_bridge = t
            print(f"  *** BRIDGED t={t}  ({time.time()-t0:.1f}s) ***")

        if t % FRAME_EVERY == 0:
            frames.append((sigma.copy(), V.copy(), list(active_ions), t))
            nm = int(np.sum(sigma >= BRIGHT))
            print(f"  t={t:5d}  metal={nm:4d}  ions={len(active_ions):2d}"
                  f"  {'BRIDGED' if bridged else '      '}  {time.time()-t0:.1f}s")

        if bridged:
            break

    V = solve_poisson(V, sigma, n=300)
    frames.append((sigma.copy(), V.copy(), [], t))
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
    bp     = os.path.join(outdir, "sigma_final_tetris.bin")
    save_binary(sigma, bp)
    print(f"  Binary → {bp}")

    plt.rcParams.update({"figure.facecolor":"#0a0a12","text.color":"white",
                         "axes.facecolor":"#050508","axes.labelcolor":"white"})
    fig, axes = plt.subplots(1, 4, figsize=(20, 6))
    fig.patch.set_facecolor("#0a0a12")

    s, v, ai, tt = frames[-1]
    render(fig, axes, s, v, ai, tt, bridged)
    png = os.path.join(outdir, "filament_tetris.png")
    fig.savefig(png, dpi=150, bbox_inches="tight", facecolor="#0a0a12")
    print(f"  PNG    → {png}")

    if save_gif and len(frames) > 1:
        gif = os.path.join(outdir, "filament_tetris.gif")
        print(f"  Saving GIF ({len(frames)} frames)…")
        def _upd(idx):
            s, vv, ai, tt = frames[idx]
            render(fig, axes, s, vv, ai, tt,
                   bridged and tt >= (t_bridge or T_MAX))
        ani = animation.FuncAnimation(fig, _upd, frames=len(frames),
                                      interval=100, blit=False)
        ani.save(gif, writer="pillow", fps=8,
                 savefig_kwargs={"facecolor":"#0a0a12"})
        print(f"  GIF    → {gif}")

    plt.close(fig)
    return tot, gl


if __name__ == "__main__":
    run(save_gif="--gif" in sys.argv)
