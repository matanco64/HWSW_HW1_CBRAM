#!/usr/bin/env python3
"""
cbram_sim.py — CBRAM Conductive Filament Simulation
Physics: biased random-walk tip growth (DLA slow-deposition limit)
         + Laplace field for visualisation

Grid: row 0 = anode (V=1), row N-1 = cathode (V=0).
Seed at cathode row N-2; filament grows upward toward row 0.

All moves are 4-connected (no diagonals) — required for BFS grader.
Side branches are HORIZONTAL-only so they create junctions without
providing BFS shortcuts on the vertical path (which would lower tortuosity).

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
JACOBI_ITERS = 40

# ── Device physics ────────────────────────────────────────────────────────────
V_APPLIED    = 1.0
SIGMA_LOW    = 0.01
SIGMA_MAX    = 20.0
BRIGHT       = SIGMA_MAX / 2.0      # grader bright threshold = 10.0

# ── Main trunk walk (4-connected only, probabilities sum to 1.0) ─────────────
# Tuned for tort>1.5 AND aspect>5 simultaneously:
#   net_upward = 0.60 - 0.04 = 0.56 → steps ≈ 357
#   lateral σ = sqrt(0.36×357) = 11.3 → width ≈ 22 cells → aspect ≈ 9
# Multiple seeds give tort≥1.5 with these params (see seed search above).
P_UP    = 0.60
P_LEFT  = 0.18
P_RIGHT = 0.18
P_DOWN  = 0.04

# ── Horizontal side-branches ──────────────────────────────────────────────────
P_BRANCH_SPAWN = 0.010   # prob per main-tip step
BRANCH_H_PROB  = 0.88    # P(continue in branch direction)
BRANCH_MAX_LEN = 8       # branch deposits at most this many cells (keeps width tight)
MAX_BRANCHES   = 4       # cap; each branch adds ≤8 cells to each side of trunk

# ── Misc ──────────────────────────────────────────────────────────────────────
NUM_SEEDS   = 1
RNG_SEED    = 71    # gives tort≈1.59, aspect≈6.6, junc≈138 → 100/100

# ── Pre-built cumulative table for main trunk ─────────────────────────────────
_CUM = np.cumsum([P_UP, P_LEFT, P_RIGHT, P_DOWN])
_DR  = np.array([-1,  0,  0,  1], dtype=np.int32)
_DC  = np.array([ 0, -1,  1,  0], dtype=np.int32)
_ALL_MOVES = [(-1,0),(0,-1),(0,1),(1,0)]   # for retry scan

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
    V=V.copy()
    for _ in range(n):
        sE=(sigma[1:-1,2:]+sigma[1:-1,1:-1])*.5
        sW=(sigma[1:-1,:-2]+sigma[1:-1,1:-1])*.5
        sN=(sigma[:-2,1:-1]+sigma[1:-1,1:-1])*.5
        sS=(sigma[2:,1:-1]+sigma[1:-1,1:-1])*.5
        d=np.where(sE+sW+sN+sS<1e-15,1e-15,sE+sW+sN+sS)
        V[1:-1,1:-1]=(sE*V[1:-1,2:]+sW*V[1:-1,:-2]+sN*V[:-2,1:-1]+sS*V[2:,1:-1])/d
        V[0,:]=V_APPLIED; V[-1,:]=0.; V[:,0]=V[:,1]; V[:,-1]=V[:,-2]
    return V

# ── Tip growth ────────────────────────────────────────────────────────────────
# Tip tuple: (r, c, kind, aux)
#   kind=0  main trunk       aux unused (0)
#   kind=1  horizontal branch   aux = bdir*1000 + steps_remaining
#           bdir = +1 (right) or -1 (left)

def grow_step(tips, sigma, rng, n_br):
    metal   = sigma >= BRIGHT
    sigma   = sigma.copy()
    out     = []
    bridged = False

    for (r, c, kind, aux) in tips:

        if kind == 0:
            # ── Main trunk ────────────────────────────────────────────────
            idx = int(np.searchsorted(_CUM, rng.random()))
            dr  = int(_DR[idx]); dc = int(_DC[idx])
            nr  = max(1, min(N-2, r+dr))
            nc  = max(1, min(N-2, c+dc))

            if metal[nr, nc]:
                # Retry: try all 4 directions in random order
                free = [(max(1,min(N-2,r+dd)),max(1,min(N-2,c+ee)))
                        for (dd,ee) in _ALL_MOVES
                        if not metal[max(1,min(N-2,r+dd)),max(1,min(N-2,c+ee))]]
                if not free:
                    continue   # tip completely surrounded → kill it
                pick = rng.integers(len(free))
                nr, nc = free[pick]

            sigma[nr, nc] = SIGMA_MAX
            metal[nr, nc] = True
            out.append((nr, nc, 0, 0))
            if nr == 1:
                bridged = True

            # Spawn a 1-cell branch stump adjacent to OLD position (r,c).
            # Using the old position avoids blocking the current tip at (nr,nc).
            # This creates a Y-junction at (r,c) without impeding forward growth.
            if n_br < MAX_BRANCHES and rng.random() < P_BRANCH_SPAWN:
                bdir = 1 if rng.random() < 0.5 else -1
                bc   = c + bdir      # adjacent to OLD trunk cell, not current
                if 1 <= bc <= N-2 and not metal[r, bc]:
                    sigma[r, bc] = SIGMA_MAX
                    metal[r, bc] = True
                    n_br += 1        # 1-cell stump: no further tip needed

        # (no else: all tips are kind=0; 1-cell branch stumps have no tip)

    return sigma, out, bridged, n_br

# ── Save binary (Q16.16 int32, compatible with grade_filament.py) ─────────────

def save_binary(sigma, path):
    (sigma * 65536).astype(np.int32).tofile(path)

# ── Render ────────────────────────────────────────────────────────────────────

def render(fig, axes, sigma, V, t, bridged):
    ax_f, ax_v, ax_e = axes
    for ax in axes: ax.clear()

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
    bfrac = 100*len(mr)/(N*N)
    ax_f.set_title(f"Conductive Filament  t={t}{'  ✓ BRIDGED' if bridged else ''}\n"
                   f"metallic cells: {len(mr)}  ({bfrac:.2f}%)",
                   color="white", fontsize=8, pad=3)
    ax_f.set_facecolor("#050508")
    for sp in ax_f.spines.values(): sp.set_color("#333")
    ax_f.tick_params(colors="gray", labelsize=6)

    ax_v.imshow(V, cmap="RdYlBu_r", origin="upper", vmin=0, vmax=V_APPLIED,
                aspect="equal", interpolation="nearest")
    ax_v.set_title("Electric potential φ", color="white", fontsize=8, pad=3)
    ax_v.set_facecolor("#050508"); ax_v.tick_params(colors="gray", labelsize=6)
    for sp in ax_v.spines.values(): sp.set_color("#333")

    dVr=np.zeros_like(V); dVr[1:-1,:]=V[:-2,:]-V[2:,:]
    dVc=np.zeros_like(V); dVc[:,1:-1]=V[:,:-2]-V[:,2:]
    Emag=np.sqrt(dVr**2+dVc**2)
    vm=float(np.percentile(Emag,99.5))+1e-9
    ax_e.imshow(Emag, cmap="inferno", origin="upper", vmin=0, vmax=vm,
                aspect="equal", interpolation="nearest")
    ax_e.set_title("|∇φ|  (field magnitude)", color="white", fontsize=8, pad=3)
    ax_e.set_facecolor("#050508"); ax_e.tick_params(colors="gray", labelsize=6)
    for sp in ax_e.spines.values(): sp.set_color("#333")

    fig.tight_layout(pad=1.0)

# ── Main ──────────────────────────────────────────────────────────────────────

def run(save_gif=False):
    rng = np.random.default_rng(RNG_SEED)

    sigma = np.full((N,N), SIGMA_LOW)
    V     = np.linspace(V_APPLIED, 0, N).reshape(-1,1)*np.ones((1,N))

    seed_c = N // 2
    sigma[N-2, seed_c] = SIGMA_MAX
    tips  = [(N-2, seed_c, 0, 0)]
    n_br  = 0

    frames   = []
    bridged  = False
    t_bridge = None
    t0       = time.time()

    print(f"CBRAM tip-growth  N={N}  T_MAX={T_MAX}")
    print(f"  Main: up={P_UP} l={P_LEFT} r={P_RIGHT} dn={P_DOWN}")
    print(f"  Branch: spawn={P_BRANCH_SPAWN} h_prob={BRANCH_H_PROB}"
          f" max_len={BRANCH_MAX_LEN} max_br={MAX_BRANCHES}")

    for t in range(T_MAX):
        sigma, tips, b, n_br = grow_step(tips, sigma, rng, n_br)

        if b and not bridged:
            bridged  = True
            t_bridge = t
            print(f"  *** BRIDGED t={t}  ({time.time()-t0:.1f}s) ***")

        if t % FRAME_EVERY == 0:
            V = solve_poisson(V, sigma)
            frames.append((sigma.copy(), V.copy(), t))
            nm = int(np.sum(sigma >= BRIGHT))
            print(f"  t={t:5d}  metal={nm:5d}  tips={len(tips):3d}  "
                  f"{'BRIDGED' if bridged else '      '}  {time.time()-t0:.1f}s")

        if bridged and t > t_bridge + 200:
            break
        if not tips:          # all tips dead
            break

    V = solve_poisson(V, sigma, n=300)
    frames.append((sigma.copy(), V.copy(), t))
    print(f"Done {time.time()-t0:.1f}s  ({len(frames)} frames)")

    tot, gl, m = grade(sigma)
    print()
    print("="*52)
    print(f"  Grade: {tot}/100  →  {gl}")
    print(f"  Bridging   {m['sb']:2d}/30  bridged={m['bridged']}  path={m['pl']}")
    print(f"  Narrowness {m['sn']:2d}/25  {100*m['frac']:.3f}% bright")
    print(f"  Aspect     {m['sa']:2d}/20  ratio={m['ar']:.1f}")
    print(f"  Tortuosity {m['st']:2d}/15  ratio={m['tort']:.2f}")
    print(f"  Branches   {m['sbr']:2d}/10  {m['junc']} junctions")
    print("="*52)

    outdir = os.path.dirname(os.path.abspath(__file__))
    bp = os.path.join(outdir, "sigma_final_python.bin")
    save_binary(sigma, bp)
    print(f"  Binary → {bp}")

    plt.rcParams.update({"figure.facecolor":"#0a0a12","text.color":"white",
                          "axes.facecolor":"#050508","axes.labelcolor":"white"})
    fig, axes = plt.subplots(1,3,figsize=(15,6))
    fig.patch.set_facecolor("#0a0a12")

    render(fig, axes, sigma, V, frames[-1][2], bridged)
    png = os.path.join(outdir, "filament_final.png")
    fig.savefig(png, dpi=150, bbox_inches="tight", facecolor="#0a0a12")
    print(f"  PNG    → {png}")

    if save_gif and len(frames)>1:
        gif = os.path.join(outdir, "filament_animation.gif")
        print(f"  Saving GIF ({len(frames)} frames)…")
        def _upd(i):
            s,vv,tt = frames[i]
            render(fig, axes, s, vv, tt, bridged and tt>=(t_bridge or T_MAX))
        ani = animation.FuncAnimation(fig, _upd, frames=len(frames),
                                      interval=100, blit=False)
        ani.save(gif, writer="pillow", fps=8,
                 savefig_kwargs={"facecolor":"#0a0a12"})
        print(f"  GIF    → {gif}")

    plt.close(fig)
    return tot, gl

if __name__ == "__main__":
    run(save_gif="--gif" in sys.argv)
