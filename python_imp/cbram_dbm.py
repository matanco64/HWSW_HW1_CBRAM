#!/usr/bin/env python3
"""
cbram_dbm.py — CBRAM filament via the Dielectric Breakdown Model.

The filament is a connected metallic cluster grown one cell at a time. Each step:
  1. solve Laplace for the potential, with the cluster pinned to the cathode
     potential (V=0) and the anode held at V_APPLIED;
  2. add one empty cell adjacent to the cluster, chosen with probability ∝ V^ETA.

ETA is the branchiness dial:
  ETA → 0   dense bush (DLA)
  ETA ≈ 1   branched lightning bolt
  ETA large single straight needle

The Jacobi Laplace solve is the computational hot loop — the target for the
C++ optimization stages — and is load-bearing: every growth decision reads the
freshly-solved field.

Usage:  python3 cbram_dbm.py [N] [ETA]
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import sys, os, time

# ── Parameters ────────────────────────────────────────────────────────────────
N            = int(sys.argv[1]) if len(sys.argv) > 1 else 200
ETA          = float(sys.argv[2]) if len(sys.argv) > 2 else 3.0   # 95/100 sweet spot
V_APPLIED    = 1.0
JACOBI_ITERS = 30          # Jacobi sweeps per growth step (warm-started)
SIGMA_MAX    = 20.0        # metallic conductivity written to the .bin
SIGMA_LOW    = 0.01
RNG_SEED     = 42

_cmap = LinearSegmentedColormap.from_list("cbram", [
    (0.00, (0.02, 0.02, 0.06)), (0.55, (0.75, 0.25, 0.00)),
    (0.80, (1.00, 0.75, 0.10)), (1.00, (1.00, 1.00, 1.00))])

# ── Laplace solve (Jacobi), cluster pinned to cathode potential ───────────────
def relax(V, metal, n):
    for _ in range(n):
        V[1:-1, 1:-1] = 0.25 * (V[:-2, 1:-1] + V[2:, 1:-1] +
                                V[1:-1, :-2] + V[1:-1, 2:])
        V[0, :]  = V_APPLIED                  # anode
        V[-1, :] = 0.0                         # cathode
        V[:, 0]  = V[:, 1]; V[:, -1] = V[:, -2]   # insulating side walls
        V[metal] = 0.0                         # cluster ≡ cathode potential
    return V

# ── Grow the filament ─────────────────────────────────────────────────────────
def run():
    rng = np.random.default_rng(RNG_SEED)
    metal = np.zeros((N, N), dtype=bool)
    metal[N - 2, N // 2] = True               # single seed at the cathode
    V = np.linspace(V_APPLIED, 0., N)[:, None] * np.ones((1, N))

    t0, step, bridged = time.time(), 0, False
    while True:
        relax(V, metal, JACOBI_ITERS)

        # candidate sites: empty cells 4-adjacent to the cluster (not electrodes)
        nbr = np.zeros((N, N), dtype=bool)
        nbr[1:-1, 1:-1] = (metal[:-2, 1:-1] | metal[2:, 1:-1] |
                           metal[1:-1, :-2] | metal[1:-1, 2:])
        cand = nbr & ~metal
        cand[0, :] = cand[-1, :] = False
        rows, cols = np.where(cand)
        if len(rows) == 0:
            break

        # DBM rule: choose one candidate with probability ∝ V^ETA
        w = np.maximum(V[rows, cols], 0.0) ** ETA
        if w.sum() <= 0:
            break
        i = rng.choice(len(rows), p=w / w.sum())
        r, c = int(rows[i]), int(cols[i])
        metal[r, c] = True
        step += 1

        if r <= 1:                            # reached the anode → bridged
            bridged = True
            break
        if step % 200 == 0:
            print(f"  step={step:5d}  cells={int(metal.sum()):5d}  "
                  f"tip_row={int(np.where(metal.any(1))[0].min())}  "
                  f"{time.time()-t0:.1f}s")

    relax(V, metal, 200)                       # clean field for the figure
    print(f"Done: {step} cells, bridged={bridged}, {time.time()-t0:.1f}s")
    return metal, V

# ── Output ────────────────────────────────────────────────────────────────────
def save(metal, V):
    outdir = os.path.dirname(os.path.abspath(__file__))
    sigma = np.where(metal, SIGMA_MAX, SIGMA_LOW)
    (sigma * 65536).astype(np.int32).tofile(os.path.join(outdir, "sigma_final_dbm.bin"))

    # field magnitude |∇φ| — what concentrates at the leading tip and drives growth
    dVr = np.zeros_like(V); dVr[1:-1, :] = V[:-2, :] - V[2:, :]
    dVc = np.zeros_like(V); dVc[:, 1:-1] = V[:, :-2] - V[:, 2:]
    Emag = np.hypot(dVr, dVc)

    plt.rcParams.update({"figure.facecolor": "#0a0a12"})
    fig, (a, b, c) = plt.subplots(1, 3, figsize=(16, 6))
    a.imshow(np.log1p(sigma), cmap=_cmap, vmin=0, vmax=np.log1p(SIGMA_MAX),
             interpolation="nearest")
    a.set_title(f"σ  (N={N}, η={ETA}, {int(metal.sum())} cells)",
                color="white", fontsize=9)
    b.imshow(V, cmap="RdYlBu_r", vmin=0, vmax=V_APPLIED, interpolation="nearest")
    b.set_title("potential φ  (Jacobi solve)", color="white", fontsize=9)
    c.imshow(Emag, cmap="inferno", vmin=0,
             vmax=float(np.percentile(Emag, 99.5)) + 1e-9, interpolation="nearest")
    c.set_title("|∇φ|  (drives growth)", color="white", fontsize=9)
    for ax in (a, b, c):
        ax.set_facecolor("#050508"); ax.tick_params(colors="gray", labelsize=6)
    fig.tight_layout()
    png = os.path.join(outdir, "filament_dbm.png")
    fig.savefig(png, dpi=150, facecolor="#0a0a12"); plt.close(fig)
    print(f"  → {png}")
    print(f"  → sigma_final_dbm.bin")

if __name__ == "__main__":
    save(*run())
