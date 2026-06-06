#!/usr/bin/env python3
"""
render_dbm_video.py — 3-panel growth video from dbm_stage0 frame dumps.

The C++ stage 0 writes per-frame V and σ binaries into frames_stage0/
(frame_NNNNNN_V.bin, frame_NNNNNN_S.bin, N×N int32, Q16.16). This reads each
pair and renders σ | φ | |∇φ| side by side, then encodes to MP4 via ffmpeg.

Usage:  python3 render_dbm_video.py [N] [frames_dir] [out.mp4]
"""

import numpy as np, glob, os, sys
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import LinearSegmentedColormap

N          = int(sys.argv[1]) if len(sys.argv) > 1 else 200
frames_dir = sys.argv[2] if len(sys.argv) > 2 else "build/frames_stage0"
out        = sys.argv[3] if len(sys.argv) > 3 else "filament_stage0.mp4"
SIGMA_MAX  = 20.0
Q          = 65536.0

_cmap = LinearSegmentedColormap.from_list("cbram", [
    (0.00, (0.02, 0.02, 0.06)), (0.55, (0.75, 0.25, 0.00)),
    (0.80, (1.00, 0.75, 0.10)), (1.00, (1.00, 1.00, 1.00))])

vfiles = sorted(glob.glob(os.path.join(frames_dir, "frame_*_V.bin")))
if not vfiles:
    sys.exit(f"no frames in {frames_dir} (run cbram_stage0 first)")
print(f"{len(vfiles)} frames from {frames_dir}")

def load(vf):
    V = np.fromfile(vf, dtype=np.int32).reshape(N, N) / Q
    S = np.fromfile(vf.replace("_V.bin", "_S.bin"), dtype=np.int32).reshape(N, N) / Q
    dVr = np.zeros_like(V); dVr[1:-1, :] = V[:-2, :] - V[2:, :]
    dVc = np.zeros_like(V); dVc[:, 1:-1] = V[:, :-2] - V[:, 2:]
    return V, S, np.hypot(dVr, dVc)

# Fix |∇φ| colour scale from the final (most developed) frame so it doesn't flicker.
e_vmax = float(np.percentile(load(vfiles[-1])[2], 99.5)) + 1e-9

plt.rcParams.update({"figure.facecolor": "#0a0a12"})
fig, (ax_s, ax_v, ax_e) = plt.subplots(1, 3, figsize=(16, 6))

def render(i):
    V, S, E = load(vfiles[i])
    for ax in (ax_s, ax_v, ax_e):
        ax.clear(); ax.set_facecolor("#050508"); ax.tick_params(colors="gray", labelsize=6)
    ax_s.imshow(np.log1p(S), cmap=_cmap, vmin=0, vmax=np.log1p(SIGMA_MAX), interpolation="nearest")
    ax_s.set_title(f"σ  filament  (frame {i+1}/{len(vfiles)}, {int((S>=SIGMA_MAX/2).sum())} cells)",
                   color="white", fontsize=9)
    ax_v.imshow(V, cmap="RdYlBu_r", vmin=0, vmax=1.0, interpolation="nearest")
    ax_v.set_title("potential φ  (Jacobi solve)", color="white", fontsize=9)
    ax_e.imshow(E, cmap="inferno", vmin=0, vmax=e_vmax, interpolation="nearest")
    ax_e.set_title("|∇φ|  (drives growth)", color="white", fontsize=9)
    fig.tight_layout()

ani = animation.FuncAnimation(fig, render, frames=len(vfiles), interval=50)
ani.save(out, writer="ffmpeg", fps=30, dpi=120,
         savefig_kwargs={"facecolor": "#0a0a12"})
plt.close(fig)
print(f"  → {os.path.abspath(out)}")
