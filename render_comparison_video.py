#!/usr/bin/env python3
"""
render_comparison_video.py — time-synchronized comparison of any number of CBRAM stages.

Video time is proportional to real wall-clock time (scaled by TIME_SCALE), so the
faster stages visibly bridge first on screen. Stages are laid out horizontally.

Usage:
  python3 render_comparison_video.py [N [stages [out.mp4 [time_scale]]]]

  N           grid size (default 200)
  stages      comma-separated stage numbers (default "0,1")
  out.mp4     output path (default comparison_<stages>.mp4)
  time_scale  slow-down factor vs real time (default 3.0)

Examples:
  python3 render_comparison_video.py 200 0,1
  python3 render_comparison_video.py 200 0,1,2 comparison_all.mp4 4.0
"""

import csv, glob, os, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import LinearSegmentedColormap

try:
    import imageio_ffmpeg
    matplotlib.rcParams["animation.ffmpeg_path"] = imageio_ffmpeg.get_ffmpeg_exe()
except ImportError:
    pass

N          = int(sys.argv[1])             if len(sys.argv) > 1 else 200
stages_str = sys.argv[2]                  if len(sys.argv) > 2 else "0,1"
stages     = [int(s) for s in stages_str.split(",")]
_default_out = f"comparison_{'vs'.join(str(s) for s in stages)}.mp4"
out        = sys.argv[3]                  if len(sys.argv) > 3 else _default_out
TIME_SCALE = float(sys.argv[4])           if len(sys.argv) > 4 else 3.0

FPS       = 30
SIGMA_MAX = 20.0
Q         = 65536.0
HOLD_SECS = 2.0  # hold after all stages finish

_cmap = LinearSegmentedColormap.from_list("cbram", [
    (0.00, (0.02, 0.02, 0.06)), (0.55, (0.75, 0.25, 0.00)),
    (0.80, (1.00, 0.75, 0.10)), (1.00, (1.00, 1.00, 1.00))])

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BUILD_DIR  = os.path.join(SCRIPT_DIR, "build")


def load_stage(stage):
    frames_dir = os.path.join(BUILD_DIR, f"frames_stage{stage}")
    ts_path    = os.path.join(frames_dir, "frame_timestamps.csv")
    vfiles     = sorted(glob.glob(os.path.join(frames_dir, "frame_*_V.bin")))

    if not vfiles:
        sys.exit(f"No frames in {frames_dir}.\n"
                 f"Run: bash make_video.sh {stage} {N}")
    if not os.path.exists(ts_path):
        sys.exit(f"No timestamps at {ts_path}.\n"
                 f"Rebuild stages and re-run make_video.sh {stage} {N}")

    ts = []
    with open(ts_path) as f:
        for row in csv.DictReader(f):
            ts.append((int(row["step"]), float(row["elapsed_ms"])))

    while len(ts) < len(vfiles):   # closing frame may add one extra file
        ts.append(ts[-1])

    return vfiles, ts


def frame_idx_at(ts, t_ms):
    """Binary search: index of the latest frame with elapsed_ms <= t_ms."""
    lo, hi = 0, len(ts) - 1
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if ts[mid][1] <= t_ms:
            lo = mid
        else:
            hi = mid - 1
    return lo


def load_sigma(vf):
    return np.fromfile(vf.replace("_V.bin", "_S.bin"), dtype=np.int32).reshape(N, N) / Q


# ── Load all stages ───────────────────────────────────────────────────────────
data = [(s, *load_stage(s)) for s in stages]   # [(stage, vfiles, ts), ...]
finish_ms  = {s: ts[-1][1] for s, vfiles, ts in data}
total_ms   = max(finish_ms.values())
ref_ms     = finish_ms[stages[0]]              # stage 0 (or first listed) is the reference
speedups   = {s: ref_ms / finish_ms[s] for s in stages}

print(f"{'Stage':<8} {'Time (ms)':<12} {'Speedup vs stage ' + str(stages[0])}")
for s, vfiles, ts in data:
    marker = "  ← reference" if s == stages[0] else f"  {speedups[s]:.2f}×"
    print(f"  {s:<6} {finish_ms[s]:<12.0f}{marker}")

n_frames = int((total_ms * TIME_SCALE + HOLD_SECS * 1000) / 1000 * FPS) + 1

# ── Layout: 1 row × N cols, width scales with stage count ────────────────────
n = len(stages)
fig, axes = plt.subplots(1, n, figsize=(6 * n, 6.5))
if n == 1:
    axes = [axes]
fig.patch.set_facecolor("#0a0a12")
plt.rcParams.update({"figure.facecolor": "#0a0a12"})


def render(vf_idx):
    real_ms = (vf_idx / FPS * 1000) / TIME_SCALE

    for ax, (stage, vfiles, ts) in zip(axes, data):
        ax.clear()
        ax.set_facecolor("#050508")
        ax.tick_params(colors="gray", labelsize=6)
        for sp in ax.spines.values():
            sp.set_color("#333")

        eff_ms        = min(real_ms, finish_ms[stage])
        idx           = frame_idx_at(ts, eff_ms)
        step, elapsed = ts[idx]
        bridged       = real_ms >= finish_ms[stage]

        ax.imshow(np.log1p(load_sigma(vfiles[idx])), cmap=_cmap,
                  vmin=0, vmax=np.log1p(SIGMA_MAX),
                  interpolation="nearest", origin="upper")

        color  = "#44ff88" if bridged else "white"
        suffix = f"\n✓ BRIDGED @ {finish_ms[stage]:.0f} ms" if bridged else ""
        spdup  = f"  ({speedups[stage]:.2f}×)" if stage != stages[0] else "  (ref)"
        ax.set_title(f"Stage {stage}{spdup}\nstep {step}  {elapsed:.0f} ms{suffix}",
                     color=color, fontsize=9, pad=4)

    speedup_str = "  |  ".join(
        f"S{s}: {speedups[s]:.2f}×" for s in stages if s != stages[0]
    )
    fig.suptitle(
        f"t = {real_ms:.0f} ms   (real-time × 1/{TIME_SCALE:.0f})    {speedup_str}",
        color="#aaaaaa", fontsize=9, y=1.01,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.97])


print(f"Encoding {n_frames} video frames → {out}")
ani = animation.FuncAnimation(fig, render, frames=n_frames, interval=1000 / FPS)
ani.save(out, writer="ffmpeg", fps=FPS, dpi=130,
         savefig_kwargs={"facecolor": "#0a0a12"})
plt.close(fig)
print(f"  → {os.path.abspath(out)}")
