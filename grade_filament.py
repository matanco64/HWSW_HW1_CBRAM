#!/usr/bin/env python3
"""
grade_filament.py — CBRAM Filament Quality Grader
Usage: python3 grade_filament.py <sigma_final.bin> [N]
  N is optional; inferred from sqrt(filesize/4) if omitted.

Reads an N×N grid of int32_t Q16.16 sigma values (row-major binary).
Scores the filament shape against 5 lightning-filament criteria (100 pts total).
Prints a human-readable report and a machine-parseable __GRADE_SUMMARY__ line.
"""
import sys, os, math, struct
from collections import deque

# Must match physics.h constants (in Q16.16 units)
SHIFT          = 16
ONE            = 1 << SHIFT
SIGMA_MAX_Q    = 20 * ONE          # SIGMA_MAX = 20.0
SIGMA_HIGH_Q   = 10 * ONE          # SIGMA_HIGH = 10.0
BRIGHT_THRESH  = SIGMA_MAX_Q // 2  # > 10.0 → "bright" (metallic filament cell)


def load_sigma(path, N=None):
    with open(path, "rb") as f:
        data = f.read()
    n_cells = len(data) // 4
    if N is None:
        N = int(math.isqrt(n_cells))
        if N * N != n_cells:
            sys.exit(f"ERROR: file size {len(data)} bytes is not (N*N*4). "
                     f"Closest N={N} gives {N*N*4} bytes.")
    else:
        if N * N * 4 != len(data):
            sys.exit(f"ERROR: expected {N*N*4} bytes for N={N}, got {len(data)}.")
    sigma = struct.unpack(f"<{n_cells}i", data)
    return sigma, N


def is_bright(s):
    return s > BRIGHT_THRESH


# ---------------------------------------------------------------------------
# Metric 1 — Bridging (0–30 pts)
# BFS: connected path of bright cells from row 1 to row N-2.
# ---------------------------------------------------------------------------
def bfs_bridge(sigma, N):
    dist = {}
    prev = {}
    queue = deque()

    for c in range(1, N - 1):
        idx = 1 * N + c
        if is_bright(sigma[idx]):
            queue.append((1, c))
            dist[(1, c)] = 0
            prev[(1, c)] = None

    found, end_cell = False, None
    while queue:
        r, c = queue.popleft()
        if r == N - 2:
            found, end_cell = True, (r, c)
            break
        for dr, dc in ((-1, 0), (1, 0), (0, -1), (0, 1)):
            nr, nc = r + dr, c + dc
            if 1 <= nr <= N - 2 and 1 <= nc <= N - 2 and (nr, nc) not in dist:
                if is_bright(sigma[nr * N + nc]):
                    dist[(nr, nc)] = dist[(r, c)] + 1
                    prev[(nr, nc)] = (r, c)
                    queue.append((nr, nc))

    if not found:
        return False, 0, []

    path = []
    cur = end_cell
    while cur is not None:
        path.append(cur)
        cur = prev[cur]
    path.reverse()
    return True, len(path), path


# ---------------------------------------------------------------------------
# Metric 2 — Narrowness (0–25 pts)
# Fraction of interior cells that are bright.
# ---------------------------------------------------------------------------
def compute_narrowness(sigma, N):
    interior = (N - 2) * (N - 2)
    bright = sum(
        1 for r in range(1, N - 1)
          for c in range(1, N - 1)
          if is_bright(sigma[r * N + c])
    )
    return bright / interior if interior > 0 else 1.0


# ---------------------------------------------------------------------------
# Metric 3 — Aspect ratio (0–20 pts)
# Height / width of bounding box of all bright cells.
# ---------------------------------------------------------------------------
def compute_aspect_ratio(sigma, N):
    bright_rows, bright_cols = [], []
    for r in range(1, N - 1):
        for c in range(1, N - 1):
            if is_bright(sigma[r * N + c]):
                bright_rows.append(r)
                bright_cols.append(c)
    if not bright_rows:
        return 0.0
    h = max(bright_rows) - min(bright_rows) + 1
    w = max(bright_cols) - min(bright_cols) + 1
    return h / max(w, 1)


# ---------------------------------------------------------------------------
# Metric 4 — Tortuosity (0–15 pts)
# BFS shortest bridge path length / straight-line distance (N-3 cells).
# ---------------------------------------------------------------------------
def compute_tortuosity(path_len, N):
    straight = N - 3  # row 1 → row N-2
    if straight <= 0:
        return 1.0
    return path_len / straight


# ---------------------------------------------------------------------------
# Metric 5 — Branch count (0–10 pts)
# Number of bright cells with 3+ bright neighbors (junction points).
# ---------------------------------------------------------------------------
def count_junctions(sigma, N):
    junctions = 0
    for r in range(1, N - 1):
        for c in range(1, N - 1):
            if not is_bright(sigma[r * N + c]):
                continue
            bn = sum(
                1 for dr, dc in ((-1, 0), (1, 0), (0, -1), (0, 1))
                if 1 <= r + dr <= N - 2 and 1 <= c + dc <= N - 2
                   and is_bright(sigma[(r + dr) * N + (c + dc)])
            )
            if bn >= 3:
                junctions += 1
    return junctions


# ---------------------------------------------------------------------------
# Scoring
# ---------------------------------------------------------------------------
def score_bridging(bridged, path_len):
    if bridged:
        return 30, f"path found, length {path_len}"
    return 0, "no bridge"


def score_narrowness(frac):
    pct = frac * 100
    if frac < 0.02:
        return 25, f"{pct:.1f}% bright"
    if frac < 0.05:
        return 20, f"{pct:.1f}% bright"
    if frac < 0.10:
        return 12, f"{pct:.1f}% bright"
    if frac < 0.20:
        return 5,  f"{pct:.1f}% bright"
    return 0, f"{pct:.1f}% bright"


def score_aspect(ar):
    if ar > 5:
        return 20, f"ratio {ar:.1f}"
    if ar > 3:
        return 15, f"ratio {ar:.1f}"
    if ar > 2:
        return 8,  f"ratio {ar:.1f}"
    return 0, f"ratio {ar:.1f}"


def score_tortuosity(bridged, tort):
    if not bridged:
        return 0, "n/a (no bridge)"
    if tort > 1.5:
        return 15, f"ratio {tort:.2f}"
    if tort > 1.2:
        return 10, f"ratio {tort:.2f}"
    if tort > 1.05:
        return 5,  f"ratio {tort:.2f}"
    return 0, f"ratio {tort:.2f}"


def score_branches(bridged, n_junctions, frac):
    # A blob floods the entire grid with junctions — only meaningful when filament is narrow.
    # Require bridging AND < 10% bright cells before awarding branch points.
    if not bridged or frac >= 0.10:
        return 0, f"{n_junctions} junctions (too wide to count)"
    if n_junctions >= 3:
        return 10, f"{n_junctions} junctions"
    if n_junctions >= 2:
        return 6,  f"{n_junctions} junctions"
    if n_junctions >= 1:
        return 3,  f"{n_junctions} junction"
    return 2, "0 junctions (straight bridge)"


def letter_grade(total):
    if total >= 90: return "A"
    if total >= 75: return "B"
    if total >= 55: return "C"
    if total >= 35: return "D"
    return "F"


VERDICTS = {
    "A": "Excellent — narrow, tortuous, branching lightning filament",
    "B": "Good — clear filament with measurable branching",
    "C": "Moderate — filament visible but broad or poorly branched",
    "D": "Poor — weak filament structure, mostly diffuse",
    "F": "Fail — no lightning shape (stripe, blob, or no bridge)",
}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    if len(sys.argv) < 2:
        print("Usage: python3 grade_filament.py <sigma_final.bin> [N]")
        sys.exit(1)

    sigma_file = sys.argv[1]
    N_arg = int(sys.argv[2]) if len(sys.argv) > 2 else None

    sigma, N = load_sigma(sigma_file, N_arg)

    bridged, path_len, bridge_path = bfs_bridge(sigma, N)
    frac = compute_narrowness(sigma, N)
    ar   = compute_aspect_ratio(sigma, N)
    tort = compute_tortuosity(path_len, N) if bridged else 1.0
    njunc = count_junctions(sigma, N)

    s_bridge,    n_bridge    = score_bridging(bridged, path_len)
    s_narrow,    n_narrow    = score_narrowness(frac)
    s_aspect,    n_aspect    = score_aspect(ar)
    s_tort,      n_tort      = score_tortuosity(bridged, tort)
    s_branch,    n_branch    = score_branches(bridged, njunc, frac)

    total  = s_bridge + s_narrow + s_aspect + s_tort + s_branch
    grade  = letter_grade(total)
    verdict = VERDICTS[grade]

    basename = os.path.basename(sigma_file)
    W = 48
    print("CBRAM Filament Grade Report")
    print("=" * W)
    print(f"Simulation : {basename}  (N={N})")
    print(f"Bridging   : {s_bridge:2d}/30  ({n_bridge})")
    print(f"Narrowness : {s_narrow:2d}/25  ({n_narrow})")
    print(f"Aspect     : {s_aspect:2d}/20  ({n_aspect})")
    print(f"Tortuosity : {s_tort:2d}/15  ({n_tort})")
    print(f"Branches   : {s_branch:2d}/10  ({n_branch})")
    print("-" * W)
    print(f"Total      : {total}/100  → Grade: {grade}")
    print(f"Verdict    : {verdict}")
    print()
    # Machine-readable line for the skill to parse
    print(f"__GRADE_SUMMARY__ {total} {grade} {basename}")


if __name__ == "__main__":
    main()
