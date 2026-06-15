#!/usr/bin/env bash
# Regenerate derived figure assets for report/cbram_report.typ.
#
#   report/assets/filament_0vs1.png   — Figure 1 (hero): the final bridged frame
#                                        extracted from comparison_0vs1.mp4.
#
# The flame graphs (results/stage*_flamegraph.svg, Figure 2) are produced by the
# main profiling pipeline (./run.sh), not here.
#
# Usage:  report/regen_assets.sh [VIDEO]
#   VIDEO defaults to comparison_0vs1.mp4 in the repo root.
set -euo pipefail

REPO="$(cd "$(dirname "$0")/.." && pwd)"
VIDEO="${1:-${REPO}/comparison_0vs1.mp4}"
OUT="${REPO}/report/assets/filament_0vs1.png"

mkdir -p "$(dirname "$OUT")"
# -sseof -0.2 seeks ~0.2 s before the end → the final, fully-bridged frame.
# -update 1 overwrites a single image; -q:v 2 is near-lossless JPEG-quality PNG.
ffmpeg -y -sseof -0.2 -i "$VIDEO" -update 1 -q:v 2 "$OUT"
echo "wrote $OUT"
