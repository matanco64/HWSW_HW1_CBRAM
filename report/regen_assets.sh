#!/usr/bin/env bash
# Regenerate derived figure assets for report/cbram_report.typ.
#
#   report/assets/filament_0vs1.png   — Figure 1 (hero): a real-time race frame
#     from comparison_0vs1.mp4, grabbed at the instant SoA (bottom row) has bridged
#     while AoS (top row) is still ~halfway. That gap on screen IS the speed-up.
#
# The flame graphs (results/stage*_flamegraph.svg, Figure 2) are produced by the
# main profiling pipeline (./run.sh), not here.
#
# Usage:  report/regen_assets.sh [TS] [VIDEO]
#   TS     timestamp (seconds) to grab. Default 5.7 — for the N=1024 video rendered
#          with TIME_SCALE=0.067, this is just after the SoA row turns green/BRIDGED.
#          If you re-render the video (different N / TIME_SCALE), eyeball the moment
#          the bottom row reads "BRIDGED" while the top still shows a step count and
#          pass that time here.
#   VIDEO  defaults to comparison_0vs1.mp4 in the repo root.
set -euo pipefail

REPO="$(cd "$(dirname "$0")/.." && pwd)"
TS="${1:-5.7}"
VIDEO="${2:-${REPO}/comparison_0vs1.mp4}"
OUT="${REPO}/report/assets/filament_0vs1.png"

mkdir -p "$(dirname "$OUT")"
# -update 1 overwrites a single image; -q:v 2 is near-lossless JPEG-quality PNG.
ffmpeg -y -ss "$TS" -i "$VIDEO" -frames:v 1 -q:v 2 "$OUT"
echo "wrote $OUT  (from $VIDEO @ ${TS}s)"
