#!/usr/bin/env bash
# make_comparison.sh — frames run for any number of stages, then render a
# time-synchronized comparison video where faster stages bridge first on screen.
#
# Usage: ./make_comparison.sh [N] [frame_interval] [time_scale] [stage ...]
#   N              grid size (default 200)
#   frame_interval dump a frame every K growth steps (default 10)
#   time_scale     slow-down factor vs real time (default 3.0)
#   stage ...      stages to compare (default: 0 1)
#
# Examples:
#   bash make_comparison.sh                      # stage 0 vs 1, N=200
#   bash make_comparison.sh 200 10 3.0 0 1 2    # all three stages, N=200

set -euo pipefail

SRCDIR="$(cd "$(dirname "$0")" && pwd)"

N="${1:-200}";          shift || true
FRAME_INT="${1:-10}";   shift || true
TIME_SCALE="${1:-3.0}"; shift || true
STAGES=("${@:-0 1}")    # remaining args are stage numbers; default 0 1
if [[ ${#STAGES[@]} -eq 0 ]]; then STAGES=(0 1); fi

STAGES_CSV=$(IFS=,; echo "${STAGES[*]}")
OUT="${SRCDIR}/comparison_${STAGES_CSV//,/vs}.mp4"

# ── 1. Build once ─────────────────────────────────────────────────────────────
echo "=== Building ==="
bash "${SRCDIR}/build.sh"
echo ""

# ── 2. Frames run for each stage ──────────────────────────────────────────────
for stage in "${STAGES[@]}"; do
    echo "=== Stage ${stage}: frames run ==="
    bash "${SRCDIR}/make_video.sh" --no-build "${stage}" "${N}" 30 "${FRAME_INT}"
    echo ""
done

# ── 3. Render comparison ──────────────────────────────────────────────────────
echo "=== Rendering comparison video ==="
python3 "${SRCDIR}/render_comparison_video.py" \
    "${N}" "${STAGES_CSV}" "${OUT}" "${TIME_SCALE}"

echo ""
echo "=== Done: ${OUT} ==="
