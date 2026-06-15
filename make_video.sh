#!/usr/bin/env bash
# make_video.sh — build, frames run, and render an MP4 for one CBRAM stage.
# Usage: ./make_video.sh [stage] [N] [fps] [frame_interval]
#   stage           which stage to run (default 0)
#   N               grid size (default 200; smaller = faster run + tiny grid)
#   fps             output video frame rate (default 30)
#   frame_interval  dump a frame every K growth steps (default 10)
#
# The binary runs until bridged (-s 0), dumping frames into build/frames_stageN/,
# then render_dbm_video.py renders them to filament_stageN.mp4.

set -euo pipefail

SRCDIR="$(cd "$(dirname "$0")" && pwd)"
BUILD_DIR="${SRCDIR}/build"

STAGE="${1:-0}"
N="${2:-200}"
FPS="${3:-30}"
FRAME_INT="${4:-10}"

BINARY="${BUILD_DIR}/cbram_stage${STAGE}"
FRAMES_DIR="${BUILD_DIR}/frames_stage${STAGE}"
OUT="${SRCDIR}/filament_stage${STAGE}.mp4"

# ── 1. Build ──────────────────────────────────────────────────────────────────
echo "=== Building ==="
bash "${SRCDIR}/build.sh"
echo ""

if [[ ! -f "${BINARY}" ]]; then
    echo "ERROR: ${BINARY} not found after build (stage ${STAGE} may not exist yet)."
    exit 1
fi

# ── 2. Frames run ─────────────────────────────────────────────────────────────
echo "=== Frames run: stage ${STAGE}, N=${N}, frame every ${FRAME_INT} steps ==="
mkdir -p "${FRAMES_DIR}"
rm -f "${FRAMES_DIR}"/frame_*.bin   # clear stale frames

(cd "${BUILD_DIR}" && \
    "${BINARY}" "${N}" -f -v -F "${FRAME_INT}")

frame_count=$(ls "${FRAMES_DIR}"/frame_*_V.bin 2>/dev/null | wc -l)
echo ""
echo "  Captured ${frame_count} frames in ${FRAMES_DIR}"

if [[ ${frame_count} -eq 0 ]]; then
    echo "ERROR: no frames were written. Check that the binary supports -f."
    exit 1
fi

# ── 3. Render ─────────────────────────────────────────────────────────────────
echo ""
echo "=== Rendering video ==="
python3 "${SRCDIR}/render_dbm_video.py" "${N}" "${FRAMES_DIR}" "${OUT}"

echo ""
echo "=== Done: ${OUT} ==="
