#!/usr/bin/env bash
# make_video.sh — convert PPM frames into MP4 videos for each stage.
# Usage: ./make_video.sh [framerate]
#   framerate  frames per second (default 30)

set -euo pipefail

SRCDIR="$(cd "$(dirname "$0")" && pwd)"
FPS="${1:-30}"

for stage_dir in "${SRCDIR}"/frames_stage*/; do
    stage=$(basename "${stage_dir}")          # e.g. frames_stage0
    out="${SRCDIR}/${stage}.mp4"

    frame_count=$(ls "${stage_dir}"frame_*.ppm 2>/dev/null | wc -l)
    if [[ ${frame_count} -eq 0 ]]; then
        echo "  Skipping ${stage} (no frames found)"
        continue
    fi

    echo "  Encoding ${stage} (${frame_count} frames @ ${FPS} fps) → ${stage}.mp4"
    ffmpeg -y \
        -framerate "${FPS}" \
        -pattern_type glob \
        -i "${stage_dir}frame_*.ppm" \
        -c:v libx264 \
        -pix_fmt yuv420p \
        -crf 18 \
        "${out}"
    echo "  Done: ${out}"
done
