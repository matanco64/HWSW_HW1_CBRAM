#!/usr/bin/env bash
# make_submission.sh — assemble the slim HW1 submission zip for the course staff.
# Rebuilds the report PDFs, copies only the required + supporting files into
# submission/, and zips it. Excludes build artifacts, raw perf .data, the laptop
# results, .idea, old_cpp, and the planning docs.
set -euo pipefail

REPO="$(cd "$(dirname "$0")" && pwd)"
OUT="${REPO}/submission"
ZIP="${REPO}/cbram_hw1_submission.zip"

# ── 1. (Re)build the report + cover PDFs ──────────────────────────────────────
# Real IDs are injected at compile time from the gitignored report/ids.local so
# they never live in the repo. Without that file the PDFs show "<ID>".
IDS="${REPO}/report/ids.local"
if [[ -f "${IDS}" ]]; then
  # shellcheck source=/dev/null
  source "${IDS}"
else
  echo "WARNING: ${IDS} not found — compiling with <ID> placeholders." >&2
fi
ID_INPUTS=(--input "matan-id=${MATAN_ID:-<ID>}" --input "yuval-id=${YUVAL_ID:-<ID>}")

typst compile "${ID_INPUTS[@]}" "${REPO}/report/cbram_report.typ"   "${REPO}/report/cbram_report.pdf"
typst compile "${ID_INPUTS[@]}" "${REPO}/report/names_and_ids.typ"  "${REPO}/report/names_and_ids.pdf"

# ── 2. Stage the files ────────────────────────────────────────────────────────
rm -rf "${OUT}" "${ZIP}"
mkdir -p "${OUT}/results"

# Component 1: the report (≤3 pp).            Component 5: names + IDs PDF.
cp "${REPO}/report/cbram_report.pdf"  "${OUT}/"
cp "${REPO}/report/names_and_ids.pdf" "${OUT}/"

# Components 2 & 3: unoptimized + optimized source (+ stage2 = the "didn't pay" one).
cp "${REPO}"/dbm_stage0.cpp "${REPO}"/dbm_stage1.cpp "${REPO}"/dbm_stage2.cpp "${OUT}/"
# Shared sources needed to compile the stages.
cp "${REPO}"/physics_dbm.h "${REPO}"/io.h "${REPO}"/io.cpp "${OUT}/"
# flush_cache (cold-cache helper built by build.sh); perf_metrics.py (run.sh step 9).
cp "${REPO}"/flush_cache.cpp "${REPO}"/perf_metrics.py "${OUT}/"

# Component 4: compile / run / profile scripts.
cp "${REPO}"/build.sh "${REPO}"/run.sh "${REPO}"/CMakeLists.txt "${OUT}/"

# Required by the brief: the AI prompt log.
cp "${REPO}/MDs/prompts.md" "${OUT}/prompts.md"

# Supporting evidence the report cites (small text + svg only — NOT the raw .data).
cp "${REPO}/results/metrics.md" "${REPO}/results/env.txt" "${OUT}/results/"
cp "${REPO}"/results/stage{0,1,2}.perf            "${OUT}/results/"
cp "${REPO}"/results/stage{0,1,2}_flamegraph.svg  "${OUT}/results/"
# Folded stacks: tiny text that regenerates the flame graphs with any options.
cp "${REPO}"/results/stage{0,1,2}.folded          "${OUT}/results/"

# Visual: the real-time race video (small).
cp "${REPO}/comparison_0vs1.mp4" "${OUT}/"

cp "${REPO}/SUBMISSION_README.txt" "${OUT}/README.txt"

# ── 3. Zip ────────────────────────────────────────────────────────────────────
( cd "${OUT}" && zip -r -q "${ZIP}" . )
echo "Built ${ZIP}"
( cd "${OUT}" && find . -type f | sort )
echo "--- size ---"; du -h "${ZIP}" | cut -f1
