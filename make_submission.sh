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
typst compile "${REPO}/report/cbram_report.typ"   "${REPO}/report/cbram_report.pdf"
typst compile "${REPO}/report/names_and_ids.typ"  "${REPO}/report/names_and_ids.pdf"

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

# Component 4: compile / run / profile scripts.
cp "${REPO}"/build.sh "${REPO}"/run.sh "${REPO}"/CMakeLists.txt "${OUT}/"

# Required by the brief: the AI prompt log.
cp "${REPO}/MDs/prompts.md" "${OUT}/prompts.md"

# Supporting evidence the report cites (small text + svg only — NOT the raw .data).
cp "${REPO}/results/metrics.md" "${REPO}/results/env.txt" "${OUT}/results/"
cp "${REPO}"/results/stage{0,1,2}.perf            "${OUT}/results/"
cp "${REPO}"/results/stage{0,1,2}_flamegraph.svg  "${OUT}/results/"

# Visual: the real-time race video (small).
cp "${REPO}/comparison_0vs1.mp4" "${OUT}/"

cp "${REPO}/SUBMISSION_README.txt" "${OUT}/README.txt"

# ── 3. Zip ────────────────────────────────────────────────────────────────────
( cd "${OUT}" && zip -r -q "${ZIP}" . )
echo "Built ${ZIP}"
( cd "${OUT}" && find . -type f | sort )
echo "--- size ---"; du -h "${ZIP}" | cut -f1
