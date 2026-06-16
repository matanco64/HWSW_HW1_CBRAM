#!/usr/bin/env bash
# run.sh — build, verify correctness, and profile all CBRAM stages.
# Usage: ./run.sh [N] [STEPS] [--no-build|-B] [--no-run|-R] [--topdown]
#                  [--perf-runs N] [--stages 0,1,2,3] [--steps K]
#   N      grid size (default 6144 — working set ≈ 300 MB, far past the 24 MB
#          L3, so the memory hierarchy is genuinely exercised)
#   STEPS  cap on growth steps (default 80; 0 = run until bridged).
#          A full bridge at large N takes thousands of steps; a fixed cap
#          keeps every run at an identical, affordable amount of work.

set -euo pipefail

SRCDIR="$(cd "$(dirname "$0")" && pwd)"
BUILD_DIR="${BUILD_DIR:-${SRCDIR}/build}"
RESULTS_DIR="${SRCDIR}/results"

N=6144
STEPS=80
NO_BUILD=0
NO_RUN=0
TOPDOWN=0
PERF_RUNS=3
STAGES="0 1 2 3"

NPOS=0   # bare numbers: first is N, second is STEPS
while [[ $# -gt 0 ]]; do
    case "$1" in
        --no-build|-B)  NO_BUILD=1 ;;
        --no-run|-R)    NO_RUN=1 ;;
        --topdown)      TOPDOWN=1 ;;
        --perf-runs)    PERF_RUNS="$2"; shift ;;
        --stages)       STAGES="${2//,/ }"; shift ;;
        --steps)        STEPS="$2"; shift ;;
        [0-9]*)         if [[ $NPOS -eq 0 ]]; then N="$1"; else STEPS="$1"; fi
                        NPOS=$((NPOS + 1)) ;;
        *)              echo "Unknown flag: $1" >&2; exit 1 ;;
    esac
    shift
done

mkdir -p "${RESULTS_DIR}"

# Helper: check if a stage binary exists
has_stage() { [[ -f "${BUILD_DIR}/cbram_stage${1}" ]]; }

# ── 1. Build ──────────────────────────────────────────────────────────────────
if [[ $NO_BUILD -eq 0 ]]; then
    echo "=== Building all stages ==="
    bash "${SRCDIR}/build.sh"
    echo ""
fi

# ── 2. Log hardware environment ───────────────────────────────────────────────
echo "=== Hardware environment ==="
lscpu | tee "${RESULTS_DIR}/env.txt"
echo ""

# ── 3. Run available stages ───────────────────────────────────────────────────
run_stage() {
    local stage=$1
    local label=$2
    local extra_env="${3:-}"
    if ! has_stage "${stage}"; then
        echo "  Skipping stage ${stage} (not built yet)"
        return
    fi
    echo "=== Running stage ${stage} (${label}, N=${N}, steps=${STEPS}) ==="
    mkdir -p "${BUILD_DIR}/frames_stage${stage}"
    # Flush caches before timing so the run starts from a cold-cache state.
    [[ -x "${BUILD_DIR}/flush_cache" ]] && "${BUILD_DIR}/flush_cache" > /dev/null
    # Run from BUILD_DIR so output files land there. Frame dumps are opt-in
    # (-f) and deliberately NOT enabled here; `time` logs wall/user per run.
    (cd "${BUILD_DIR}" && \
        time env ${extra_env} "${BUILD_DIR}/cbram_stage${stage}" "${N}" -s "${STEPS}")
    cp "${BUILD_DIR}/V_final_stage${stage}.bin" \
       "${RESULTS_DIR}/V_final_stage${stage}.bin"
    echo ""
}

if [[ $NO_RUN -eq 0 ]]; then
    run_stage 0 "naive AoS"
    run_stage 1 "SoA"
    run_stage 2 "SoA + time skewing"
    run_stage 3 "SoA + skewing + SIMD"
fi

# ── 4. Correctness gate ───────────────────────────────────────────────────────
echo "=== Correctness verification (cmp against stage 0 reference) ==="
REF="${RESULTS_DIR}/V_final_stage0.bin"
all_ok=1
for stage in 1 2 3; do
    if ! has_stage "${stage}"; then continue; fi
    FILE="${RESULTS_DIR}/V_final_stage${stage}.bin"
    if cmp -s "${REF}" "${FILE}"; then
        echo "  Stage ${stage}: BIT-IDENTICAL OK"
    else
        echo "  Stage ${stage}: MISMATCH — investigate before trusting perf numbers!"
        all_ok=0
    fi
done
if [[ $all_ok -eq 0 ]]; then
    echo ""
    echo "ERROR: correctness gate failed. Fix before profiling."
    exit 1
fi
echo ""

# ── 5. Perf profiling ─────────────────────────────────────────────────────────
echo "=== Profiling (perf stat -r ${PERF_RUNS}) ==="

PERF_EVENTS="cycles,instructions,cache-references,cache-misses,\
L1-dcache-loads,L1-dcache-load-misses,LLC-loads,LLC-load-misses,\
dTLB-loads,dTLB-load-misses"

for stage in ${STAGES}; do
    if ! has_stage "${stage}"; then continue; fi
    echo ""
    echo "--- Stage ${stage} ---"
    [[ -x "${BUILD_DIR}/flush_cache" ]] && "${BUILD_DIR}/flush_cache" > /dev/null
    (cd "${BUILD_DIR}" && \
        perf stat -r "${PERF_RUNS}" \
        -e "${PERF_EVENTS}" \
        "${BUILD_DIR}/cbram_stage${stage}" "${N}" -n -s "${STEPS}" \
        2> "${RESULTS_DIR}/stage${stage}.perf")
    cat "${RESULTS_DIR}/stage${stage}.perf"
done

# ── 5b. Topdown microarchitecture analysis ───────────────────────────────────
# Where do the cycles actually go? Level-2 topdown splits Backend Bound into
# Memory Bound vs Core Bound — the Memory-Bound % is the number that should
# shrink after the SoA layout change, directly demonstrating the HW/SW insight.
topdown_one() {
    local stage=$1
    has_stage "${stage}" || return 0
    local out="${RESULTS_DIR}/stage${stage}.topdown"
    echo ""
    echo "--- Stage ${stage} ---"
    if (cd "${BUILD_DIR}" && perf stat --topdown --td-level 2 \
            -- "${BUILD_DIR}/cbram_stage${stage}" "${N}" -n -s "${STEPS}" >/dev/null) 2> "${out}"; then
        :
    elif (cd "${BUILD_DIR}" && perf stat --topdown \
            -- "${BUILD_DIR}/cbram_stage${stage}" "${N}" -n -s "${STEPS}" >/dev/null) 2> "${out}"; then
        :
    else
        echo "  (topdown not supported by this perf/CPU — skipping)" | tee "${out}"
        return 0
    fi
    cat "${out}"
}
if [[ $TOPDOWN -eq 1 ]]; then
    echo ""
    echo "=== Topdown analysis (frontend / backend / memory bound) ==="
    for stage in ${STAGES}; do
        topdown_one "${stage}"
    done
fi

# ── 7. Hotspot confirmation & Flamegraph (every built stage) ─────────────────
echo ""
echo "=== Hotspot check & Flamegraphs (perf record, per stage) ==="

# Ensure FlameGraph tools are available
FLAMEGRAPH_DIR="${SRCDIR}/FlameGraph"
if [[ ! -d "${FLAMEGRAPH_DIR}" ]]; then
    echo "Cloning FlameGraph repository..."
    git clone https://github.com/brendangregg/FlameGraph.git "${FLAMEGRAPH_DIR}"
fi

# Record one stage with user-space (:u) call-graphs and emit its flamegraph.
# -F 999 samples at 999 Hz; --call-graph dwarf gives correct stacks at -O2.
# -n keeps the profile I/O-free so only the compute hot path shows.
flamegraph_one() {
    local stage=$1
    has_stage "${stage}" || return 0
    local data="${RESULTS_DIR}/stage${stage}.data"
    local folded="${RESULTS_DIR}/stage${stage}.folded"
    local svg="${RESULTS_DIR}/stage${stage}_flamegraph.svg"
    echo ""
    echo "--- Stage ${stage} ---"
    [[ -x "${BUILD_DIR}/flush_cache" ]] && "${BUILD_DIR}/flush_cache" > /dev/null
    (cd "${BUILD_DIR}" && \
        perf record -e cpu-clock:u -F 999 --call-graph dwarf \
        -o "${data}" \
        -- "${BUILD_DIR}/cbram_stage${stage}" "${N}" -n -s "${STEPS}")
    # Keep the folded stacks: a tiny (KB) plain-text artifact that regenerates the
    # flame graph anywhere with any options, needing neither the .data nor the
    # binary. We omit --inline: at -O2 it synthesizes a phantom 'main' leaf frame
    # above solve_potential (an inline-attribution artifact, not a real call).
    perf script -i "${data}" | \
        "${FLAMEGRAPH_DIR}/stackcollapse-perf.pl" > "${folded}"
    "${FLAMEGRAPH_DIR}/flamegraph.pl" \
        --title "cbram_stage${stage} (N=${N})" "${folded}" > "${svg}"
    echo "  Flamegraph: ${svg}  (folded stacks: ${folded})"
    # Quick text hotspot summary (top self-time symbols).
    perf report -i "${data}" --stdio --no-children 2>/dev/null | \
        grep -E '^\s+[0-9]+\.[0-9]+%' | head -8
}
pids=(); logs=()
for stage in ${STAGES}; do
    log=$(mktemp); logs+=("${stage}:${log}")
    flamegraph_one "${stage}" > "${log}" 2>&1 &
    pids+=($!)
done
wait "${pids[@]}"
for entry in "${logs[@]}"; do cat "${entry#*:}"; rm -f "${entry#*:}"; done

# ── 8. Diff summary ───────────────────────────────────────────────────────────
echo ""
echo "=== Source diffs between stages ==="
for pair in "dbm_stage0.cpp dbm_stage1.cpp" "dbm_stage1.cpp dbm_stage2.cpp" "dbm_stage2.cpp dbm_stage3.cpp"; do
    a="${pair%% *}"
    b="${pair##* }"
    [[ -f "${SRCDIR}/${a}" && -f "${SRCDIR}/${b}" ]] || continue
    echo "--- ${a} → ${b} ---"
    diff -u --label "${a}" --label "${b}" \
        "${SRCDIR}/${a}" "${SRCDIR}/${b}" | head -60 || true
    echo ""
done

# ── 9. Derived metrics summary (IPC, miss rates, MPKI, speedup vs stage 0) ────
echo ""
echo "=== Derived metrics summary ==="
python3 "${SRCDIR}/perf_metrics.py" "${RESULTS_DIR}" || \
    echo "  (perf_metrics.py failed — check the stage*.perf files)"
echo ""

echo "=== All done. Results in ${RESULTS_DIR}/ ==="
