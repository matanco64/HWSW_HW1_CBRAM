#!/usr/bin/env bash
# run.sh — build, verify correctness, and profile all CBRAM stages.
# Usage: bash run.sh [N]
#   N  grid size (default 1024; tune so stage 0 takes 10-30 s)
#
# Must be invoked with `bash run.sh` (not ./run.sh) because the OneDrive
# rclone mount is noexec — shebangs cannot be executed directly from it.

set -euo pipefail

SRCDIR="$(cd "$(dirname "$0")" && pwd)"
BUILD_DIR="${BUILD_DIR:-/tmp/cbram_build}"
N="${1:-1024}"
RESULTS_DIR="${BUILD_DIR}/results"

mkdir -p "${RESULTS_DIR}"

# Helper: check if a stage binary exists
has_stage() { [[ -f "${BUILD_DIR}/cbram_stage${1}" ]]; }

# ── 1. Build ──────────────────────────────────────────────────────────────────
echo "=== Building all stages ==="
bash "${SRCDIR}/build.sh"
echo ""

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
    echo "=== Running stage ${stage} (${label}, N=${N}) ==="
    mkdir -p "${BUILD_DIR}/frames_stage${stage}"
    # Run from BUILD_DIR so output files land there
    (cd "${BUILD_DIR}" && \
        env ${extra_env} "${BUILD_DIR}/cbram_stage${stage}" "${N}")
    cp "${BUILD_DIR}/V_final_stage${stage}.bin" \
       "${RESULTS_DIR}/V_final_stage${stage}.bin"
    echo ""
}

run_stage 0 "naive AoS"
run_stage 1 "SoA"
run_stage 2 "SoA + time skewing"
run_stage 3 "SoA + skewing + OpenMP" "OMP_PROC_BIND=close OMP_PLACES=cores"

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
echo "=== Profiling (perf stat -r 3) ==="

PERF_EVENTS="cycles,instructions,cache-references,cache-misses,\
L1-dcache-loads,L1-dcache-load-misses,LLC-loads,LLC-load-misses,\
dTLB-loads,dTLB-load-misses"

for stage in 0 1 2; do
    if ! has_stage "${stage}"; then continue; fi
    echo ""
    echo "--- Stage ${stage} ---"
    (cd "${BUILD_DIR}" && \
        perf stat -r 3 \
        -e "${PERF_EVENTS}" \
        "${BUILD_DIR}/cbram_stage${stage}" "${N}" \
        2> "${RESULTS_DIR}/stage${stage}.perf")
    cat "${RESULTS_DIR}/stage${stage}.perf"
done

if has_stage 3; then
    echo ""
    echo "--- Stage 3 (all physical cores) ---"
    (cd "${BUILD_DIR}" && \
        OMP_PROC_BIND=close OMP_PLACES=cores \
        perf stat -r 3 \
        -e "${PERF_EVENTS}" \
        "${BUILD_DIR}/cbram_stage3" "${N}" \
        2> "${RESULTS_DIR}/stage3.perf")
    cat "${RESULTS_DIR}/stage3.perf"
fi

# ── 6. Thread scaling sweep (stage 3) ────────────────────────────────────────
if has_stage 3; then
    PHYS_CORES=$(lscpu | awk '/^Core\(s\) per socket/ {cores=$NF}
                              /^Socket\(s\)/           {sockets=$NF}
                              END {print cores*sockets}')
    echo ""
    echo "=== Thread scaling sweep (stage 3, physical cores: ${PHYS_CORES}) ==="
    for T in 1 2 4 8; do
        [[ $T -gt $PHYS_CORES ]] && continue
        echo -n "  OMP_NUM_THREADS=${T} ... "
        (cd "${BUILD_DIR}" && \
            OMP_NUM_THREADS=${T} OMP_PROC_BIND=close OMP_PLACES=cores \
            perf stat -r 3 -e cycles,instructions,LLC-load-misses \
            "${BUILD_DIR}/cbram_stage3" "${N}" \
            2> "${RESULTS_DIR}/stage3_t${T}.perf")
        grep "seconds time elapsed" "${RESULTS_DIR}/stage3_t${T}.perf" | head -1
    done
fi

# ── 7. Hotspot confirmation (stage 0) ────────────────────────────────────────
echo ""
echo "=== Hotspot check (perf record on stage 0) ==="
(cd "${BUILD_DIR}" && \
    perf record -g \
    -o "${RESULTS_DIR}/stage0.data" \
    "${BUILD_DIR}/cbram_stage0" "${N}")
perf report -i "${RESULTS_DIR}/stage0.data" --stdio 2>/dev/null | head -30

# ── 8. Diff summary ───────────────────────────────────────────────────────────
echo ""
echo "=== Source diffs between stages ==="
for pair in "naive.cpp opt1_soa.cpp" "opt1_soa.cpp opt2_blocked.cpp" "opt2_blocked.cpp opt3_omp.cpp"; do
    a="${pair%% *}"
    b="${pair##* }"
    [[ -f "${SRCDIR}/${a}" && -f "${SRCDIR}/${b}" ]] || continue
    echo "--- ${a} → ${b} ---"
    diff -u --label "${a}" --label "${b}" \
        "${SRCDIR}/${a}" "${SRCDIR}/${b}" | head -60 || true
    echo ""
done

echo "=== All done. Results in ${RESULTS_DIR}/ ==="
