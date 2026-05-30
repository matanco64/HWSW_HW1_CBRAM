#!/usr/bin/env bash
# build.sh — compile all CBRAM simulation stages
# Outputs go to BUILD_DIR (default /tmp/cbram_build) because the OneDrive
# rclone mount is noexec; source files stay in the repo.

set -euo pipefail

SRCDIR="$(cd "$(dirname "$0")" && pwd)"
BUILD_DIR="${BUILD_DIR:-/tmp/cbram_build}"

CXX="${CXX:-g++}"
CXXFLAGS="-O2 -march=native -std=c++17 -Wall -Wno-unused-result -I${SRCDIR}"

mkdir -p "${BUILD_DIR}"

echo "=== CBRAM build ==="
echo "  SRCDIR    = ${SRCDIR}"
echo "  BUILD_DIR = ${BUILD_DIR}"
echo "  CXX       = ${CXX}"
echo "  CXXFLAGS  = ${CXXFLAGS}"
echo ""

build_stage() {
    local src="$1"
    local out="$2"
    local extra="${3:-}"
    if [[ ! -f "${SRCDIR}/${src}" ]]; then
        echo "  Skipping ${out} (${src} not yet written)"
        return
    fi
    echo -n "  Building ${out} ... "
    # shellcheck disable=SC2086
    ${CXX} ${CXXFLAGS} ${extra} "${SRCDIR}/${src}" "${SRCDIR}/io.cpp" \
        -o "${BUILD_DIR}/${out}"
    echo "OK"
}

build_stage naive.cpp        cbram_stage0
build_stage opt1_soa.cpp     cbram_stage1
build_stage opt2_blocked.cpp cbram_stage2
build_stage opt3_omp.cpp     cbram_stage3 "-fopenmp"

echo ""
echo "All binaries in ${BUILD_DIR}/"
