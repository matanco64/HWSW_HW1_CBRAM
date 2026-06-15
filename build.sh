#!/usr/bin/env bash
# build.sh — compile all CBRAM simulation stages

set -euo pipefail

SRCDIR="$(cd "$(dirname "$0")" && pwd)"
BUILD_DIR="${BUILD_DIR:-${SRCDIR}/build}"

CXX="${CXX:-g++}"
CXXFLAGS="-O2 -g -march=native -std=c++17 -Wall -Wno-unused-result -I${SRCDIR}"

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

build_stage dbm_stage0.cpp   cbram_stage0
build_stage dbm_stage1.cpp   cbram_stage1
build_stage dbm_stage2.cpp   cbram_stage2

echo ""
echo "All binaries in ${BUILD_DIR}/"
