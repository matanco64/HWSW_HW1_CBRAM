cmake_minimum_required(VERSION 3.16)
project(CBRAM_Simulation CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Common flags for all stages — must be identical across stages
set(COMMON_FLAGS "-O2 -g -march=native -Wall -Wno-unused-result")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${COMMON_FLAGS}")

# Shared I/O library (compiled once, linked into every stage)
add_library(cbram_io STATIC io.cpp)
target_include_directories(cbram_io PUBLIC ${CMAKE_SOURCE_DIR})

# ── Stage 0: DBM AoS naive ────────────────────────────────────────────────────
add_executable(cbram_stage0 dbm_stage0.cpp)
target_link_libraries(cbram_stage0 PRIVATE cbram_io)

# ── Stage 1: Struct of Arrays ─────────────────────────────────────────────────
add_executable(cbram_stage1 dbm_stage1.cpp)
target_link_libraries(cbram_stage1 PRIVATE cbram_io)

# ── Stage 2: SoA + time skewing ───────────────────────────────────────────────
add_executable(cbram_stage2 dbm_stage2.cpp)
target_link_libraries(cbram_stage2 PRIVATE cbram_io)

# ── Convenience target: build all stages ──────────────────────────────────────
add_custom_target(all_stages
    DEPENDS cbram_stage0 cbram_stage1 cbram_stage2
    COMMENT "Building all available CBRAM simulation stages"
)

# ── Print build summary ───────────────────────────────────────────────────────
message(STATUS "")
message(STATUS "CBRAM Simulation — Build Configuration")
message(STATUS "  Compiler : ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "  Flags    : ${CMAKE_CXX_FLAGS}")
message(STATUS "")
message(STATUS "Run:  cmake -B build -S . && cmake --build build")
message(STATUS "")
