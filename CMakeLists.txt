cmake_minimum_required(VERSION 3.16)
project(CBRAM_Simulation CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Common flags for all stages — must be identical (except OpenMP on stage 3)
set(COMMON_FLAGS "-O2 -march=native -Wall -Wno-unused-result")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${COMMON_FLAGS}")

# TODO: uncomment when cbram_stage3 is added
# find_package(OpenMP REQUIRED)

# Shared I/O library (compiled once, linked into every stage)
add_library(cbram_io STATIC io.cpp)
target_include_directories(cbram_io PUBLIC ${CMAKE_SOURCE_DIR})

# ── Stage 0: AoS naive ────────────────────────────────────────────────────────
add_executable(cbram_stage0 naive.cpp)
target_link_libraries(cbram_stage0 PRIVATE cbram_io)

# ── Stage 1: Struct of Arrays ─────────────────────────────────────────────────
# TODO: uncomment when opt1_soa.cpp is written
# add_executable(cbram_stage1 opt1_soa.cpp)
# target_link_libraries(cbram_stage1 PRIVATE cbram_io)

# ── Stage 2: SoA + time skewing ───────────────────────────────────────────────
# TODO: uncomment when opt2_blocked.cpp is written
# add_executable(cbram_stage2 opt2_blocked.cpp)
# target_link_libraries(cbram_stage2 PRIVATE cbram_io)

# ── Stage 3: SoA + time skewing + OpenMP ─────────────────────────────────────
# TODO: uncomment when opt3_omp.cpp is written
# add_executable(cbram_stage3 opt3_omp.cpp)
# target_link_libraries(cbram_stage3 PRIVATE cbram_io OpenMP::OpenMP_CXX)

# ── Convenience target: build all stages ──────────────────────────────────────
# TODO: add cbram_stage1 cbram_stage2 cbram_stage3 once their sources exist
add_custom_target(all_stages
    DEPENDS cbram_stage0
    COMMENT "Building all available CBRAM simulation stages"
)

# ── Print build summary ───────────────────────────────────────────────────────
message(STATUS "")
message(STATUS "CBRAM Simulation — Build Configuration")
message(STATUS "  Compiler : ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "  Flags    : ${CMAKE_CXX_FLAGS}")
message(STATUS "  OpenMP   : ${OpenMP_CXX_VERSION}")
message(STATUS "")
message(STATUS "Run:  cmake -B build -S . && cmake --build build")
message(STATUS "")
