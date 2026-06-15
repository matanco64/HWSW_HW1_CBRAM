HW1 — Code Profiling, Optimization, and HW/SW Understanding
CBRAM Filament-Growth Simulation (Dielectric Breakdown Model)
Submitters: Matan Cohen (<ID>), Yuval Kogan (<ID>)

================================================================================
MAPPING TO THE REQUIRED SUBMISSION COMPONENTS
================================================================================
1. PDF document (<=3 pages) ............ cbram_report.pdf
2. Unoptimized C++ source .............. dbm_stage0.cpp   (AoS baseline)
3. Optimized C++ source ................ dbm_stage1.cpp   (Struct-of-Arrays)
4. Shell script (compile/run/profile) .. build.sh, run.sh
5. PDF with names + IDs ................ names_and_ids.pdf

Also included:
  dbm_stage2.cpp ......... a further optimization we explored (time-skewing /
                           temporal blocking) that did NOT pay off — analysed in
                           the report as our "approach that didn't work".
  physics_dbm.h, io.h, io.cpp ... shared sources needed to compile the stages.
  CMakeLists.txt ......... alternative build (cmake -B build -S . && cmake --build build).
  prompts.md ............. log of the AI prompts used (required by the brief).
  comparison_0vs1.mp4 .... real-time race of stage 0 vs stage 1 (the speed-up, visualised).
  results/ ............... the perf evidence the report cites:
                             metrics.md ............ derived IPC / miss-rate table
                             env.txt ............... lscpu of the profiling machine
                             stage{0,1,2}.perf ..... raw `perf stat -r 3` output
                             stage{0,1,2}_flamegraph.svg ... `perf record` flame graphs

================================================================================
HOW TO BUILD AND RUN
================================================================================
  ./build.sh                     # compiles cbram_stage0/1/2 into ./build/
  ./build/cbram_stage0 6144 -s 80   # unoptimized run (N=6144, 80 growth steps)
  ./build/cbram_stage1 6144 -s 80   # optimized run
  cmp build/V_final_stage1.bin build/V_final_stage0.bin   # bit-identical check

Full pipeline (build + correctness gate + perf stat + flame graphs):
  ./run.sh 6144 80

Profiling machine for the numbers in the report:
  single-core Intel Xeon E5-2630 v3 @ 2.40 GHz (L1d 32 KiB, L2 4 MiB, L3 16 MiB).
