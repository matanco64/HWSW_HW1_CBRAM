# prompts.md — AI Prompt Log

Log of prompts used throughout this project, as required by the assignment.

---

## Planning Phase

### Prompt 1 — Initial plan review and merge
> "Hello, Look at the 2 Plans for creating a CBRAM simulation for a HW perf task, both plans sit at the current directory CBRAM_PLAN{,_V2}.md the V2 Have improvements over the first version, Let me know what you think about them, HW1.pdf is the task, Let's merge them into a singular clear plan with steps, one change that I don't want from V2 is that we should create a different file for each optimization, but let's make the diff clear in some way."

**Tool:** Claude Code (claude-sonnet-4-6)
**Output:** `CBRAM_PLAN_FINAL.md` — merged plan combining V1's detailed physics spec with V2's integer arithmetic, 3-stage optimization layering (SoA → time skewing → OpenMP), and ablation study methodology. Separate source files per stage retained from V1 preference; per-file diff header blocks added to make changes explicit.

---

## Implementation Phase

### Prompt 2 — Root-cause analysis of naive.cpp (Yuval)
> "Look at the naive.cpp code (and its supporting files) and tell me what's wrong with it. I want to understand what is producing those 3 vertical lines instead of a lightning. give me a deep analysis of the physics. Use the cbram_memristor_guide.md as a reference."

**Tool:** Claude Code (claude-sonnet-4-6)
**Output:** `cbram_analysis.md` — full root-cause analysis identifying 7 physics defects: no stochastic nucleation (primary), no lateral diffusion (primary), linear flux instead of sinh, sigma-everywhere growth, seeds at anode not cathode, no Butler-Volmer injection, and unused Joule heating. Defined 7-step cumulative fix plan with code snippets.

---

### Prompt 3 — Grade naive baseline (Yuval)
> "Let's add a grade script, start from step 0 and see where we are."

**Tool:** Claude Code (claude-sonnet-4-6)
**Output:** `grade_filament.py` — 5-metric grader (bridging 30 pts, narrowness 25 pts, aspect ratio 20 pts, tortuosity 15 pts, branches 10 pts). Step 0 grade: **0/100 F** (34.7% bright, no bridge, aspect 1.0).

---

### Prompt 4 — Implement steps 1–5 (Yuval)
> (Series of prompts to implement each fix step cumulatively in improved.cpp)

**Tool:** Claude Code (claude-sonnet-4-6)
**Output:** `improved.cpp` + `physics_improved.h` updated through step 5 (sinh field-enhanced hopping). `cbram_progress.md` tracking results per step.

| Step | Change | Score |
|---|---|---|
| 1 | Lateral diffusion | 5/100 F |
| 2 | Stochastic deposition | 20/100 F |
| 3 | Cathode seeds + sigma-decoupled drift | 25/100 F |
| 4 | Tip-only deposition (+ ION_HIGH prefill) | 30/100 F |
| 5 | Sinh field-enhanced hopping | 30/100 F |

---

### Prompt 5 — Debug grader hang, implement steps 6 and 7 (Yuval)
> "I think we have a problem when running: 'python3 grade_filament.py sigma_final_improved.bin' - claude stucks around this point. now implement last fixed 6 and 7 and run the analysis - let's see the lightning working."

**Tool:** Claude Code (claude-sonnet-4-6)
**Diagnosis:** Grader itself runs fine (<1 s for N=200). The "stuck" issue was a previous session running the simulation with DEFAULT_N=4096 (takes hours). Current sigma_final_improved.bin is correctly N=200.

**Step 6 — Butler-Volmer injection:** Added `q_exp_approx()` to `physics_improved.h`. Replaced constant anode injection with `ION_INJECT × exp(BV_ALPHA × η)`. Initial version using `ION_INJECT × sigma × exp(...)` broke the bridge (sigma at row 1 is SIGMA_LOW since seeds are at cathode after step 3 — catch-22). Fixed to voltage-only: `inj = q_mul(ION_INJECT, bv)`.

**Step 7 — Joule heating:** Added `update_temperature()` after Jacobi sweeps: computes P = σ·|∇V|² per cell, sets temp = TEMP_INIT + P·R_TH. In drift, replaced fixed ION_MOBILITY with Arrhenius-scaled `mob_t = ION_MOBILITY × exp(EA_OVER_K × ΔT)`.

**Result:** Grade **30/100 F** — bridge restored, but blob persists (90.9% bright, aspect 1.0). See `summary.md` for diagnosis and next steps.

---

### Prompt 6 — File hygiene, prompts log, commit prep, and summary (Yuval)
> "OK, now we are going to do the following: 1. make sure all files are up-do-date. 2. update @HWSW_HW1_CBRAM/prompts.md with all the prompts I gave you today. make sure you mark it was given by Yuval. 3. make sure we are ready for commiting, we need to ignore artifacts and only keep real files - DON'T COMMIT YOURSELF. 4. create a summary.md file that will include: what we did, where we stand, future work to make it bocome lightning (based on @HWSW_HW1_CBRAM/cbram_memristor_guide.md) and a ramp-up prompt for the next session."

**Tool:** Claude Code (claude-sonnet-4-6)
**Output:** This prompts.md update; `.gitignore` verified (frames*, *.bin, build/ already excluded); `summary.md` created.

---

## Session 2 — Fix re-implementation after context compaction

### Prompt 7 — Apply Fix 1: lateral diffusion (Yuval, via /compact instruction)
> "then continue with Option A and apply fix 1"

**Tool:** Claude Code (claude-sonnet-4-6)
**Output:** Added explicit-Euler lateral diffusion pass in `drift_diffusion` using `D_LAT = ONE/20`. Narrowness improved (34.7% → 15.8%) but aspect ratio regressed (1.0 → 0.2) — diffusion homogenised ions into a horizontal band at the anode. Grade: **5/100 F**.

---

### Prompt 8 — Apply Fix 2: stochastic deposition (Yuval)
> "proceed to Fix 2"

**Tool:** Claude Code (claude-sonnet-4-6)
**Output:** Replaced deterministic `σ += SIGMA_GROWTH × ion` with per-cell xorshift32 probabilistic deposition (threshold = `SIGMA_GROWTH × ion` mapped to `[0, 2³²)`). Narrowness jumped to 2.3% (20/25 pts). Aspect still 0.0 — bright cells all near anode seeds. Grade: **20/100 F**.

---

### Prompt 9 — Apply Fix 3: cathode nucleation (Yuval)
> "proceed to Fix 3"

**Tool:** Claude Code (claude-sonnet-4-6)
**Output:** Three combined changes: (1) seeds relocated from row 1 → row N-2; (2) anode injection made constant (not σ-dependent, since σ at row 1 is now SIGMA_LOW); (3) ion drift decoupled from σ — `flux = ION_MOBILITY × E × ion` instead of `× σ × E × ion` (ionic and electronic conductivities are independent). Perfect narrowness 1.8% (25/25 pts). Aspect still 0.0 — ions don't reach the cathode seeds within 3000 steps without tip-only growth. Grade: **25/100 F**.

---

### Prompt 10 — Apply Fix 4: tip-only deposition (Yuval)
> "proceed to Fix 4"

**Tool:** Claude Code (claude-sonnet-4-6)
**Output:** Four combined changes: (1) cathode seeds initialised at `SIGMA_MAX` (fully metallic so tip-only fires immediately); (2) electrolyte pre-filled with `ion = ION_HIGH` (physically: dissolved metal in equilibrium before SET); (3) `update_sigma` gated on adjacency to metallic cell (`σ > SIGMA_HIGH`); (4) `SIGMA_GROWTH_IMP = ONE` (10× baseline, ~10 deposits to metallic per cell → bridge in ~2000 steps); (5) early-exit check moved from row N-2 to row 1 (SET complete when filament reaches anode). Filament bridged at **t=2033**. However blob persists (91.2% bright, aspect 1.0) — tip-only gate failed to suppress the blob because the entire electrolyte becomes adjacent-to-metal once the filament spans the grid. Grade: **30/100 F** (bridging 30/30, all other metrics 0).

---

### Prompt 11 — Status check and command-line help (Yuval)
> "what is the status? it seemed like you stuck" / "how do you run the simulation? what is the command line?"

**Tool:** Claude Code (claude-sonnet-4-6)
**Output:** Confirmed grader runs fine (<1 s for N=200). Session appeared stuck due to large verbose output being paginated. Commands: `./build/cbram_improved 200 -v` to run; `python3 grade_filament.py sigma_final_improved.bin` to grade. Current best: 30/100 F.

---

### Prompt 12 — Update prompts.md (Yuval)
> "read @HWSW_HW1_CBRAM/prompts.md and update with what we added"

**Tool:** Claude Code (claude-sonnet-4-6)
**Output:** This update.

---

## Session 3 — Python ground truth and physics model decision (Matan)

### Prompt 13 — Python implementation analysis and C++ port decision (Matan)
> "look at the folder python_imp, there is a python implemenation created by Yuval, Unlike our CPP implementation, it works and outputs a good result, For the project we would like to move it to C/C++, what would you prefer on doing (C imp or C++ imp) also, please tell me about the implementation is it still cache optmizable? and again we should start with a naive implemenation"

**Tool:** Claude Code (claude-sonnet-4-6)
**Output:** Analyzed `python_imp/cbram_sim.py` — identified it uses biased random-walk tip-growth (DLA-style), not continuum drift-diffusion. Jacobi solver is visualization-only; tip probabilities are hardcoded constants (P_UP=0.60 etc.), V never passed to `grow_step()`. Recommended C++ to reuse existing infrastructure. Confirmed Jacobi kernel is still the cache-optimization target.

---

### Prompt 14 — Questioning the role of Jacobi (Matan)
> "If it is just a random walk why the Jacobi kernel is important for the simulation? is it not just a random walk sim?"

**Tool:** Claude Code (claude-sonnet-4-6)
**Output:** Confirmed V is not used inside `grow_step()` at all — Jacobi output is purely cosmetic. Concluded the pure random walk has no meaningful bottleneck to optimize, making it unsuitable for the HW optimization story.

---

### Prompt 15 — Physics model selection (Matan)
> "I want you to read the @HW1.pdf I think doing B will be a better option for doing the task at hand (a computationally intensive sim), but please suggest more options, I think the original Continuum Nernst-Planck drift was a bit hard (it didnt work well for the branching) but maybe there is a middle ground I am missing"

**Tool:** Claude Code (claude-sonnet-4-6)
**Output:** Read HW1.pdf. Presented 4 options: (A) pure random walk — Jacobi cosmetic, nothing to optimize; (B) continuum Nernst-Planck — correct bottleneck, poor branching; (C) physical transport + stochastic tip-only deposition — sweet spot, Jacobi load-bearing, branching from stochastic deposition; (D) interface DLA with Butler-Volmer probabilities. Recommended Option C.

---

### Prompt 16 — Python as ground truth (Matan)
> "I agree with C, but first I want to fix the Current python simulation, it is easier to iterate and it will be our ground truth, what do you think?"

**Tool:** Claude Code (claude-sonnet-4-6)
**Output:** Agreed. Clarified "fix" = replace pure random walk with Option C physics (real ion transport + stochastic tip-adjacent deposition). Noted Python float ≠ Q16.16, so it is a qualitative ground truth, not bit-exact.

---

### Prompt 17 — Plan Python rewrite, visual-first (Matan)
> "Let's start with the plan, currently DONT give the grader too much meaning, let me look visually first and then we will fix the grader together, but first, add the last couple of promts I gave you to promts.md"

**Tool:** Claude Code (claude-sonnet-4-6)
**Output:** Updated `prompts.md` with prompts 13–17. Next: plan the Python rewrite of `cbram_sim.py` with Option C physics, visual iteration before grader tuning.

---

## Session 4 — DBM pivot and naive C++ port (Matan)

### Prompt 18 — Critique the Tetris approach (Matan)
> "currently I am working on a Tetris inspired implementation of the sim cbram_tetris.py, I don't think it works well enough … what do you think I should change?"

**Tool:** Claude Code (claude-opus-4-8)
**Output:** Ran it (12/100 F) — diagnosed it as plain DLA producing a wide bush, widest at the feeding (anode) side, never bridging. Root cause: ions deposit on first contact with the broad top surface before reaching the high-field tips. Recommended the Dielectric Breakdown Model ingredient — sticking probability ∝ (local field)^η — as the bush→bolt knob.

### Prompt 19 — Implement field-gated tip-only sticking (Matan)
> "Yes"

**Tool:** Claude Code (claude-opus-4-8)
**Output:** Replaced the additive `field_factor` / neighbor-count sticking in `cbram_tetris.py` with a field-gated power-law rule (`P ∝ (|∇V|/E_REF)^ETA`, hard threshold below `E_STICK_MIN`). The fan collapsed (narrowness 0.1%) but growth stalled at 41 cells (25/100 F) — a documented dead-end. Concluded particle DLA is the wrong substrate.

### Prompt 20 — Switch to direct DBM, minimal and clean (Matan)
> "Maybe direct DBM is the correct approach here and let's drop the particles entirely, Let's maybe make it in a new file entirely, with as little as code as possible, and clean as possible"

**Tool:** Claude Code (claude-opus-4-8)
**Output:** Wrote `python_imp/cbram_dbm.py` (~95 lines): grow the metallic cluster one cell at a time; each step solve Laplace (cluster pinned to V=0, anode at V_APPLIED), add one cluster-adjacent empty cell with probability ∝ V^η. η-sweep: 1→62, 2→85, **3→95 (A)**, 4→95, 6→90. Locked η=3. The Jacobi field is genuinely load-bearing (drives every growth decision).

### Prompt 21 — Restore the field panels (Matan)
> "what about the third panel … I think it is nice to see them since this is what Our jacobi actually solves for"

**Tool:** Claude Code (claude-opus-4-8)
**Output:** Final figure is 3 panels — σ (filament), φ (Jacobi solve), |∇φ| (drives growth, glows at the leading tip).

### Prompt 22 — Log prompts + plan the naive C++ port (Matan)
> "Let's log the prompts and Let's plan the move to the naive.cpp implementation … as naive as possible in terms of caches … plan for future improvements like cache optimization and multiprocessing and I also read about Red-Black Gauss-Seidel solver"

**Tool:** Claude Code (claude-opus-4-8)
**Decisions:** Q16.16 fixed-point everywhere → optimization stages stay bit-identical (`run.sh` `cmp` gate preserved). RB-GS **dropped** once it became clear it uses updated values mid-sweep and so cannot be bit-identical to the Jacobi baseline.
**Output:** Archived old continuum C++ to `old_cpp/`. New `physics_dbm.h` + `dbm_stage0.cpp` (naive AoS `struct Cell{int32 V,metal}`, ~50% cache-line use; warm-started Jacobi hot loop; full-`__int128` V^3 weighted pick to avoid fixed-point underflow near the cathode). Wired `cbram_stage0` to it in `build.sh`/`CMakeLists.txt`. Result: **90/100 A** at N=200 (bridged step 588, 1.5% bright), byte-identical `V_final` across runs. Future ladder (documented): SoA → cache blocking/time-skewing → OpenMP, all bit-identical Jacobi.

### Prompt 23 — 3-panel growth video (Matan)
> "I would like to maybe create a video from it with 3 panels so generate accordingly"

**Tool:** Claude Code (claude-opus-4-8)
**Output:** Changed `dbm_stage0.cpp` per-frame output from a σ-only PPM to raw V+σ binaries (`frames_stage0/frame_NNNNNN_{V,S}.bin`), since the φ and |∇φ| panels need the potential field; always dumps a final frame on the bridged state. New `python_imp/render_dbm_video.py` reads the frame pairs, renders σ | φ | |∇φ|, and encodes to MP4 via ffmpeg (|∇φ| scale fixed from the final frame to avoid flicker). Produced `python_imp/filament_stage0.mp4` (1920×720, 295 frames, ~10 s @ 30 fps). Note: old `make_video.sh` (globs `*.ppm`) no longer applies to the DBM stages.

### Prompt 24 — Verify bit-identical runs (Matan)
> "I don't think 2 runs are bit identical, can you make sure of that?"

**Tool:** Claude Code (claude-opus-4-8)
**Output:** Verified determinism three ways across independent runs at N=200: `V_final` identical (md5 `1aa2bb34…`), `sigma_final` identical (md5 `693987aa…`), and all 590 frame files identical (`diff -rq` clean). Root cause confirmed deterministic: fixed-seed `xorshift64`, row-major candidate order, pure `__int128` integer pick, integer Jacobi, single-threaded, no float in the compute path. Clarified scope: this is same-stage/same-N reproducibility; the stage0-vs-stageN bit-identity gate is separate (`run.sh` `cmp`), and the C++ is intentionally not bit-identical to the Python float ground truth.
