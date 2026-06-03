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
