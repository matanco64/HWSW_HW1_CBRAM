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

*(prompts will be added here as implementation proceeds)*
