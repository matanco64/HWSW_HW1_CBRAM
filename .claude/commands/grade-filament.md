# Grade CBRAM Filament

Grade a CBRAM filament simulation output against theoretical lightning-filament criteria,
then record the result in the progress tracker.

## Usage

`/grade-filament [name]`

`name` is the simulation tag — one of: `stage0`, `improved`, or any custom string.
If omitted, auto-detect the most recently modified `sigma_final_*.bin` in `HWSW_HW1_CBRAM/`.

## Steps to Follow

1. **Locate the binary file.**
   - Full path: `HWSW_HW1_CBRAM/sigma_final_$ARGUMENTS.bin`
   - If `$ARGUMENTS` is empty, run: `ls -t HWSW_HW1_CBRAM/sigma_final_*.bin | head -1` to find the newest file.
   - Verify it exists; if not, tell the user to run the simulation first.

2. **Run the grading script.**
   ```bash
   python3 HWSW_HW1_CBRAM/grade_filament.py HWSW_HW1_CBRAM/sigma_final_$ARGUMENTS.bin
   ```
   Capture the full output.

3. **Display the grade report** to the user verbatim (everything up to and including the `Verdict` line).

4. **Parse the `__GRADE_SUMMARY__` line** at the end of the output.
   Format: `__GRADE_SUMMARY__ <score> <letter> <filename>`
   Extract: score (integer), letter grade (A/B/C/D/F), filename.

5. **Determine the step number** by counting the data rows already in `HWSW_HW1_CBRAM/cbram_progress.md`
   (rows that start with `|` and contain a digit in the first cell). Next step = count + 1.
   Step 0 is the baseline (naive / stage0).

6. **Determine the key change description.**
   - If grading `stage0`: "baseline — no improvements"
   - Otherwise: look at what was most recently changed in `improved.cpp` and summarize in ≤8 words.
     If you can't determine this, ask the user to describe the change briefly.

7. **Append a row to `HWSW_HW1_CBRAM/cbram_progress.md`:**
   ```
   | <step> | <name> | <key change> | <score>/100 | <grade> | <verdict (short)> |
   ```
   Use the `Verdict` line from the report, shortened to ≤10 words.

8. **Show a one-line summary** after appending:
   `Step <N> recorded: <score>/100 (Grade <letter>) — <verdict short>`

## Notes
- The grading script infers N automatically from the file size.
- Run `make_video.sh` manually to generate animations from frames (not part of this skill).
- To inspect the frames visually: `frames_stage0/` or `frames_improved/` contain `.ppm` files.
