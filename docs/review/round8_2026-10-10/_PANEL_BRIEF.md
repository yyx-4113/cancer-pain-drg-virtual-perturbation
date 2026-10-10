## Independence discipline (mandatory)

This is **Round 8** of an independent review panel. The manuscript has survived several prior rounds. You MUST treat it as a **first submission** and form every judgement from the text and source data you read yourself.

**Forbidden to read (do not open, do not grep, do not ask another reviewer about them):**
- `docs/review/REVIEW_round*.md`, `docs/review/RESPONSE_*.md`, `docs/review/REVISION_*.md`
- `docs/review/round1_*/`, `round2_*/`, `round3_*/`, `round4_*/`, `round5_*/`, `round6_*/`, `round7_*/` (any prior round dir)
- `author_verification_statement.md`, `GITHUB_DEPOSIT_SOP.md`, `manuscript_v0*.md`
- any `.workbuddy/` memory files, any task-status / project-overview files
- **other reviewers' output files in this `round8_2026-10-10/` directory** (read only your own)

Do NOT assume the manuscript is mature or has passed previous review. Any claim you CAN verify, you MUST verify yourself.

## Output contract

For every item, four mandatory parts:
- 【Problem】 one sentence
- 【Evidence】 pinned to file:line, table/section, or exact numbers; numbers you cite must be ones you recomputed yourself from the artefacts
- 【Why it matters】 concrete effect on conclusions / credibility / acceptance
- 【Specific fix】 a paste-ready English replacement sentence, or an explicit spec for a new analysis

"Consider strengthening the discussion" is banned.

## Also required
- § Stands up (≥3, with evidence) — explicitly mark things you suspected but found correct.
- § Questions for the authors — state what you need to know; do not guess answers.
- § What I actually checked — files read, commands run, values recomputed vs the manuscript, with discrepancy stated.

## Environment notes (for recomputation)
- Managed Python: `C:/Users/Administrator/.workbuddy/binaries/python/envs/default/Scripts/python.exe`
- Use **PowerShell** (not Bash) for shell commands; set `$env:PYTHONIOENCODING="utf-8"` before running Python (the console is GBK and will otherwise crash on the `−` glyph).
- Runnability gate: `scripts/98_check_manuscript_numbers.py manuscript_final.md` (expect 291 checks, 0 failed); unit tests: `scripts/99_unit_tests.py` (expect ALL TESTS PASSED).
- Artefacts live in `results/tables/` (CSV/JSON); the gate re-derives tables from them at run time.
- Do not mention what tools you use; write review comments only.
