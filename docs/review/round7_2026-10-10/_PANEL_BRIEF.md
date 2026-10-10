# Round 7 Independent Panel Brief — fresh re-review after Round 6 revisions

**Manuscript:** `docs/manuscript/manuscript_final.md` (single-author computational reanalysis; bone-cancer-pain vs CIPN DRG; exact cross-study permutation + calibrated in-silico perturbation + B16F10 single-cell validation).
**Context:** A Round 6 four-reviewer panel produced a consolidated verdict "Revise (moderate)" with P0–P2 findings. The author has since applied every P0, P1 and the substantive P2 edits and re-ran the reproducibility gate, which now reports **291/0 checks passing** and unit tests all pass. This Round 7 panel is a FRESH, independent re-review of the *revised* manuscript. Treat it as a first submission of the new draft.

**Your independence is enforced.** Do NOT open or grep any of:
- `docs/review/round1..6_*/`, `docs/review/REVIEW_round*`, `docs/review/R3_editorial_review.md`
- `docs/review/round6_2026-10-10/*` (the Round 6 reports and its own _PANEL_BRIEF.md)
- `author_verification_statement.md`, `GITHUB_DEPOSIT_SOP.md`, `manuscript_v0.*.md`
- any other reviewer's file under `docs/review/round7_2026-10-10/`

Assume the manuscript may still have flaws; do not assume it passed review.

## What was revised (verify these landed correctly — your job is to confirm the fixes are sound AND find anything new)
1. **§3.3** heading/lead reframed from "interferon event, not an ion-channel event" to "prominent interferon-stimulated-gene programme and myeloid marker up-regulation, on a background of neuronal-transcript dilution"; added a sentence engaging the source study [10] conclusion and a probe of type-I specificity in §3.8; added citations [17] (STING/TBK1/NF-κB in Walker-256 DRG) and [18] (DRG macrophage infiltration in CIBP).
2. **§3.8 / ref [14] (Baruch 2025 Nature 646:462–473):** replaced the "conflict our design cannot adjudicate" sentence with a narrowed, honest statement — neuronal arm partly adjudicated (Ifit2 down in neurons q=0.046), myeloid arm untestable (CD45 depletion), and notes Baruch's IFN-I is an ATF3-driven nerve-healing programme consistent with the paper's neuronal-injury reading.
3. **§4 Discussion:** dropped "(Gfap up)" from the *secure* conclusions with a caution that the companion mouse dataset shows satellite glia carry almost no Gfap; reframed the clinical-pharmacology sentence from assertion to hypothesis-generating (now cites [16] for anti-resorptive, [1,2] for duloxetine/CIPN).
4. **§2.3:** clarified that §3.3 reports muscle-adjusted limma-trend *moderated* logFC while Table 3 reports *unadjusted per-sample mean* logFC on the same 13,601 genes (agree within ~0.02 for Cd68/Aif1 but ~0.1–0.2 for Gfap/Ccl2).
5. **§3.6:** softened "132 independent tests" to "132 tests (correlated within contrast…)"; softened T2 binomial to an "approximate descriptive gauge" while keeping the exact p = 0.0067.
6. **§3.5:** added citation [19] (Peters 2007, paclitaxel DRG ATF3/macrophage).
7. **§5:** unified MDE "≈0.052" (was ≈0.05) and "6.1×10⁻⁴" (was 6×10⁻⁴).
8. **AI-disclosure (§):** strengthened — names an LLM assistant, attests no AI in design/analysis/interpretation, affirms human authorship of all conclusions.
9. **References:** added [17],[18],[19] (3 real DOIs). Gate confirms all 19 refs carry a DOI.
10. **Table 1** footnote: notes glial/immune module n=5 (one gene dropped at expression filter); **README** GSE292522 control label unified to "non-tumorigenic keratinocytes (MPEK-BL6)".

## Each reviewer's mandatory focus (see individual prompts). Deliverables: a dated review file in `docs/review/round7_2026-10-10/` with the four-part item format (Problem / Evidence / Why-it-matters / Specific-fix), plus §Stands-up, §Questions-for-authors, §What-I-checked. Recompute at least 3–5 key numbers from the deposited artefacts; state explicitly what you verified vs the gate.
