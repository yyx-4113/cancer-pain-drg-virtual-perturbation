# R5-A5 Review — Venue Editor (independent journal editor, submission-readiness)

Reader posture: first submission. Sources read: `_PANEL_BRIEF.md`, `manuscript_final.md`,
the deposited CSVs/JSON in `results/tables/`, `scripts/98_check_manuscript_numbers.py`,
and live DOI resolution. No earlier manuscript version, no other file under `docs/review/`,
no `.workbuddy/` file, no manifest/deposit document was opened.

## 1. Verdict

**Minor Revision.** No blocking issue. The manuscript is a coherent, internally consistent,
fully declared soundness study; it is **not yet publishable as-is** only because three
cosmetic house-style items remain (abstract length under a strict counter, table numbering
order, a forward-looking data-availability promise). All are fixable in one editing pass
without new analysis.

## 2. Summary of the manuscript

Single-author computational reanalysis (Yongxin Yang, one affiliation, no post-nominals).
Two bulk rat DRG RNA-seq studies on the same platform (GSE149648 Walker-256 bone cancer
pain; GSE160543 paclitaxel/oxaliplatin) are compared by exhaustive exact permutation, and
an in-silico perturbation layer (17 target modules vs 5 negative controls, degree-matched
nulls, an oracle anchor, additive/ridge baselines) is calibrated against its own controls.
A mouse single-cell dataset (GSE292522, B16F10, 4+4 animals, 55,967 cells) tests whether
the Walker-256 interferon programme generalises. Two claims: a large, direction-specific
biological divergence, and a methodological negative — the perturbation layer is
anti-selective and adds nothing over an additive baseline. 16 references, 5 tables, 3 figures.

## 3. Major issues

None. No stale number, misstatement, unsupported claim, or missing artefact found.
The numbers gate ends `291 checks, 0 failed`.

## 4. Minor issues

1. **Abstract length is convention-dependent and sits on the 300-word line.**
   `manuscript_final.md` §Abstract. Plain whitespace count = **310 words**; with the four
   section labels counted = **314**. The venue checker's own counter
   (`98_check_manuscript_numbers.py:493–499`) strips section labels *and* standalone
   glyphs/superscripts and counts only alnum-bearing tokens, giving **291** ("[PASS]
   abstract <= 300 words: True"). A submission system or editor counting whitespace tokens
   would read > 300. The brief's "295" matches the lenient count only.
   *Fix:* trim ~15 words (e.g. compress the Methods sentence listing calibrations) so the
   abstract is ≤300 under **every** counter, labels included.

2. **Tables are not numbered in order of first citation.** First-mention order in the body
   is **Table 1, Table 4, Table 2, Table 3, Table 5** (Table 4 appears at §3.4, before
   Tables 2–3 at §3.5). All five are cited exactly once each — the earlier "uncited table"
   defect is fixed — but most venues require ascending citation order.
   *Fix:* either renumber the tables to match appearance (perturbation table becomes
   Table 2, etc.) or move the Table 3 reference so the sequence is ascending.

3. **Data/code availability is a forward promise, not a live resource.**
   `## Data and code availability`: code "**will be made available** at the project
   repository (https://github.com/yyx-4113/cancer-pain-drg-virtual-perturbation) and
   archived at Zenodo **upon acceptance**." The wording is honest about a not-yet-created
   repository and does not overclaim — good — but many SCIE venues ask for a resolvable
   URL at review. *Fix:* confirm the GitHub URL is live (even as a skeleton with
   `requirements.txt` and the test suite) or state plainly "repository created and
   anonymised/timestamped on acceptance" and drop the implied live link.

## 5. Numbers I independently recomputed (claim | manuscript | recomputed | match?)

| Claim | Manuscript | Recomputed | Match |
|---|---|---|---|
| Abstract ≤300 words | "295" (brief) / gate PASS | 310 whitespace, 314 w/ labels, 291 gate-counter | **convention-dependent — see M1** |
| Keywords | 6 | 6 (`cancer pain; dorsal root ganglion; …; negative control`) | yes |
| Running head | — | "Divergent DRG programs: bone cancer pain vs CIPN" = 47 chars / 9 words | fine |
| Author post-nominals | single author | "Yongxin Yang¹", no MD/PhD suffix | yes |
| Tables 1–5 each cited | — | 1,1,1,1,1 in body; order 1,4,2,3,5 | order off (M2) |
| Figures cited in order | — | Fig 1, 2, 3; files fig1–3.{png,pdf} present | yes |
| References sequential, DOIs | 16 | first-appearance 1..16; 16 refs, 16 unique DOIs | yes |
| DOI resolution | "all resolve" | 16/16 registered; 5 returned HTTP 200, 11 publisher-403 but confirmed via Crossref API with matching titles | yes |
| Cd68 Δ (D7 vs PTX) | +3.42, q 0.012 | +3.4183, q 0.01249 (`etiology_marker_cross.csv`) | yes |
| Aif1 Δ (D7 vs PTX) | +3.08, q 0.012 | +3.0757, q 0.01249 | yes |
| Gal Δ | −2.59, q 0.012 | −2.5921, q 0.01249 | yes |
| Median \|z\| ctrl vs target | 2.47 vs 1.07, p 0.0019 | 2.4702 vs 1.0729, p 0.001914 | yes |
| Dedup universe / r | 132 tests, r = −1.000 | n=132 (102 tgt/30 ctrl), r=−1.000 | yes |
| T5 MDE | ≈0.052 (normal), 0.10 (t) | 0.0518 / 0.1046; 2.8×SE | yes |
| Single-cell power | d 1.27 vs MDE 2.38 | 1.267 vs 2.38 | yes |
| Gate | 291 / 0 failed | `291 checks, 0 failed` reproduced | yes |

*Gate audit:* the checker's `chk(label, got, want)` compares **runtime-derived artefact
values** (`got`) against **hardcoded literal expectations** (`want`, e.g. `0.0019`,
`3/4,901`), and `in_text()` greps the manuscript for the same literals. It therefore
catches drift between text and artefacts and between artefacts and the frozen constants;
it would not catch a constant that is wrong in the artefact *and* baked into both files,
but no such case was found and every value I recomputed independently agreed.

## 6. Publishability assessment

**Ready?** No — publishable after a Minor Revision pass; no re-analysis, no new artefact,
no reviewer-requested computation.

**Ordered minimum to reach acceptance-ready:**
1. Trim the abstract to ≤300 words under the strictest (whitespace, labels-included) count.
2. Renumber or reorder tables so first-citation order is ascending 1→5.
3. Resolve the data/code-availability wording to a live-or-timestamped resource.
4. (Optional, cosmetic) confirm the running head against the target venue's char limit.

**Desk-reject risk at a soundness-oriented SCIE venue** (negative-result/reanalysis):
low, but not zero. The paper is not a bare "method failed" report — it pairs a robust,
exactly-testable positive (myeloid/interferon vs neuronal-injury divergence, effect sizes
3.0–3.4 log₂, exact permutation) with the negative. The framing is strong enough; the
residual risk is a desk editor reading the title's tail ("negative-control calibration of
in-silico perturbation") as incremental tool-criticism. **Single change that most reduces
risk:** foreground the calibrations as a *transferable, reusable deliverable* — state in
the abstract's Conclusions and the cover letter that the paper hands the field a concrete
negative-control battery (control modules inside the perturbation + degree-matched nulls +
additive baseline + forced target inclusion) that any virtual-perturbation pain study can
adopt, so the contribution reads as a method standard plus a biological finding, not a
single failed experiment. Journal-tier metrics are deliberately omitted: unverified and
venue-specific, and not needed for this judgement.
