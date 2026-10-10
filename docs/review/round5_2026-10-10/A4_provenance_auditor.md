# R5-A4 Review — Data provenance / reproducibility auditor

## 1. Verdict

**Minor Revision.** The manuscript is **publishable at a soundness-oriented SCIE venue** once
one concrete, easily-fixed defect is corrected (the data-availability URL currently 404s, §4.1).
No scientific, statistical or provenance claim failed verification: the gate passes
`291 checks, 0 failed`, I independently recomputed **>35 numbers** from the deposited CSVs/JSON
(all match), Tables 1–5 agree cell-for-cell with the artefacts, references run in
first-appearance order with DOIs that resolve, and the figures are generated from the CSVs their
captions cite. Nothing else blocks acceptance.

## 2. Summary of the manuscript

A single-author reanalysis of two public rat DRG bulk RNA-seq studies (GSE149648 bone-cancer pain,
GSE160543 paclitaxel/oxaliplatin CIPN, same platform GPL24782, n=4/group) plus an independent
mouse single-cell validation (GSE292522, B16F10, 4+4 animals, 55,967 cells). Two questions are
kept separate: (Q1) do the two aetiologies drive qualitatively different DRG programmes, tested by
exhaustive label permutation (exact p); and (Q2) can in-silico perturbation of bulk transcriptomes
rank analgesic targets, tested against negative-control modules, degree-matched nulls, an oracle
anchor and additive/ridge baselines. Q1 is positive and large (myeloid/interferon in Walker-256;
neuronal-injury in the others); Q2 is a calibrated methodological negative (anti-selective: controls
out-rank targets). The interferon programme is shown to be model-specific by the single-cell layer.

## 3. Major issues

None. No number was untraceable or contradicted; no table/text/figure conflict; no fabricated
reference. The negative layer is internally consistent across Abstract, §3.6, §5 Limitations and
§6 Conclusions (median |z| 2.47 vs 1.07, p=0.0019; min q=0.165; MDE ≈0.052; 41/102 targets above
MDE; 31/132 exact-p<0.05 — all reproduced).

## 4. Minor issues

1. **Data/code availability URL is dead (the one item to fix).** §"Data and code availability"
   names `https://github.com/yyx-4113/cancer-pain-drg-virtual-perturbation` and says code "will be
   made available at the project repository … upon acceptance". I queried the URL: **HTTP 404**
   (unreachable/non-existent or private). The future-tense "upon acceptance" is defensible, but the
   wording "**the** project repository (URL)" implies an existing repo and presents a link that
   currently fails — a reproducibility wart an editor/copyeditor will hit. Fix: either create the
   public repo, or reword to avoid a specific dead URL (e.g. "code, artefacts and tests will be
   deposited in a public repository and archived at Zenodo upon acceptance"). This is the only
   reason the verdict is Minor rather than Accept.

2. **Gate script (98_check_manuscript_numbers.py): scope of the Table 2/3 regeneration.**
   The script *does* regenerate Table 2/3 cells from the CSVs, but only the four per-aetiology
   logFC columns (`logFC_A`/`logFC_B`); the "Largest Δ" and "q" columns are **not** regenerated
   cell-by-cell — they are covered only by a hardcoded subset of (module/gene × pair) checks
   (lines 185–220). I verified all remaining Δ/q cells by hand (all match), but for the round-5
   "every cell vs the CSV" requirement the gate should regenerate Δ and q for every marker/module
   row, not a subset.

3. **Gate script: several checks are string-presence only.** `in_text(...)` and the final keyword
   loop (lines 390–491) assert substrings such as `"4,600"`, `"0.991"`, `"anti-selectiv"`. These
   are weak guards (a wrong-but-superstring value, e.g. `"124.95"` containing `"24.9"`, would
   pass) and cannot detect a numerically wrong statement that happens to embed the expected token.
   Not tautological (they can and do fail), but the numeric `chk()` re-derivations — not the string
   checks — carry the real assurance. Consider tightening the numeric-bearing strings to regex with
   word boundaries.

4. Cosmetic: Table 1 rows are ordered by mean t descending except Purinergic P2X (+0.69) sits
   above Prostanoid/COX (+1.07); Table 2/3 use a different logFC scale (unadjusted per-sample mean)
   than §3.3 (muscle-adjusted limma) — this dual scale is explicitly declared in §2.3, so it is a
   reader-clarity point, not an inconsistency.

## 5. Numbers you independently recomputed (claim | manuscript | recomputed | match?)

All read at run time from `results/tables/*`; python 3.13.14, pandas 2.3.3.

| # | Claim (source) | Manuscript | Recomputed | Match |
|---|---|---|---|---|
| 1 | Cd68 Δ Bone_D7 vs PTX (`etiology_marker_cross.csv`) | +3.42 | 3.4183 | ✓ |
| 2 | Aif1 Δ Bone_D7 vs PTX | +3.08 | 3.0757 | ✓ |
| 3 | Gal Δ Bone_D7 vs PTX | −2.59 | −2.5921 | ✓ |
| 4 | Cd68 Δ Bone_D14 vs PTX | +3.07 | 3.0710 | ✓ |
| 5 | Ccl2 Δ Bone_D7 vs OXA | +1.39 | 1.3907 | ✓ |
| 6 | Il1b Δ Bone_D7 vs PTX | +1.23 | 1.2305 | ✓ |
| 7 | Calca Δ / Scn10a Δ | −0.77 / −0.70 | −0.7698 / −0.7022 | ✓ |
| 8 | Glial/immune Δ (Table 3) | +2.70 | 2.6954 | ✓ |
| 9 | Neuropeptides Δ | −1.35 | −1.3453 | ✓ |
| 10 | Erythrocyte module Δ Bone_D7 vs OXA | −2.48 | −2.4831 | ✓ |
| 11 | Interferon module mean t D7 (`module_level_signal.csv`) | +24.9 | 24.923 | ✓ |
| 12 | CYP q=6.8e-5 / Glial t=+11.4 | ✓ | 6.8e-5 / 11.412 | ✓ |
| 13 | Interferon KO reversal / z / %oracle (`virtual_perturbation_scores.csv`) | +0.228 / +4.61 / 75 | 0.22825 / 4.611 / 75.43 | ✓ |
| 14 | Spearman(network, additive) (`head_to_head_baselines.json`) | 0.991 | 0.99115 | ✓ |
| 15 | mean |rev| network/additive | 0.057 / 0.033 | 0.0568 / 0.0328 | ✓ |
| 16 | Dedup independent tests (`etiology_bootstrap_differences.csv`) | 132 | 132 | ✓ |
| 17 | min q_permutation (dedup) | 0.165 | 0.16545 | ✓ |
| 18 | randset-significant 32 (14 ctrl / 18 tgt) | 32/14/18 | 32/14/18 | ✓ |
| 19 | largest |z| (erythrocyte, D14 vs D7) | 7.39 | 7.3878 | ✓ |
| 20 | z erythrocyte PTX vs OXA | −5.78 | −5.7822 | ✓ |
| 21 | Mann–Whitney |z| controls vs targets | p=0.0019 | 0.001914 | ✓ |
| 22 | median |z| control / target | 2.47 / 1.07 | 2.4702 / 1.0729 | ✓ |
| 23 | binomial over-representation | 0.0067 | 0.00671 | ✓ |
| 24 | anti-selectivity excl. interferon | 2.33 vs 1.07, p=0.019 | 2.33 / 1.07 / 0.0187 | ✓ |
| 25 | T5 MDE (normal) = 2.802×median SE | ≈0.052 | 0.0518 | ✓ |
| 26 | T5 MDE (t-based) | ≈0.10 | 0.1046 | ✓ |
| 27 | glial Δ / CI / perm p | 0.274 / 0.203–0.287 / 0.0014 | 0.27429 / 0.20307–0.28662 / 0.001428 | ✓ |
| 28 | target tests above MDE (41) / exact-p<0.05 (31/132) | 41 / 31 | 41 / 31 | ✓ |
| 29 | KO–OE mirror r | −1.000 | −1.000 | ✓ |
| 30 | myeloid fractions Con/B16F10 (`sc_GSE292522_summary.json`) | 0.12% / 0.05% | 0.00117 / 0.000518 (55.8% drop) | ✓ |
| 31 | myeloid Cohen's d (pooled, n−2) | 1.27 | 1.2679 | ✓ |
| 32 | MDE Cohen's d (n=4/grp, t, df6) | 2.38 | 2.371 | ✓ |
| 33 | sc pseudobulk DE up/down | 1517 / 703 | 1517 / 703 | ✓ |
| 34 | Sprr1a (+1.55, t=54.5); Ifit2 q=0.042; Cd68/Aif1 Δ≈0 | ✓ | 1.5539/54.49; 0.0417; −0.0005/0.0009 | ✓ |
| 35 | concordance ρ D7/D1/D14 | 0.098/0.112/0.030 | 0.0985/0.1115/0.0299 | ✓ |
| 36 | QC: PCA 46.8/17.4/7.2/4.8; corr 0.913–0.994 | ✓ | 46.85/17.45/7.23/4.76; 0.9133–0.9935 | ✓ |
| 37 | exact-p floors 3/4901, 3/71, 2/71 | ✓ | 6.12e-4, 0.04225, 0.02817 | ✓ |
| 38 | erythrocyte bin 2-control/4-tumour (`sc_GSE292522_composition.csv`) | ✓ | counts 5,4,22,3 / 7,5,26,15 | ✓ |

Gate: `python scripts/98_check_manuscript_numbers.py` → **291 checks, 0 failed**. Unit tests:
`python scripts/99_unit_tests.py` → ALL TESTS PASSED, **14 test groups** (T1–T14) as §2.9 claims.
§2.9 versions match `requirements.txt` exactly (numpy 2.4.6, scipy 1.18.1, pandas 2.3.3,
statsmodels 0.15.0, matplotlib 3.11.1; Python 3.13).

**References.** Numbering runs in strict first-appearance order 1→16 (verified by parsing the
body text); all 16 carry a DOI. 8 DOIs spot-checked: 6 resolve HTTP 200 via `doi.org` and the two
that returned 403 (Wiley, ASCO/JNCCN — anti-bot blocks) resolve to the correct titles via the
Crossref API. Refs [9], [14], [15] titles/journal/year verified against Crossref. **No fabricated
reference.**

**Figures.** `fig1..3.{png,pdf}` exist. `35_figures.py` shows fig2 reads `etiology_marker_cross.csv`
(7 markers = caption) and fig3 reads `etiology_bootstrap_differences.csv` — captions cite the same
CSVs. Consistent.

## 6. Publishability assessment

**Publishable, contingent on one fix.** The science is sound, the provenance is unusually clean
(every headline number traces to a CSV/JSON; the gate is genuine, not cosmetic), and the negative
result is honestly reported. Ordered minimum to publish:

1. Fix the data-availability statement — create the named GitHub repo or reword to a non-dead
   URL/"deposited upon acceptance" (§4 item 1). **Blocking for the Minor-revision decision only.**
2. (Recommended, non-blocking) Extend the gate to regenerate the "Largest Δ" and "q" columns of
   Tables 2/3 for every row, and tighten the string-presence checks to bounded numerics.

With item 1 done, I would raise this to **Accept**.
