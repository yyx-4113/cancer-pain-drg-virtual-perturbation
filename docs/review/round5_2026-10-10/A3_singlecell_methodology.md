# R5-A3 Review — single-cell / pseudobulk methodology

## 1. Verdict (Accept / Minor Revision / Major Revision / Reject)

**Minor Revision** — publishable at a soundness-oriented SCIE venue once the minor
items below are addressed.

**No blocking issue was found.** I was instructed to treat this as a first
submission and to attack §2.8 / §3.8 (single-cell layer) against the deposited
artefacts. Every quantity in my remit reproduces from the CSVs/JSON to the
precision quoted, the inference unit is the animal everywhere a group test is
made, and the five targeted caveats are all present and accurate. The gate script
`scripts/98_check_manuscript_numbers.py` ends `291 checks, 0 failed`. The two
issues I raise are (i) a small disagreement between two deposited "muscle-adjusted"
artefacts that touches one §3.3 gene value and (ii) a framing nit; neither changes
a conclusion.

## 2. Summary of the manuscript

A single-author reanalysis of public rat bulk DRG RNA-seq (GSE149648 Walker-256
bone-cancer pain, GSE160543 paclitaxel/oxaliplatin, n = 4/group, same platform)
compared by exhaustive label permutation, plus a calibrated in-silico perturbation
layer tested against negative controls, degree-matched nulls, an oracle anchor and
an additive/ridge baseline, plus an independent mouse single-cell validation
(GSE292522, B16F10 cancer-nerve-injury, 4+4 animals, 55,967 cells). Two claims:
(1) biological — bone-cancer pain and CIPN drive qualitatively different DRG
programmes (myeloid/interferon vs neuronal-injury); (2) methodological negative —
in-silico perturbation cannot rank analgesic targets (anti-selective and
power-limited).

## 3. Major issues

**None.** No claim in §2.8/§3.8 was found to be unsupported by, or in conflict
with, the deposited artefacts. Specifically:

- **Pseudoreplication (candidate Major — cleared).** The unit of inference is the
  animal everywhere a group test is made in the single-cell layer. The whole-DRG
  pseudobulk (script 31, `31_process_sc.py:122-134`) is gene × *sample* (8 samples).
  The per-compartment extension (script 32, `32_percelltype_sc.py:108-159`) forms a
  *per-animal* pseudobulk per compartment (`agg[ctype][prefix]`, one vector per
  animal) and runs limma-trend across the 4 tumour vs 4 control animals; no
  test touches individual cells. The composition test is a Mann–Whitney U on the 4
  vs 4 *per-animal* fractions. The manuscript's statement at
  `manuscript_final.md:97` ("the unit of inference is the animal, not the cell, so
  the extension introduces no single-cell pseudoreplication") is therefore correct.
  No double-dipping was found in the test statistics; the one shared-marker concern
  (the same neuronal markers define a bin *and* are tested in it) is explicitly
  flagged by the author (`manuscript_final.md:159`, "the genes and the bins are
  defined by the same neuronal markers").

## 4. Minor issues

1. **Two deposited "muscle-adjusted" D7 artefacts disagree on Gfap (and mildly on
   Ccl2).** `DE_D7_vs_Sham_adj.csv` gives Gfap logFC = **1.294**, Ccl2 = **0.787**;
   `etiology_marker_de.csv` (script 22's muscle-adjusted limma, `22_etiology_comparison.py:400-424`)
   gives Gfap = **1.328**, Ccl2 = **0.806**. §3.3 (`manuscript_final.md:115`) quotes
   "+1.33 … Ccl2 +0.81", i.e. the `etiology_marker_de.csv` values, while §2.2/§3.2
   counts and the Table-5 concordance use `DE_D7_vs_Sham_adj.csv`. The gap is ~2.6%
   and changes no conclusion, but a reader cross-checking §3.3 against the canonical
   adjusted file would find 1.29. *Fix:* reconcile the two fits or state which file
   is canonical for the §3.3 marker values (one sentence in §2.3).

2. **"Genuine null for large effects" for Cd68/Aif1 sits in mild tension with the
   CD45 blind spot.** `manuscript_final.md:155` infers the whole-DRG Cd68/Aif1 null
   is "genuine" because the design "detected Ifit2 at log₂ fold change −0.048".
   Power is gene- and cell-type-specific: with myeloid cells at 0.05–0.12% of the
   preparation (0.015–0.2% per sample; 46 cells total), a whole-DRG pseudobulk has
   little sensitivity to a *myeloid-restricted* shift, so an ISG power benchmark does
   not transfer. The same paragraph correctly calls the *fraction* result
   "inconclusive, not evidence of equality". *Fix:* qualify the marker statement as
   "no shift at whole-DRG resolution" to keep the two nulls on the same footing.

3. **Gate expectations are partly hardcoded.** `98_check_manuscript_numbers.py`
   compares many artefacts against literal constants (module-level t at lines
   132-139; Table-2/3 marker dicts at 185-197, 207-215, 290-295). That still catches
   artefact drift (constant vs artefact), but the strongest guarantees are the
   genuinely runtime-derived sections — the T1/T2/T5 recomputation (lines 320-402)
   and the cell-by-cell regeneration of Tables 2/3 from the CSVs (lines 426-480).
   Worth a comment in §2.9 distinguishing the two kinds of check. Cosmetic.

4. **"U = 2/14" notation** (`manuscript_final.md:155` and Table 5) is unusual.
   Recomputation confirms the two U statistics are 2 (tumour) and 14 (control) out of
   a maximum of 16; a parenthetical clarification would help.

## 5. Numbers you independently recomputed

All recomputed at run time from `results/tables/` with the managed Python.

| Claim | Manuscript | Recomputed | Match? |
|---|---|---|---|
| Neuron share, Con → B16F10 | 79% → 67% | 78.59% → 66.95% | yes |
| Schwann share, Con → B16F10 | 19% → 29% | 18.85% → 29.28% | yes |
| Neuron & Schwann MWU two-sided p | 0.029 (both) | 0.0286 (both) | yes |
| Myeloid fraction, Con vs B16F10 | ~0.12% vs ~0.05% | 0.117% vs 0.0517% | yes |
| Myeloid drop | 55% | 55.8% | yes |
| Myeloid MWU | U = 2/14; p₂ = 0.114; p₁ = 0.057 | U_B16=2, U_Con=14; 0.1143; 0.0571 | yes |
| Myeloid cells total / ≥5-cell animals | 46 / 5 of 8 | 46 / 5 (S2,S3,S4,S5,S7) | yes |
| T cells total | 1 | 1 (Sample7) | yes |
| Testable compartments | 8 | 8 (from `sc_GSE292522_celltype_DE.csv`) | yes |
| Erythrocyte bin group counts | 2 control, 4 tumour | 2 Con (S1,S3) / 4 B16 (S5–S8) | yes |
| Satellite-glia Gfap UMI | 7 across 8 animals, 0 in five | 0,0,0,3,0,2,2,0 = 7; five zeros | yes |
| Pseudobulk scale | log₂(counts per 10,000 + 1) | script 31:126 & 32:154, ×1e4 then +1 then log2 | yes |
| Concordance ρ (D7 / D1 / D14) | 0.098 / 0.112 / 0.030 | 0.0985 / 0.1115 / 0.0299 | yes |
| Sprr1a | +1.55 log₂, t = 54.5 | +1.554, t = 54.49 | yes |
| Ifit2 (pseudobulk) | −0.048, q = 0.042 | −0.0476, q = 0.0417 | yes |
| Oas1c (pseudobulk) | q = 0.014 | q = 0.0136 | yes |
| Gfap (pseudobulk) | +0.004, q = 0.006 | +0.0036, q = 0.0061 | yes |
| Cd68 / Aif1 (pseudobulk) | unchanged | logFC −0.0005 / +0.0009 | yes |
| Gal in neurons / NF neurons | +1.47 (q<0.001) / +2.54 | 1.472 (q=0.0001) / 2.537 (q=0) | yes |
| Gal whole-DRG average | +1.42 | +1.425 | yes |
| Scn10a in neurons | q = 0.012 | q = 0.0117 | yes |
| Ifit1 / Isg15 in satellite glia | q = 0.004 / 0.007 | 0.0037 / 0.0071 | yes |
| Neuron ISG composite (Welch p) | 0.93 | 0.9255 | yes |
| Satellite-glia ISG composite p | 0.14 | 0.1418 | yes |
| LOO sign consistency (sg / Schwann) | 8/8, sig ≤ 2 | 8/8; 1 and 2 | yes |
| Myeloid Cohen's d (n−2 denom) | 1.27 | 1.267 | yes |
| MDE Cohen's d (df = 6, 80%) | 2.38 | 2.381 | yes |
| T1 median \|z\| control vs target | 2.47 vs 1.07, p = 0.0019 | 2.4702 vs 1.0729, p = 0.001914 | yes |
| Erythrocyte \|z\| (max) | 7.39 | 7.39 (Bone_D14 vs Bone_D7) | yes |
| Erythrocyte z (PTX vs OXA) | −5.78 | −5.782 | yes |
| T4 target percentile | 44.7 | 44.7 | yes |
| DE counts q<0.05 (D1/D7/D14, adj) | 0 / 4,600 / 235 | 0 / 4,600 / 235 | yes |
| DE counts q<0.05 (unadj) | 1 / 2,563 / 958 | 1 / 2,563 / 958 | yes |
| Gfap D7 bulk (adj, §3.3) | +1.33 | 1.328 (`etiology_marker_de.csv`) vs 1.294 (`DE_D7_vs_Sham_adj.csv`) | see Minor 1 |

Also verified: the abstract is ≤300 words; references run in first-appearance
order 1→16; all 16 DOIs returned HTTP 200 from `api.crossref.org`; `fig1–3.png`
exist; Figures 2 and 3 visually match their captions (Fig 2: seven markers,
bone-d7 red vs paclitaxel blue with Δ connectors; Fig 3: 132 deduplicated tests,
median \|z\| bars at 1.07 and 2.47, MWU p = 0.0019, negative controls red).

## 6. Publishability assessment

**Publishable after minor revision.** There is no blocking issue in the
single-cell/pseudobulk layer: the inference unit is correctly the animal, the
composition, Gfap-purity, compartment-count, concordance and power numbers all
recompute exactly, and the caveat list (argmax-over-marker-sets coarse, overlapping
panels, CD45 blind spot, Baruch disagreement, composition as a marker-mean proxy)
is complete and accurate. Ordered minimum to clear the remaining items:

1. Reconcile or explicitly label the two "muscle-adjusted" D7 artefacts so §3.3's
   Gfap/Ccl2 values trace to a single canonical file (Minor 1).
2. Re-word the Cd68/Aif1 marker null to whole-DRG resolution (Minor 2).
3. Optional: distinguish runtime-derived from hardcoded gate checks in §2.9, and
   gloss the "U = 2/14" notation (Minors 3–4).

None of these requires new analysis or changes any reported conclusion.
