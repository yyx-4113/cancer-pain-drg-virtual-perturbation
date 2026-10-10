# R5-A1 Review — cancer-pain / DRG neuroimmunology

**Manuscript:** `docs/manuscript/manuscript_final.md` (first submission)
**Reviewer:** R5-A1, independent. Read only the panel brief, the manuscript, `results/tables/`, `scripts/98_check_manuscript_numbers.py`, and the figures. No round-1–4 material, no `manuscript_v0.*`, no `.workbuddy/` file, no other reviewer output was opened.

---

## 1. Verdict

**Accept.**

The manuscript is publishable, as it stands, at a soundness-oriented SCIE venue. I recomputed the domain-critical numbers from the deposited CSVs/JSON (marker Δ and q, module mean *t* and logFC, dedup mirror *r* = −1.000, T1/T4/T5 statistics, permutation floors, Cohen's *d*) and every one matches the text. The gate `scripts/98_check_manuscript_numbers.py` ends `291 checks, 0 failed`. The biological claims are supported by marker numbers, the glial/myeloid narrative is internally consistent across §3.3/§3.5/§3.8/Discussion, the interferon finding is scoped correctly to the Walker-256 day-7 peak, and §3.8's non-replication is stated fairly. Every DOI I relied on resolves and says what is claimed. I found **no blocking issue**; the two items below are cosmetic and do not require another round.

---

## 2. Summary of the manuscript

A single-author computational reanalysis of two bulk rat DRG RNA-seq studies on the same platform (GSE149648 Walker-256 tibia bone cancer pain; GSE160543 paclitaxel/oxaliplatin CIPN; n = 4/group), compared by exhaustive label permutation (C(8,4)² = 4,900 pairs, exact *p*), plus an in-silico perturbation layer calibrated against five negative-control modules, degree-matched nulls, an oracle and an additive/ridge baseline, plus an independent mouse single-cell validation (GSE292522, B16F10 cancer-nerve-injury, n = 4 + 4, 55,967 cells). Two claims: (1) the two aetiologies drive qualitatively different, exactly testable DRG programmes — myeloid/interferon-driven in Walker-256 at day 7, neuronal-injury-driven in CIPN; (2) in-silico perturbation of bulk transcriptomes cannot rank analgesic targets — it fails its own negative controls and is anti-selective (negative controls carry larger standardised effects than targets). The interferon programme is shown to be model-specific, not reproduced in B16F10.

---

## 3. Major issues

**None.** I could not construct a blocking scientific or numerical objection. The points the panel brief asks to re-check all hold:

- **Negative layer internal consistency (Abstract / §3.6 / Limitations / Conclusions).** The anti-selectivity triad is stated identically in all four places and matches `negative_result_quantification.json`: median |z| 2.47 (control) vs 1.07 (target), MWU *p* = 0.0019 (T1); MDE ≈ 0.052, median target |Δ| = 0.031, 41/102 above MDE, 31/132 exact *p* < 0.05, minimum *q* = 0.165 (T5); 14/32 randset-significant are controls, binomial *p* = 0.0067 (T2); target percentile 44.7 (T4). All verified.
- **§3.8 single-cell claims vs artefacts.** 1,517 up / 703 down; Sprr1a +1.55 (*t* = 54.5); Ifit2 ↓ (−0.048, *q* = 0.042), Oas1c ↓ (*q* = 0.014); Gfap +0.004 (*q* = 0.006); Cd68/Aif1 ≈ 0; myeloid 0.117% vs 0.052%, *d* = 1.267 vs MDE 2.38, U = 2/14, *p* = 0.114/0.057; neuronal 79%→67%, Schwann 19%→29% (both *p* = 0.029 = the *n*=4+4 MWU floor 2/70); erythrocyte bin 2 control / 4 tumour; 46 myeloid + 1 T cell; 7 Gfap UMI across 8 animals. All verified.
- **Tables 2/3/5 vs CSVs.** Regenerated cell-by-cell; all match (incl. Table 1, checked independently for both mean *t* and mean logFC).
- **Reference order and DOIs.** §References numbering runs strictly 1→16 in first-appearance order; DOIs [8,9,10,11,14] resolve and say what is claimed.
- **Figures vs captions.** fig1 (design), fig2 (seven markers, values = Table 2), fig3 (132 points, median lines 1.07/2.47, "p = 0.0019") all match their captions.
- **Declarations.** Funding, COI, CRediT, Acknowledgements, Ethics, generative-AI declaration, and data/code availability are all present and coherent for a secondary reanalysis.

---

## 4. Minor issues

1. **§2.3 (line 59): "the adjusted value is the more conservative and is the one quoted."** For the myeloid headline marker the muscle-adjusted Cd68 at day 7 is **+2.86** (`etiology_marker_de.csv`, Bone_D7 Cd68 = 2.858) while the unadjusted Table 2 value is **+2.84** — i.e. the adjusted value is very slightly *larger*, not more conservative (Gfap, by contrast, is conservative: 1.33 vs 1.50; Aif1 2.18 vs 2.19; Ccl2 0.81 vs 0.90). The sentence is a general rule; for Cd68 the direction is inverted by 0.02 log₂. Either quote the unadjusted value in §3.3 or soften "more conservative" to "the adjusted fit is quoted for consistency". Cosmetic; no number in the manuscript changes.

2. **Artefact label, not a manuscript error.** `results/tables/sc_GSE292522_celltype_summary.json` names the annotation source "Usoskin 2015 DRG markers, ref [16]", whereas the manuscript consistently cites Usoskin as **[15]** (Luger is [16]). The manuscript is internally correct; the artefact's internal index is stale and could confuse a replicator. Optional tidy-up.

3. **§5 vs §3.6 rounding of the MDE.** §3.6 gives "≈0.052 log-units" and the Limitations give "≈0.05 log-units" for the same quantity (T5 = 2.802 × 0.0185 = 0.0518). Intentional rounding, but a reader cross-checking may pause. No action required.

---

## 5. Numbers you independently recomputed

Recomputed with `C:/Users/Administrator/.workbuddy/binaries/python/envs/default/Scripts/python.exe` from `results/tables/`.

| # | Claim (manuscript) | Manuscript | Recomputed | Match? |
|---|---|---|---|---|
| 1 | Cd68 Δ, Bone-D7 vs PTX (Table 2; Abstract) | +3.42, *q* = 0.012 | 2.8390 − (−0.5793) = **+3.4183**, *q* = **0.01249** | ✔ |
| 2 | Aif1 Δ, Bone-D7 vs PTX | +3.08, *q* = 0.012 | 2.1922 − (−0.8835) = **+3.0757**, *q* = **0.01249** | ✔ |
| 3 | Gal Δ, Bone-D7 vs PTX | −2.59, *q* = 0.012 | −0.6724 − 1.9197 = **−2.5921**, *q* = **0.01249** | ✔ |
| 4 | Cd68 / Aif1 Δ, Bone-D14 vs PTX | +3.07 / +2.84 | **+3.0710 / +2.8378** | ✔ |
| 5 | Exact-permutation floor, unpaired logFC | 3/4,901 ≈ 6.1×10⁻⁴ | 3/4901 = **0.00061212** | ✔ |
| 6 | Exact-permutation floor, paired reversal | 2/71 ≈ 0.028 | paired reversal min *p* = **0.028169** = 2/71 | ✔ |
| 7 | Dedup universe | 132 (102 target / 30 control) | per_contrast & KO: **132 / 102 / 30** | ✔ |
| 8 | KO/OE mirror | *r* = −1.000 | Pearson *r* = **−1.000**, max|KO+OE| = **0.0** | ✔ |
| 9 | T5 MDE = 2.802 × median bootstrap SE | ≈0.052 | SE = **0.01851**, 2.8016·SE = **0.05184** | ✔ |
| 10 | T4 target percentile | 44.7th | median target |z| 1.0729 = **44.70th** | ✔ |
| 11 | T1 anti-selectivity | 2.47 vs 1.07, *p* = 0.0019 | **2.4702 vs 1.0729**, MWU *p* = **0.001914** | ✔ |
| 12 | T1 without interferon control | 2.33 vs 1.07, *p* = 0.019 | **2.3327 vs 1.0729**, *p* = **0.0187** | ✔ |
| 13 | Randset-significant controls | 14/32, base 30/132, binom *p* = 0.0067 | **14/32**, **0.2273**, binomial-overrep *p* = **0.00671** | ✔ |
| 14 | Glial/immune reversal Δ (min *q*) | +0.274 (0.203–0.287), *p* = 0.0014, *q* = 0.165 | **+0.27429**, CI **0.20308–0.28662**, *p* = **0.00143**, *q* = **0.16545** | ✔ |
| 15 | Erythrocyte largest |z| | 7.39 (D14 vs D7); −5.78 (PTX vs OXA) | **7.3878**; **−5.7822** | ✔ |
| 16 | Pan-neuronal control z (§3.6) | −3.74 | min KO |z| = **−3.736** | ✔ |
| 17 | Myeloid-fraction Cohen's *d* (n−2 pooled) | 1.27 | **1.2670** (means 0.1170% / 0.0517%) | ✔ |
| 18 | Interferon module mean *t*, day 7 | +24.9, *q* = 6.8×10⁻⁵ | **24.923**, *q* = **0.000068** | ✔ |
| 19 | Table 1 mean logFC (22 rows) | e.g. +2.61 / +2.01 / +0.87 | all 22 match `module_level_signal.csv` | ✔ |
| 20 | Module ΔlogFC glial/immune D7 vs PTX | +2.70, *q* = 0.0073 | **2.6954**, *q* = **0.007345** | ✔ |

Additional corroborations: PCA 46.8/17.4/7.2/4.8% (recomputed 46.85/17.45/7.23/4.76); sample correlation 0.913–0.994; muscle ANOVA F = 0.637, *p* = 0.605; 13,893 genes; §3.3 muscle-adjusted Cd68 +2.86 and Gfap +1.33 (from `etiology_marker_de.csv` 2.858 / 1.328); sc whole-DRG Gal +1.42; the 3.8 non-replication ISG q values (Ifit1 sat-glia 0.0037, Ifit2 neuron 0.0459, Isg15 sat-glia 0.0071); concordance ρ = 0.098 / 0.112 / 0.030.

Reference verification (WebSearch): [14] Baruch et al., *Nature* 2025;646:462–473, doi:10.1038/s41586-025-09370-8 — confirmed, and the abstract confirms the "injured neurons respond by **autonomously** initiating IL-6- and type I interferon-mediated inflammation", i.e. the manuscript's "neuron-intrinsic type-I interferon response" attribution is accurate. [9] Ahlmann-Eltze, Huber & Anders, *Nat Methods* 2025;22(8):1657–1661, doi:10.1038/s41592-025-02772-6 — confirmed ("None outperformed the baselines"). [10] Zhai et al., *Front Mol Neurosci* 2021;14:665085, doi:10.3389/fnmol.2021.665085 — confirmed (rat tibia TCI, DRG RNA-seq). [11] Sun et al., *Neuroscience* 2023;516:1–14, doi:10.1016/j.neuroscience.2023.02.007 — confirmed (GSE160543, paclitaxel/oxaliplatin). [8] Roohani, Huang & Leskovec, *Nat Biotechnol* 2024;42(6):927–935, doi:10.1038/s41587-023-01905-6 — confirmed. No fabricated or mischaracterised citation found.

**Gate audit.** `scripts/98_check_manuscript_numbers.py`: 291 checks, 0 failed. Expectations are **mixed** — the discovery/module/marker blocks hard-code literals (e.g. lines 132–139, 185–220) that are *compared* against the CSVs (so drift would still fail), while the negative-result block (lines 320–402) recomputes MWU/binomial/MDE/percentile at run time from the deposited CSV and requires the derived strings to appear in the manuscript. The design is sound: the runtime-derived block is the one the brief warned about (stale hard-coded expectations), and it is genuinely runtime-derived.

---

## 6. Publishability assessment

**Ready — accept.** The manuscript is internally consistent, quantitatively faithful to its artefacts, appropriately hedged on composition/confound and on the CD45-depleted single-cell blind spot, and honest about its own two recorded bugs. The negative (methodological) result is a legitimate contribution, and its calibration protocol (negative-control modules through the perturbation itself, degree-matched nulls, additive baseline, forced target inclusion) is exactly what a soundness-oriented venue should want. The biological claim rests on large, direction-specific, exactly-tested marker effects.

No further round is required. If the editorial office wishes to fold the two cosmetic minor issues (the §2.3 "more conservative" phrasing and the stale `[16]` label in one JSON) into copy-editing, that is sufficient; neither touches a reported value.
