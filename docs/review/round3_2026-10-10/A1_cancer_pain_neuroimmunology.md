# R3-A1 Review — cancer-pain neuroimmunology expert

Manuscript: `docs/manuscript/manuscript_final.md` (318 lines; line numbers below refer to this file).
All recomputations used `C:/Users/Administrator/.workbuddy/binaries/python/envs/default/Scripts/python.exe`
from the project root, against the deposited CSV/JSON artefacts only.

---

## 1. Verdict

**Minor Revision.**

The two headline claims — (i) Walker-256 bone cancer pain produces a myeloid/type-I-interferon DRG programme at day 7 whereas paclitaxel neuropathy produces a neuronal-injury programme, demonstrable on canonical markers with exact permutation tests; (ii) in-silico perturbation of these bulk datasets is anti-selective and must not be used to rank analgesic targets — are both supported by the artefacts, and I could reproduce essentially every number I checked. The permutation implementation (scripts/22_etiology_comparison.py, scripts/_limma.py:386–426) is correct, including the paired 8-case exchange null and a true-assignment self-check (scripts/22_etiology_comparison.py:265–272). Section 3.8 tells the single-cell compartment story with unusual honesty (myeloid blind spot quantified, localisation circularity named, disagreement with ref [15] disclosed and not adjudicated).

What holds it back from accept is a set of consistency and provenance defects, none of which overturns a conclusion but each of which a careful reader will find: the prose and Table 2 draw the same quantities from two mutually inconsistent deposited tables without disclosure (Major 1); one §3.2 sentence materially misdescribes DE effect sizes (Major 2); the §3.8 compartment count is wrong and one of its QC arguments is a logical slip (Major 3); and several artefact-level bugs (T3 power number, 26-vs-31 contradiction, nonstandard Cohen's d) survive because the checker does not cover them (Minor 1–3).

---

## 2. Summary of the manuscript

A single-author computational reanalysis of public rat DRG bulk RNA-seq (GSE149648, Walker-256 tibial bone cancer pain, sham/D1/D7/D14, n=4; GSE160543, paclitaxel/oxaliplatin/vehicle, n=4; same platform GPL24782), plus an independent mouse single-cell validation (GSE292522, B16F10 cancer-nerve-injury, 4+4 animals, 55,967 cells). Methods: median-of-expressed-genes scaling, log₂ FPKM, skeletal-muscle contamination covariate (20 sarcomere genes; one D7 sample at 6.26% muscle FPKM carried through leave-one-out), limma-trend moderated t, 22 predefined modules (17 analgesic-target + 5 negative-control), exhaustive label permutation over C(8,4)²=4,900 null pairs (floor 3/4,901 unpaired, 3/71 shared-control), and a diffusion-based in-silico perturbation layer calibrated against negative controls, degree-matched nulls, an oracle anchor and additive/ridge baselines.

Results: the D7 bone peak is an interferon/myeloid event (ISG module mean moderated t = +24.9; Cd68 +2.86, Aif1 +2.18 adjusted) with every neuronal module falling in lockstep with the pan-neuronal control (a composition effect). Cross-aetiology, myeloid markers diverge strongly (Cd68 Δ = +3.42, Aif1 Δ = +3.08, both q = 0.012) while the injury neuropeptide galanin diverges oppositely (Δ = −2.59, q = 0.012). The perturbation layer correlates 0.991 with an additive no-network baseline, ranks an erythrocyte control above every analgesic target (|z| = 7.39), and shows anti-selectivity (median |z| 2.47 controls vs 1.07 targets, MWU p = 0.0019) with adequate per-test power (MDE ≈ 0.105; 31/132 tests at exact p < 0.05). The B16F10 single-cell dataset shows a neuronal-injury programme (Sprr1a +1.55, Atf3 +0.84), no ISG up-regulation in any testable compartment, and weak concordance with Walker-256 (ρ = 0.098), so the interferon signature is declared model-specific.

---

## 3. Major issues

### M1. The prose and Table 2 quote the same quantities from two mutually inconsistent deposited tables, and the headline cross-aetiology tests are computed on *unadjusted* logFC without disclosure

**Citation.** §3.3 (L115), §3.5 (L131), Table 2 (L263–273); artefacts `etiology_marker_cross.csv` vs `etiology_marker_de.csv`.

**Evidence (recomputed).** The two artefacts disagree for the same contrast (Bone_D7 vs Sham):

| Gene | `etiology_marker_de.csv` (muscle-adjusted limma) | `etiology_marker_cross.csv` (raw group means) | §3.3 / Table 2 quote |
|---|---|---|---|
| Cd68 | 2.858 | 2.839 | 2.86 (L115) / +2.84 (Table 2) |
| Aif1 | 2.180 | 2.192 | +2.18 / +2.19 |
| Gfap | 1.328 | 1.498 | **+1.33 (L115) / +1.50 (Table 2)** |
| Ccl2 | 0.806 | 0.897 | **+0.81 (L115) / +0.90 (Table 2)** |
| Atf3 (D7) | 0.222 | 0.313 | **+0.22 (L131) / +0.31 (Table 2)** |
| Atf3 (D14) | 0.511 | 0.317 | **a 61% relative discrepancy, never surfaced** |

A reader comparing §3.3's "Gfap +1.33" with Table 2's "Gfap +1.50" for the identical quantity finds two numbers in the same paper; §3.5's "Atf3 … +0.22 at day 7 of bone cancer pain" (L131) directly contradicts Table 2's own Atf3 row (+0.31). Code inspection shows why: `cross_marker_logfc` (scripts/22_etiology_comparison.py:516–563) computes its logFC and the exact permutation Δ/q on **raw log₂FPKM group means** (values identical to `DE_D7_vs_Sham_unadj.csv`, e.g. Cd68 2.8390), while `de_positive_controls` (L400–423) fits a muscle-covariate limma-trend model (source of §3.3/§3.5 prose).

This has a substantive dimension beyond cosmetics: the paper's flagship QC finding is that GSM4506604 (a D7 sample) carries 6.26% sarcomeric FPKM (§3.1, L107), and D7 is one side of the headline contrasts (Bone_D7 vs PTX). Yet the muscle covariate — the paper's own confounding solution — is **not** carried into the cross-aetiology exact-permutation tests that produce the abstract's q = 0.012 claims, and the manuscript nowhere states this. No leave-one-out sensitivity for the cross-aetiology marker/module endpoint is deposited (the LOO columns in `virtual_perturbation_scores.csv` cover only the discovery day-7 reversal scores; the Limitations claim at L177 that "leave-one-out ranges are deposited with the artefacts" is only true for that subset).

**Required fix.** (a) State explicitly, in §2.7 and Table 2/3 footnotes, that the cross-aetiology permutation tests use unadjusted logFC and why. (b) Reconcile §3.3/§3.5 prose with Table 2 — either quote one variant consistently or label the variants. (c) Run the exact-permutation test once on the muscle-adjusted logFC (the infrastructure exists in `de_positive_controls`) or deposit a leave-one-animal-out range for the cross-aetiology Δs; report that the Cd68/Aif1/Gal q-values survive.

### M2. §3.2 misattributes the median |logFC|: the quoted values are medians over all 13,893 tested genes, not over the significant sets, and understate the D14 effect size ten-fold

**Citation.** L111: "0, 4,600 and 235 genes reached q < 0.05 on days 1, 7 and 14 … with median |logFC| of 0.073, 0.158 and 0.119."

**Evidence (recomputed).** From `DE_D*_vs_Sham_adj.csv`:

| Day | median |logFC| among q<0.05 genes | median |logFC| over ALL tested genes | manuscript |
|---|---|---|---|
| D1 | undefined (0 genes) | 0.0731 | 0.073 |
| D7 | 0.3398 | 0.1576 | 0.158 |
| D14 | 1.1994 | 0.1187 | 0.119 |

The sentence structure binds the medians to the significant sets. Read that way, the reader concludes the D7 peak hits have median |logFC| ≈ 0.16 (actually 0.34) and the D14 significant genes have median |logFC| ≈ 0.12 (actually 1.20 — a ten-fold misstatement). The quoted numbers are real but answer a different question than the sentence asks.

**Required fix.** Rewrite the sentence to state the population explicitly (e.g. "median |logFC| across all tested genes was X; among the q < 0.05 hits it was Y").

### M3. §3.8/Table 5: the compartment count is wrong, and "contains no detectable Gfap" is argued from a logFC that cannot distinguish absence from no-change

**Citation.** L157 ("Among the seven testable compartments…"), Table 5 row (L316, "7 testable compartments"); L157 ("contains no detectable Gfap (log₂ fold change 0.000)"); L97 ("the assigned 'neuron' compartment is 11.9% pure").

**Evidence.**
1. **Count.** §2.8 declares 12 candidate compartments and states four were excluded (T cells; NP and TH neurons; myeloid). 12 − 4 = **8**. The deposited `sc_GSE292522_celltype_DE.csv` contains exactly 8 compartments (Neuron, Neuron.NF, Neuron.PE, SatelliteGlia, Schwann, Endothelial, Fibroblast, Erythrocyte). "Seven testable compartments" (L157, Table 5) contradicts the manuscript's own arithmetic and the artefact. I verified the exclusion arithmetic itself is right: myeloid has ≥5 cells in exactly 5 of 8 animals (per-sample counts 3, 10, 12, 5 | 5, 3, 7, 1 = 46 total, matching the "46 cells" claim), T cells 1 cell total, consistent with L97/L157.
2. **Gfap logic.** The evidence offered that the satellite-glia bin "contains no detectable Gfap" is a B16F10-vs-control log₂FC of 0.000 (q = 0.9997, `sc_GSE292522_celltype_DE.csv`). A *differential* expression of exactly zero is also what a genuinely Gfap-expressing bin with no group difference would produce; logFC = 0 supports "no change", not "no expression". The claim may well be true, but the quoted statistic cannot establish it; a mean-expression/mean-count QC number per bin is needed (and is not deposited).
3. **Purity provenance.** The "11.9% pure" neuron-compartment QC figure (L97) has no deposited artefact and no implementing script: no script in `scripts/` computes the Snap25/Tubb3/Syp-positive-and-myelin-negative purity criterion, and no JSON/CSV in `results/tables/` contains it. For a paper whose selling point is that every number is machine-checked, an unverifiable QC number in Methods is a hole — and it is also the number that does the most work in justifying the "compartments are marker-score bins" caveat.

**Required fix.** Correct "seven" → eight (and Table 5); replace or supplement the Gfap-logFC argument with an expression-level QC statistic; deposit the purity computation (script + artefact) or remove the 11.9% figure.

---

## 4. Minor issues

1. **T3 artefact: gene-level power number is junk (near-zero-count genes dominate the median SE).** `negative_result_quantification.json` T3 reports `median_se_stable_genes: 0.0003` → `detectable_abs_log2FC_genes: 0.001`. Recomputed: of 21,369 "stable" (q>0.5) genes in `sc_GSE292522_pseudobulk_DE.csv`, 13,810 have implied SE < 0.001 — these are near-zero-count genes where both logFC and t are degenerate. A realistic gene-level MDE at n=4+4 is ~0.01–0.07 (e.g. Sprr1a's SE = 1.554/54.49 = 0.0285). The manuscript wisely argues power empirically (Ifit2 detected at −0.048, q = 0.042, L155) and never quotes the 0.001, but the deposited JSON asserts an absurdity and the checker does not test it. Fix the computation (filter to expressed genes) or delete the field. *(scripts/33_negative_result_quantification.py, T3 block.)*
2. **26-vs-31 contradiction inside the T5 artefact.** The JSON field `n_exact_permutation_p_lt_0.05` = 31 (manuscript follows this, correctly), but the same JSON's `interpretation` string says "26/132 exact permutation p<0.05". Recomputed from `etiology_bootstrap_differences.csv`: 31 is correct. The stale interpretation string should be fixed; it is the kind of drift the project's own documentation warns about.
3. **Nonstandard Cohen's d for the myeloid fraction.** Manuscript and Table 5 quote d = 1.37 (L155, L315). `scripts/33` computes pooled SD as sqrt(SS₁/7 + SS₂/7); the standard pooled estimator sqrt((SS₁+SS₂)/(n₁+n₂−2)) gives d = 1.267 (I recomputed both; the script's divisor of 7 per group is not n₁+n₂−2 = 6 pooled nor either group's df = 3). Either quote 1.27 or justify the estimator. The conclusion (1.27 or 1.37 both < MDE d = 2.38) is unaffected.
4. **Myeloid fraction range slightly overstated.** L157 says "0.02–0.2% of cells per sample"; recomputed range from `sc_GSE292522_composition.csv` is 0.0147%–0.1994%, i.e. 0.01–0.2%.
5. **"Aif1 and Cd68 are if anything slightly down" (L165).** PTX Aif1 = −0.88 log₂ (Table 2: −0.88); −0.88 is not "slight". Either quantify or soften differently.
6. **"Myeloid infiltration" stated flatly in the Discussion (L165).** Bulk Cd68/Aif1 cannot distinguish infiltration from activation of resident cells; the Limitations (L177) concede exactly this. The Discussion should hedge once ("myeloid-marker accumulation") rather than assert anatomy the data cannot show.
7. **Table 3 unexplained dashes.** The OXA column for Neurotrophin, Ca²⁺ α2δ, Neuropeptides, pan-neuronal (and D14/PTX for erythrocyte) are dashed although the values exist in `etiology_module_logfc_cross.csv` (e.g. Neurotrophin Bone_D7 vs OXA = −0.051). A dash convention ("ns not shown") should be stated in the caption; as printed they read as missing data.
8. **Anti-selectivity statistics lump the interferon module in with the negatives.** T1/T2 count CTRL interferon — which is the paper's own headline *positive* biology at D7 — among the five "control" modules. I reran both excluding it: T1 median |z| 2.33 (controls) vs 1.07 (targets), MWU p = 0.019; T2 controls 10/28 significant vs base rate 24/126 = 19%, binomial p = 0.029. The anti-selectivity conclusion survives but weakens (p = 0.0019 → 0.019), and the driver is clearly the erythrocyte module (5 of the 14 significant controls) — which is fine, because erythrocyte content is pure technique. The manuscript should report this stratification (technical vs compositional vs pathway controls) so the anti-selectivity claim does not look propped up by a module whose "control" status is a labelling convention.
9. **Cross-study exchangeability is assumed, not argued.** The unpaired 4,900-pair permutation null treats rat DRG samples from two laboratories as exchangeable under H₀. Platform/species/tissue matching makes this defensible, but the exact p-values for the borderline markers (Scn10a q = 0.026, Calca q = 0.027) rest on it. One sentence in §2.7 ("the permutation null assumes between-lab offsets are absent under H₀; residual lab effects would…") plus the module-level robustness (glial/immune, cytokine, TLR all q = 0.0073) would close this.
10. **§2.9 overstates the checker.** "A separate checker asserts every number reported below against the deposited CSV and JSON artefacts" (L101). In reality `98_check_manuscript_numbers.py` regenerates Tables 2 and 3 cell-by-cell and runtime-derives T1/T2/T5 (genuinely good — I verified its T1 MWU, T2 binomial and T5 MDE are recomputed from the CSV, not hardcoded, L329–368), but most other expectations are hardcoded constants compared against artefacts, only ~25 prose literals are string-checked, Tables 4 and 5 are not regenerated, and `chk(want=None)` auto-passes (L45–46). Notably, the §3.2 median |logFC| numbers (Major 2) are checked nowhere — which is why that error survived 276/276. Tighten the claim or the checker.
11. **Journal readiness.** The abstract is ~440 words (Background→Conclusions), over the 250–350 cap of most SCIE venues; the data-availability GitHub URL is paired with "to be archived at Zenodo upon acceptance", which is a placeholder a submission system may bounce; no ARRIVE checklist is cited although the brief positions ARRIVE-style reporting as the standard in play (largely moot for secondary data, but say so). Declarations (funding, COI, CRediT, AI-use) are present and adequate.
12. **Marker-panel choices are defensible (checked, not just asserted).** Cd68/Aif1 (myeloid), Ccl2/Il1b (inflammatory), Gal/Calca/Tac1/Npy (neuropeptidergic injury axis), Scn10a/Scn9a (Nav1.8/1.7), Atf3 (universal injury marker), Gfap (satellite glia), Mpz/Nefl (myelin/axonal) — all canonical for DRG pain biology, and the two-aetiology direction pattern in Table 2 matches the clinical pharmacology framing. One addition worth considering: *Atf3* is the field's default injury marker and the manuscript buries its non-significant cross-aetiology q; a sentence on why Atf3 fails the exact test while Gal survives (n=4 noise vs marker composition) would preempt reviewer questions.

---

## 5. Numbers I independently recomputed

| # | Claim (location) | Manuscript | Recomputed | Match |
|---|---|---|---|---|
| 1 | Cd68 Δ Bone_D7 vs PTX (Abstract L21, Table 2) | +3.42, q = 0.012 | +3.4183, q = 0.01249 | ✓ |
| 2 | Aif1 Δ Bone_D7 vs PTX (Abstract, Table 2) | +3.08, q = 0.012 | +3.0757, q = 0.01249 | ✓ |
| 3 | Gal Δ Bone_D7 vs PTX (Abstract, Table 2) | −2.59, q = 0.012 | −2.5921, q = 0.01249 | ✓ |
| 4 | Ccl2 Δ vs OXA / Il1b Δ vs PTX (L129, Table 2) | +1.39 (q 0.017) / +1.23 (q 0.018) | +1.3907 (0.01734) / +1.2305 (0.01821) | ✓ |
| 5 | Calca / Scn10a Δ (L131, Table 2) | −0.77 (q 0.027) / −0.70 (q 0.026) | −0.7698 (0.02722) / −0.7022 (0.02602) | ✓ |
| 6 | Cd68/Aif1 D14 Δ (L129) | +3.07 / +2.84 | +3.0710 / +2.8378 | ✓ |
| 7 | Table 2 marker logFC cells (all 4 aetiology columns × 9 genes) | see L263–273 | all match `etiology_marker_cross.csv` to ±0.005 | ✓ (but see M1: unadjusted provenance) |
| 8 | Table 3 module logFC/Δ cells | see L277–289 | all non-dash cells match `etiology_module_logfc_cross.csv` (glial/immune Δ +2.6954 q 0.0073; erythrocyte −2.4831 q 0.0073; neuropeptides −1.3453 q 0.0145; etc.) | ✓ |
| 9 | Glial/immune reversal Δ (L137, L143) | +0.274, CI 0.203–0.287, p 0.0014, q 0.165 | 0.2743 [0.2031, 0.2866], p 0.00143, q 0.1654; min q over 132 tests = 0.1654, argmin = Glial/immune | ✓ |
| 10 | Erythrocyte extremes (L137) | \|z\| = 7.39 (Bone D14 vs D7); −5.78 (PTX vs OXA); pan-neuronal −3.74 | 7.3878; −5.7822; −3.7363; max target \|z\| = 4.724 < 7.39 ("above every analgesic target" holds) | ✓ |
| 11 | T1 anti-selectivity (L141, Abstract) | median \|z\| 2.47 vs 1.07, MWU p = 0.0019 | 2.4702 vs 1.0729, p = 0.001914 (dedup KO universe, 102 target + 30 control) | ✓ |
| 12 | T2 control over-representation (L139) | 32 sig, 14 controls (44%) vs base 30/132 (23%), binomial p = 0.0067 | 32, 14 (43.8%), base 22.73%, p = 0.00671 | ✓ |
| 13 | T5 power (L143) | MDE ≈ 0.105; 0.274 exceeds 2.6×; 31/132 at p<0.05; median SE 0.019 | SE from CI width/3.92: median 0.01851; MDE d(2,2) = 5.653 × 0.01851 = 0.1046; 0.2743/0.1046 = 2.62; 31 | ✓ |
| 14 | KO/OE mirror (L137, L143) | r = −1.000; identical permutation p | r(KO Δ, OE Δ) = −1.000 (n=132); max\|p_KO − p_OE\| = 0 | ✓ |
| 15 | Table 4 (L293–304) | interferon +0.228/z +4.61/75.4%/p 0.002/q 0.011; glial +0.155/3.27/61.5%; pan-neuronal −0.109/−3.56; Na⁺ −0.090/−2.63; myelin −0.090/−2.71; erythrocyte −0.074/−2.57 | all match `virtual_perturbation_scores.csv` D7 KO rows | ✓ |
| 16 | Network vs additive (L125, L304) | ρ = 0.991 (n=132); mean \|rev\| 0.057 vs 0.033; ridge 0.712; 20.5% oracle | ρ = 0.9911; 0.0568 / 0.0328 / 0.7119; 20.51% | ✓ |
| 17 | sc pseudobulk DE (L155, Table 5) | 1,517 up / 703 down at q<0.05 | 1517 / 703 | ✓ |
| 18 | Sprr1a / Atf3 / Gadd45a / Serpinb1a / Ecel1 (L155) | +1.55 (t 54.5) / +0.84 / listed | 1.5539 (54.49) / 0.8355 / 0.5267 / 0.7716 / 0.1020 | ✓ |
| 19 | ISGs (L155, Table 5) | Ifit2 q 0.042, Oas1c q 0.014 down; none up; Gfap +0.004 q 0.006 | −0.0476 q 0.04165; −0.0291 q 0.01356; 0 up-regulated ISG rows; +0.0036 q 0.00615 | ✓ |
| 20 | Cd68/Aif1 in sc (L155) | unchanged | −0.0005 (q 0.78) / +0.0009 (q 0.92) | ✓ |
| 21 | Myeloid fraction (L155, Table 5) | 0.12% vs 0.05%; 55% drop; d 1.37 vs MDE 2.38; U = 2/14, p 0.114/0.057 | means 0.00117 / 0.00052; ratio 0.442; d = 1.267 (standard) / 1.369 (script 33's estimator); MDE d = 2.386; U = 2, p = 0.1143 / 0.0571 | ✓ except d estimator (Minor 3) |
| 22 | Compartment ISGs (L157, Table 5) | Ifit1 sat q 0.004; Ifit2 neuron q 0.046; Isg15 sat q 0.007; none up anywhere | −0.5230 q 0.0037; −0.0490 q 0.0459; −0.6300 q 0.0071; 0 up rows | ✓ (but "seven" compartments ≠ 8 in file, M3) |
| 23 | Localisation circularity numbers (L157) | Gal +1.47 neuron (q<0.001), +2.54 NF, Scn10a − neuron q 0.012, whole-DRG Gal +1.42 | 1.4720 (q 0.0001); 2.5370 (q≈0); −0.2470 (q 0.0117); 1.4247 | ✓ |
| 24 | Concordance (L155, Table 5) | ρ = 0.098 (D7), 0.112 (D1), 0.030 (D14) | 0.0985 / 0.1115 / 0.0299 (`sc_GSE292522_concordance_per_day.json`) | ✓ |
| 25 | 46 myeloid cells, 1 T cell (L97) | 46 / 1 | Σ(fraction × cells/sample) = 46.0 / 1.0 | ✓ |
| 26 | QC (L107) | corr 0.913–0.994; PCA 46.8/17.4/7.2/4.8; muscle ANOVA F 0.637 p 0.605; GSM4506604 6.26%, next 0.16% | 0.913/0.994; 46.85/17.45/7.23/4.76; matches; 6.2573/0.1604; 1 sample >1% | ✓ |
| 27 | §3.2 median \|logFC\| (L111) | 0.073 / 0.158 / 0.119 | over ALL genes: 0.0731/0.1576/0.1187; among q<0.05 hits: —/0.3398/1.1994 | ✗ (Major 2: wrong population implied) |
| 28 | §3.3 gene list (L115) | Cd68 +2.86, Aif1 +2.18, Gfap +1.33, Ccl2 +0.81, Il1b +0.80 | matches `etiology_marker_de.csv` (2.858/2.180/1.328/0.806/0.804) but *not* `etiology_marker_cross.csv`/Table 2 | ✗ (Major 1: two sources mixed) |
| 29 | §3.5 Atf3 (L131) | +1.13 PTX, +0.22 D7 | marker_de: 1.126/0.222 ✓; marker_cross/Table 2: 1.1219/0.3128 | ✗ (contradicts Table 2's own +0.31) |
| 30 | Ref [15] (L227) | Baruch et al., Nature 646:462–473, 2025; "reports a neuron-intrinsic type-I interferon response" | citation verified (PMID 40836096, DOI correct); source paper's abstract indeed states injured neurons "autonomously initiat[e] IL-6- and type I interferon-mediated inflammation" — the disagreement is described accurately and the refusal to adjudicate is appropriate | ✓ |

Checker run: `python scripts/98_check_manuscript_numbers.py` → **276 checks, 0 failed** (reproduced). The gate is meaningful where it runtime-derives expectations (T1/T2/T5 recomputed from the CSV, Tables 2/3 regenerated cell-by-cell against the manuscript text), but most remaining expectations are hardcoded constants, only ~25 prose literals are string-checked, and Tables 4/5 plus the §3.2/§3.3/§3.5 prose numbers are not covered — which is exactly where issues M1–M2 and Minor 1–2 live (see Minor 10).

---

## 6. What would make this publishable (ordered)

1. **Reconcile the two marker tables and disclose the unadjusted basis of the cross-aetiology tests** (fixes M1). One sentence in §2.7 + a Table 2 footnote; then either run the permutation on muscle-adjusted logFC or deposit its LOO sensitivity. Given `de_positive_controls` already fits the adjusted model, this is an afternoon of work and would remove the only substantive methodological objection I have.
2. **Fix the §3.2 median-|logFC| sentence** (M2) — state both populations explicitly.
3. **Correct §3.8/Table 5 compartment count (8, not 7), fix the Gfap-absence argument, deposit the purity computation** (M3).
4. **Repair the artefact layer**: T3 median-SE computation (filter expressed genes), the 26/132 stale string in T5's interpretation, and the Cohen's d estimator (standard pooled → 1.27, or justify). Extend `98_check_manuscript_numbers.py` to Tables 4/5 and to the §3.2/§3.3/§3.5 prose numbers so these classes of error cannot recur.
5. **Report the control-module stratification** for T1/T2 (erythrocyte vs compositional vs pathway controls; include the exclude-IFN sensitivity: p = 0.019 / 0.029). This makes the anti-selectivity claim stronger, not weaker, because the surviving driver is the purely technical erythrocyte set.
6. **Small language fixes**: "slightly down" for Aif1 −0.88 (L165); hedge "myeloid infiltration" once in the Discussion; explain Table 3 dashes; state the cross-lab exchangeability assumption behind the exact p (Minor 6–9).
7. **Submission hygiene**: cut the abstract to the target venue's cap; replace the Zenodo placeholder with a DOI (or state "will be deposited upon acceptance" in the venue's accepted formula); cite the ARRIVE checklist and note its applicability to secondary data.

None of these requires new data collection or changes a conclusion. The negative-result layer (T1/T2/T5) and the exact-permutation machinery are, in my experience, unusually well built for this manuscript class, and the biological story told in §3.3/§3.5/§3.8 — myeloid/interferon for the Walker-256 tibia model, neuronal injury for CIPN, model-specificity of the interferon signature — is calibrated at exactly the right level of caution.

*— Reviewer A1 (cancer-pain neuroimmunology), Round 3, 2026-10-10.*
