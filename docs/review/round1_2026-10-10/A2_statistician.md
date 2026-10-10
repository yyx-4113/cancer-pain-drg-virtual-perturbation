# A2 — Statistician / Study-design review

**Manuscript:** *Divergent dorsal root ganglion programs in bone cancer pain and chemotherapy-induced peripheral neuropathy: a calibrated in-silico perturbation analysis* (single author)
**Review layer:** Design & statistical validity
**Date:** 2026-10-10
**Independence note:** All judgements below were derived from `manuscript_final.md` and the deposited CSV/JSON artefacts in `results/tables/`; every quantitative claim was recomputed by me. I did not read any other reviewer's output or the author's response material.

---

## Summary verdict

The statistical core of this manuscript is, on the whole, *better executed than most computational-biology submissions*: the exact-permutation engine is correctly built, the BH denominator for the mirrored KO/OE contrasts is correctly reduced to 132, and the marker-level max-Δ selection is honestly folded into a single 102-test BH step. The manuscript also earns credit for killing its own headline result with its own controls. My concerns are therefore not about fraud or sloppiness but about (a) a uniform floor claim that is false for 2 of 6 contrasts, (b) structural under-power that is never disclosed, (c) a validation arm whose weakest statistic (myeloid MWU) is both mislabeled and meaningless at the stated prevalence, and (d) a cross-model concordance that is over-read. Each is detailed below with recomputed evidence.

---

## 【Problem 1】 The "floor 1/4,901" claim is false for 2 of the 6 aetiology pairs; those pairs use a 70-null (floor 1/71).

**【Evidence】** `manuscript_final.md:19` (abstract) and `:85` (methods) and `:165` (limitations) all state "C(8,4)² = 4,900 null pairs … its floor is 1/4,901." But in `etiology_bootstrap_differences.csv`, `etiology_module_logfc_cross.csv` and `etiology_marker_cross.csv` the `paired_control=True` rows for `PTX|OXA` and `Bone_D14|Bone_D7` carry `n_null = 70` with `p_permutation_exact = 0.042254 = 3/71`, not 4,900/1/4,901. The paired enumeration is the *correct* one (permuting 8 shared-control samples gives C(8,4)=70), so the implementation is right — only the prose is wrong and over-states resolution.

**【Why it matters】** The abstract and limitations sell "exact p with floor 1/4,901" as a uniform property. For the two within-study, shared-control contrasts the true floor is ~22× coarser (1/71 ≈ 0.0141). A reader replicating the "exact" claim, or trusting the stated resolution, will be misled about how fine the paired tests can discriminate.

**【Specific fix】** Replace every occurrence of "C(8,4)² = 4,900 … floor 1/4,901" with a conditional statement, e.g.: *"For the four contrasts with independent control groups the null is the exhaustive 70×70 = 4,900 (A,B) split pairs, floor 1/4,901. For the two contrasts that share a control group (PTX vs OXA, Bone_D14 vs Bone_D7) the control is fixed and only the 8 case samples are permuted, giving C(8,4) = 70 null pairs, floor 1/71."*

---

## 【Problem 2】 The two paired contrasts are structurally incapable of reaching q < 0.05, so "no significant divergence" there is a power artefact, not a finding.

**【Problem】** (one sentence — see above.) The paired pairs have a BH floor far above 0.05 by construction, yet are presented alongside the four independent pairs as if symmetrically informative.

**【Evidence】** For modules, BH is over 132 independent tests (`etiology_summary.json:47`); for the paired pairs the smallest achievable exact p is 1/71 = 0.01408, so the smallest possible q is ≈ 0.01408 × 132 = **1.86** (and the reported paired module q-values are 0.099, 0.268, 0.485, 0.659, 0.944, 1.000 — none below 0.05; see `etiology_module_logfc_cross.csv`). For markers, BH over 102 (`etiology_marker_cross.csv`, verified: recomputed min q = 0.012487 with denom 102, not 17) gives a paired floor of ≈ 0.01408 × 102 = **1.44**. So `PTX|OXA` and `Bone_D14|Bone_D7` can *never* be significant regardless of effect size. In `Table 3` every `PTX|OXA` and `Bone_D14|Bone_D7` row sits at q ≥ 0.099 — consistent with this ceiling, not with biological equivalence.

**【Why it matters】** The reader is invited to conclude the aetiologies "diverge" from the significant rows, but the two within-study, shared-control comparisons are censored by design. Conflating structurally-insignificant with biologically-null over-states the uniformity of the divergence and hides that two of six contrasts could not have spoken even with perfect separation.

**【Specific fix】** Add a sentence to §3.5/§5: *"The two shared-control pairs (PTX vs OXA, Bone_D14 vs Bone_D7) use a 70-null and therefore have a BH ceiling of q ≈ 1 for modules and q ≈ 1.4 for markers; they cannot reach significance at n = 4 and should be read as under-powered rather than as evidence of equivalence."*

---

## 【Problem 3】 The degree-matched randset null is itself confounded by cell composition, which the authors acknowledge only partially.

**【Evidence】** `etiology_bootstrap_differences.csv` (per_contrast, deduplicated to 132 KO): 32/132 tests are "significant" against the randset null (q_randset < 0.05), of which 14 are negative-control modules and 18 are targets (`etiology_summary.json:53-56`). The erythrocyte set — a technical control with no pain biology — has z = −5.78 for `PTX|OXA` (`manuscript_final.md:133`) and is in fact the *globally* largest |z| in the family: |z| = **7.39** for `Bone_D14|Bone_D7` (recomputed from the CSV). I verified 32/132 and the 14/18 split exactly.

**【Why it matters】** The manuscript correctly reads this as "the signal is non-specific." But the deeper point is under-stated: a degree-matched null that ranks a red-cell gene set *above every analgesic target* is not merely a neutral negative control — it shows that degree-matching on network degree **fails to remove the compositional confound** the paper's own §3.3/§3.6 identifies as the dominant effect. Erythrocyte genes inherit a large |Δ| because the two DRG preparations differ in perfusion/blood content (ΔlogFC = −2.48, `Table 3`), exactly the same mechanism that moves neuronal modules. So the randset null is calibrated against the wrong nuisance axis.

**【Specific fix】** Either (i) state explicitly that degree-matching is *insufficient* because network degree correlates with cell-type composition, and that the erythrocyte extreme-z is the composition confound leaking through the null; or (ii) replace the degree-matched null with a null matched on a compositional proxy (e.g., module mean logFC vs the erythrocyte/pan-neuronal control, or a null drawn from modules of equal "cell-compartment load"). Sentence to add: *"Because a hub-rich module and a perfusion-driven control share the same topological footprint, degree-matching alone cannot separate target pharmacology from preparation composition; the erythrocyte z = −5.78 … −7.39 demonstrates this residual confound."*

---

## 【Problem 4】 500 degree-matched draws is thin for a z that enters a 132-test BH; the randset p-floor (2/501) is coarse relative to the permutation floor (1/4901).

**【Evidence】** `etiology_summary.json:4` `n_randset = 500`. The randset z = (Δ_obs − randset_mean)/randset_sd with `p_randset` floored at 2/501 ≈ 0.0040 (observed `p_randset = 0.003992` in the CSV). After BH over 132, `min_q_randset_per_contrast = 0.02196` (`etiology_summary.json:52`). With only 500 draws the randset_sd (hence z) carries non-trivial Monte-Carlo error, and the 2/501 floor is ~20× coarser than the exact-permutation 1/4901 floor used for the primary endpoint.

**【Why it matters】** The randset null is used to make the strong claim that the perturbation hierarchy is non-specific. That claim rests on a *secondary* null whose resolution is an order of magnitude worse than the primary test and whose normality/extremity is Monte-Carlo-noisy. The conclusion is robust (the erythrocyte control is decisive either way), but the z-values and their q's should not be quoted to 2 decimals as if exact.

**【Specific fix】** Increase `n_randset` to ≥ 2,000 (cheap) and report z and q to 1 significant figure, or replace the Gaussian z with an exact empirical rank within the 500 draws plus the observed (i.e., floor (rank+1)/(501)). State the floor explicitly.

---

## 【Problem 5】 The pseudobulk DE on the 10x validation is statistically acceptable but omits the composition covariate the authors correctly used in the bulk arm, given that cell composition differs between groups.

**【Evidence】** `manuscript_final.md:93` collapses 10x counts to pseudobulk (sum → log2(CPM+1)) and runs limma-trend with n = 4/group. This is the standard, pseudoreplication-avoiding approach and is valid for *sample-level* inference. However `sc_GSE292522_composition.csv` shows the two groups do **not** have the same composition: neuron fraction 0.792→0.684, Schwann 0.189→0.267, fibroblast 0.004→0.033 (B16F10 higher in the latter two). The bulk pipeline explicitly added a muscle-contamination covariate (`manuscript_final.md:55-57`); the sc pseudobulk adds **none**. The bulk pipeline's own central thesis (`§3.6`, `§4`) is that composition moves neuronal markers — yet that control is absent here.

**【Why it matters】** The B16F10 "neuronal-injury programme" (Sprr1a +1.55, Atf3 +0.84) and the "interferon flat" call could both be partly a composition artefact (e.g., Schwann/fibroblast influx diluting neuronal CPM). The validation therefore does not fully honour the methodological standard the paper demands of itself, and weakens the "interferon is model-specific" inference it is meant to support.

**【Specific fix】** Add a cell-composition covariate (e.g., the Schwann+fibroblast fraction, or the neuronal fraction) to the pseudobulk limma model, exactly as the muscle score was added in the discovery arm; or explicitly show the top hits survive after including it. Report the covariate coefficient.

---

## 【Problem 6】 The Mann–Whitney on myeloid fraction (p = 0.97) is mislabeled as "no difference," is two-sided p = 0.114, and is meaningless at ~0.1 % prevalence with n = 4/group.

**【Evidence】** `sc_GSE292522_summary.json:29` reports `myeloid_MWU_p_greater = 0.9714`. Recomputed from `sc_GSE292522_composition.csv`: Con myeloid fraction = {0.00046, 0.00137, 0.00199, 0.00086} (mean 0.00117), B16F10 = {0.00081, 0.00040, 0.00071, 0.00015} (mean 0.00052). Two-sided Mann–Whitney U gives **p = 0.114** (recomputed). The reported 0.97 is a one-sided tail in the *wrong* direction: the data actually show B16F10 **lower** (0.00052 < 0.00117), so a "greater"-tail p near 1.0 is the opposite of the observed effect. Fractions are ~0.05–0.12 %, i.e., ~15–60 cells out of 55,967; cell-type assignment by "max marker-set" makes a fraction this small exquisitely sensitive to a handful of misclassified cells.

**【Why it matters】** p = 0.97 is presented (`manuscript_final.md:145`, `Table 5`) as evidence "myeloid unchanged," but (a) it is the wrong tail — the correct two-sided p is 0.114, which is *not* "no difference," it is simply underpowered; (b) at 0.1 % prevalence the metric is dominated by annotation noise, so even p = 0.114 is not interpretable as a clean null; (c) n = 4/group cannot detect a 2.3× relative change (0.12 % vs 0.05 %) at any conventional power. The "myeloid unchanged in B16F10" conclusion therefore rests on the weakest, noisiest, misreported statistic in the paper, while the stronger evidence (ISG flat in pseudobulk DE) is separate.

**【Specific fix】** Report the two-sided p (0.114), not a one-sided 0.97; state the test direction correctly; and replace the fraction-based MWU with a count-based test (e.g., a binomial/negative-binomial test on myeloid cell counts per sample, or a hurdle model), or simply demote the claim to "myeloid fraction was not higher in B16F10 (two-sided p = 0.11, n = 4/group; note the metric is unreliable below ~0.2 % prevalence)." Do not let p = 0.97 stand as "unchanged."

---

## 【Problem 7】 A near-zero cross-model Spearman ρ (0.03–0.11) is over-read as "model-specific"; it licenses only "uncorrelated," not "opposed," and is confounded by species + platform + quantification.

**【Evidence】** Recomputed from `sc_GSE292522_pseudobulk_DE.csv` vs `DE_D{1,7,14}_vs_Sham_adj.csv` over n = 12,158 shared genes: ρ(D1) = 0.1115, ρ(D7) = 0.0985, ρ(D14) = 0.0299 — matching `manuscript_final.md:145`/`Table 5` (0.112 / 0.098 / 0.030) exactly. But the accompanying permutation p is **highly significant** (p = 5.8e-35, 1.4e-27, 9.6e-4) purely because n = 12,158 genes. The comparison is mouse scRNA-seq (10x) vs rat bulk FPKM — different species, platform, and quantification (counts+CPM vs FPKM), so even a "true" shared program would be attenuated by technical noise.

**【Why it matters】** ρ ≈ 0.1 is *statistically distinguishable from zero* (so "no relationship at all" is also wrong) but *biologically negligible*. "Uncorrelated across all genes" does not license "the interferon program is model-specific / opposed" — absence of global correlation is consistent with (i) genuinely different programs, (ii) one program being a technical subset of the other, or (iii) species/platform noise swamping a real shared signal. The authors do hedge ("model-specific rather than universal"), which is defensible, but the stronger wording in the abstract ("the interferon programme is model-specific") and the use of ρ ≈ 0.1 as positive evidence of *divergence* slightly over-reaches. The decisive, clean evidence is the ISG-flat finding in the sc DE, not the global ρ.

**【Specific fix】** Reframe: *"Global logFC concordance between the B16F10 pseudobulk and the Walker-256 bulk is weak (ρ = 0.03–0.11; significant only because of n = 12,158 genes) and is confounded by species, platform and quantification, so it argues against over-generalising a single model rather than proving the programs are opposed. The specific claim that the type-I interferon program does not generalise rests on the ISG-flat pseudobulk DE, not on the global ρ."*

---

## 【Problem 8】 The within-study bootstrap CI at n = 4 is coarse and its "0 degenerate resamples" claim is implausible; the two null families are not jointly corrected.

**【Evidence】** `manuscript_final.md:87` supplies a 95 % CI from 2,000 bootstrap resamples with cases/controls drawn separately; `etiology_summary.json:45` `max_degenerate_boot = 0`. With n = 4 per group, a bootstrap of 4 with replacement has a ≈ 1.6 % chance of drawing all four samples identical per group (4/4⁴), i.e., ~30+ of 2,000 resamples degenerate — `max_degenerate_boot = 0` is therefore surprising and should be reconciled with the code. The CIs themselves (e.g., glial/immune Δ 0.274, 95 % CI 0.203–0.287, `etiology_bootstrap_differences.csv`) are plausibly reported but rest on a discrete, low-resolution null: at n = 4 the bootstrap of a difference of means has very few attainable values, so nominal 95 % coverage is optimistic. Separately, the permutation family (BH over 132) and the randset family (BH over 132) are kept separate and never jointly corrected — acceptable, but the marker family (BH over 102) and module family (BH over 132) and reversal family (BH over 132) are three distinct FDR controls presented as one narrative; no cross-family error budget is stated.

**【Why it matters】** Bootstrap CIs are presented as if they carry the same exactness as the permutation p; at n = 4 they do not, and a "0 degenerate" flag that contradicts the sampling maths undermines confidence in the interval code. The absence of a joint family correction is minor (each family is internally valid) but should be stated so the reader does not infer a single coherent FDR across all three endpoint types.

**【Specific fix】** Either show the bootstrap resamples exclude all-identical draws by construction (and document it), or report the expected degeneracy; widen the stated interpretation of the CIs ("approximate, n = 4"). Add one sentence: *"Three FDR families are reported (marker logFC, module logFC, reversal), each BH-corrected within family; no cross-family correction is applied."*

---

## 【Problem 9】 The max-Δ selection across the 6 pairs is correctly folded into the BH (102 tests) — this is done right and is noted as a positive, but the paired-pair treatment inside that BH mixes two permutation spaces.

**【Evidence】** `etiology_marker_cross.csv` contains 102 rows (17 genes × 6 pairs); I recomputed BH over 102 and reproduced `q_permutation` to 2e-16 (min q = 0.012487). So selecting the "largest Δ" per gene for display does **not** inflate Type I error, because all 102 (gene,pair) tests enter one BH step. Good. However that single BH mixes paired rows (null size 70, floor 1/71) with independent rows (null size 4,900, floor 1/4,901); both are valid exact p-values so BH is technically admissible, but the mixed floors mean the effective null distribution of the p-values is heterogeneous — worth one disclosure sentence.

**【Why it matters】** This is the one place a naive reading would suspect hidden multiplicity; it is actually handled. The only refinement needed is transparency that the 102-test BH combines two different permutation spaces.

**【Specific fix】** One sentence: *"The 102 marker tests (17 genes × 6 pairs) enter a single BH step, so the displayed max-Δ per gene is already multiplicity-corrected; note the set mixes 70-null (paired) and 4,900-null (independent) p-values."*

---

## Stands up (what is sound and should stay)

1. **The exact-permutation engine and the mirrored-KO/OE BH correction are correct.** I verified `reversal(KO) = −reversal(OE)` to max|Δ| = 0 and recomputed BH over 132: min q = 0.165447 (per_contrast) and 0.116711 (pooled), matching the deposited `q_permutation` to 1e-16; 0 tests survive q < 0.05. The "no module reversal survives correction" conclusion is exactly right. (`etiology_bootstrap_differences.csv`, `etiology_summary.json`.)
2. **The negative-control design is genuinely informative and honestly interpreted.** Five a priori controls (neuronal, myelin, housekeeping, ISG, erythrocyte) expose the compartment/composition effect in §3.3 and the non-specificity of the perturbation ranking in §3.6; the erythrocyte ΔlogFC = −2.48 is a real, clever measurement of preparation difference.
3. **The paired-control permutation is implemented correctly** (70-null, not a bogus 4,900), and the marker max-Δ selection is honestly embedded in a 102-test BH. These are above the field's usual standard.
4. **The pseudobulk DE approach is the right call** for the sc validation (avoids pseudoreplication); the replication of the neuronal-injury program (Sprr1a/Atf3) and the flat ISGs are credible at the count level.

---

## Questions for the authors

1. In `etiology_summary.json` `max_degenerate_boot = 0`: with n = 4 per group, how are bootstrap resamples drawn so that no all-identical draw occurs? Please show the resampling code or the expected degeneracy.
2. For the paired pairs (`PTX|OXA`, `Bone_D14|Bone_D7`), do you agree the BH ceiling is ~1.0 (modules) / ~1.4 (markers), making them structurally incapable of significance? Should §3.5 disclose this rather than present them as null?
3. The myeloid MWU is reported as p = 0.97 ("unchanged"). My two-sided recomputation is 0.114, and the direction shows B16F10 *lower*. Was a one-sided tail reported, and in which direction? Will you switch to a count-based test or at least report the two-sided p?
4. The sc pseudobulk omits the composition covariate you correctly used in the bulk arm, even though neuron/Schwann/fibroblast fractions differ between groups. Do the top B16F10 hits survive after adding a composition covariate?
5. The cross-model ρ is significant (p ≈ 1e-27) yet trivial (ρ ≈ 0.1). Do you agree the decisive evidence for "interferon not generalisable" is the ISG-flat pseudobulk DE, not the global ρ, and will you reframe accordingly?
6. Will you correct the uniform "floor 1/4,901" claim to the paired 1/71 in the abstract, methods and limitations?

---

## What I actually checked

- **Mirror / BH (Q1):** Read `etiology_bootstrap_differences.csv` (528 rows, 2 scopes). Confirmed max|delta_KO + delta_OE| = 0 and identical exact p across all 132 module-pairs; recomputed BH over 132 → min q = 0.165447 (per_contrast) and 0.116711 (pooled_by_study), matching deposited `q_permutation` to 1e-16; 0/132 significant at q < 0.05. Matched `etiology_summary.json` (`n_independent_tests_per_scope=132`, `min_q_permutation_per_contrast=0.165447`).
- **Exact-permutation floors (Q2):** Confirmed `paired_control=True` rows for `PTX|OXA` and `Bone_D14|Bone_D7` carry `n_null=70`, `p=0.042254=3/71` in all three cross files; non-paired rows carry `n_null=4900`. Verified the prose claim "floor 1/4,901" is false for the two paired pairs.
- **Randset null (Q3/Q4):** From `etiology_bootstrap_differences.csv` (per_contrast, 132 KO) recomputed 32/132 significant (q_randset<0.05), 14 control / 18 target (matches `etiology_summary.json`). Confirmed erythrocyte z = −5.78 (PTX|OXA) and global max |z| = 7.39 (Bone_D14|Bone_D7). Noted `n_randset=500`, floor 2/501.
- **Myeloid MWU (Q5):** Recomputed from `sc_GSE292522_composition.csv`: two-sided Mann–Whitney p = 0.114 (author reported one-sided "p_greater" = 0.9714). Fractions 0.05–0.12 %.
- **Cross-model ρ (Q6):** Recomputed Spearman(B16F10 pseudobulk logFC, Walker-256 bulk logFC) over 12,158 shared genes: D1 = 0.1115, D7 = 0.0985, D14 = 0.0299 — exact match to manuscript; noted corresponding p = 5.8e-35 / 1.4e-27 / 9.6e-4.
- **Marker BH (Q7/Q9):** Read `etiology_marker_cross.csv` (102 rows). Recomputed BH over 102 → min q = 0.012487 (matches to 2e-16); BH over 17 gives 0.002081 (mismatch), confirming the correct denominator is 102 and that max-Δ selection is already corrected.
- **Module logFC paired handling:** `etiology_module_logfc_cross.csv` confirms paired pairs use n_null=70 and their q-values are ≥ 0.099 (structurally censored).
- **Reproducibility claims:** Did not execute the author's pipeline; verified the deposited artefacts are internally consistent with the prose numbers I could recompute (ρ network vs additive = 0.991, `head_to_head_baselines.json`; DE q<0.05 counts; Spearman ρ values). I did **not** read any other reviewer file, the editorial review, response/revision documents, or the deposit SOP.
