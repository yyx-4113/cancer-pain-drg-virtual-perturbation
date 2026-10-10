# Reviewer B8 — Independent statistical / causal-inference / reproducible-research review

**Manuscript:** `docs/manuscript/manuscript_final.md`
**Role:** Reviewer B8 (statistician / causal-inference / reproducible research)
**Mode:** Fresh, independent first-submission review. I did not read any prior review/response files, the verification statement, the deposit SOP, `manuscript_v0*.md`, or other reviewers' `round8_*` files. Every judgement below comes from the manuscript text and artefacts I recomputed myself.

---

## Findings (focus layer: design validity, exact-permutation, multiplicity, framing)

### F1 — Exact-permutation floors 3/4,901, 3/71, 2/71 are reproduced and internally consistent

【Problem】 None for the floor *values*; one minor looseness in the justification of the unpaired-reversal floor (see Fix).

【Evidence】 From `etiology_marker_cross.csv` and `etiology_bootstrap_differences.csv` (per_contrast knockout rows, n=132):
- Marker unpaired (paired_control=False): observed min `p_permutation_exact` = 0.00061212; ×4901 = **3.000** → equals 3/4,901.
- Marker paired (paired_control=True): min p = 0.0422535; ×71 = **3.000** → equals 3/71.
- Reversal unpaired: min p = 0.00142828; ×4901 = **7.000** → ≥ floor (floor not empirically attained; smallest *attained* reversal p is 7/4,901).
- Reversal paired: min p = 0.0281690; ×71 = **2.000** → equals 2/71.

The attained minimums are all ≥ the stated floors, and the two *attained* floors (3/71 marker-paired, 2/71 reversal-paired) match the manuscript exactly (`manuscript_final.md:87-89`, `:177`). This validates the "observed and −observed both enumerated → 3/(n+1)" rule and the "reversal complement is not an exact mirror → 2/71" rule.

【Why it matters】 Exact permutation with a stated floor is only honest if the floor is the *actual* smallest attainable value under the enumeration actually used. The data confirm the floors are real lower bounds and that the manuscript never reports a p below them.

【Specific fix】 (Tier‑3, not a blocker) The unpaired-reversal floor is asserted as 3/4,901 "because their logFC difference is antisymmetric" (`manuscript_final.md:89`, parenthetical). This is a hand-wave: the reversal statistic is not the logFC difference, and its floor is properly justified by the *full* mirror (flipping both studies gives −Δ_rev), which is always in the 70×70 grid. The sentence is defensible but should say "because the full label complement of both studies negates Δ_rev" rather than conflating it with the logFC antisymmetry. Also note explicitly that 7/4,901 is the smallest *attained* unpaired-reversal p (the 3/4,901 floor is a bound, not an observed value) so a reader does not think 3/4,901 was attained there.

---

### F2 — Multiplicity dedup (264→132) and BH correctness are exact

【Problem】 None.

【Evidence】 `etiology_bootstrap_differences.csv` per_contrast has 264 rows = 22 modules × 6 aetiology pairs × 2 directions. KO and OE are exact mirrors: across all 132 shared (module, pair) combinations, `max|KO + OE| = 0.0` and `max|p_KO − p_OE| = 0.0` (my recompute). So 264 collapses to 132 independent tests. BH recomputed from the 132 per_contrast knockout `p_permutation_exact` values (`statsmodels multipletests fdr_bh`) matches the stored `q_permutation` to `3.3e-16`. Stored `q_randset` yields exactly 32 significant (14 control, 18 target), matching `negative_result_quantification.json` and the text (`:137`, `:139`).

【Why it matters】 The headline claim "no module's reversal Δ reached q < 0.05; smallest q = 0.165" (`:137`) rests on the dedup being exact and BH correct. Both are.

【Specific fix】 None.

---

### F3 — §2.3 two-logFC mechanism is now CORRECT (covariate re-estimation, not moderation shrinkage)

【Problem】 None. The previously-wrong "moderation shrinkage" framing is gone.

【Evidence】 `manuscript_final.md:59` states the §3.3 value is the limma-trend coefficient "with empirical-Bayes moderation of the standard error only, not of the fold change," and the Table‑3 value is the raw per-sample mean log₂FC with no covariate. I verified the mechanism directly:
- For Cd68/Aif1/Gfap/Ccl2 at Bone_D7, the cross-unadjusted `logFC_A` equals the **no-covariate** limma coefficient exactly: Cd68 cross 2.8390 = `DE_D7_vs_Sham_unadj.csv` 2.8390; Aif1 2.1922 = 2.1922; Gfap 1.4984 = 1.4984; Ccl2 0.8974 = 0.8974 (my recompute vs the adj file).
- Gap (unadjusted minus muscle-adjusted, `etiology_marker_de.csv` vs `etiology_marker_cross.csv`): Gfap **+0.170**, Ccl2 **+0.091**, Cd68 **−0.019**, Aif1 **+0.012**. All small; Cd68/Aif1 essentially zero, exactly as the text says ("2.86 versus 2.84; 2.18 versus 2.19", `:59`).

Because the unadjusted endpoint equals the no-covariate limma coefficient, and limma's EB step shrinks the *variance* (hence SE/t/p) and never the OLS coefficient, the gap is genuinely the muscle-covariate re-estimation of the group coefficient. The mechanism statement is correct.

【Why it matters】 This was a known earlier error (attributing the gap to EB moderation). It is now fixed; the two logFC paths are kept apart and the within-study vs cross-aetiology distinction is honest.

【Specific fix】 (Tier‑3 wording) `:57/:59` calls the unadjusted endpoint "raw per-sample mean log₂ fold change." It is in fact the no-covariate limma OLS coefficient (which, for balanced n=4, equals the simple group-mean logFC). The word "raw" is fine but could invite a reader to think it is a naive ratio-of-means; a one-line clarification ("equals the no-covariate limma coefficient") would pre-empt that.

---

### F4 — T2 binomial is correctly demoted to a NON-independent descriptive gauge

【Problem】 None. It is already framed as a descriptive gauge, not a second formal significance result.

【Evidence】 `manuscript_final.md:139` explicitly states: "binomial p = 0.0067; **this is a non-independent descriptive gauge, not a formal test** — the per-trial success probability depends on module kind and the tests share contrasts, and it captures the same signal as T1 rather than an independent confirmation." My recompute: 32 randset-significant of 132 (base rate control = 30/132 = 0.2273); 14 of 32 are controls → `scipy.stats.binomtest(14, 32, 0.2273, 'greater').pvalue = 0.00671`, matching the text. The same 14-control/18-target split appears in T1's universe, so T2 is indeed the same signal re-expressed.

【Why it matters】 Reporting a second "p = 0.0067" as independent confirmation would be double-counting. The manuscript avoids this; the framing is appropriate.

【Specific fix】 None.

---

### F5 — T3 Ifit2 exemplar is now framed via median-SE MDE, not "well powered" on a near-zero-SE gene

【Problem】 None. The earlier "well powered" framing on a near-zero-SE gene is corrected.

【Evidence】 `manuscript_final.md:155`: "Gene-level DE at q < 0.05 resolves effects down to **~0.49 log₂ units at the median moderated SE** (n = 4 per group); a few low-variance genes such as **Ifit2 (log₂ −0.048, q = 0.042) are significant because their SE is unusually small**, so the myeloid-marker null (Cd68/Aif1 unchanged) is a genuine null for *large* effects, not a generic power claim." `negative_result_quantification.json` T3: `detectable_abs_log2FC_genes = 0.488`, `median_se_expressed_genes = 0.2048`, `mde_cohens_d_genes = 2.38` (0.2048 × 2.38 = 0.487 ≈ 0.49). The Ifit2 exemplar is explicitly used to *illustrate* the SE-driven exception, not as evidence of broad power.

【Why it matters】 A reviewer could otherwise read "Ifit2 significant at q=0.042" as the assay being powered for small effects; the text now correctly pre-empts that.

【Specific fix】 None.

---

### F6 — T5 MDE is correctly framed as a per-test effect-size detectability floor, not target-vs-control discrimination

【Problem】 None. The framing distinction (detectability vs discrimination) is explicit and correct.

【Evidence】 `manuscript_final.md:143`: "This is a **per-module power floor at the Δ = 0 boundary and speaks to effect-size detectability, not to target-vs-control discrimination, which T1 addresses directly (and negatively).**" My recompute from `etiology_bootstrap_differences.csv` (per_contrast KO, 132 rows): median bootstrap SE from CI95 = (hi−lo)/(2·1.96) = **0.01850**; × 2.802 = **0.05185 ≈ 0.052** (matches text). Counts: 41/102 target tests exceed the normal-approx MDE; 31/132 reach exact permutation p < 0.05. Both match `negative_result_quantification.json` (T5) and the text exactly.

【Why it matters】 Conflating "the layer can detect a large Δ" with "the layer discriminates targets from controls" would invert the conclusion. The manuscript keeps them separate: T5 = detectability floor (a design property), T1 = anti-selectivity (a validity verdict). Correct.

【Specific fix】 None. (Note the conservative t-based MDE 0.105 is reported only as an upper bound, consistent with the normal-approx 0.052 being the operative floor — appropriate.)

---

### F7 — T1 anti-selectivity MWU is reproduced; clustered non-independence is acknowledged but mild

【Problem】 None blocking. The MWU p treats 132 effects as independent although they are correlated within the 6 contrasts; the text already flags this.

【Evidence】 `etiology_bootstrap_differences.csv` per_contrast KO, 102 target vs 30 control |z|: `mannwhitneyu(target, control, 'two-sided').pvalue = 0.001914`; medians 1.0729 (target) vs 2.4702 (control). Matches text (`:141`, p=0.0019, 1.07, 2.47) and JSON T1 exactly. Excluding the interferon control still gives median |z| 2.33 vs 1.07, p=0.019 (`:141`), so anti-selectivity is not carried by the single strongest module.

【Why it matters】 The 22 modules within a contrast share the same pain axis and network, so the 132 |z| are not 132 independent draws; a clustered test would be more exact. However the median gap is ~2.3× and survives even after dropping the strongest control, so the conclusion is robust to clustering.

【Specific fix】 (Tier‑3) Optionally add a within-contrast permutation (permute module labels, keep contrast structure) as a sensitivity check; the result will almost certainly stay significant given the effect size, and it would close the only residual caveat on T1.

---

### F8 — Design-level residual a green gate cannot catch

【Problem】 One genuine design assumption worth stating explicitly; not a blocker.

【Evidence】 The cross-aetiology exact permutation treats the 8 samples (4 bone-cancer + 4 CIPN) as exchangeable labels under H₀ (`manuscript_final.md:83-89`). This is a valid exact test of *label exchangeability*, but its interpretation as "aetiology divergence" assumes group assignment is the dominant axis. The manuscript already quantifies the competing batch axis: erythrocyte ΔlogFC = −2.48 (`:145`, `:177`) and the §5 limitation that "platform, species and tissue are matched but batch is not." This is handled honestly — the erythrocyte control is precisely the batch gauge and it is reported, not hidden.

【Why it matters】 A green gate (number-match, BH, mirror symmetry) cannot adjudicate whether "exact p" answers the causal question the authors intend. The authors' own negative-control design (erythrocyte/module compositional controls) is the correct guard, and it is used.

【Specific fix】 (Tier‑3) Consider one sentence in §5 stating that the exact-permutation p tests *label exchangeability* under the stated H₀, and that the erythrocyte Δ is the explicit batch-axis confound check — i.e., make the causal-interpretation boundary explicit. Not required for acceptance.

---

## Stands up

1. **Exact-permutation validity.** All three floors (3/4,901, 3/71, 2/71) are reproduced from the artefacts and the attained minimums never fall below them; the 2/71-vs-3/71 distinction for reversal vs marker paired contrasts is empirically confirmed. This is a genuine strength — no Monte-Carlo error, no seed dependence, and the floor accounting is correct.
2. **Multiplicity discipline.** KO/OE exact mirrors (max|KO+OE|=0, max|p diff|=0) justify 264→132; BH recomputed to machine precision; "no q<0.05" at min q=0.165 is exact. The headline null result on the perturbation layer is arithmetically sound.
3. **Framing of the negative-result quantification (T1/T2/T3/T5).** Every one of the four focus framings the panel asked about is now correct: §2.3 gap = muscle-covariate re-estimation (not EB shrinkage, verified by cross = no-covariate limma coefficient); T2 = non-independent descriptive gauge; T3 Ifit2 = median-SE MDE exemplar (not "well powered"); T5 MDE = per-test detectability floor (not target-vs-control discrimination). The numbers behind all four (p=0.0019, 0.0067, 0.052, 41/102, 31/132) reproduce exactly.
4. **Reproducibility gate.** `scripts/98_check_manuscript_numbers.py` (291 checks, 0 failed) and `scripts/99_unit_tests.py` (14 groups, all passed) both run clean; they recompute T1/T2/T5 and regenerate the table cells from deposited CSV/JSON rather than storing expected literals.

---

## Questions

1. For the unpaired-reversal contrasts, the smallest *attained* p is 7/4,901, not the 3/4,901 floor. Is the 3/4,901 floor for unpaired reversal derived from the full-both-studies mirror (which does negate Δ_rev), or is it carried over from the logFC-endpoint argument? The value is a correct lower bound either way, but the stated justification could be tightened (F1).
2. T1 treats 132 |z| as independent. Have the authors considered a within-contrast clustered permutation as a sensitivity check? (Effect size suggests it will hold; F7.)
3. The cross-aetiology H₀ is label exchangeability. Is the erythrocyte ΔlogFC = −2.48 intended as the explicit batch-confound check, and should that causal-interpretation boundary be stated once in §5? (F8.)

---

## What I actually checked

**Commands run (managed Python `C:/Users/Administrator/.workbuddy/binaries/python/envs/default/Scripts/python.exe`, `$env:PYTHONIOENCODING="utf-8"`):**

1. `_b8_recompute.py` — recomputed from `results/tables/`:
   - Floors: marker_unpaired min p ×4901 = 3.000 (→3/4,901); marker_paired ×71 = 3.000 (→3/71); reversal_unpaired min p ×4901 = 7.000 (bound, not attained); reversal_paired ×71 = 2.000 (→2/71). `3/4901=6.121e-4`, `3/71=4.225e-2`, `2/71=2.817e-2`.
   - T1: n_target=102, n_control=30; median|z| target=1.0729, control=2.4702; MWU two-sided p=0.001914.
   - T2: 32 significant (14 control, 18 target); base 30/132=0.2273; binomtest(14,32,0.2273,'greater')=0.006708.
   - T5: median boot SE=0.018505; ×2.802=0.05185; 41/102 targets > MDE; 31/132 exact p<0.05.
   - §2.3 gaps (muscle-adj `etiology_marker_de.csv` vs cross `etiology_marker_cross.csv`, Bone_D7): Gfap +0.170, Ccl2 +0.091, Cd68 −0.019, Aif1 +0.012.
2. `_b8_dedup.py` — confirmed KO/OE exact mirrors (max|KO+OE|=0.0, max|p diff|=0.0), 132 KO rows, BH recompute matches stored q to 3.3e-16, 32 randset-significant among KO.
3. Cross-vs-no-covariate limma check — confirmed cross-unadjusted `logFC_A` equals `DE_D7_vs_Sham_unadj.csv` coefficient exactly for Cd68/Aif1/Gfap/Ccl2, so the §2.3 gap is purely the muscle covariate.
4. `scripts/98_check_manuscript_numbers.py manuscript_final.md` — **291 checks, 0 failed**; all T1/T2/T5 recomputations and Table 2/3 cell regenerations PASS.
5. `scripts/99_unit_tests.py` — **14 groups, ALL PASSED** (moderated-t recovery, BH, diffusion operator, network structure, KO/OE mirror, identifier-matrix rejection, D14≠D1, degree-matched bins, exhaustive splits C(8,4)=70, two-sided p, ridge LOSO, cosine-batch 2-D, muscle score).

**Discrepancies found:** None between my recomputations and the manuscript/JSON. Minor: `etiology_marker_de.csv` adj Cd68 = 2.858 vs `DE_D7_vs_Sham_adj.csv` adj = 2.862 (0.004, attributable to shared-gene-set vs full-gene fit); immaterial to the 2.86-vs-2.84 narrative. The 98-checker's first run reported "manuscript missing" for a few T1/T2/T5 text strings only because its regex did not match the spaced formatting ("p = 0.0019"); a corrected path run then showed those strings present (PASS). No substantive checker failure.

---

## VERDICT

**ACCEPTABLE — no remaining design/statistics blocker.**

All four panel-mandated focus items (§2.3 two-logFC mechanism, T2 demotion, T3 Ifit2 framing, T5 MDE framing) are correct; the exact-permutation floors and multiplicity dedup/BH are reproduced exactly; T1/T2/T5 numbers match to the printed precision; and the manuscript passes its own 291-check number verifier and 14-group test suite. The only outstanding items are Tier‑3 clarifications (tighten the unpaired-reversal-floor justification, optionally add a clustered T1 sensitivity, state the label-exchangeability/blood-batch boundary in §5), none of which affects the conclusions.
