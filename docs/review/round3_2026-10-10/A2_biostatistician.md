# R3-A2 Review — Biostatistician (multiple testing, permutation theory, power)

**Manuscript:** `docs/manuscript/manuscript_final.md` (first-submission posture)
**Recompute script:** `docs/review/round3_2026-10-10/A2_recompute.py` (all numbers below re-derivable by running it from the project root)
**Gate:** `scripts/98_check_manuscript_numbers.py` run by this reviewer: **276 checks, 0 failed**.

---

## 1. Verdict

**Minor Revision.**

The quantitative skeleton of the paper — the exact permutation layer, the KO/OE deduplication, the negative-result quantification T1/T2/T4/T5 — reproduces exactly from the deposited artefacts, and the permutation code in `scripts/22_etiology_comparison.py` is correct in its indexing and orientation (both self-checks run and are enforced at runtime; no bug found). The revision need is concentrated in three statistical claims that are wrong or under-derived as *text*, not as results: (i) the stated floor of the paired permutation test (3/71) is contradicted by the paper's own deposited table, which contains six tests at the true floor 2/71; (ii) the MDE of 0.105 mixes a df-2 t-design factor with a bootstrap SE and is roughly twice the coherent value; (iii) "the failure cannot be attributed to statistical power" overstates what the power numbers support. All are fixable in text plus two small artefact/string repairs.

---

## 2. Summary of the manuscript

The paper reanalyses two public rat DRG bulk RNA-seq series (GSE149648 Walker-256 tibial bone cancer; GSE160543 paclitaxel/oxaliplatin neuropathy; n = 4/group, same platform) and one mouse single-cell DRG series (GSE292522, B16F10, 4+4 animals). It makes two claims. **Biological:** bone cancer pain and CIPN drive qualitatively different DRG programmes — myeloid/type-I-interferon in Walker-256 at day 7, neuronal-injury in CIPN — demonstrated on canonical markers and 22 predefined modules with exhaustive label-permutation tests (C(8,4)² = 4,900 null pairs; two-sided floor 3/4,901). **Methodological (negative):** in-silico perturbation of the same data cannot rank analgesic targets. The perturbation layer correlates 0.991 with an additive no-network baseline, ranks an erythrocyte control above every analgesic target (|z| = 7.39), and is *anti-selective*: negative-control modules carry larger standardised effects than targets (median |z| 2.47 vs 1.07, two-sided Mann–Whitney p = 0.0019; controls over-represented among "significant" hits, binomial p = 0.0067), while remaining powered for large effects (MDE ≈ 0.105; 31/132 exact p < 0.05). The single-cell section shows the interferon programme is model-specific. The paper is internally consistent about this story across Abstract, §3.6, Limitations and Conclusions — the numbers agree everywhere — and carries an unusually strong reproducibility apparatus (self-checking permutation code, a 276-check artefact gate).

---

## 3. Major issues

### M1. The stated floor of the paired (shared-control) permutation test is wrong: it is 2/71, not 3/71 — and the paper's own table proves it

**Citation:** §2.7 (manuscript_final.md:87): "the enumeration contains only 70 distinct split pairs, and the floor is therefore 3/71 ≈ 0.042"; repeated in Limitations (line 177). Gate enforcement of the error: `98_check_manuscript_numbers.py:462` requires the literal string `"3/71"` to be present.

**Evidence:**
- The implementation `exact_two_sided_p` (`_limma.py:420-426`) computes p = (ge + 1)/(n + 1) where ge = #{|null| ≥ |obs|} over the enumerated null. For the **unpaired** 70×70 = 4,900 enumeration the floor claim **is** correct: the mirror pair (ā, b̄) is a member of the enumerated grid and yields exactly −Δ (cosine is odd in its second argument after mean-centering), so ge ≥ 2 and floor = 3/4,901 ≈ 6.12×10⁻⁴. I verified the mirror presence synthetically (ge ≥ 2 in 50/50 random trials; `A2_recompute.py` §9).
- For the **paired** 8-case exchange null (`22_etiology_comparison.py:301-324`) there is **no mirror**: the complement split gives cos(δ_A, c0a − mean(B)) − cos(δ_B, c0b − mean(A)), which is not −Δ of the true split (synthetic check: complement-split null −0.23136 vs −obs +0.33404; ge = 1 occurs in 4/300 synthetic trials). The true attainable floor is therefore (1+1)/71 = **2/71 = 0.028169**.
- The deposited artefact already contains values at the true floor: **6 of the 44 paired knockout rows** have p_permutation_exact = 0.028169 exactly — CTRL erythrocyte (both Bone_D14 vs Bone_D7 and PTX vs OXA), CTRL pan-neuronal (PTX vs OXA), Ca²⁺ α2δ (PTX vs OXA), Glial/immune (Bone_D14 vs Bone_D7), Injury/regeneration (PTX vs OXA).

**Required fix:** change "3/71 ≈ 0.042" to "2/71 ≈ 0.028" in §2.7 and Limitations, and state explicitly why the two-sided floor differs between the two enumeration types (mirror present in the 4,900 grid, absent from the 70-value exchange null). Update the gate's `in_text("3/71")` assertion (98:462) accordingly — as written, the reproducibility gate protects the error. No p-value in the results changes; the reported values are correct outputs of the (correct) implementation; only the Methods description is wrong.

### M2. The T5 power derivation is framework-incoherent, and "the failure cannot be attributed to statistical power" overstates the evidence

**Citation:** Abstract (line 21), §3.6 (lines 143), Limitations (line 177); derivation in `33_negative_result_quantification.py:144-152`.

**Evidence (recomputed):**
- d = 5.6535 = `TTestIndPower().solve_power(nobs1=2, alpha=0.05, power=0.8, ratio=1)` — the standardised MDE of a two-sample *t* with **2 animals per group** (the huge value is driven by the df = 2 critical value, t₀.₉₇₅,₂ = 4.303). This design factor is then multiplied by a **bootstrap SE** ((ci95_hi − ci95_lo)/3.92; median 0.01850 across the 132 dedup tests) which is not a df-2 quantity — it is a 2,000-resample SE of a complex statistic. Multiplying the two has no principled reading. The coherent normal-approximation MDE at 80% power is (1.960 + 0.842) × 0.01850 = 2.802 × SE = **0.0518**, i.e. roughly half the claimed 0.105. The error is *conservative* for the paper's own conclusion (a smaller true MDE makes the layer better powered than claimed), but "the minimum detectable Δ at 80% power is ≈0.105" is presented in three places as a design constant, and the phrase "at the unpaired n = 2 + 2 design" is nowhere derived — each cross-aetiology Δ consumes 12 (paired contrasts) to 16 (unpaired) animals; if a small-df correction is intended, its df must be stated and justified.
- The prose overreach: "Nor can the failure be attributed to statistical power" (§3.6:143). Recomputed: the **median target |Δ| is 0.0313** — 3.3× below the claimed MDE and 1.7× below the corrected MDE. Only **23/102** target tests exceed 0.1046 (controls: 15/30). So for the *typical* target, power IS low; what the data support is narrower and still damning: the test detects large effects (31/132 exact p < 0.05; the largest target Δ 0.2743 exceeds even the corrected MDE ~5.3×), and the hits it does produce are *control-enriched* (T2). The failure is anti-selectivity **plus** low per-test power for typical target effects — not "not power".
- Related consistency point: §3.6:143 says "each cross-aetiology Δ is a single unpaired estimate", but 44 of the 132 rows (Bone_D14 vs Bone_D7; PTX vs OXA) are shared-control contrasts evaluated under the 8-case exchange null, not an unpaired estimate. The T5 machinery treats all 132 alike.

**Required fix:** (a) either derive the n = 2+2/d = 5.653 choice explicitly with its df justification, or replace 0.105 with the coherent 2.802 × SE = 0.052 (and update "2.6-fold" → ~5.3-fold, or keep the conservative version but label it as such) in Abstract, §3.6, Limitations; (b) reword the power sentence to "powered for large effects; the typical target effect (median |Δ| = 0.031) lies below the detectable threshold, so the ranking failure is one of anti-selectivity, not merely of insensitivity"; (c) align Conclusions (line 181) with the same wording.

### M3. The randset q-values are tie-dominated at the resolution floor of the 500-draw null

**Citation:** §2.7 (line 89: "over 500 degree-matched draws"); §3.6 (lines 137-139); `22_etiology_comparison.py:77, 359-372`.

**Evidence (recomputed):** with N_RANDSET = 500 the two-sided empirical p has floor 2/501 = 0.003992, and **24 of the 132** deduplicated tests sit exactly at that floor (i.e. the observed Δ exceeded all 500 degree-matched draws). The BH q-values that define "randset-significant at q < 0.05" (32 tests) are computed over this heavily tied p vector, so the count and the q values (0.0220 for the whole floor block) are artefacts of the tie structure at this granularity — a rerun with a different N_RANDSET could move the boundary of the "32". The T2 conclusion itself (14/32 controls vs base rate 30/132, binomial p = 0.0067) is computed on the counts and is very likely robust, but the reported q layer is fragile.

**Required fix:** recompute the randset layer with N_RANDSET ≥ 5,000 (cheap — it is a batched matmul), or report the resolution floor next to every randset q and soften "significant at q < 0.05" to "all 500 degree-matched draws exceeded (empirical p at resolution 2/501)".

---

## 4. Minor issues

1. **Stale strings in `33_negative_result_quantification.py`** — the module docstring (line 15) says "MDE = 2.381 × median bootstrap SE" (2.3808 is the n=4+4 d used in T3, not T5's 5.653), and the T5 `interpretation` string hardcodes "**26**/132 exact permutation p<0.05" (line 169) while the computed field in the same JSON object is **31**. The deposited `negative_result_quantification.json` thus contradicts itself (field 31, prose 26). The manuscript correctly uses 31. Fix both strings; note the gate cannot catch prose inside JSON interpretation fields.
2. **Stale sentence in §2.7 (line 87):** "the split index is paired so that both sides draw the same control subset" describes the *deprecated* implementation (see the code's own comments at `22_etiology_comparison.py:283-291`); the current code fixes the shared control and exchanges the 8 case animals. Delete or rewrite the sentence.
3. **Latent degenerate-bootstrap trim bug** (`22_etiology_comparison.py:341-357`): `dB[b]` is assigned even for degenerate resamples, and the trim `dB = dB[:N_BOOT - ndeg]` drops the *last* ndeg entries rather than the degenerate ones. Currently inert — `n_degenerate_boot` = 0 across all 528 rows — but fix before it can fire.
4. **Exactness claim needs a conditional qualifier** (§2.7, line 87-89): for the reversal endpoint the diffused seed δ = P·seed and the network operator P are held fixed while only the pain axis is permuted (`22:252, 269-332`), so the p is exact *conditional on the per-contrast network operator*. The module-logFC endpoint (the primary endpoint) is exactly exchangeable with no such qualifier. One sentence in Methods suffices.
5. **SE from percentile-CI width** (`33:146`): (hi − lo)/3.92 assumes a symmetric bootstrap distribution; the deposited CIs are visibly asymmetric (glial/immune: 0.0712 below vs 0.0123 above the point estimate). Store the bootstrap SD directly in the CSV to remove the approximation.
6. **Gate takes T5's d from the JSON** (`98:360`) instead of recomputing `solve_power` at check time, unlike MWU/binomial which are recomputed (98:330-357). Recompute d in the gate too; otherwise a corrupted JSON field would self-verify.
7. **Rounding drift:** §3.6 (line 143) says "median 0.019" for the bootstrap SE; the artefact value is 0.01850 (rounds to 0.019 only half-up). Use 0.0185 consistently.
8. **"Minimum detectable Δ" wording** (Abstract, line 21): the MDE is quoted from the *median* SE; per-test MDE varies with the SE (IQR 0.0148–0.0250, i.e. ~0.041–0.070 under the z-approximation). Say "detectable at the median test" rather than implying a single design constant.

---

## 5. Numbers I independently recomputed

| # | Claim | Manuscript | Recomputed (from `etiology_bootstrap_differences.csv`, dedup KO per_contrast unless noted) | Match? |
|---|-------|-----------|------------------------------|--------|
| 1 | Universe size per scope | 264 tests (§3.6:137) | 264 rows per_contrast (528 total incl. pooled scope) | ✓ |
| 2 | Deduplicated hypotheses | 132 = 22 modules × 6 pairs (§3.6:137) | 132 (102 target + 30 control) | ✓ |
| 3 | KO/OE mirror | reversal(OE) = −reversal(KO), r = −1.000; p identical to 8 decimals, max\|Δ\| = 0 (§3.6:137) | corr(δ_KO, δ_OE) = −1.000000; max\|δ_KO + δ_OE\| = 0.000e+00; max\|p_KO − p_OE\| = 0.000e+00 | ✓ |
| 4 | n_null | 4,900 unpaired / 70 paired (§2.7:87) | n_null = 4900 (88 unpaired rows) and 70 (44 paired rows) | ✓ |
| 5 | Unpaired two-sided floor | 3/4,901 ≈ 6.1×10⁻⁴ (§2.7:87) | Correct as a floor: mirror pair always in grid (synthetic ge ≥ 2 in 50/50); smallest *observed* p = 7/4,901 = 0.0014283 | ✓ |
| 6 | **Paired two-sided floor** | **3/71 ≈ 0.042** (§2.7:87, §5:177) | **2/71 = 0.028169 — attained by 6 of 44 paired rows** (erythrocyte ×2, pan-neuronal, Ca²⁺ α2δ, glial/immune, injury/regeneration); complement split is not a mirror (synthetic + verified) | **✗ (M1)** |
| 7 | T1 median \|z\| target / control | 1.07 / 2.47 (§3.6:141, Abstract, Conclusions) | 1.0729 / 2.4702 | ✓ |
| 8 | T1 Mann–Whitney two-sided p | 0.0019 | 0.001914 | ✓ |
| 9 | T1 means \|z\| (JSON) | — | target 1.4635, control 2.5768 | ✓ (vs JSON) |
| 10 | T2 randset-significant | 32 total / 18 target / 14 control (§3.6:139) | 32 / 18 / 14 (q_randset < 0.05) | ✓ |
| 11 | T2 base rate | 30/132 (23%) | 30/132 = 0.2273 | ✓ |
| 12 | T2 binomial p (one-sided greater) | 0.0067 | 0.00671 | ✓ |
| 13 | 32 randset-sig all with CI excluding zero | "every one of them" (§3.6:139) | 32/32 | ✓ |
| 14 | 64/264 randset-sig raw | §3.6:137 | 64 (32 KO + 32 mirror OE) | ✓ |
| 15 | T4 target percentile | 55.7th (§3.6:141) | 55.7 | ✓ |
| 16 | Largest \|z\| overall | erythrocyte, \|z\| = 7.39, bone D14 vs D7 (§3.6:137) | CTRL erythrocyte, Bone_D14 vs Bone_D7, z = +7.3878, p_randset at floor 2/501 | ✓ |
| 17 | Erythrocyte PTX vs OXA | z = −5.78 (§3.6:137) | −5.7822 | ✓ |
| 18 | Pan-neuronal z | −3.74 (§3.6:137) | −3.7362 (Bone_D7 vs PTX) | ✓ |
| 19 | T5 median bootstrap SE | 0.019 (§3.6:143); 0.0185 (JSON) | 0.01850 (from CI width/3.92) | ✓ (rounding note, m7) |
| 20 | T5 d (n = 2+2) | 5.653 (JSON) | 5.6535 = solve_power(nobs1=2) | ✓ arithmetic, ✗ derivation (M2) |
| 21 | T5 MDE | ≈0.105 (Abstract:21, §3.6:143, §5:177) | 5.6535 × 0.01850 = 0.1046 (reproduces); coherent z-approx MDE = 2.802 × 0.01850 = **0.0518** | ✗ framework (M2) |
| 22 | Largest target Δ | 0.274, glial/immune, Bone_D7 vs PTX (§3.6:143) | 0.2743 (Glial/immune, Bone_D7 vs PTX) | ✓ |
| 23 | Ratio to MDE | 2.6× | 0.2743/0.1046 = 2.622 (5.30× vs corrected MDE) | ✓ arithmetic |
| 24 | Exact p < 0.05 count | 31 of 132 (Abstract, §3.6:143, §5:177) | 31 (25 target + 6 control; 10 paired + 21 unpaired). JSON prose says 26 — stale (m1) | ✓ |
| 25 | Glial/immune row | Δ = +0.274, CI 0.203–0.287, p = 0.0014, q = 0.165 (§3.6:137) | Δ = 0.2743, CI [0.2031, 0.2866], p = 7/4,901 = 0.0014283, q = 0.1654 | ✓ |
| 26 | Min q (exact permutation) | 0.165 (Abstract:21) | 0.1654; 0 tests q < 0.05 | ✓ |
| 27 | BH family level | 132 per scope, KO dedup (§3.6:137; `22:687-695`) | BH applied within scope over the 132 dedup KO p's, separately for permutation and randset families | ✓ (sound) |
| 28 | 8-case exchange null validity | §2.7:87 | Genuinely exchangeable under H0 (labels of the 8 case animals; shared control fixed); true-assignment self-checks run and are enforced (`22:272-282, 311-324`); no index/orientation bug found | ✓ |
| 29 | Bootstrap covers observed | gate-checked | 132/132 (132 dedup); among 31 exact-p<0.05, 30 have CI excluding zero (1 discrepancy — different inference targets, not a bug) | ✓ |
| 30 | ρ(network, additive) = 0.991 | "across all 132 module × contrast × perturbation conditions" (§3.4:125) | head_to_head_baselines.json: 0.99115, n_conditions = 132 = 22 modules × 3 discovery contrasts × 2 perturbations (discovery space, not cross-aetiology) | ✓ (phrasing could name the space) |
| 31 | Cross-study ρ = 0.929 | §3.7:151 | sensitivity_summary.json `spearman_network_vs_additive` = 0.929 (gate-checked) | ✓ |
| 32 | Randset resolution | 500 draws (§2.7:89) | empirical p floor 2/501 = 0.003992; **24/132 tests at the floor** → q tie-dominated | ✗ granularity (M3) |
| 33 | Degenerate bootstrap resamples | — | max n_degenerate_boot = 0 across 528 rows (trim bug latent only) | ✓ (m3) |
| 34 | Targets exceeding MDE | implied by "adequately powered" (§3.6:143) | median target \|Δ\| = 0.0313; 23/102 > 0.1046 (controls 15/30); 41/102 > 0.0518 | ✗ wording (M2b) |

---

## 6. What would make this publishable (concrete, ordered)

1. **Fix the paired-floor claim** (M1): "2/71 ≈ 0.028" in §2.7 and Limitations, with the one-line explanation of why the two-sided floor differs between the 4,900-grid and the 70-value exchange null; update the gate's `"3/71"` string assertion. No results change.
2. **Re-derive or replace the MDE** (M2a): either justify d = 5.653 with an explicit df argument, or report the coherent 2.802 × SE ≈ 0.052 (making the largest target effect ~5.3× MDE). Propagate consistently to Abstract, §3.6, Limitations, Conclusions.
3. **Recalibrate the power sentence** (M2b): state that the median target |Δ| (0.031) lies below the MDE and that the argument against "underpowering" rests on the control-enrichment of the significant set (T2), not on blanket adequacy.
4. **Re-run the randset layer with N_RANDSET ≥ 5,000** (M3) so the q_randset values are not tie artefacts; the T1/T2 conclusions should be re-confirmed on the finer grid.
5. **Artefact hygiene** (m1): remove the stale "2.381" and "26/132" strings from script 33 so the deposited JSON is internally consistent.
6. **Remove the stale §2.7 sentence** about paired split indices (m2) and add the conditional-exactness qualifier for the reversal endpoint (m4).
7. **Gate hardening** (m6): recompute T5's d at check time; add a check that no paired p falls below the true floor 2/71 and no unpaired p below 3/4,901.
8. Optionally (m3, m5): fix the degenerate-resample trim and store bootstrap SDs directly in the CSV.

All of the above are text/code repairs; none requires changing any scientific conclusion. The negative-result layer as quantified (T1 anti-selectivity p = 0.0019; T2 control over-representation p = 0.0067; no reversal q < 0.05 with min q = 0.165) is exactly reproducible from the deposited artefacts and is, in this reviewer's judgement, the strongest and most honest part of the paper.
