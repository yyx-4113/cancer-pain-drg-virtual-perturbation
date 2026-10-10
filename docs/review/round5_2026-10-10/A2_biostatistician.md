# R5-A2 Review — biostatistics / multiple testing / power

Reviewer role: independent biostatistician (multiple testing, exact inference, power).
Posture: first submission, read-only, everything re-derived from source artefacts.
Manuscript: `docs/manuscript/manuscript_final.md`.

## 1. Verdict

**Minor Revision** (acceptance-level; publishable after minor statistical caveats).

No blocking biostatistical error. Every number I was asked to re-derive reproduces from
the deposited artefacts to the precision the manuscript quotes. The gate
`scripts/98_check_manuscript_numbers.py` runs to **291 checks, 0 failed**. The endpoint-
dependent permutation floor is real, the KO/OE mirror deduplication is exact, and the
negative-layer power/quotation (T1–T5) is internally consistent. The remaining items are
caveats about test independence and the nature of the gate, not corrections.

## 2. Summary of the manuscript

Single-author reanalysis of two bulk rat DRG RNA-seq studies (GSE149648 Walker-256 tibia,
GSE160543 paclitaxel/oxaliplatin, same platform GPL24782, n=4/group) plus an independent
mouse scRNA-seq validation (GSE292522, B16F10, 4+4 animals, 55,967 cells). Two claims:
(1) bone cancer pain and CIPN drive qualitatively different DRG programmes, testable on
canonical markers by an exhaustive label permutation (exact p); (2) in-silico perturbation
of 17 analgesic-target modules cannot rank targets—it fails negative controls, is
anti-selective, and is only powered for large effects. Statistics: moderated-t (limma-trend,
muscle covariate), BH q, exact label permutation over C(8,4)^2 = 4,900 split pairs,
degree-matched random-set nulls, within-study bootstrap CIs, and a binomial/MWU
anti-selectivity analysis.

## 3. Major issues

None. I looked specifically for a wrong denominator, an inflated p, or a power claim that
the effect sizes do not support, and found none. The three places a biostatistical review
usually breaks a paper of this kind all hold up:

- **Endpoint-dependent floor is genuine, not rhetorical.** The floor follows from the
  *antisymmetry* of each endpoint, not from a coding artefact. Unpaired logFC Δ over the
  4,900 split pairs is antisymmetric (complement split = exact negation), so |null| ≥ |obs|
  for at least the observed pair and its mirror → minimum attainable p = 3/4,901; the
  attained minimum equals 3/4,901 exactly. Paired logFC Δ is likewise antisymmetric under
  the exchange-8-case null (C=70) → floor 3/71, attained exactly. Paired reversal Δ is *not*
  antisymmetric (the projection at the complement split is not the exact negation), so only
  the observed split itself is as extreme → floor 2/71, attained exactly (`min paired KO
  reversal p × 71 = 2.00`). This is a correct and non-trivial distinction.

- **KO/OE deduplication is exact.** `reversal(OE) = −reversal(KO)` for all 132
  (module, pair) cells (corr(rev_KO, −rev_OE) = 1.000, max|rev_KO + rev_OE| = 0,
  max|Δp| = 0.0), so 132 independent tests per network scope is the right denominator and
  BH is applied to the deduplicated knockout set. I reproduced the small formal
  inconsistency that the *randset* p is not bit-identical across mirrors (different
  degree-matched draws) but is treated as one hypothesis—this is the conservative choice
  and the manuscript states it.

- **The power arithmetic is right.** MDE = 2.802 × median bootstrap SE = 2.802 × 0.01850 =
  0.0518 log-units is the correct two-sided 80%-power normal reference (z0.975 + z0.80 =
  1.96 + 0.8416). The t-based reference (0.1046) is labelled conservative; the largest
  target Δ (0.274) clears both, the median target |Δ| (0.031) misses both, and 41/102
  (23/102) targets clear the normal (t) threshold with 31/132 exact p<0.05. The threshold
  logic and the "resolves large effects, not typical ones" wording match the numbers.

## 4. Minor issues

1. **Independence assumption behind T1 (Mann–Whitney) is not caveated, though T2's is.**
   The 132 deduplicated tests are *not* independent: they share 4 contrasts and 22 modules,
   so both the MWU (T1, §3.6, p = 0.0019) and the binomial (T2, p = 0.0067) sit on
   correlated units. The manuscript explicitly flags the shared-contrast problem for T2
   ("the tests share contrasts, so the p is approximate") but quotes T1's p=0.0019 without
   the same caveat. Required fix: one clause extending the same caveat to T1 (or, better, a
   contrast-blocked / module-stratified re-permutation of the anti-selectivity statistic).
   Not blocking: the anti-selectivity conclusion rests on the effect size (median |z| 2.47
   vs 1.07, and 2.33 vs 1.07 with the interferon control removed), which a blocked null would
   not overturn.

2. **The 98-gate is a drift detector, not an independent re-derivation.** I instrumented the
   checker: of the 291 executed assertions, **≈86 are runtime-derived** (the expected value
   is recomputed at run time from the CSV and compared to another artefact/JSON—the whole
   T1/T2/T5 block, plus the 68 Table 2/3 cells regenerated from `etiology_marker_cross.csv`
   and `etiology_module_logfc_cross.csv`, plus 8 derived `in_text` checks) and **≈205 assert
   an artefact value against a hardcoded literal** (e.g. `chk("network genes", s21[...], 3042)`,
   the `want_t`/`want_mod`/`want_mk`/`want_lf` dictionaries, the 29 literal `in_text`
   strings). This is a strong guard against manuscript↔artefact drift but cannot detect an
   error frozen *into* an artefact; e.g. the "0/4,600/235" DE counts, the module means, and
   §3.1–3.4 are trusted, not re-derived, by the gate. This is disclosed in the script's own
   header ("Hardcoding a stale number here would 'protect' the error"), so I record it as a
   reproducibility limitation, not a defect.

3. **§3.6 percentile wording.** "Targets sit just below the median … (44.7th percentile)"
   is `scipy.stats.percentileofscore(allz.absz, target_median_absz, kind="rank")` over the
   132 deduplicated |z| values (script 33 T4), i.e. the rank-percentile of the target
   *median* value, not the fraction of the 102 targets below the pooled median. The value is
   correct and consistent with the script; the sentence could be misread as a target-count
   share. Cosmetic rewording only.

4. **Table 4 shows 8 of 22 modules without saying it is the top-8 by p.** The caption reads
   "Calibration of the in-silico perturbation layer (day 7, discovery)" but the table is the
   eight lowest-p modules; the other 14 (q up to 0.759) are in the CSV. Add "top 8 by
   p (of 22)" to the caption to avoid implying completeness.

## 5. Numbers independently recomputed

All comparisons were made with the manuscript text and the deposited CSV/JSON; `p` is the
exact label-permutation p unless stated.

| Claim | Manuscript | Recomputed | Match? |
|---|---|---|---|
| Unpaired logFC Δ floor | 3/4,901 | min attained = 3/4,901 | yes |
| Paired logFC Δ floor | 3/71 | min attained = 3/71 | yes |
| Paired reversal Δ floor | 2/71 | min attained = 2/71 | yes |
| Unpaired reversal (theoretical floor 3/4,901) | n/a | min attained 7/4,901 (≥ floor) | consistent |
| KO/OE mirror | r = −1.000, max\|Δp\| = 0 | corr = 1.000, max\|rev_KO+rev_OE\| = 0, max\|Δp\| = 0.0 | yes |
| Independent tests per scope | 132 | 132 (KO dedup, per_contrast) | yes |
| T1 median \|z\| control vs target | 2.47 / 1.07 | 2.4702 / 1.0729 | yes |
| T1 MWU two-sided p | 0.0019 | 0.001914 | yes |
| T1 excl. interferon control | 2.33 / 1.07, p=0.019 | 2.3327 / 1.0729, p=0.0187 | yes |
| T2 controls among randset-sig | 14/32 vs base 30/132, p=0.0067 | 14/32, base 0.2273, binom p=0.00671 | yes |
| Randset-sig total (dedup) | 32 (18 target / 14 control) | 32 (18 / 14) | yes |
| T5 median bootstrap SE | 0.0185 (≈0.019) | 0.01850 | yes |
| T5 MDE (normal) | ≈0.052 (2.802×SE) | 0.0518 | yes |
| T5 MDE (t-based) | ≈0.10 | 0.1046 | yes |
| T5 max target Δ | 0.274 | 0.2743 (Glial/immune, D7 vs PTX) | yes |
| T5 median target \|Δ\| | 0.031 | 0.0313 | yes |
| T5 targets above MDE | 41 (23 conservative) | 41 (23) | yes |
| T5 exact p<0.05 | 31/132 | 31/132 | yes |
| §3.6 target percentile | 44.7th | percentileofscore kind="rank" = 44.7 | yes |
| Smallest permutation q | 0.165 | min over 132 (BH re-run) = 0.16545 | yes |
| §3.8 / Table 5 Cohen's d | 1.27 (n−2 denom) | 1.267 (pooled SD, df=6) | yes |
| §3.8 / Table 5 MDE d | 2.38 | 2.38 (n=4, 80% power) | yes |
| Myeloid fraction | 0.12% vs 0.05%; U=2/14, p=0.114 (0.057) | 0.00117 vs 0.00052; U=2/14, p=0.1143 (0.0571) | yes |
| Erythrocyte bin imbalance | 2 control / 4 tumour | `celltype_scores`: Erythrocyte present in 6 animals {Con:2, B16F10:4} | yes |
| Table 2 q (spot ≥3) | 0.012 / 0.017 / 0.018 / 0.026 / 0.027 | 0.012487 / 0.017343 / 0.018211 / 0.026015 / 0.027216 | yes |
| Table 3 q (spot ≥3) | 0.0073 / 0.0145 / 0.0162 / 0.0310 | 0.007345 / 0.014503 / 0.016160 / 0.030973 | yes |
| Table 4 q (spot ≥3) | 0.011 / 0.011 / 0.020 / 0.020 | BH re-run over 22 D7 modules = 0.010989 / 0.010989 / 0.019536 | yes |
| sc pseudobulk DE | 1,517 up / 703 down | 1,517 / 703 | yes |
| sc Sprr1a / Ifit2 | +1.55 (t=54.5) / −0.048 (q=0.042) | 1.5539 (t=54.49) / −0.0476 (q=0.0417) | yes |
| sc Oas1c down q | 0.014 | 0.0136 | yes |
| Normalisation scale | log2(counts-per-10,000 + 1) | `32_percelltype_sc.py:154` `cp10k =…*1e4+1; log2` | yes |
| Gate | 291 checks, 0 failed | 291 checks, 0 failed | yes |
| Gate expectations runtime-derived | — | ≈86 executed (28%) derived; ≈205 (72%) hardcoded literals | reported |

## 6. Publishability assessment

**Publishable after minor revision** at a soundness-oriented SCIE venue. The biostatistical
core—exact permutation inference with a correctly derived endpoint-dependent floor, exact
KO/OE mirror deduplication, a properly referenced MDE, and an anti-selectivity result that
survives removal of its strongest control—is sound and reproducible from the artefacts. There
is no blocking issue.

Ordered minimum to clear Minor Revision:
1. Add the shared-contrast / non-independence caveat to T1 (§3.6) as it already appears for
   T2, so p=0.0019 is not read as a calibrated independent-trial p.
2. Reword the §3.6 percentile sentence to make clear it is the rank-percentile of the target
   *median* |z| (not the share of targets below the pooled median).
3. State in Table 4's caption that it is the top-8-by-p subset of 22 modules.
4. (Optional, reproducibility) Note in §2.9 that the 98-checker's expectations are
   predominantly hardcoded literals protecting against manuscript/artefact drift, with the
   negative-result block and Table 2/3 regenerated at run time.
