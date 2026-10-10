# R4-A2 Review — Biostatistics, multiple-testing correction and statistical power

Reviewer role: independent biostatistician (multiplicity / power / permutation inference).
Manuscript: `docs/manuscript/manuscript_final.md` (read as a first submission).
Independence: only `_PANEL_BRIEF.md` read under `docs/review/`; no earlier manuscript
version, no `.workbuddy` file, no other reviewer output.

---

## 1. Verdict

**Minor Revision.** The manuscript is **publishable** at a soundness-oriented SCIE
journal after minor revision.

Every statistical quantity I was asked to re-derive reproduces the manuscript to the
quoted precision, including the three exact-permutation p-value floors, the KO/OE mirror
identity, the deduplicated multiplicity denominator, T1/T2/T5, the MDE chain, Cohen's d,
and the BH q values in Tables 2/3/4. The reproducibility gate runs clean
(`287 checks, 0 failed`). The negative-result argument is internally coherent across the
Abstract, §3.6, §5 and §6. Nothing I found blocks acceptance. The issues below are real
but non-blocking: they concern the *interpretation* of correlated test statistics and
wording, not the correctness of any number.

---

## 2. Summary of the manuscript

The paper is a single-author computational reanalysis of public bulk rat DRG RNA-seq
comparing two cancer-pain aetiologies. GSE149648 (Walker-256 tibial bone cancer, sham /
d1 / d7 / d14, n = 4) and GSE160543 (paclitaxel / oxaliplatin / vehicle, n = 4) share
species, tissue and platform (GPL24782), so a cross-aetiology comparison is interpretable.
Differential expression uses a moderated t with a muscle-contamination covariate (§2.3–2.4).

Two claims are kept apart (lines 15–25, 37):

1. **Biological (positive).** The two aetiologies drive qualitatively different DRG
   programmes. Bone cancer pain at its d7 peak is a type-I interferon / myeloid event
   (interferon module t = +24.9, q = 6.8e-5; Cd68 Δ = +3.42, Aif1 Δ = +3.08, galanin
   Δ = −2.59, all q ≈ 0.012 by exact permutation), whereas paclitaxel neuropathy is a
   neuronal-injury programme. The divergence is confirmed by exhaustive label permutation
   (§2.7, §3.5, Tables 2–3).

2. **Methodological (negative).** The same data do **not** support ranking analgesic
   targets by in-silico perturbation. The layer is calibrated against five negative-control
   modules, degree-matched nulls, an oracle anchor, an additive baseline and a ridge
   baseline (§2.6). It correlates 0.991 with the additive model, ranks an erythrocyte
   control above every target (|z| = 7.39), assigns controls larger standardised effects
   than targets (median |z| 2.47 vs 1.07, MWU two-sided p = 0.0019), and clears the median
   target effect only weakly (§3.6). An independent mouse single-cell dataset (GSE292522,
   B16F10, n = 4 + 4) reproduces the neuronal-injury programme but not the interferon
   programme, so the interferon signature is model-specific (§3.8).

The statistics I scrutinised: (i) the exact permutation scheme and its three declared
p-value floors; (ii) the KO/OE mirror deduplication; (iii) T1 (Mann–Whitney on |z|), T2
(binomial over-representation of controls); (iv) T5 (MDE from the median bootstrap SE);
(v) Cohen's d and MDE-d in §3.8; (vi) BH q in Tables 2/3/4.

---

## 3. Major issues

None. I could not construct a statistical defect that would stop me recommending
acceptance. The issues that follow are Minor.

---

## 4. Minor issues

### M1. T1/T2 treat 132 correlated test statistics as independent samples
*Citation:* §3.6 lines 139, 141, 185; `scripts/33_negative_result_quantification.py:40–89`.

T1 is a Mann–Whitney U over the 132 deduplicated tests (102 target vs 30 control |z|), and
T2 is a binomial over the 32 randset-significant tests. Both tests assume the analysed
units are independent, but the 132 units are **22 modules × 6 aetiology pairs** and the
six pairs reuse the same four contrasts, so module |z| values are strongly dependent
across pairs. The nominal p (0.0019, 0.0067) is therefore not an exact error rate.

*Required fix (text only):* state that the MWU/binomial p-values are descriptive summaries
under an independence assumption that does not strictly hold, and note that the
anti-selectivity conclusion does not depend on them (it already survives dropping the
interferon control, p = 0.019, and is visible directly in the median |z| gap and in
Figure 3). No renumbering needed.

### M2. The "t-based reference" MDE uses an unexplained n = 2 per group
*Citation:* §3.6 line 143, §5 line 179; `scripts/33_negative_result_quantification.py:153–159`.

The normal-approximation MDE (2.802 × 0.0185 = 0.052) is well defined. The "more
conservative t-based reference ≈ 0.10" is computed with `solve_power(nobs1=2, ...)` →
d = 5.653. Nothing in the design has 2 animals per group; §5's phrase "at these degrees of
freedom" is vague. The arbitrariness of n = 2 makes the conservative bound look chosen
rather than derived.

*Required fix:* give the df/n behind d = 5.653 explicitly (or replace it with a clearly
labelled "illustrative" bound). The conclusion is unchanged either way.

### M3. Bootstrap SE is back-derived from a percentile interval under a normality assumption that one row visibly violates
*Citation:* `scripts/33_negative_result_quantification.py:156`;
`results/tables/etiology_bootstrap_differences.csv`.

T5 defines SE = (ci95_hi − ci95_lo)/(2 × 1.96). For the glial/immune (Bone_D7 vs PTX) row
the interval is [0.2031, 0.2866] around Δ = 0.2743 — the point estimate sits near the
*upper* bound, i.e. the bootstrap distribution is markedly left-skewed, so the
normal-theory SE is a poor local summary. At the median across 132 tests the effect is
small and the MDE conclusion is unaffected, but the method should be named honestly.

*Required fix:* one sentence stating the SE is a normal-approximation summary of the
bootstrap interval and that individual CIs may be asymmetric (the glial/immune CI is the
clearest case). Cosmetic.

### M4. "design-independent" is stronger than the evidence supports
*Citation:* Abstract line 23 ("its failure is design-independent ...").

The failure is demonstrated for two network scopes (`per_contrast`, `pooled_by_study`)
and one hyper-parameter grid on two datasets. "Design-independent" reads as a universal
claim; "robust across the network scopes and hyper-parameters tested" is what the
artefacts support.

*Required fix:* soften to the tested scope.

### M5. Misleading key name in an artefact (`myeloid_MWU_p_less`)
*Citation:* `results/tables/sc_GSE292522_summary.json`, key `myeloid_MWU_p_less`;
`scripts/98_check_manuscript_numbers.py:280`.

The stored value 0.057143 is the **one-sided *greater*** p (control myeloid fraction >
B16F10): I recomputed `mannwhitneyu(Con, B16F10, alternative='greater') = 0.057143`,
whereas the actual `alternative='less'` gives 0.971429. The manuscript does not mis-state
the number (it reports "one-sided p = 0.057" without direction), but the key name is wrong
and would mislead a future reader of the artefact.

*Required fix:* rename the key (e.g. `myeloid_MWU_p_one_sided_greater`). Artefact-only.

### M6. Marker set in Table 2 is not explicitly declared pre-specified
*Citation:* §2.5 (modules declared a priori, line 67) vs Table 2 (lines 269–279).

The 17 target modules are declared a priori, but the canonical markers of Table 2 are not
stated as pre-selected. With 19 markers × 6 pairs = 102 tests, an undeclared marker set
invites a garden-of-forking-paths reading. The BH correction (across 102 tests) and the
exact permutation floor mitigate this, and the markers are conventional injury/immune
genes, so this is a reporting point only.

*Required fix:* state in §2.5 or §2.7 that the marker panel was fixed before analysis.

---

## 5. Numbers you independently recomputed

All recomputations were done from `results/tables/*.csv`, `negative_result_quantification.json`
and (for the endpoint floors) a direct reconstruction from `data/raw/` with the managed
Python. I did **not** trust the text.

### 5a. Exact permutation p — the three floors, and whether endpoint dependence is real

The single most important thing the brief asked me to adjudicate: **the endpoint
dependence is real and structural, not a tie-tolerance or one-sidedness artefact.**

Code path: `exact_two_sided_p` (`scripts/_limma.py:420–426`) uses
`ge = Σ(|null| ≥ |obs| − 1e-12)`, `p = (ge+1)/(n+1)`. The tolerance is 1e-12; the test is
two-sided for **both** endpoints. So neither tolerance nor convention can produce the
difference.

I reconstructed the paired nulls directly (module *Glial/immune*, Bone_D7 vs Bone_D14,
shared sham; per-contrast diffusion operator, α = 0.55, k = 12, 3,038 nodes):

| endpoint | observed Δ | complement-split value | obs + complement | ge | p |
|---|---|---|---|---|---|
| logFC (paired) | +0.750362 | −0.750362 | **0.00e+00** | 4 | 5/71 |
| reversal (paired) | +0.026717 | −0.019268 | **+0.007449** | 1 | **2/71** |

Mechanism: for the logFC endpoint the complement split is the *exact algebraic negation*
of the observed value (`Δ = mean(A_case) − mean(B_case)`; swapping the 4-vs-4 assignment
negates it), so both the observation and its mirror are in the null → attainable floor
3/71. For the reversal endpoint the two axes at the complement split are
`c0 − mean(B_case)` and `c0 − mean(A_case)`, which are *not* mutual negations; the mirror
therefore does not reproduce `−obs`, and only the observed split is as extreme → attainable
floor 2/71. `scripts/22_etiology_comparison.py:294–327` implements exactly this, and the
comment at lines 283–291 states the reasoning correctly.

Artefact check of the attained minima:

| claim | manuscript | recomputed (artefact) | match |
|---|---|---|---|
| unpaired floor = 3/4,901 | 3/4,901 = 6.1e-4 | min unpaired p, module & marker = 0.0006121 = **3/4,901** | ✓ |
| paired logFC floor = 3/71 | 3/71 ≈ 0.042 | min paired p, module & marker = 0.0422535 = **3/71** | ✓ |
| paired reversal floor = 2/71 | 2/71 ≈ 0.028 | min paired reversal p = 0.0281690 = **2/71** | ✓ |
| strongest reversal Δ = glial/immune | p = 0.0014 | p = 0.001428 = 7/4,901 | ✓ |

(Note: the unpaired reversal endpoint also has attainable floor 3/4,901 by the same
antisymmetry; its *observed* minimum is 7/4,901 — the strongest reversal. The manuscript's
"3/4,901" statement is scoped to the logFC endpoint and is correct.)

### 5b. Dedup: KO/OE are exact mirrors

| claim | manuscript | recomputed | match |
|---|---|---|---|
| reversal(KO) = −reversal(OE) | r = −1.000 | corr = −1.000000; max&#124;rev_A(KO)+rev_A(OE)&#124; = 0.0; max&#124;rev_B(KO)+rev_B(OE)&#124; = 0.0 | ✓ |
| exact-permutation p identical | max&#124;Δ&#124; = 0 (8 dp) | max&#124;p_KO − p_OE&#124; = **0.0**; max&#124;q_KO − q_OE&#124; = 0.0 | ✓ |
| independent tests per scope | 132 | 132 (KO), target 102 / control 30 | ✓ |
| raw tests | 264 | 264 (22 modules × 6 pairs × 2) | ✓ |

### 5c. T1, T2, T5 / MDE

| claim | manuscript | recomputed | match |
|---|---|---|---|
| T1 n | 102 target / 30 control | 102 / 30 | ✓ |
| T1 median &#124;z&#124; control | 2.47 | 2.4702 | ✓ |
| T1 median &#124;z&#124; target | 1.07 | 1.0729 | ✓ |
| T1 Mann–Whitney two-sided p | 0.0019 | 0.0019139 | ✓ |
| T1 excl. interferon, p | 0.019 | 0.0186968 (median 2.33 vs 1.07) | ✓ |
| T2 randset-significant | 32/132 (14 control) | 32 (14 control, 18 target) | ✓ |
| T2 base rate | 30/132 (23%) | 0.22727 | ✓ |
| T2 binomial p (greater) | 0.0067 | 0.0067079 | ✓ |
| T5 median bootstrap SE | 0.0185 | 0.018500 | ✓ |
| T5 MDE = 2.802 × SE | 0.052 | 2.801585 × 0.0185 = 0.05183 | ✓ |
| T5 t-based reference | ≈0.10 | 5.653 × 0.0185 = 0.10458 | ✓ |
| T5 largest target Δ | 0.274 | 0.27429 (Glial/immune, Bone_D7 vs PTX) | ✓ |
| T5 median target &#124;Δ&#124; | 0.031 | 0.03130 | ✓ |
| T5 targets > MDE (normal) | 41/102 | 41/102 | ✓ |
| T5 targets > MDE (t-based) | 23/102 | 23/102 | ✓ |
| T5 exact p < 0.05 | 31/132 | 31/132 | ✓ |
| T5 min q | 0.165 | 0.16545 | ✓ |
| largest &#124;z vs randset&#124; = control | 7.39 | 7.388 (CTRL erythrocyte, Bone_D14 vs Bone_D7) | ✓ |

### 5d. Cohen's d / Table 5 (single-cell myeloid)

| claim | manuscript | recomputed | match |
|---|---|---|---|
| myeloid fraction control | 0.12% | 0.001170 (0.1170%) | ✓ |
| myeloid fraction B16F10 | 0.05% | 0.000517 (0.0517%) | ✓ |
| drop | 55% | 55.8% | ✓ |
| Cohen's d (n−2 denom = 6) | 1.27 | 1.26700 | ✓ |
| MDE d (n = 4 vs 4, 80%) | 2.38 | 2.38080 | ✓ |
| myeloid MWU two-sided | 0.114 | 0.114286 | ✓ |
| myeloid one-sided | 0.057 | 0.057143 (greater) | ✓ |
| pseudobulk up/down q<0.05 | 1,517 / 703 | 1,517 / 703 | ✓ |
| Sprr1a | +1.55 (t = 54.5) | +1.554 (t = 54.49) | ✓ |
| Ifit2 | −0.048, q = 0.042 | −0.048, q = 0.0417 | ✓ |
| Oas1c q | 0.014 | 0.0136 | ✓ |
| Gfap | +0.004, q = 0.006 | +0.004, q = 0.0061 | ✓ |
| concordance ρ (D7/D1/D14) | 0.098 / 0.112 / 0.030 | 0.0985 / 0.1115 / 0.0299 | ✓ |

### 5e. BH q spot-checks (≥ 3 required; I checked 20)

Recomputed BH independently and compared to the deposited `q_*` column:

| table | rows checked | result |
|---|---|---|
| Table 2 (markers, BH over 102) | Cd68, Aif1, Gal (0.0125); Ccl2 (0.0173); Il1b (0.0182); Calca (0.0272); Scn10a (0.0260) | all match to 4 dp |
| Table 3 (modules, BH over 132) | Glial/immune, Cytokine, TLR, CTRL erythrocyte (0.0073); Neuropeptides (0.0145) | all match to 4 dp |
| Table 4 (discovery D7 KO, n = 22) | interferon, glial/immune, pan-neuronal, K+/HCN (0.011); Na+, myelin, erythrocyte, neuropeptides (0.0195) | all match to 4 dp |

### 5f. Gate audit (`scripts/98_check_manuscript_numbers.py`)

Ran it myself: **`287 checks, 0 failed`** (exit 0), including `68 table cells regenerated`.

Runtime-derived vs hardcoded (instrumented the module and counted executed assertions):
- **≈ 80 of 287 are genuinely runtime-derived:** ~68 parse the manuscript's Table 2/3
  cells at check time and compare them to `etiology_marker_cross.csv` /
  `etiology_module_logfc_cross.csv` (`98_...:447–476`); ~12 recompute a statistic from the
  bootstrap CSV at run time and compare it to `negative_result_quantification.json`
  (MWU, binomial, MDE; `98_...:330–388`).
- **≈ 207 of 287 have hardcoded expectations:** ~164 compare an artefact value against a
  bare numeric literal or a hardcoded dict of manuscript numbers (e.g. `want_t`,
  `want_mod`, `want_mk`, lines 132–220), ~6 are string-literal `chk`s, and ~37 are
  `in_text` substring-presence guards (lines 390–399, 481–488).

Interpretation: the gate is a **manuscript↔artefact drift detector**, not an independent
statistical re-derivation. It cannot catch an error that is present identically in the
manuscript and the artefact it was copied from; conversely, its three most valuable
features — the manuscript-parsed Table 2/3 cell loop, the runtime MWU/binomial/MDE
recomputation, and the `in_text` guards — are exactly the three I recommend it keep. My
independent recomputation (5a–5e) shows the underlying numbers are nonetheless correct, so
the gate's hardcoded literals are not masking any stale value.

---

## 6. Publishability assessment

**Ready after Minor Revision.** The statistical core is verified end-to-end: the exact
permutation scheme is mathematically sound, its endpoint-dependent floors are real and
correctly derived in the code, the multiplicity denominator (132, not 264) is justified by
an identity I confirmed to machine precision, and every headline number in §3.5–3.8 and
Tables 2/3/4/5 reproduces. The negative-result narrative is consistent between the
Abstract, §3.6, the Limitations and the Conclusions, and does not overstate the power
claim (it explicitly separates an anti-selectivity failure from a median-effect power
floor).

Ordered minimum to acceptance:

1. **M1** — one sentence in §3.6 acknowledging that the T1/T2 p-values assume
   independence of the 132 correlated tests; keep the conclusion.
2. **M2** — state the n/df behind the t-based MDE reference (or relabel it
   "illustrative").
3. **M3** — one sentence that the T5 SE is a normal-approximation summary of a possibly
   skewed bootstrap interval.
4. **M4** — soften "design-independent" to the tested scopes.
5. **M5** — rename `myeloid_MWU_p_less` in the artefact and the gate (no manuscript
   number changes).
6. **M6** — record that the Table 2 marker panel was fixed before analysis.

None of these requires re-running the analysis or changes any reported number. Items
1–4 are wording; item 5 is a non-breaking artefact rename; item 6 is one clause. With
these, I would support acceptance.
