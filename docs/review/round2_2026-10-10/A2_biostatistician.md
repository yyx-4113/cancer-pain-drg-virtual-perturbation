# A2 — Biostatistics (design) review

**Manuscript:** `docs/manuscript/manuscript_final.md`
**Round:** 2 (independent panel)
**Reviewer role:** statistical validity and inferential honesty
**Date:** 2026-10-10

---

## Scope note

I was asked to judge p-values, tests, power statements and multiplicity
corrections. I read `manuscript_final.md` once, then
`scripts/22_etiology_comparison.py`, `scripts/_limma.py`,
`scripts/33_negative_result_quantification.py`, `scripts/32_percelltype_sc.py`,
`scripts/20_multistage_pipeline.py`, `scripts/23_sensitivity.py`, and every
table named in the shared brief. I wrote six independent check scripts
(`review_A2/A2_check1..6`) that re-parse the raw GEO member files and rebuild
the statistics without importing the project's own primitives.

**Bottom line.** The exact-permutation machinery is real and the headline
marker results are reproducible to six decimal places from raw counts. But
four things in §3.6 and §3.8 are statistically wrong, not merely imprecise:
the two **p-value floors**, the **shared-control ("paired") null**, the
**power/MDE calculation**, and the **direction of the T1 test**. The power
error is the most serious: it is the single new statistical claim of this
round, it is built on a test that does not exist for this design, and
correcting it **reverses** the conclusion the paper draws from it. Three
smaller items (the Table 1 `q` column, the myeloid `d`, and the `U = 2/14`
denominator) are mislabelled rather than wrong.

---

# PART 1 — FINDINGS

Each finding follows the mandated four-part contract.

---

## F1. The two p-value floors are wrong by a factor of three

**【Problem】** The stated minimum attainable p-values — 1/4,901 for the four
unpaired contrasts and 1/71 for the two shared-control contrasts — are not
attainable. The true floors are 3/4,901 and 3/71.

**【Evidence】** The enumeration is exhaustive, so the observed labelling is
itself a member of the null, and the 70 splits come in 35 complement pairs
whose axes are exact negations. I verified this directly:
`all_splits(4,4)` gives 70 rows, each with exactly 4 `True`, all distinct;
the complement map is an involution with **0 fixed points**; and
`max |logFC_r + logFC_complement(r)| = 0.000e+00` across the transcriptome
(`review_A2/A2_check1_permutation.py`, section (a)). Therefore for any
observed Δ both `+Δ` and `−Δ` are always present in the null, so
`#{null : |null| ≥ |obs|} ≥ 2`, and the repo's rule
`p = (ge+1)/(n+1)` (`scripts/_limma.py:420-426`) can never return its own
nominal minimum of `1/(n+1)`.

| contrast type | n_null | paper's stated floor | repo formula's best case | **attainable floor** |
|---|---|---|---|---|
| unpaired (70×70) | 4900 | 1/4901 = 2.04 × 10⁻⁴ | 2/4901 = 4.08 × 10⁻⁴ | **3/4901 = 6.12 × 10⁻⁴** |
| shared-control (index-paired) | 70 | 1/71 = 1.41 × 10⁻² | 2/71 = 2.82 × 10⁻² | **3/71 = 4.23 × 10⁻²** |

This is not pedantic. It is visible in the deposited table: the smallest
p-value anywhere in `etiology_marker_cross.csv` is **0.000612 = 3/4901**, and
the smallest on the shared-control pairs is **0.042254 = 3/71**. Neither can
be expressed as a multiple of 1/4901 or 1/71. The manuscript
(§2.7 line 87, §5 line 173) and the code comment at
`scripts/22_etiology_comparison.py:33-36` both state the floor as 1/4901 and
1/71.

Note also that the floors imply a hard resolution limit that the paper does
not acknowledge: with 22 modules × 6 pairs, no permutation result in this
study can ever be reported at p < 6.12 × 10⁻⁴, so "exact p with no resampling
error" is true but "exact" here means exact *within a coarse lattice*.

**【Why it matters】** The floors are the paper's single explicit argument for
why its permutation scheme is superior to a bootstrap ("a hard floor on the
p-value of 1/4901 rather than an artificial 1/2,000", `22_etiology_comparison.py:35`).
The realised floor of 6.12 × 10⁻⁴ is barely finer than the 5 × 10⁻⁴ Monte-Carlo
resolution the paper is arguing against, and the shared-control floor of
0.0423 is worse than a conventional bootstrap would give. A reviewer who
checks the arithmetic will conclude the Methods section was not reconciled
with the code that produced the table — which puts the credibility of §2.7
as a whole in question, including the parts that are correct.

**【Specific fix】** Replace the sentence at §2.7 line 87 (and its twin at
§5 line 173) with:

> Because the enumeration is exhaustive, the observed labelling is itself a
> member of the null, and the 70 splits form 35 complement pairs whose axes
> are exact negations. Every observed Δ therefore has its mirror image −Δ in
> the null, so at least two of the 4,900 (unpaired) or 70 (shared-control)
> null values are always as extreme as the observed one. The minimum
> attainable p under the (ge+1)/(n+1) rule used here is consequently
> **3/4,901 = 6.12 × 10⁻⁴** for the four unpaired contrasts and **3/71 =
> 0.0423** for the two shared-control contrasts. These are lattice limits of
> the design, not properties of the effect.

---

## F2. The "paired" treatment of the shared-control contrasts uses the wrong null

**【Problem】** For the two contrasts that share a control group (Bone_D14 vs
Bone_D7; PTX vs OXA), the code pairs the split *index* so that both sides draw
the same control subset. This is not the exchangeable-label null for the
estimand, and it makes the shared-control p-values too small.

**【Evidence】** Read `scripts/22_etiology_comparison.py:286-292`:

```python
paired = (list(perm[ka]["shares_control"]) == list(perm[kb]["shares_control"]))
if paired:
    m = min(len(sa), len(sb))
    null = sa[:m] - sb[:m]
```

For a shared control the estimand collapses:

```
Δ = logFC(A) − logFC(B) = [mean(A_case) − mean(C0)] − [mean(B_case) − mean(C0)]
                          = mean(A_case) − mean(B_case)
```

The control term cancels identically, so the only exchangeable labelling is
over the **eight case samples** (4 A-case + 4 B-case), giving C(8,4) = 70
splits of those eight. The repo instead indexes two *independent* 70-split
enumerations of `[A_case(4) | C0(4)]` and `[B_case(4) | C0(4)]` by position
and subtracts them, which forces the shared control into the same
case/control role on both sides while leaving the A-case and B-case samples
free to be permuted against each other in lockstep — a scheme that
corresponds to no randomisation of the eight case labels.

The two nulls are numerically different objects. For Cd68, Bone_D14 vs
Bone_D7: SD of the repo's 70-long index-paired null = **0.2635**, SD of the
exchange-8-case null = **0.4130**; only **16 of 70** values in the repo null
are distinct (`A2_check1_permutation.py`, section (d)). Resulting p-values:

| gene | Δ | repo index-paired p | exchange-8-case p |
|---|---|---|---|
| Cd68 | −0.347 | 0.1549 | 0.4648 |
| Gfap | −0.904 | 0.0423 | 0.0986 |
| Gal | +0.686 | 0.0423 | 0.0704 |
| Scn10a | +0.254 | 0.4366 | 0.1831 |
| Il1b | −0.060 | 1.0000 | 0.9437 |

For Cd68 the repo null is **3× too narrow** and inflates apparent evidence by
a factor of three.

**【Why it matters】** §2.7 presents the index-pairing as a methodological
refinement ("the split index is paired so that both sides draw the same
control subset"). It is not a refinement; it is a different and incorrect
null, and it is wrong in the direction that makes results look more
significant. It does not currently flip any conclusion — I checked: **0 of 44**
module-level shared-control tests reach q < 0.05, and the smallest p attained
on those pairs is already the lattice floor 3/71 — so this is a Methods error
rather than a Results error. But it is presented as a point of rigour, and
pairing the wrong index is exactly the class of error that makes a
"calibrated" analysis claim look undeserved. If a future revision promotes any
shared-control contrast to a headline (Bone_D14 vs Bone_D7 is the obvious
candidate, being the paper's own primary timepoint), the error becomes
material.

**【Specific fix】** Either delete the paired branch and enumerate the eight
case labels for shared-control pairs, or state explicitly what the index
pairing does and does not claim. Paste-ready Methods text:

> For the two contrasts that share a control group, the control term cancels
> identically from Δ, so the estimand is the difference in means of the eight
> case samples. We therefore enumerate the C(8,4) = 70 assignments of those
> eight samples and compute Δ under each; the control samples are held fixed
> because they do not enter Δ. Pairing the split index of the two
> case-versus-control enumerations is *not* equivalent, because it permutes
> the case labels in lockstep with the control labels and yields a null with
> approximately two-thirds the spread of the exchangeable-label null.

---

## F3. The power analysis tests a comparison that does not exist, and its conclusion is backwards

**【Problem】** §3.6 states that "the perturbation contrast is a 2-versus-2
comparison (knockout versus overexpression on the same four animals), for
which the minimum detectable Cohen's d at 80% power is 5.65", and concludes
the layer "is only powered to detect reversal-score differences exceeding
~0.42 units, whereas the largest observed reversal difference is ≈0.27", so
the screen is "hypothesis-generating, not a target-ranking tool". The design
premise is wrong, the effect-size scale is wrong, and the conclusion does not
follow.

**【Evidence】** Three independent problems, each verified.

**(a) There is no 2-vs-2 comparison.** The reversal score is
`cos(P_k @ seed, axis_k)` where `axis_k = mean(control) − mean(case)` over 4+4
animals and `P_k` is built on those same 8 animals
(`22_etiology_comparison.py:205-210`). It is a **single scalar per
(module, contrast, perturbation)** — a functional of 8 animals, not a
per-animal measurement. KO and OE are the *same* animals with the seed sign
flipped, so reversal(KO) = −reversal(OE) identically. I rebuilt all 22
modules × 4 contrasts from raw FPKM with a separate network/diffusion
implementation: `corr(KO, OE) = −1.000000` over 176 scores and
`max |OE − (−KO)| = 0.000e+00`. A two-sample t-test requires two groups of
independent observations; here the two "groups" are related by an exact sign
flip, their correlation is −1, and the within-group difference *is* the
quantity of interest. `TTestIndPower(nobs1=2, ratio=1)` is a test of a
contrast that cannot be formed. `TTestIndPower(nobs1=2, …)` does return
d = 5.6535, so the number is arithmetically reproducible — it is the
*design assignment* that is wrong.

**(b) The effect-size scale is the wrong dispersion.** The MDE is computed as
`d × rev_sd` where `rev_sd` is the mean of the per-column standard deviations
of the 22 reversal scores (`33_negative_result_quantification.py:124-126`) —
that is the **between-module spread**, i.e. real biological heterogeneity
plus noise, not the sampling SD of one module's score. I recomputed it:
SD(rev_Bone_D14) = 0.0792, SD(rev_Bone_D7) = 0.0845, SD(rev_PTX) = 0.0708,
SD(rev_OXA) = 0.0618; mean = 0.0741, matching the deposited 0.075.

The dispersion a power analysis needs is the SE of Δ under resampling of the
animals. Two independent estimates, both from the paper's own artefacts:

- **Bootstrap:** the deposited `ci95_lo/hi` imply
  SE(Δ) = (ci95_hi − ci95_lo)/3.92; median over the 132 tests = **0.0185**
  (mean 0.0205, range 0.0022–0.0509). This is 4× smaller than the 0.074
  used.
- **Leave-one-animal-out**, which I ran myself for Glial/immune at Bone_D7,
  rebuilding both the network and the axis each time: reversal ranges
  +0.1421 to +0.1601, **LOO SD = 0.0055**, and dropping GSM4506605 alone
  moves the score by −0.0158.

**(c) The corrected MDE is ~10× smaller than the reported one, and the
design is not underpowered.**

| quantity | d | × SE | = MDE (reversal units) |
|---|---|---|---|
| paper: `nobs1=2` × between-module SD | 5.653 | 0.0741 | **0.419** |
| corrected: unpaired 4v4 × bootstrap SE(Δ) | 2.381 | 0.0185 | **0.044** |
| corrected: paired n=4 × bootstrap SE(Δ) | 2.128 | 0.0185 | **0.039** |
| LOO-based: unpaired 4v4 × LOO SD | 2.381 | 0.0055 | **0.013** |

The largest observed target Δ is 0.2743 — **6.2× the corrected unpaired MDE**
and 7.0× the corrected paired MDE. The paper's own exact permutation test
agrees and is decisive on this point: **26 of 132** tests have
`p_permutation_exact < 0.05` and the smallest is **p = 0.00143**. A design
that could only see 0.42-unit differences cannot produce p = 0.0014.

**(d) Is the reversal score paired?** Partly, and the paper does not say so.
For the two shared-control pairs Δ *is* paired (common control animals); for
the four cross-study pairs it is not. The correct analysis is therefore
split: paired *t*/bootstrap on the two shared-control contrasts, unpaired on
the four cross-study ones. The paper uses a single unpaired two-sample
calculation for everything.

**(e) The claim about the largest observed difference is internally
inconsistent across the paper's own artefacts.** §3.6 says the largest
observed reversal difference is ≈0.27. I recomputed: max |Δ| over the 132
per_contrast KO tests = **0.3073**, attained by CTRL erythrocyte (PTX vs OXA);
max |Δ| among target modules = 0.2743. The `T5` block of
`negative_result_quantification.json` says "observed cross-aetiology reversal
differences are <0.2". Three artefacts, three different numbers
(0.20 / 0.27 / 0.307).

**【Why it matters】** This is the strongest new statistical claim of the round
and it is the load-bearing sentence for the paper's methodological
conclusion, repeated in §3.6, §5 and §6. As written it says the study is too
underpowered to detect what it is trying to detect. Recomputed on the correct
design and the correct dispersion, the study is *well* powered for the
reversal contrast and produces 26 nominally significant permutation results
whose failure to survive BH is a statement about multiplicity across 132
correlated hypotheses, not about sample size. The paper therefore mis-attributes
its own null result. Two consequences: (i) a reader who accepts the MDE will
reasonably discount the whole perturbation layer as uninterpretable for
reasons of design, when the real problem is that the layer's z-scores do not
separate targets from controls — which is a *different*, and more interesting,
criticism; (ii) the claim is repeated in the Limitations as a reason the
conclusions "can only exclude very large target effects", which is false and
will be checked.

Note that this does **not** rescue the perturbation layer. The layer's
failure is established by other evidence (ρ = 0.99 with the additive model,
controls outranking targets, erythrocyte |z| = 7.39). The problem is that the
paper reaches that conclusion partly through a power argument that, if the
authors had checked it, would have pointed the other way.

**【Specific fix】** Delete the "2 versus 2" sentence and replace the whole
power passage in §3.6 (line 139) with:

> The reversal score is a single cosine computed from the 4 case and 4 control
> animals of one contrast, so knockout and overexpression are the same eight
> animals with the seed sign flipped and are exact mirror images
> (r = −1.000); there is no 2-versus-2 comparison and a two-sample t power
> calculation does not apply. What can be powered is the cross-aetiology
> difference Δ = reversal(A) − reversal(B), which is a paired quantity for
> the two contrasts that share a control group and unpaired for the four that
> do not. Taking the standard error of Δ from our own bootstrap
> (median SE = 0.0185 over the 132 tests) and the exact 4-versus-4 or paired
> n = 4 requirements (d = 2.38 and d = 2.13 at 80% power, two-sided
> α = 0.05), the design resolves Δ of approximately 0.044 (unpaired) or 0.039
> (paired) reversal units. The largest observed target Δ is 0.274, six times
> that threshold, and 26 of the 132 exact permutation tests return p < 0.05.
> The reversal layer is therefore adequately powered for the effect sizes it
> produces; what limits it is not sample size but the absence of any
> separation between analgesic-target and negative-control modules (see
> above), together with the fact that 132 hypotheses share 22 modules and six
> linearly dependent aetiology contrasts.

Replace the Limitations sentence at §5 line 173 ("with a minimum detectable
Cohen's d of 5.65 at 80% power, so the layer can only exclude very large
target effects") with the same corrected figures. Also reconcile the three
inconsistent "largest observed reversal difference" values to 0.307 (all
modules) / 0.274 (targets only), and fix the `T5` interpretation string in
`scripts/33_negative_result_quantification.py` to use the SE of Δ from the
bootstrap rather than the between-module SD of reversal.

---

## F4. T1 is a one-sided test oriented against the data, and "statistically indistinguishable" is the wrong reading

**【Problem】** §3.6 concludes that the 17 analgesic-target modules are
"statistically indistinguishable from the five negative controls" on the
basis of p = 0.50 (signed) and p = 1.00 (|z|). Both p-values come from
`alternative="greater"` with target in the first argument — i.e. the
alternative is "targets larger than controls", which is the direction the data
reject.

**【Evidence】** `scripts/33_negative_result_quantification.py:41-42`:

```python
u1 = mannwhitneyu(tgt.z_vs_randset, ctl.z_vs_randset, alternative="greater")
u2 = mannwhitneyu(tgt.absz,        ctl.absz,        alternative="greater")
```

I recomputed all four combinations on the deduplicated KO universe
(102 target rows, 30 control rows):

| | target > control | control > target | two-sided |
|---|---|---|---|
| signed z | 0.3529 | 0.6491 | 0.706 |
| \|z\| | **0.9991** | **9.57 × 10⁻⁴** | **1.91 × 10⁻³** |

(p = 0.9991/0.9991 in the 264-row universe matches the deposited 0.4973/1.00
after deduplication.) The mean |z| figures are as the paper reports —
target 1.462, control 2.589 — and Welch *t* on |z| gives *t* = −4.49,
p = 2.6 × 10⁻⁵.

Two further problems with T1 as specified:

- **Pseudoreplication.** The test treats 264 rows (or 132 after dedup) as
  independent observations. They are not: 22 modules are scored on the same 4
  contrasts and the same 6 pair-differences. Pairwise Pearson *r* among the
  132 Δ values: median 0.032, max 0.980, **60% of pairs exceed |r| = 0.5**. The
  effective rank of the 22×22 Δ correlation matrix is **1.44 of 22** (top
  eigenvalue 18.02, next 3.14). A rank test on 132 such rows has an unknown
  and certainly inflated type-I error.
- **"Indistinguishable" over-claims.** A non-significant difference test does
  not license a claim of indistinguishability; it licenses "no evidence of a
  difference at this n". Given that the point estimates go the *other* way
  and the two-sided test rejects at p ≈ 0.002, "indistinguishable" is not
  merely unlicensed — it is contrary to the data.

Correct analysis, respecting the block structure: permuting the 22
**module→class labels** (the exchangeable unit, since modules are what is
being classified) over 20,000 draws gives two-sided **p = 0.0331** for
(mean |z|_target − mean |z|_control) = −1.113. Effect size: Cliff's δ =
−0.374, rank-biserial r = −0.374; median |z| 1.073 (targets) vs 2.470
(controls), difference −1.397, bootstrap 95% CI [−2.130, −0.272].

**【Why it matters】** The paper's methodological conclusion — that the layer
does not prioritise analgesics — is *understated* by its own statistic.
Reporting the one-sided p in the anti-hypothesis direction converts a
significant result (controls have larger perturbation effects than targets,
p ≈ 0.002 two-sided, p = 0.033 by the correct permutation) into an apparent
null. This is the opposite of the usual over-claiming failure mode, but it is
still a defect: the manuscript asserts a *specific* statistical conclusion
("indistinguishable") that its own numbers contradict, and a reader who
recomputes will find the two-sided p immediately. It also weakens T4, whose
"55.7th percentile" is computed over a universe that is 77.3% targets by
construction — controls sit at the **24.1st** percentile, i.e. they have
*larger* effects, and the paper reports only the target median.

**【Specific fix】** Replace the sentence at §3.6 line 139 beginning "Comparing
the two module classes directly" with:

> Comparing the two module classes directly, negative-control modules have
> **larger** standardised perturbation effects than analgesic targets, not
> smaller: mean |z| 2.59 versus 1.46, median |z| 2.47 versus 1.07, Cliff's δ =
> −0.374 (two-sided Mann–Whitney p = 1.9 × 10⁻³ on the 132 deduplicated
> tests). Because the 132 tests share 22 modules scored on the same six
> aetiology contrasts (effective rank 1.44 of 22), we additionally permute the
> 22 module-to-class labels, which is the exchangeable unit: p = 0.033 for the
> difference in mean |z| of 1.11. Analgesic-target modules therefore do not
> merely fail to separate from negative controls — the negative controls are
> preferentially flagged, which is a stronger failure of specificity than we
> originally stated.

Also report the control percentile (24.1st) alongside the target median in
T4, and state that the T1/T4 universe is 77.3% targets by construction.

---

## F5. T2's conclusion is stated as the opposite of what the numbers show

**【Problem】** §3.6 concludes that "the randset null is non-selective" /
"flags analgesic targets and controls equally", on the basis that 28 of 64
randset-significant tests are controls.

**【Evidence】** The count is right (I confirm 64 = 32 KO + 32 OE; 28 controls
= 43.8%). But no base rate is compared against it. The universe is 17 target
+ 5 control modules, so the control base rate is **5/22 = 22.7%**. Observing
28/64 = 43.8% controls is close to **twice** the base rate:

- exact binomial, H0: p_control = 0.227, n = 64, k = 28 → **two-sided p = 2 × 10⁻⁴**.

On the deduplicated universe: 14/32 = 43.8% controls, same direction.

The bootstrap-CI part of the argument does not carry the weight placed on it.
The CI is on Δ versus 0; the randset p is on Δ versus a degree-matched random
gene set. **They test different nulls.** Of the 264 rows, **141** have a
bootstrap CI excluding zero, of which only 64 are randset-significant —
so **77 non-significant rows also carry a CI excluding zero**. "All 64
significant tests have a CI excluding zero" is therefore uninformative about
selectivity: it is what you expect whenever Δ ≠ 0 with n = 4 and a bounded
cosine, and it holds for the 77 non-significant rows too. A percentile
bootstrap on 4 + 4 samples is also a weak interval — it cannot express
skewness or the boundedness of a cosine, and the code silently drops
degenerate resamples (`dB = dB[:N_BOOT - ndeg]`, line 323) rather than
resampling them.

**【Why it matters】** "Non-selective" and "flags controls at twice the base
rate, p = 2 × 10⁻⁴" are different claims with different implications. The first
says the null is uninformative; the second says the null is *anti*-informative
— it preferentially selects against analgesic targets. The second is the
stronger and better-supported statement, and it is the one that matches the
paper's own erythrocyte result (|z| = 7.39, the largest in the study) and the
T1 finding in F4. As written, the paper understates its own central
methodological criticism and simultaneously attaches an inferential
non-sequitur (the CI) to it.

**【Specific fix】** Replace the second sentence of §3.6 line 139 ("All 64
significant tests carry bootstrap 95% confidence intervals that exclude zero,
yet 28 of the 64 (44%) are negative-control modules, so the null flags
analgesic targets and controls equally") with:

> Of the 32 deduplicated randset-significant tests, 14 (43.8%) are negative-
> control modules against a base rate of 5/22 = 22.7% (exact binomial
> two-sided p = 2 × 10⁻⁴). The randset null therefore does not flag targets
> and controls equally: it flags controls at approximately twice their base
> rate, i.e. it is anti-selective for analgesic targets. Separately, 141 of
> the 264 tests have a bootstrap 95% CI for Δ excluding zero while only 64 are
> randset-significant; a CI excluding zero tests Δ ≠ 0 and carries no
> information about whether this module is more extreme than a degree-matched
> random gene set, so it does not bear on the selectivity question.

---

## F6. Table 1's `q` column is not a p-value for the `mean t` printed beside it

**【Problem】** Table 1 presents `mean t` and `q` in the same row, which reads as
"this mean t has this q". The two columns answer different questions.

**【Evidence】** `scripts/20_multistage_pipeline.py:391`:

```python
_, pu = stats.mannwhitneyu(sub["t"].values, bg, alternative="two-sided")
```

where `sub` is the module's 3–11 genes and `bg` is `r["t"].values` for the
**whole transcriptome** (~13,893 genes). So `q_comp` tests "are this module's
genes' moderated *t* larger than the transcriptome's", a rank test with
*n* = 5 versus 13,900. It is **not** a test of H0: mean moderated *t* = 0, and
there is no interval or test for `mean t` anywhere in
`module_level_signal.csv` (columns: `contrast, module, kind, n, mean_t,
mean_logFC, p_comp, genes_used, q_comp`).

The consequence is visible in the table. K⁺/HCN channels has mean *t* = −5.64
with q = 0.0053; opioid receptors has mean *t* = −3.76 with q = 0.077 — the
ordering by q is not the ordering by |mean t|. If a reader takes the naive
route (mean *t* ± SD/√n over genes), CTRL interferon gives z = 15.3,
pan-neuronal −7.4, K⁺/HCN −8.5, GABAergic −7.9 — none of which is what
`q_comp` computes, and none of which is reported. Both routes also treat the
module's genes as independent replicates, which they are not by construction
(a module *is* a co-expression set).

I note the paper's substantive use of this table is sound and I endorse it:
CTRL pan-neuronal, a negative control, has the largest |mean t| of any module
at day 7 (−5.20, with 8 of 10 genes individually at q < 0.05), which is
exactly the compartment-effect argument the paper makes. The defect is that
the inferential object attached to it is mismatched.

**【Why it matters】** Table 1 is the paper's headline biological result
("day 7 was dominated by a type-I interferon programme, module mean moderated
t = +24.9, q = 6.8 × 10⁻⁵", §3.3, repeated in the Abstract). A reader who
takes 6.8 × 10⁻⁵ as the significance of +24.9 will be misled about what was
tested. Given that the paper's own negative-control module has the second
largest |mean t|, the honest inferential statement is about the *comparison*
between modules, not about each module against zero — and that comparison is
what the exact permutation test in §3.5/Table 3 provides, not `q_comp`.

**【Specific fix】** Rename the Table 1 column from `q` to
`q vs transcriptome (rank test)`, and add a footnote:

> The q column is a Benjamini–Hochberg value over 66 module × stage
> combinations for a Mann–Whitney rank test of the module's 3–11 genes'
> moderated *t* against the *t* of all 13,893 genes. It is **not** a p-value
> for the mean moderated *t* in the adjacent column, and no test of H0: mean
> moderated *t* = 0 is reported. Module genes are co-regulated by construction
> and are not independent replicates, so we do not attach a p-value to the
> module mean. Cross-aetiology comparisons of module mean logFC are tested
> exactly in Table 3.

---

## F7. The single-cell power split uses a circular yardstick and a non-standard effect size

**【Problem】** §3.8 argues gene-level DE was "well powered" from the median SE
of "stable" genes, and calls the myeloid-fraction result "underpowered" using
Cohen's *d* = 2.38 and *d* = 1.37. The first argument is circular; the second
mis-states the effect size and mis-attributes the non-significance.

**【Evidence】**

**(a) The median-SE argument is circular.** `33_negative_result_quantification.py:95-97`
computes `se = stable.logFC.abs()/stable.t.abs()` over genes with q > 0.5 and
takes the median = **0.000268**, then multiplies by d = 2.38 to claim the
design detects |logFC| ≥ 0.001. I recomputed the SE distribution: 5th
percentile 4 × 10⁻⁵, 25th 6 × 10⁻⁵, 50th 2.7 × 10⁻⁴, 75th 4.2 × 10⁻³, 95th
2.7 × 10⁻². The median is small because the "stable" set is dominated by genes
with essentially zero variance, not because the assay is uniformly precise.
Meanwhile the design detected **2,220 genes** at q < 0.05 with median
|logFC| = **0.0448** — two orders of magnitude above the claimed 0.001
threshold. The stated MDE cannot be right when the design exceeds it 40×.

**(b) The specific inference drawn from Ifit2 does not follow.** §3.8 says the
design "detected Ifit2 at log₂ fold change −0.048 (q = 0.042) — so the
unchanged myeloid markers Cd68 and Aif1 are a genuine null for large
effects." I verified Ifit2: logFC = −0.0476, t = −3.86, q = 0.0417. But
|logFC| = 0.048 is a **3.3% change** — a tenths-of-a-percent effect. Detecting
a 3% change says the assay is sensitive to tiny effects; it does not license
"a genuine null for large effects" in a *different* gene. The correct object
is the interval on Cd68/Aif1 themselves, which is available and which I
recomputed: Cd68 logFC = −0.0005 (q = 0.783), Aif1 logFC = +0.0009
(q = 0.920). Those values do support "null" — but they support it directly,
and the Ifit2 detour is what a reader must be misled by.

**(c) The myeloid Cohen's *d* uses a non-standard pooled SD.**
`33_negative_result_quantification.py:101-103` computes

```python
pooled_sd = sqrt(((con-con.mean())**2).sum()/7 + ((b16-b16.mean())**2).sum()/7)
```

which equals `sqrt((SS_total)/7)`. The conventional pooled SD is
`sqrt(SS_total/(n1+n2-2)) = sqrt(SS_total/6)`. Recomputed: SS_Con =
1.3165 × 10⁻⁶, SS_B16 = 2.7438 × 10⁻⁷, mean difference 6.524 × 10⁻⁴;
repo SD = 4.767 × 10⁻⁴ → *d* = **1.3685**; conventional SD = 5.149 × 10⁻⁴
→ *d* = **1.2670**. So the reported 1.37 is inflated by 8% relative to the
standard convention. The direction of the paper's conclusion is unchanged, but
the number is not the standard *d*.

**(d) The "underpowered" label is a property of the rank test, not of the
effect size.** Cohen's *d* on a proportion bounded in [0, 1] with observed
values ~5 × 10⁻⁴ and skewness 0.87 is fragile: on the log scale the same data
give *d* = **0.083** instead of 1.267 — a 15-fold change from a monotone
reparameterisation. Cohen's *d* is not scale-invariant, and reporting it as
*the* effect size for a bounded fraction with a log link available is a poor
choice. More importantly, a Student *t*-test on the same eight numbers gives
*t* = 1.792, df = 6, **p = 0.123** — essentially the MWU's p = 0.114. The
non-significance is the *n* = 4 vs 4 resolution of a rank test (its own floor
is 2/70 = 0.0286), not a failure to detect a Cohen's *d* = 1.37.

**(e) The one-sided p is reported after the direction was seen.** §3.8 gives
both "two-sided p = 0.114" and "one-sided p = 0.057" and calls the result
"inconclusive". The direction (Con > B16) was established by looking at the
data. I verified U(B16 > Con) = **2.0** and both p-values exactly. Reporting
the one-sided value alongside the two-sided one, in a sentence that draws a
conclusion from it, is a choice the reader cannot audit.

**【Why it matters】** This passage is the paper's defence of a null it does not
like. It is also the passage where the paper models good practice for a
negative result — and then uses a circular yardstick, a non-standard
denominator, a scale-dependent effect size on a bounded fraction, and a
post-hoc one-sided p. The substantive conclusion (the myeloid drop is
inconclusive, not evidence of equality) is defensible and I would keep it, but
it currently rests on reasoning that would not survive a methods reviewer.
Giving the direct Cd68/Aif1 intervals is both simpler and stronger.

**【Specific fix】** Replace the gene-level and myeloid sentences in §3.8
(line 151) with:

> Gene-level pseudobulk DE in this design resolves small effects: 2,220 genes
> reach q < 0.05 with a median |logFC| of 0.045, and the interferon-
> stimulated gene Ifit2 is recovered at log₂FC −0.048 (q = 0.042). We
> therefore assess the myeloid markers directly rather than by inferring
> power from the median standard error of non-significant genes, which is
> dominated by near-zero-variance genes and is not a valid yardstick: Cd68
> (log₂FC −0.0005, q = 0.78) and Aif1 (+0.0009, q = 0.92) are unchanged with
> no gene-level signal to explain away. The myeloid *fraction*, by contrast,
> is a different quantity with a different resolution: 0.117% in controls
> versus 0.052% in B16F10, a 55.8% relative drop (Mann–Whitney U = 2,
> two-sided p = 0.114; Student t on the same eight values p = 0.123). At
> n = 4 per group a rank test cannot resolve below p = 0.029, so this is an
> inconclusive null rather than evidence of equality. We report the effect on
> the log scale (log₂ ratio = −1.18) rather than as Cohen's *d*, which is not
> invariant to monotone rescaling of a bounded fraction.

---

## F8. Table 2's "Largest Δ" column reports a selection the q-value does not account for

**【Problem】** Table 2 reports, per gene, the largest |Δ| over the six aetiology
pairs, together with the q-value of that particular pair.

**【Evidence】** I recomputed the per-gene maxima from
`etiology_marker_cross.csv`: Cd68 D7|PTX (|Δ| = 3.418, p = 0.000612,
q = 0.012487), Aif1 D7|PTX (3.076, q = 0.012487), Gal D7|PTX (2.592,
q = 0.012487), Npy D7|PTX (1.961, q = 0.105), Gfap D7|OXA (1.450,
q = 0.091). The q-values are BH over the 102 deposited marker tests
(17 genes × 6 pairs — the table has 102 rows, not the 19 × 6 = 114 the code's
marker list implies, because two markers are absent from the shared gene
space). Choosing the maximum over six pairs and then quoting a q computed
over all 102 rows does not correct for the selection step: the relevant family
for a "largest Δ per gene" claim is 17 × 6 selected maxima, not 102 rows.

To be fair to the authors: the BH family (102 rows) and the set of per-gene
maxima (17 rows) are different objects, and the *selected* q's are the
uncorrected-for-selection ones. I verified BH is internally consistent
everywhere I checked (`np.allclose` to 1e-9 in all six families), so this is
a presentation/selection issue, not an arithmetic error.

**【Why it matters】** Table 2 is the paper's positive headline. The selection
is benign here — the top genes are significant by a wide margin under any
reasonable correction — but reporting a selected maximum with a
selection-unadjusted q is the kind of thing a statistical reviewer flags, and
it costs nothing to fix by stating the family explicitly.

**【Specific fix】** Add to the Table 2 caption:

> For each gene we report the aetiology pair with the largest |Δ|. The q
> values are Benjamini–Hochberg over all 102 gene × pair tests (17 markers ×
> 6 pairs; 2 of the 19 candidate markers are absent from the shared gene
> space). Because the reported pair is selected as the maximum over the six
> pairs, the q should be read as an unadjusted-for-selection value; the
> selection does not affect the top-ranked genes, whose p is at the lattice
> floor of the permutation scheme.

---

## F9. The §3.7 lowess-span robustness sweep is not on the same gene universe as the number it is meant to protect

**【Problem】** §3.7 reports "316–325 genes at q < 0.05 across spans 0.2–0.6"
as evidence that differential expression is robust to the lowess span. §3.2
reports **235** genes at q < 0.05 for day 14. These are not on the same gene
universe, and the reported reference value (235) falls outside the swept range.

**【Evidence】** `scripts/23_sensitivity.py:75-80` builds the gene universe as
the **intersection with the chemo dataset** (13,601 genes shared), whereas
`scripts/20_multistage_pipeline.py` reports the bone-only universe
(`DE_D14_vs_Sham_adj.csv`, 13,893 genes). The sensitivity sweep is on D14
(`23_sensitivity.py:164`, `Xc = Lb[ca14 + sh]`) — the same contrast as §3.2's
"235". Within the sensitivity artefact itself the reference span (0.4) gives
**320**, and the sweep spans 316–325. The main-pipeline value 235 sits outside
that range and is never reconciled. I confirmed the main-table counts:
D1 = 0, D7 = 4,600, D14 = 235 at q < 0.05, matching §3.2 exactly, and
0 / 5,867 / 860 at q < 0.10, also matching.

The other §3.7 numbers check out: `sensitivity_summary.json` gives
`grid_discovery_min_spearman` = 0.9593 (paper: min 0.959),
`grid_discovery_median_spearman` = 0.9876 (paper: median 0.988),
`grid_delta_min_spearman` = 0.9887 (paper: min 0.989),
`grid_discovery_min_top5_overlap` = 5 (paper: "top five modules identical in
all 60 settings"). The grid has 75 points (60 α×k×hv + 15 α×k), matching
"60 settings" then "15 settings". So the robustness claim itself is sound;
only the DE-count bridge is missing.

**【Why it matters】** Minor in isolation, but §3.7 is the paper's argument
that "sections 3.2, 3.3 and 3.5 are therefore not artefacts of tuning", and
it currently cannot be checked for §3.2 without re-running the pipeline. A
reader comparing 235 to 316–325 sees a 40% discrepancy with no explanation.

**【Specific fix】** Replace the lowess sentence in §3.7 (line 145) with:

> Differential expression barely moved with the lowess span: on the
> 13,601-gene universe shared with the chemotherapy dataset, day 14 yields
> 316–325 genes at q < 0.05 across spans 0.2–0.6 and 323 with the trend
> switched off, and Spearman ρ of the moderated *t* against the reference is
> ≥ 0.9999. (The 235 genes quoted in Section 3.2 are on the 13,893-gene
> bone-only universe used for the primary analysis; the two differ by the 292
> genes absent from the chemotherapy dataset.)

---

## F10. A significant compartment-level ISG result is omitted from §3.8

**【Problem】** §3.8 states "the only significant ISG changes were
down-regulation of Ifit2 in neurons (q = 0.046) and Isg15 in satellite glia
(q = 0.007)". There are three such changes in the deposited table, and the
omitted one has the **smallest** q of the three.

**【Evidence】** From `sc_GSE292522_celltype_DE.csv`, rows with an ISG at
q < 0.05:

| compartment | gene | logFC | t | q |
|---|---|---|---|---|
| SatelliteGlia | Isg15 | −0.630 | −4.54 | **0.0071** |
| Neuron.NF | Ifit2 | −0.060 | −3.85 | **0.0412** |
| Neuron | Ifit2 | −0.049 | −3.89 | 0.0459 |

The manuscript reports the first and the third and omits Neuron.NF
(q = 0.0412). Both reported values match the deposited table, so this is
omission, not error. I also checked the whole-DRG claim: the 38 ISGs in
`scripts/32_percelltype_sc.py` are all present in the pseudobulk table, and
**zero are up-regulated at q < 0.05** (one, Ifit2, is down at q = 0.042) —
so the substantive claim "no ISG was up-regulated" is correct and I verified
it.

**【Why it matters】** Small in effect (a third instance of the same direction
in a neuronal subtype), but the omission is in the paper's favour
elsewhere — no, actually it is *against* the paper's own "the interferon
programme is absent at the compartment level" framing only in the sense that
it removes a data point from a list. Selective omission of the strongest
member of a three-item list, in a paper that elsewhere makes a point of
negative controls and full disclosure, is the kind of omission that invites
suspicion. Also note the multiplicity family: `32_percelltype_sc.py:160`
applies BH **within each compartment separately** (8 independent families),
while §3.8 reads across all 8 compartments at once. For a claim of the form
"no compartment showed X", the family is 8 × n_ISG.

**【Specific fix】** Replace the relevant sentence in §3.8 (line 153) with:

> No compartment showed interferon-stimulated-gene up-regulation. The
> significant ISG changes were all down-regulation: Isg15 in satellite glia
> (log₂FC −0.63, q = 0.007), Ifit2 in neurofilament-positive neurons (−0.060,
> q = 0.041) and Ifit2 in neurons (−0.049, q = 0.046). Benjamini–Hochberg was
> applied within each compartment separately (eight families); because the
> claim of interest is "in no compartment", we additionally report the
> per-compartment interferon-stimulated-gene score, which is not higher in
> B16F10 than control in any of the eight compartments (and is lower in
> satellite glia, MWU p = 0.029).

---

## F11. `U = 2/14` is a malformed statistic

**【Problem】** §3.8 reports "the Mann–Whitney U test gives U = 2/14,
two-sided p = 0.114".

**【Evidence】** I recomputed U(B16F10 > Con) = **2.0** exactly, and both
p-values exactly (two-sided 0.1143, one-sided 0.0571). The U value is right.
But for n₁ = n₂ = 4 the maximum attainable U is n₁n₂ = **16**, not 14, and
the "14" appears to be the number of remaining animals (28 total − 14
myeloid-bearing?) rather than any property of the test. As written "U = 2/14"
is not interpretable by a reader.

**【Why it matters】** Trivial in isolation. Included because it is in a Table
5 cell that a reader will check, and because it is the kind of notation error
that suggests the numbers were transcribed rather than recomputed at the
proof stage.

**【Specific fix】** Replace "U = 2/14" with "U = 2 of 16 possible".

---

# PART 2 — § STANDS UP

Findings I specifically suspected and then verified as **correct**. These are
deliverables, not filler.

### S1. The exhaustive enumeration is genuinely exhaustive

I expected to find a non-exhaustive or Monte-Carlo shortcut. There is none.
`all_splits(4,4)` produces exactly C(8,4) = 70 rows, each with exactly 4 case
labels, all 70 subsets distinct, and the complement map is a fixed-point-free
involution. The 70 × 70 = 4,900 unpaired enumeration is the complete product
orbit. All 70 axes are distinct and exactly antisymmetric under complement
(`max |axis_r + axis_comp(r)| = 0.000e+00`). **The claim "enumerated
exhaustively … exact p with no resampling error and no dependence on a seed"
(§2.7) is true.** Only the stated floor is wrong (F1).

### S2. Cd68 day-7 vs paclitaxel reproduces exactly from raw counts

I re-parsed the GEO member files with my own parser (no shared code),
re-applied the median-of-expressed-genes scaling, and rebuilt the exhaustive
null. Bone_D7 vs PTX for Cd68: logFC(D7) = +2.839023, logFC(PTX) = −0.579270,
**Δ = +3.418293**, p = 0.000612, q = 0.012487 — identical to
`etiology_marker_cross.csv` to six decimals. The same holds for Aif1
(Δ = +3.075712), Gal (Δ = −2.592142), Ccl2 (Δ = +1.390698), Il1b
(Δ = +1.230465), Calca (Δ = −0.769814), Scn10a (Δ = −0.702180). Table 2 is
solid. My gene universe also reproduces exactly: 13,601 shared genes.

### S3. The knockout/overexpression mirror property is exact, not approximate

The paper says "identical exact-permutation p to eight decimals"; the truth is
stronger. Over all 132 matched (module, pair) cells I recomputed
`max |rev_A(KO) + rev_A(OE)| = 0`, `max |rev_B(KO) + rev_B(OE)| = 0`,
`max |Δ(KO) + Δ(OE)| = 0`, and `max |p_perm(KO) − p_perm(OE)| = 0` — **exact
zeros, not eight-decimal agreement**. I then rebuilt all 22 modules × 4
contrasts from raw FPKM with an independent network/diffusion
implementation: `max |rebuilt − deposited| = 9.52 × 10⁻¹⁵` and
`max |OE − (−KO)| = 0.000e+00`. The mirror property is a theorem here (cosine
is odd in its first argument and the seed is the only thing that changes), and
the code realises it exactly.

### S4. Collapsing 264 → 132 by deduplicating KO is the right call, with one caveat

The manuscript's reasoning is sound: KO and OE carry one piece of information,
so BH over the raw 2N rows double-counts each hypothesis and is
anti-conservative. Deduplicating to the 132 KO rows is correct and I verified
the deposited q's are BH over exactly those 132. The `max|Δ| = 0` result above
is what licenses it.

**The caveat** (a caveat, not an error): the 132 remaining hypotheses are *not*
independent, so "132 independent hypotheses" overstates. 22 modules are scored
on the same 4 contrasts; the 6 aetiology pairs are C(4,2) = 6 comparisons among
4 aetiologies, of which only **3 are independent** (I verified the algebraic
identity numerically: for Glial/immune, d(D14,PTX) − d(D7,PTX) = −0.0267 =
d(D14,D7) exactly; same for CTRL erythrocyte, +0.1215 = +0.1215). Pairwise
Pearson *r* among the 132 Δ values exceeds |r| = 0.5 for 60% of pairs;
effective rank of the 22 × 22 Δ correlation matrix is 1.44 of 22. BH's
guarantee requires independence or PRDS, so the reported q's are not
formally valid FDR controls. In practice this does not harm the paper's
conclusion (which is a null: min q = 0.165) and it is conservative in the
direction that matters here, but the word "independent" in the Abstract should
be "distinct".

### S5. min q = 0.165, p = 0.0014, and the glial/immune row all reproduce

§3.6's specific quoted row — glial/immune, Δ = +0.274, 95% interval
0.203–0.287, p = 0.0014, q = 0.165 — matches the deposited table exactly
(Bone_D7 vs PTX: Δ = +0.2743, CI [0.2034, 0.2872], p = 0.00143,
q = 0.1654). I recomputed BH over the 132 p-values and reproduced 0.1654 to
4 dp. "No module's reversal Δ reached q < 0.05" is correct: **0 of 132**.

### S6. The erythrocyte result is the paper's best negative control and it holds

The single strongest claim in §3.6 — that the erythrocyte set carries the
largest standardised effect in the study — reproduces exactly. Over the 132
per_contrast KO tests, max |z_vs_randset| = **7.3878**, attained by CTRL
erythrocyte Bone_D14 vs Bone_D7, matching the quoted 7.39. Same module:
PTX vs OXA z = **−5.7822**, matching the quoted −5.78. CTRL pan-neuronal
z = −3.74 as quoted. The gene set carries no pain biology, so this is a clean
demonstration that the layer reports composition, not pharmacology.

### S7. The 32/132 randset-significant count and its 14-control breakdown reproduce

The deduplicated count is exact: 32 KO tests with q_randset < 0.05 (18 target,
14 control), doubling to 64 across both perturbations with 36 target / 28
control — matching the T2 JSON. The claim "64 of the 264 (or 32/132
deduplicated) … 14 of those 32 were negative-control modules" is arithmetically
correct. (The *interpretation* is what fails — F5.)

### S8. "No ISG was up-regulated" is correct

I read the 38-gene ISG list out of `scripts/32_percelltype_sc.py`, confirmed
all 38 are present in the pseudobulk table, and verified **zero** are
up-regulated at q < 0.05 (and zero at q < 0.10). Only Ifit2 is down
(q = 0.042). The per-compartment ISG *score* is not higher in B16F10 than
control in any of the 8 compartments either. This is a genuine negative,
correctly obtained.

### S9. The whole-DRG and compartment-level single-cell counts reproduce exactly

55,967 cells; 1,517 up and 703 down at q < 0.05 (2,220 total) — matches
Table 5. Myeloid fractions 0.117% (control) and 0.052% (B16F10), ratio 0.4423,
i.e. a 55.8% drop, matching "0.12% / 0.05% / 55% drop". MWU two-sided
p = 0.1143 and one-sided p = 0.0571, both exact. Gfap q = 0.0061, Ifit2
q = 0.0417, Oas1c q = 0.0136, Atf3 q = 0.0001, Sprr1a +1.5539 with
t = 54.49, Gal +1.4247, Scn10a −0.2781 (q = 0.0064). BH reproduces the
deposited genome-wide q to 5.9 × 10⁻¹⁵ over all 34,603 genes.

### S10. The §3.2 and §3.7 discovery numbers reproduce

Day 1/7/14 at q < 0.05: **0 / 4,600 / 235** — exact. At q < 0.10:
**5 / 5,867 / 860** — exact. Median |logFC| **0.0731 / 0.1576 / 0.1187**,
matching 0.073 / 0.158 / 0.119. The contaminated sample GSM4506604 carries
6.257% muscle FPKM with a robust z of 120.9, an order of magnitude above the
next highest (GSM4506613 at 0.160%) — exactly as described. The §3.7
robustness numbers all reproduce: min ρ 0.9593, median 0.9876, Δ min
0.9887, top-5 overlap 5/5 in all 75 grid points, 316–325 genes across spans.

### S11. The design is honestly described elsewhere

§2.4 correctly states that at n = 4 the ordinary *t* has 2–6 residual df and
is poorly calibrated, and that moderation borrows strength across ~14,000
genes. §2.7 correctly states that a permutation test gives no interval and
that a bootstrap is needed for the CI. §2.8 correctly states that the unit of
inference in the compartment extension is the animal, not the cell — and the
deposited `sc_GSE292522_celltype_DE.csv` confirms per-animal pseudobulk
limma-trend across 4 + 4 animals. The manuscript's methodological
self-awareness is real; F3, F4 and F7 are exceptions to it, not the rule.

---

# PART 3 — § QUESTIONS FOR THE AUTHORS

I do not guess answers. Each of these is something I could not resolve from
the deposited artefacts.

**Q1 (F3, blocking).** What was the intended comparison behind "2 versus 2
(knockout versus overexpression on the same four animals)"? As implemented,
KO and OE are the same eight animals with the seed sign flipped
(r = −1.000000). If a different estimand was intended — for example the
per-animal reversal score under a leave-one-out scheme — please state it and
report the MDE for that estimand.

**Q2 (F3, blocking).** Which dispersion did you intend as the effect-size
scale for the MDE? The between-module SD of reversal (0.0741) and the
bootstrap SE of Δ (0.0185) differ by 4×, and the leave-one-out SD of the
reversal score is 0.0055. The reported MDE of 0.42 rests on the first.

**Q3 (F3).** Given that 26 of 132 exact permutation tests return p < 0.05 and
the smallest is 0.00143, how do you reconcile this with a design "only powered
to detect reversal-score differences exceeding ~0.42 units"? One of the two
statements must go.

**Q4 (F3).** Three artefacts give three different "largest observed reversal
difference" values: §3.6 says ≈0.27, the T5 JSON says "<0.2", and the
deposited table's maximum is 0.3073. Which is intended?

**Q5 (F2).** Was the index-pairing for shared-control contrasts chosen
deliberately, and if so on what theoretical grounds? My reading is that the
exchange-8-case-label null is the correct one because the control cancels from
Δ; I would like to know whether you intended a different randomisation.

**Q6 (F4).** Why was `alternative="greater"` with targets in the first
argument chosen for T1? The data reject that direction
(p = 9.6 × 10⁻⁴ one-sided the other way, 1.9 × 10⁻³ two-sided). A two-sided
test would have supported a *stronger* version of your criticism.

**Q7 (F4).** Do you accept the target/control row-level MWU at all, given the
22 modules share 4 contrasts and 6 (algebraically 3) aetiology differences?
A module-label permutation gives p = 0.033; a mixed model would be another
route.

**Q8 (F5).** Was any base rate considered when interpreting 28/64 controls as
"flags targets and controls equally"? The base rate is 5/22 = 22.7%, so 43.8%
is enrichment, two-sided binomial p = 2 × 10⁻⁴.

**Q9 (F7).** What justified the median-SE-of-stable-genes yardstick for
"gene-level DE was well powered", given the design detected 2,220 genes at
q < 0.05 with median |logFC| = 0.045, far above the derived 0.001 threshold?

**Q10 (F7).** Was the one-sided myeloid p = 0.057 pre-specified before the
direction (Con > B16) was observed? If not, please drop it from §3.8.

**Q11 (F7).** Was the denominator 7 in the pooled SD
(`33_negative_result_quantification.py:101-103`) intended? It is neither the
pooled SD (6) nor the within-group SD.

**Q12 (F6).** Is there an inferential object you intended for the module mean
moderated *t*, separate from `q_comp`? If not, please relabel the Table 1
column so it cannot be read as a p-value for the mean.

**Q13 (F9).** Should the §3.7 lowess sweep be re-run on the 13,893-gene
bone-only universe so that it actually brackets the reported 235?

**Q14 (F10).** Was the Neuron.NF Ifit2 result (q = 0.0412) excluded
deliberately? It is the smallest q among the three significant compartment
ISG changes.

**Q15 (housekeeping).** `etiology_summary.json` records
`n_sig_randset_q05_percontrast = 68`, but the deposited CSV gives 64 (32 KO +
32 OE) and the same JSON's `_indep` field and `randset_sig_breakdown` both say
32. Is the 68 a stale field?

**Q16 (unverifiable).** §3.8 reports concordance with Walker-256 at day 1
(ρ = 0.112) and day 14 (ρ = 0.030), but `sc_GSE292522_summary.json` contains
only the day-7 value (ρ = 0.09845, n = 12,158, p = 1.42 × 10⁻²⁷). Where are
the day-1 and day-14 concordance statistics computed, and can they be added to
the deposited summary?

**Q17 (unverifiable).** Table 4's headline "Spearman ρ = 0.991 (n = 132)"
between network and additive reversal cannot be recomputed from
`etiology_bootstrap_differences.csv`, which carries `rev_A`/`rev_B` but no
additive columns. The value appears only in `head_to_head_baselines.json`.
Should the per-condition network and additive reversal vectors be deposited
so the paper's central "the graph is decoration" claim is auditable?

---

# PART 4 — § WHAT I ACTUALLY CHECKED

## Files read
`docs/manuscript/manuscript_final.md` (full);
`docs/review/round2_2026-10-10/_PANEL_BRIEF.md`;
`scripts/22_etiology_comparison.py`, `scripts/_limma.py`,
`scripts/33_negative_result_quantification.py`,
`scripts/32_percelltype_sc.py`, `scripts/20_multistage_pipeline.py`
(lines 370–500), `scripts/23_sensitivity.py` (lines 75–200);
`results/tables/`: `etiology_marker_cross.csv`,
`etiology_marker_de.csv`, `etiology_module_logfc_cross.csv`,
`etiology_bootstrap_differences.csv`,
`etiology_reversal_point_estimates.csv`, `etiology_summary.json`,
`negative_result_quantification.json`, `module_level_signal.csv`,
`perturbation_bone.csv`, `head_to_head_baselines.json`,
`sensitivity_summary.json`, `sensitivity_limma_trend.csv`,
`qc_muscle_scores.csv`, `DE_D1_vs_Sham_adj.csv`, `DE_D7_vs_Sham_adj.csv`,
`DE_D14_vs_Sham_adj.csv`, `sc_GSE292522_pseudobulk_DE.csv`,
`sc_GSE292522_composition.csv`, `sc_GSE292522_celltype_DE.csv`,
`sc_GSE292522_celltype_scores.csv`, `sc_GSE292522_summary.json`.
Raw data: all 16 GSE149648 and all 12 GSE160543 member files.

I did **not** read any Round-1 file, `R3_editorial_review.md`, any
`manuscript_v0.*`, or anything under `.workbuddy/`, and no other reviewer's
output.

## Scripts I wrote and ran
All in `review_A2/`, all independent of the project's own primitives except
where explicitly stated:

| script | what it does |
|---|---|
| `A2_check1_permutation.py` | Own GEO parser + normalisation + split enumeration; Cd68/Aif1/Gal/Ccl2/Il1b/Calca/Scn10a/Atf3/Gfap over all 6 pairs; exhaustiveness proof; floor derivation; repo vs exchange-8 nulls for shared controls; BH recomputation |
| `A2_check2_power_and_tests.py` | KO/OE mirror over 132 cells; 264→132 dependence structure (correlations, eigenvalues, effective rank, algebraic redundancy); MDE arithmetic; T1/T2/T3/T4 recomputation |
| `A2_check3_corrections.py` | Corrected myeloid pooled SD; paired-vs-unpaired MDE structure; 64-vs-68 discrepancy; §3.6 quoted numbers; BH family audit across all six families; Table 2 max-selection audit |
| `A2_check4_singlecell.py` | Genome-wide BH check; ISG list extracted from source and tested (38 genes); per-compartment ISG q and ISG-score MWU; Table 5 row-by-row; U statistic |
| `A2_check5_reversal_from_scratch.py` | **Complete independent rebuild** of the network (kNN, no self-loops), diffusion operator (closed form), and reversal cosine for all 22 modules × 4 contrasts × 2 perturbations from raw FPKM; leave-one-animal-out with network + axis rebuilt each drop |
| `A2_check6_discovery_and_T1.py` | Table 1 `q_comp` semantics and per-gene t values; §3.2 counts; muscle-score table; T1 in all four tail directions; module-label permutation; Cliff's δ and bootstrap CI |

## Values recomputed vs the manuscript

| manuscript value | my recomputation | verdict |
|---|---|---|
| Cd68 Δ = +3.42, q = 0.012 | +3.418293, q = 0.012487 (from raw FPKM) | **match** |
| Aif1 Δ = +3.08, q = 0.012 | +3.075712, q = 0.012487 | **match** |
| Gal Δ = −2.59, q = 0.012 | −2.592142, q = 0.012487 | **match** |
| Ccl2 Δ = +1.39, q = 0.017 | +1.390698, q = 0.017343 | **match** |
| p floor 1/4,901 | attainable floor 3/4,901 = 6.12 × 10⁻⁴ | **discrepancy** |
| p floor 1/71 | attainable floor 3/71 = 0.0423 | **discrepancy** |
| 70×70 = 4,900 exhaustive | 70 rows, 4 True each, all distinct, complement involution free | **match** |
| KO/OE p identical to 8 dp | max difference exactly 0 | **match (stronger)** |
| reversal scores (Table 4, §3.4) | max deviation from deposited 9.52 × 10⁻¹⁵ | **match** |
| 264 → 132 hypotheses | 528 rows deposited = 264/scope; 132 distinct per scope | **match, wording imprecise** |
| min q_permutation = 0.165 | 0.1654, BH reproduced | **match** |
| glial/immune p = 0.0014, q = 0.165 | 0.00143, 0.1654 | **match** |
| erythrocyte \|z\| = 7.39, −5.78, pan-neuronal −3.74 | 7.3878, −5.7822, −3.74 | **match** |
| 64/264 randset-sig, 32/132, 14 controls | 64 = 32+32; 32 = 18 target + 14 control | **match** |
| T1 MWU p = 0.497 / 1.00 | 0.4973 / 0.9991 as reported; two-sided \|z\| p = 1.9 × 10⁻³ | **discrepancy (sidedness)** |
| target median pct = 55.7 | 55.7 (control median 24.1, unreported) | **match, incomplete** |
| MDE d = 5.65 (2 vs 2) | 5.6535 arithmetically; but no such contrast exists | **wrong design** |
| MDE ≈ 0.42 reversal units | 0.044 (unpaired) / 0.039 (paired) on the paper's own SE | **discrepancy** |
| largest reversal Δ ≈ 0.27 | 0.2743 targets / 0.3073 all modules; T5 JSON says "<0.2" | **internally inconsistent** |
| myeloid d = 1.37 | 1.3685 with SS/7; 1.2670 conventional | **discrepancy (8%)** |
| myeloid MDE d = 2.38 | 2.3808 | **match** |
| myeloid MWU p = 0.114 / 0.057, U = 2 | 0.1143 / 0.0571, U = 2.0; "U = 2/14" denominator wrong (max 16) | **match, notation error** |
| 55% myeloid drop | 55.8% (ratio 0.4423) | **match** |
| Ifit2 logFC −0.048, q = 0.042 | −0.047582, q = 0.0417 | **match** |
| 1,517 up / 703 down, 55,967 cells | 1,517 / 703 / 55,967 | **match** |
| Gfap q = 0.006, Oas1c q = 0.014 | 0.0061, 0.0136 | **match** |
| "only significant ISG changes: Ifit2 neurons q=0.046, Isg15 glia q=0.007" | plus **Neuron.NF Ifit2 q = 0.0412** | **omission** |
| "no ISG up-regulated" | 0 of 38 ISGs up at q < 0.05 (also 0 at q < 0.10) | **match** |
| D1/D7/D14 q<0.05 = 0 / 4,600 / 235 | 0 / 4,600 / 235 | **match** |
| q<0.10 = 5 / 5,867 / 860 | 5 / 5,867 / 860 | **match** |
| median \|logFC\| = 0.073 / 0.158 / 0.119 | 0.0731 / 0.1576 / 0.1187 | **match** |
| ρ robustness min 0.959, median 0.988, Δ min 0.989 | 0.9593 / 0.9876 / 0.9887 | **match** |
| top five identical in all 60 settings | overlap 5/5 over 75 grid points | **match** |
| lowess sweep 316–325 genes (protecting 235) | 316–325 on the 13,601-gene universe; 235 is on the 13,893-gene universe | **discrepancy (unreconciled)** |
| ρ(net, add) = 0.991, n = 132 | 0.991145 in `head_to_head_baselines.json`; not recomputable from the CSV | **unverifiable** |
| mean \|rev\| 0.057 / \|add\| 0.033 / \|ridge\| 0.712 | 0.0568 / 0.0328 / 0.7119 | **match** |
| GSM4506604 at 6.26% muscle FPKM | 6.257%, next highest 0.160%, robust z 120.9 | **match** |
| concordance ρ = 0.098 (D7), 0.112 (D1), 0.030 (D14) | 0.09845 deposited; D1/D14 not in any artefact | **partly unverifiable** |

## Errors in my own work, corrected

For transparency, since one of them affected a number I quote:

1. My first pass at the myeloid pooled SD printed *d* = 3.10 because I wrote
   `sqrt(SS_total)/6` instead of `sqrt(SS_total/6)`. Corrected in
   `A2_check3_corrections.py`: the conventional pooled SD is 5.149 × 10⁻⁴ and
   *d* = **1.2670** (not 3.10). The manuscript's 1.37 therefore *overstates*
   the conventional *d* by 8%, not understates it.
2. My first reversal rebuild used only 6 of the 22 modules, giving a 3,015-gene
   node set instead of 3,038 and a mismatched `depk` index (both network
   scopes). After fixing to all 22 modules and filtering to
   `network_scope == "per_contrast"`, agreement is 9.52 × 10⁻¹⁵.
3. An early version of the module-label permutation indexed a 22-element
   boolean array against a 132-element vector. Corrected to permute the
   22 module→class labels as blocks.

## What I could not verify

- Table 4's ρ = 0.991 (n = 132): `etiology_bootstrap_differences.csv` has no
  additive-model columns, so the correlation cannot be rebuilt from deposited
  artefacts (Q17).
- The day-1 and day-14 cross-model concordance coefficients (Q16).
- Whether the one-sided myeloid p was pre-specified (Q10).
- The `n_sig_randset_q05_percontrast = 68` field in `etiology_summary.json`
  (Q15) — the CSV and the JSON's other fields both say 64/32.
- Whether the shared-control index-pairing and the T1 tail choice were
  deliberate (Q5, Q6).

## Overall assessment

The permutation machinery is sound and the positive results reproduce exactly
from raw counts; I found no fabrication and no arithmetic error in any
headline marker number. The problems are concentrated in the inferential
layer added this round, and they run in an unusual direction: **F3, F4 and F5
each make the paper's own negative case weaker than its own data support.**
A power analysis that is wrong by a factor of ~10 in the MDE, a one-sided test
aimed at the direction the data reject, and an enrichment described as
"equality" are three instances of a systematic pattern — the calibration
results are being reported in the most self-effacing way the numbers permit.
Correcting all three would *strengthen* the paper's methodological thesis
while removing the claims that do not survive checking. That is a
substantially better paper, and it is reachable with the analyses already
deposited.
