# Reviewer B2 — Design-layer review (statistics / causal inference / reproducible research)

**Manuscript.** "Divergent dorsal root ganglion programs in bone cancer pain and
chemotherapy-induced peripheral neuropathy: an exact cross-study reanalysis and
a negative-control calibration of in-silico perturbation."

**Reviewer role.** Design-layer: statistician / causal-inference / reproducible
research. This is a fresh, independent read of the revised single-author
manuscript; no prior review round was available to me, so I treat it as a first
submission and assume it may contain flaws.

**Charge.** Verify the three claimed revisions — (i) §2.3 dual-logFC clarification,
(ii) §3.6 "132 tests per scope" + T2 "descriptive gauge" softening, (iii) §5 unified
MDE "≈0.052" and "6.1×10⁻⁴" — and audit the exact-permutation validity, the
multiplicity handling, and the T1–T5 negative-result quantification. Every
quantitative claim below was independently recomputed from the deposited CSV/JSON
artefacts and, where relevant, from the project's own `exact_two_sided_p` and the
`22_*` / `33_*` scripts.

**Headline.** The arithmetic of the revision is sound: every quoted number in
§2.3, §3.5, §3.6, §5 and Tables 3–5 recomputes to the stated precision, and the
exact-permutation floors are real. The problems I found are *interpretive*, not
numerical: one incorrect mechanism label in §2.3, one under-played dependence in
T2, one over-stated gene-level-power exemplar in T3, and one frame that is
narrower than the conclusion it is used to support in T5.

---

## Finding 1 — §2.3 dual-logFC: the numbers are correct, but the claimed
mechanism ("moderation shrinkage") is wrong

【Problem】 The revised §2.3 attributes the gap between the muscle-adjusted
limma-trend logFC and the unadjusted per-sample mean logFC (Gfap 1.33 vs 1.50,
Ccl2 0.81 vs 0.90) to "moderation shrinkage," yet empirical-Bayes moderation
shrinks the *variance* (and therefore the moderated t and SE), not the logFC
*point estimate* — the gap is produced by the muscle covariate, not by moderation.

【Evidence】 I recomputed both estimators directly from the deposited artefacts.
The muscle-adjusted limma-trend coefficient (Bone_D7, `etiology_marker_de.csv`,
produced by `de_positive_controls` which calls `fit_limma` with the 20-gene muscle
score as a covariate) is: Gfap = 1.328, Ccl2 = 0.806, Cd68 = 2.858, Aif1 = 2.180.
The unadjusted per-sample mean log₂ fold change (Bone_D7 column of
`etiology_marker_cross.csv`, which is `mean(case) − mean(control)` over the four
samples with no covariate) is: Gfap = 1.498, Ccl2 = 0.897, Cd68 = 2.839,
Aif1 = 2.192. The four pairwise gaps (1.33/1.50; 0.81/0.90; 2.86/2.84; 2.18/2.19)
all recompute exactly as the manuscript claims, so the *separation of the two
columns* is correct and useful.

However, `fit_limma` (`_limma.py:76–124`) computes the coefficient `logFC` as the
ordinary least-squares solution of the design matrix and applies empirical-Bayes
moderation only to `s2_post` (the posterior variance), the standard error, and the
moderated t. There is no shrinkage of the logFC coefficient itself — limma does
not pull fold changes toward zero. The 0.17-unit Gfap gap and 0.09-unit Ccl2 gap
are therefore the consequence of *adding the muscle-score covariate to the design*,
which re-estimates the group coefficient conditional on muscle contamination, and
of comparing a model coefficient to a raw per-sample mean. They are not
"moderation shrinkage."

The direction makes the mislabel obvious. Moderation, even if it touched the
coefficient, would not systematically *reduce* a low-precision gene's estimated
logFC by 0.17 while leaving a high-precision gene's estimate within 0.02; it only
changes the SE. Here Gfap and Ccl2 (lower expression / noisier) move *down* by
~0.1–0.2 when the muscle covariate is added, which is exactly what a covariate that
removes a shared compositional component does — not what variance moderation does.

【Why it matters】 The dual-logFC clarification is the manuscript's direct answer to
a legitimate prior criticism (two different "logFC" values appear in the same
paper, one in Table 1/§3.3 and one in Table 3, both on the same 13,601 genes).
The separation is correct and worth keeping. But if the *mechanism* is mis-stated,
a careful reader who opens `fit_limma` will conclude the author does not
understand their own pipeline, and the otherwise-correct distinction between
"muscle-adjusted model coefficient" and "raw per-sample mean" loses its anchor.
It also repeats a widespread textbook error — that limma "shrinks log fold
changes" — which is worth not perpetuating in a methods-heavy paper.

【Specific fix】 In §2.3, replace the sentence

> "they differ in estimator, not only gene set: the within-study programme
> statements of Section 3.3 report the muscle-adjusted limma-trend *moderated*
> log₂ fold change, whereas the cross-aetiology differences of Table 3 report the
> *unadjusted per-sample mean* log₂ fold change … and can differ … explained by
> moderation shrinkage"

with

> "the two differ because the §3.3 value is the *coefficient* of the
> muscle-adjusted linear model (limma-trend, with empirical-Bayes moderation of
> the standard error only, not of the fold change), whereas the Table 3 value is
> the raw per-sample mean log₂ fold change with no covariate; the gap reflects the
> muscle-covariate re-estimation of the group coefficient, not moderation of the
> fold change. For well-expressed markers (Cd68, Aif1) the covariate moves the
> estimate by only ~0.02 log₂ units; for lower-precision markers (Gfap, Ccl2) it
> moves it by ~0.1–0.2 log₂ units."

---

## Finding 2 — Exact-permutation floors hold, and the observed statistics lie
inside the enumerated nulls

【Problem】 The claimed permutation floors — 3/4,901 for the four unpaired
cross-aetiology contrasts, 3/71 for the two shared-control marker contrasts, and
2/71 for the same shared-control *reversal* contrasts — must be confirmed, and the
observed statistics must be confirmed to sit at or above those floors (i.e. inside
the null), because exact permutation is the statistical backbone of the paper's
two secure claims.

【Evidence】 I rebuilt `exact_two_sided_p` (`_limma.py:420–426`:
`ge = Σ |null| ≥ |obs|; p = (ge + 1)/(n + 1)`) and checked it against the
deposited tables.

- **Marker, unpaired** (`etiology_marker_cross.csv`, `paired_control = False`,
  n_null = 4,900): the minimum observed p over all markers is
  0.0006121199… = 3/4,901 exactly. This appears for Cd68, Aif1, Gal, Ccl2, Il1b
  etc. (the strongest markers).
- **Marker, paired** (`paired_control = True`, n_null = 71): minimum observed p =
  0.0422535211… = 3/71 exactly (Atf3/PTX|OXA, Gal/PTX|OXA, Calca/PTX|OXA,
  Cd68/PTX|OXA, Nefl/PTX|OXA).
- **Reversal, paired** (`etiology_bootstrap_differences.csv`, `paired_control =
  True`, n_null = 71; 88 such rows): minimum observed p = 0.0281690141… = 2/71
  exactly, and *every* paired-reversal p is ≥ 2/71. No paired-reversal test falls
  below the floor.
- **Reversal, unpaired** (n_null = 4,900): by the same antisymmetry argument as
  the logFC endpoint, obs and −obs are both enumerated, so its floor is also
  3/4,901; the smallest observed unpaired-reversal p is 0.001428, well above the
  floor.

I also reconstructed the counting rule on synthetic nulls: when obs and −obs are
both present, `exact_two_sided_p(obs)` returns 3/(n+1); when only obs is present
it returns 2/(n+1). This matches the manuscript's prose explanation in §2.7
("the logFC difference is antisymmetric, so the complement split is an exact
mirror … whereas the reversal projection at the complement split is not the exact
negation … so only the observed split itself is as extreme"). The paired branch in
`22_etiology_comparison.py:283–329` (the exchange-8-case null with the shared
control fixed) implements exactly this, and the self-checks at lines 313–324
confirm the observed split reproduces the point estimates.

【Why it matters】 If the floors were wrong, or if an observed statistic fell
*below* its floor (outside the enumerated null), the "exact, no Monte-Carlo
error, hard floor" selling point would be false and the p-values overstated.
Neither happens. The construction is valid: the cross-aetiology null genuinely
enumerates the entire label space, the floors are real, and the observed
statistics are contained in their own nulls. This is the part of the revision
that most fully delivers.

【Specific fix】 No correction to the floors or to the §2.7 / §5 statements.
One optional editorial tightening: state explicitly in §2.7 that the unpaired
*reversal* contrast also has floor 3/4,901 (obs and −obs both enumerated), so the
reader is not left to infer it from the logFC case. This is editorial, not
substantive, and does not change any number.

---

## Finding 3 — T1 anti-selectivity (MWU, median |z| 2.47 vs 1.07, p = 0.0019)
recomputes and is not driven solely by interferon

【Problem】 The headline anti-selectivity claim — negative-control modules carry
larger standardised effects than analgesic-target modules — must recompute, and
the "removing interferon still leaves controls larger" check (p = 0.019) must also
recompute, because both are used to argue the failure is about validity rather
than one rogue module.

【Evidence】 From `etiology_bootstrap_differences.csv`, the per_contrast ×
knockout deduplicated universe (132 tests): target n = 102, control n = 30.
Two-sided Mann–Whitney U on |z_vs_randset|:

- full universe: p = 0.0019139 (manuscript 0.0019); median |z| control =
  2.4702, target = 1.0729 (manuscript 2.47 vs 1.07).
- after dropping the interferon control module (24 control tests remain):
  p = 0.0186968 (manuscript 0.019); median |z| control = 2.3327 vs 1.0729 for
  targets.

Both recompute to the quoted precision. The direction is robust: even with the
single largest |z| (erythrocyte) and the interferon module removed, controls still
sit above targets. The magnitude also barely moves (2.47 → 2.33), so the result is
not carried by one extreme control.

【Why it matters】 T1 is the strongest, cleanest argument in §3.6. It uses a
standard, assumption-light rank test directly on the 132 effect sizes and does not
depend on the degree-matched random-gene-set null (unlike T2). That it survives
removal of the dominant interferon signal means the conclusion is not a
single-module artefact and genuinely reflects the layer ranking composition/
topology controls above pharmacology targets. This finding stands and supports
the revision.

【Specific fix】 No numerical correction. Two optional additions: (a) report the
MWU U statistic and the Hodges–Lehmann shift of |z| alongside the medians, so the
*magnitude* of the control–target gap is visible and not just its significance;
(b) state explicitly that the 132 observations are clustered within 6 aetiology
pairs, so the MWU p is approximate (the manuscript already flags this for the
exact-permutation p, but should say the same for T1).

---

## Finding 4 — T2 binomial p = 0.0067 recomputes, but the "descriptive gauge"
softening is insufficient: the test is correlated, module-labelled, and partly
circular

【Problem】 The manuscript softened T2 to "an approximate descriptive gauge rather
than a formal test" yet still prints the exact binomial p = 0.0067 in the body and
uses it as a second, independent-looking piece of evidence that the degree-matched
null "confirms red blood cells at nearly twice their expected rate." The
binomial's independence and base-rate assumptions are both violated, so presenting
it next to T1 overstates how many independent confirmations of anti-selectivity
exist.

【Evidence】 Recomputed from `etiology_bootstrap_differences.csv` (per_contrast ×
knockout, `q_randset < 0.05`): randset-significant total = 32, of which 14 are
controls and 18 are targets; control base rate in the 132-test universe =
30/132 = 0.2273; `binomtest(14, 32, 0.2273, alternative = "greater")` =
0.0067079 (manuscript 0.0067). The number is correct.

Two design problems remain despite the softening:

1. **Non-independence.** The 32 "significant" trials are not independent
   Bernoulli draws. They are 22 modules × 6 aetiology pairs that share contrasts
   and — for the degree-matched null specifically — share the same 500
   degree-matched gene-set draws. Clustering inflates the effective sample size
   relative to the binomial's n = 32, so the p-value is anti-conservative.

2. **Base rate is not the true per-trial success probability.** The binomial
   tests H₀: each trial's success probability = 30/132. But the success
   probability is *structurally higher for control modules* — erythrocyte,
   pan-neuronal and myelin effects are composition- and topology-driven and are
   therefore almost always "significant" against a degree-matched null. The
   module *kind* is a fixed label, not a random draw, so the trial success
   probability is not homogeneous. Testing "are controls over-represented among
   significant results?" against a base rate of 30/132 is therefore partly
   circular: H₀ is already implausible by construction, and the p = 0.0067
   measures an effect the design guarantees rather than one it discovered.

The manuscript's new sentence does acknowledge both issues ("the per-trial success
probability depends on module kind and the tests share contrasts") — but it still
lands the exact p = 0.0067 as a confirmatory number in the same paragraph as T1,
which lets a reader count T1 + T2 as two independent, formally-significant
replications of anti-selectivity. They are not: T1 says control |z| are larger
(rank test); T2 says controls are over-represented among the "significant" subset
of that *same* |z| distribution. Same phenomenon, two views.

【Why it matters】 This is the one revision I would push further. The softening is
necessary and honest, but insufficient: as written, T2 still functions as a second
significance result. Because its assumptions are violated in the conservative
direction (dependence + module-fixed success probability), the p is not a valid
frequentist quantity and should not be presented on equal footing with T1's MWU.

【Specific fix】 Either (a) demote T2 explicitly — e.g.

> "As a *non-independent* illustration of the same effect, 14/32 randset-
> significant tests are control modules versus a 30/132 base rate (binomial
> p = 0.0067, a descriptive gauge only; this is the same signal captured by T1,
> not an independent test)"

or (b) replace the binomial with a cluster-robust assessment: permute the
target/control *label* within each aetiology pair (or block-bootstrap the 6 pairs)
and report the permutation p for "controls over-represented among significant."
That respects the within-contrast correlation and removes the circular base-rate
assumption. Do not present T1 and T2 as two separate significant results.

---

## Finding 5 — T3 single-cell myeloid power is defensible ("underpowered, not
equality"), but the gene-level "well powered" exemplar overstates typical power

【Problem】 The myeloid-fraction Cohen's d and Mann–Whitney U recompute, and the
"inconclusive null, not evidence of equality" reading is correct — but the
companion claim that gene-level DE is "well powered (detected Ifit2 at logFC
−0.048, q = 0.042)" uses a near-zero-effect, near-zero-SE gene as the exemplar and
overstates what the n = 4 design typically resolves.

【Evidence】 From `sc_GSE292522_composition.csv`: Con myeloid fractions =
{0.00046, 0.00137, 0.00199, 0.00086} (mean 0.00117); B16F10 =
{0.00081, 0.00040, 0.00071, 0.00015} (mean 0.00052); B16/Con = 0.44 (≈56% drop).
Pooled SD with denominator n₁ + n₂ − 2 = 6 gives observed Cohen's d = 1.267
(manuscript 1.27). MDE in Cohen's d at 80% power, two-sided α = 0.05, n = 4 per
group = 2.38 (recomputed). Mann–Whitney U on the four-vs-four fractions gives
U = 2/14 (scipy reports the complementary statistic 14; the pair is 2 and 14),
two-sided p = 0.114, and the directional (Con > B16) one-sided p = 0.057 — all as
quoted. So d = 1.27 < MDE 2.38 ⇒ underpowered, while the point estimate is a large
56% drop ⇒ "not evidence of equality" is the honest reading. This part stands.

Separately, `negative_result_quantification.json` reports the gene-level MDE in
|logFC| = 2.38 × median_SE (0.2048) = 0.488. Ifit2 was called significant at logFC
−0.048, which is far below 0.488 — i.e. Ifit2's standard error is anomalously
small (a low-count gene), so it is a poor exemplar of "typical" gene-level power.
The substantive gene-level statement should be: gene-level DE resolves effects
down to ~0.49 log₂ units at the median SE, and resolves several ISGs; a handful of
low-SE genes are significant at much smaller effects. Citing Ifit2 at −0.048 as
proof the test is "well powered" invites the reply that the design only detects
effects whose SE happens to be tiny.

【Why it matters】 The myeloid conclusion ("inconclusive, not null") is correct and
well-stated — that is the right honest reading of a 56% drop that misses
significance at n = 4. The only blemish is the Ifit2 exemplar, which makes the
gene-level power claim look stronger and more uniform than it is. A reviewer or
reader who checks the median-SE MDE (0.49) against the Ifit2 call (−0.048) will
see the exemplar is a special case.

【Specific fix】 Reword the §3.8 sentence

> "Gene-level differential expression was well powered in this design — it
> detected Ifit2 at log₂ fold change −0.048 (q = 0.042)"

to

> "Gene-level DE at q < 0.05 resolves effects down to ~0.49 log₂ units at the
> median moderated SE (n = 4 per group); a few low-variance genes such as Ifit2
> (log₂ −0.048, q = 0.042) are significant because their SE is unusually small,
> so the myeloid-marker null (Cd68/Aif1 unchanged) is a genuine null for *large*
> effects, not a generic power claim for the whole transcriptome."

---

## Finding 6 — T5 reversal MDE (≈0.052, z_sum = 2.802) recomputes, but the
per-test one-sample frame is narrower than the "power floor" conclusion draws

【Problem】 The unified MDE "≈0.052" (was ≈0.05) and "6.1×10⁻⁴" revisions recompute,
and the per-test one-sample MDE is a legitimate computation — but it measures
per-module detectability of Δ = 0 at the null, which is a different question from
the paper's actual validity claim (can the layer discriminate targets from
controls?), so the "power floor" should not be read as the layer's power to rank.

【Evidence】 From `etiology_bootstrap_differences.csv` (per_contrast × knockout,
132 tests): bootstrap SE of Δ = (ci95_hi − ci95_lo)/(2 × 1.96); median SE =
0.01850 (manuscript 0.0185/0.019). Normal-approx 80%-power reference
z_sum = 1.959964 + 0.841621 = 2.801585 (manuscript 2.802). MDE = 2.801585 ×
0.01850 = 0.05184 (manuscript ≈0.052). Largest target Δ = 0.2743 (glial/immune,
Bone_D7 vs PTX), exceeds MDE by 5.29×; median target |Δ| = 0.0313 < MDE;
n_target exceeding MDE = 41/102; n exact-permutation p < 0.05 = 31/132 (manuscript
"31 of the 132"). All recompute.

The one-sample MDE = (z_{α/2} + z_{β}) × SE(Δ̂) is the correct dispersion for a
single contrast's Δ = 0 test, and KO/OE being exact mirrors (r = −1.000) genuinely
means there is no second arm, so a one-sample frame is appropriate. The arithmetic
is sound and "resolves large effects, not the median target" is defensible.

The only risk is over-reading. T5's MDE says the study can detect a per-module Δ of
~0.05 *if it exists*; it does **not** say the study can tell targets from controls
— that is T1 (controls have *larger* |z|, p = 0.0019). A reader could conflate
"median target below MDE" (an effect-size detectability floor) with "the layer
cannot rank targets" (a discrimination failure). The stronger, separate statement
is that the layer ranks controls *above* targets, which is a validity failure of a
different kind from a power failure.

【Why it matters】 Numbers are fine; the framing is the only issue. If T5 is read
as the layer's "power to rank," it under-sells the actual critique (anti-
selectivity) and over-states what the MDE measures (per-module effect-size
detectability). Keeping the two distinct protects the manuscript's main negative
conclusion.

【Specific fix】 In §3.6 / §5, append one clause to the MDE sentence, e.g.

> "this is a per-module power floor at the Δ = 0 boundary and speaks to effect-
> size detectability, not to target-vs-control discrimination, which T1 addresses
> directly (and negatively)."

No number changes are needed; the §5 "≈0.052" and "6.1×10⁻⁴" are correct.

---

## Finding 7 — Multiplicity dedup (264 → 132) and BH per scope are correct;
min q = 0.165 stands, but the 132 are themselves correlated

【Problem】 The claimed revision "132 tests per network scope (correlated within
contrast, so the effective number of independent tests is smaller)" and the
"min q = 0.165 for exact permutation" must be confirmed, and the KO/OE dedup logic
in `main()` must be checked for double-counting or dropped hypotheses.

【Evidence】 `etiology_summary.json`: n_independent_tests_per_scope = 132,
n_independent_tests_total = 264, min_q_permutation_per_contrast = 0.165447;
`negative_result_quantification.json` T5: n_independent_tests = 132,
n_exact_permutation_p_lt_0.05 = 31. KO and OE are exact mirrors
(reversal(OE) = −reversal(KO)), so `exact_two_sided_p` is identical to eight
decimals; `22_etiology_comparison.py:687–703` keeps the knockout subset, applies
BH within each scope over the 132 p-values, and merges the q back onto the OE twin
— no hypothesis is dropped and none is double-counted. The randset q is handled
the same way (lines 691–692). Because the 132 tests are positively dependent within
contrast, BH under PRDS remains valid (and is conservative), so q = 0.165 is an
honest floor and "no perturbation Δ reaches q < 0.05" is correct. The "effective
number of independent tests is smaller" caveat is present and appropriate.

【Why it matters】 This is the multiplicity handling a prior round presumably
criticised; it is now correct and the manuscript no longer inflates the denominator
by counting KO and OE as two independent hypotheses. The only residual point is
that BH is applied at nominal n = 132 while the effective count is smaller — under
positive dependence this makes q *conservative*, which is the safe direction, so no
action is required. The finding stands.

【Specific fix】 No change required. Optionally, report a correlation-adjusted
effective-number estimate (e.g. a permutation-based effective n, or Gao–Becker)
alongside 132 for completeness, but the current BH is defensible and the conclusion
holds.

---

## § Stands up

1. **Exact-permutation validity (Finding 2).** The floors 3/4,901, 3/71 and 2/71
   all recompute exactly from the deposited artefacts and from `exact_two_sided_p`;
   every observed statistic lies at or above its floor. The "exact, no Monte-Carlo
   error, hard floor" claim is genuine and is the paper's strongest methodological
   asset.

2. **T1 anti-selectivity and its interferon-robustness (Finding 3).** MWU p =
   0.0019 and median |z| 2.47 vs 1.07 recompute; removal of interferon still gives
   p = 0.019 (median |z| 2.33 vs 1.07). This is the paper's strongest,
   assumption-light result and it survives the obvious stress test. It stands.

3. **Multiplicity dedup and BH (Finding 7).** The 264 → 132 dedup is correctly
   implemented (no dropped or double-counted hypotheses), BH per scope is valid
   under positive dependence, and min q = 0.165 stands — "no perturbation Δ
   survives correction" is correct. It stands.

4. **T3 inconclusive-null framing for the myeloid fraction (Finding 5).** d = 1.27,
   MDE d = 2.38, MWU 2/14, p = 0.114 all recompute; the "underpowered, not evidence
   of equality" reading is the honest one and is well-stated. It stands (only the
   Ifit2 exemplar needs rewording).

5. **T5 MDE arithmetic (Finding 6).** 0.0518 and z_sum = 2.8016 recompute; the
   per-test one-sample frame is appropriate given the exact KO/OE mirror, and the
   "resolves large effects, not the median target" reading is defensible. It stands
   (with the one clarifying clause suggested).

---

## § Questions for the authors

1. In §2.3 you now separate "muscle-adjusted limma-trend moderated logFC"
   (Table 1 / §3.3) from "unadjusted per-sample mean logFC" (Table 3). Do you
   agree the gap is the muscle *covariate* and not moderation, and will you correct
   the "moderation shrinkage" wording (Finding 1)?

2. For T2 (Finding 4): do you accept that the binomial p = 0.0067 is not an
   independent confirmation of T1, and will you either demote it to an explicitly
   non-independent illustration or replace it with a cluster-permutation that
   respects the within-contrast correlation and the module-fixed success
   probability?

3. For T5 (Finding 6): can you add one sentence making clear the per-test MDE is an
   effect-size detectability floor and not a measure of target-vs-control
   discrimination (which T1 negates)?

4. The unpaired *reversal* contrast has floor 3/4,901 by the same obs-and-−obs
   argument as the unpaired logFC — will you state this explicitly in §2.7 so the
   reader does not assume only the logFC endpoint carries that floor?

5. In T3 (Finding 5), will you replace the Ifit2-at-−0.048 "well powered" exemplar
   with the median-SE MDE (~0.49 log₂ units) so the gene-level power claim is not
   carried by a low-SE special case?

6. The exact-permutation p for the cross-aetiology marker/module differences
   conditions on the observed lab/batch labels (different laboratories, different
   dissection/perfusion protocols; the erythrocyte-control Δ of −2.48 is itself a
   direct measure of this). You acknowledge this as a limitation. Is there any
   sensitivity you can report — e.g. a regression on a batch proxy, or a leave-lab-
   out check — that qualitatively bounds how much of the Cd68/Aif1 Δ could be batch
   rather than biology, even if only descriptively?

7. For T1 (Finding 3) and the §3.6 anti-selectivity claim more broadly: the 132
   effects are clustered within 6 aetiology pairs. Will you report the effective
   number of independent contrasts (or a cluster-robust p) so the "p = 0.0019"
   MWU is not read as 132 fully independent observations?

---

## § What I actually checked

**Manuscript / methods read in full.** `docs/manuscript/manuscript_final.md`
(Abstract through References, including Tables 1–5 and Figures 1–3 captions);
`scripts/22_etiology_comparison.py`; `scripts/_limma.py`;
`scripts/33_negative_result_quantification.py`. I did not open any
`docs/review/round*` file, `author_verification_statement.md`,
`GITHUB_DEPOSIT_SOP.md`, or `manuscript_v0.*.md`, per the review rules.

**Artefacts recomputed against (values stated to deposited precision).**

- `results/tables/etiology_marker_cross.csv` — marker unpaired floor
  p = 0.00061212 = 3/4,901; paired floor p = 0.04225352 = 3/71; Cd68/Aif1/Gal/
  Ccl2/Il1b Δ and q as printed.
- `results/tables/etiology_bootstrap_differences.csv` — T1 (MWU p = 0.001914;
  medians 2.4702 / 1.0729), REMOVE-IFN (p = 0.01870; 2.3327 / 1.0729), T2 (32
  significant, 14 controls, base 30/132, binomial p = 0.00671), T5 (median SE
  0.01850, MDE 0.05184, z_sum 2.801585, max target Δ 0.2743, 41/102 exceed MDE,
  31/132 exact p < 0.05), reversal paired floor p = 0.02816901 = 2/71 with
  n_null = 71, unpaired reversal n_null = 4,900.
- `results/tables/etiology_marker_de.csv` — muscle-adjusted Bone_D7 logFC: Gfap
  1.328, Ccl2 0.806, Cd68 2.858, Aif1 2.180 (the §2.3 dual-logFC "adjusted"
  column).
- `results/tables/etiology_summary.json` — n_independent = 132 / 264,
  min_q_permutation = 0.165447, n_sig_randset_q05_percontrast_indep = 32.
- `results/tables/negative_result_quantification.json` — T1 / T2 / T3 / T5 values
  match the body text (T3 myeloid d = 1.267, MDE d = 2.38, Con 0.00117 / B16F10
  0.00052).
- `results/tables/sc_GSE292522_composition.csv` — Con / B16F10 myeloid fractions,
  pooled-SD denominator 6, observed d = 1.267, B16/Con = 0.44.
- `results/tables/sc_GSE292522_pseudobulk_DE.csv` — median SE of expressed genes
  = 0.2048; gene-level MDE ≈ 0.49 log₂ units.

**Logic reconstructed / scripts re-run.** I re-executed the project's own T1–T5
quantification and independently recomputed each of the seven findings above from
the deposited CSV/JSON; I reconstructed `exact_two_sided_p` and the two floor
conventions (obs-and-−obs present ⇒ 3/(n+1); only obs present ⇒ 2/(n+1)) on
synthetic nulls to confirm the counting rule; and I cross-checked the KO/OE dedup
and BH merge in `22_etiology_comparison.py:687–703`.

**Discrepancies found (numerical).** None. Every quoted number in §2.3, §3.5,
§3.6, §5 and Tables 3–5 recomputes to the stated precision; the three claimed
revisions (dual-logFC, "132 tests"/T2 softening, "≈0.052"/"6.1×10⁻⁴") all hold
arithmetically.

**Discrepancies found (interpretive).** Four, all wording/framing on top of
arithmetically correct numbers:

- (i) §2.3 attributes the dual-logFC gap to "moderation shrinkage" when it is the
  muscle covariate (Finding 1).
- (ii) T2's binomial is presented as a near-independent confirmation although it is
  the same signal as T1 and violates independence / base-rate assumptions (Finding 4).
- (iii) T3's Ifit2-at-−0.048 exemplar overstates typical gene-level power (Finding 5).
- (iv) T5's per-test MDE is narrower than the "power floor" conclusion draws and
  should not be read as the layer's power to rank (Finding 6).

No arithmetic error was identified in the revision; the recommended changes are
mechanism-label, dependence-framing, and exemplar-selection corrections that
strengthen an already numerically sound design layer.
