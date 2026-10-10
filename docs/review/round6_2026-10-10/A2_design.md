# Independent peer review — Design layer (Reviewer A2)

**Manuscript:** *Divergent dorsal root ganglion programs in bone cancer pain and
chemotherapy-induced peripheral neuropathy: an exact cross-study reanalysis and a
negative-control calibration of in-silico perturbation*

**Reviewer role:** Design-layer expert — statistician / causal-inference /
reproducible-research methods.

**Status of review:** Fresh, first-submission pass. I had no access to prior review
rounds and did not open any file in the forbidden set.

**Scope of this review:** Only the design-layer questions I was assigned —
exact-permutation validity, multiplicity, the five negative-result
quantifications T1–T5, and the two-logFC-source design. I do not assess
biological interpretation, wet-lab plausibility, or the single-cell annotation
methodology beyond its statistical claims.

---

## 0. What a design reviewer looks for, and how I read this paper

A design-layer review is not a re-run of the authors' arithmetic.
It asks four questions that the authors (understandably) tend to skip:

1. Is each test *valid* — does its null hypothesis match the claim?
2. Is each correction *appropriate* — does it match the dependence structure?
3. Is each number *traceable* and *internally consistent* with its own artefacts?
4. Where a claim is about *statistical validity itself*, is the claim held to a
   higher standard than the work it critiques?

I read `scripts/_limma.py`, `scripts/22_etiology_comparison.py`,
`scripts/33_negative_result_quantification.py`, and `scripts/20_multistage_pipeline.py`
in full, then recomputed every headline number from the deposited CSVs/JSONs
(see §5). The arithmetic reproduces. The concerns below are about questions
(1)–(4), concentrated on the negative-result quantification layer.

---

## 1. Summary verdict

This is a careful, unusually honest computational manuscript.
The exact-permutation machinery is implemented correctly.
The floors it advertises (3/4,901; 3/71; 2/71) hold to the digit.
The KO/OE dedup and BH correction are correct.
The anti-selectivity conclusion (T1) is a genuine and well-executed
methodological finding.

I reproduced every headline number I was asked to check from the deposited
artefacts (see §7), and where I disagree it is on *interpretation* and
*assumption*, not on arithmetic.

My concerns are concentrated on three design-layer decisions that weaken the
*interpretation* of the negative-result quantification, plus one concrete
artefact bug in the T3 power script:

1. **T2's binomial test is statistically incoherent (F1)** because its success
   probability is not constant across trials — it depends on module kind — and
   it is effectively a restatement of T1, so its p-value is not a valid test of
   selectivity.

2. **The "132 independent tests" label overstates independence (F2).** The 132
   (module × contrast) statistics are correlated; T1's MWU and T2's binomial
   silently assume independence, and the effective sample size is never
   quantified.

3. **T3's reported gene-level MDE (0.488 log₂) is computed on a 5-gene high-SE
   tail and is mislabeled "median SE of expressed genes" (F3).** The true
   typical detectable effect is ≈0.0014 log₂; the deposited figure is
   internally inconsistent with the genes the manuscript itself flags as
   detected.

4. **T5's per-test MDE is defensible in frame but fragile (F4)** in its
   dependence on bootstrap-SE symmetry, and the t-based variant uses a
   mis-specified two-sample n=2/n=2 geometry.

None of these overturns the paper's conclusions.
The central claim — in-silico perturbation fails its own negative controls and
its additive baseline — survives all of them.
But items 1–3 should be corrected before publication, because they are claims
about *statistical validity* made by a validity study, and a validity study must
hold itself to a higher standard than the work it critiques.

---

## 2. Findings (four-part contract: Problem / Evidence / Why it matters / Specific fix)

### F1 — T2 binomial selectivity test is incoherent: the success probability is not constant across trials

**【Problem】** The binomial test that "control modules are over-represented among
randset-significant tests (p = 0.0067)" assumes 32 independent Bernoulli trials
each with a constant success (control) probability of 0.227, but significance
is itself strongly correlated with module kind, so the null it tests is the very
effect T1 already demonstrates.

**【Evidence】** `scripts/33_negative_result_quantification.py:67-74` builds
`sig = pc[pc.q_randset < 0.05]`, then calls
`binomtest(k_ctl, len(sig), base_rate, alternative="greater")` with
`base_rate = 30/132 = 0.2273`. I recomputed from
`etiology_bootstrap_differences.csv` (per_contrast, knockout): 32 significant
tests, of which 14 are controls (43.75%). The binomial returns p = 0.0067,
which I reproduced exactly. But T1 (`33_negative_result_quantification.py:47`)
has already shown that control modules carry systematically larger |z| (median
2.47 vs 1.07, MWU p = 0.0019). That means the conditional probability of being
randset-significant is *not* 0.227 — it is much higher for control modules than
for targets. A binomial with a single constant p = 0.227 therefore tests a null
that the data violate by construction. The manuscript's admission "the tests
share contrasts, so the p is approximate" (§3.6) acknowledges only the
*independence* problem and not the more damaging *non-constant-p* problem.

Worked illustration of the circularity.
Let p_c be the probability a control test is significant and p_t the same for a
target. The binomial H0 forces p_c = p_t = 0.227. But T1 establishes
p_c > p_t. So the binomial rejects its own H0 *precisely because* p_c ≠ p_t —
which is what T1 already showed. Testing "are controls over-represented among
significant?" with a flat-rate binomial is circular, and the p-value it returns
is the probability of an outcome under a model that the data demonstrably
violate.

**【Why it matters】** T2 is presented as independent corroboration of T1 ("the
null that confirms analgesic targets confirms red blood cells at nearly twice
their expected rate"). In reality it is a re-expression of the same
kind→effect-size relationship viewed through a binary threshold, and its
p-value is not interpretable as a clean test of selectivity. A reader who treats
0.0067 as evidence of a *second, independent* failure is over-reading the
statistic. Worse, because the 32 trials are not independent either (one module
can be significant across several of the 6 contrasts, contributing several
correlated trials), the effective n is well below 32; 0.0067 is therefore not
even a conservative bound — its actual calibration is unknown.

**【Specific fix】** Replace the binomial with a label-permutation test: hold the
132 (module, contrast) |z| values fixed, permute only the `kind`
(target/control) label across modules while preserving the observed degree/size
structure, and compute the null distribution of "number of controls among the
32 largest |z|". Report that permutation p in place of the binomial p.
Alternatively, and more conservatively, demote T2 from a significance test to a
descriptive cross-tabulation: state "14/32 (44%) of randset-significant tests
are controls versus 30/132 (23%) expected" and explicitly note that this is a
restatement of T1, not an independent test; drop the 0.0067 p-value as a
validity claim. In either case, do not keep a binomial p-value as a validity
claim.

---

### F2 — The "132 independent tests" label overstates independence; T1/T2 assume it and the effective n is never quantified

**【Problem】** The manuscript repeatedly calls the deduped set "132 independent
tests," but the 132 (22 modules × 6 contrasts) statistics are correlated both
across contrasts that share a module and across modules that share a contrast;
only the KO/OE mirror is genuinely removed. T1's MWU and T2's binomial consume
these correlations as if they were independent observations.

**【Evidence】** `scripts/22_etiology_comparison.py:687-695` dedups by keeping
knockout only. That is correct: KO/OE are exact mirrors — I verified max|Δ| = 0
across all 132 in `etiology_summary.json`, and directly in
`etiology_bootstrap_differences.csv` where Glial/immune Bone_D7↔PTX KO Δ =
+0.2742888 and OE Δ = −0.2742888 with identical p_permutation_exact =
0.00142828. So the 264→132 reduction is a legitimate dedup of the *mirror*, not
a claim that the 132 are mutually independent. Yet `etiology_summary.json:42-43`
and manuscript §3.6 say "132 independent tests per network scope."

Within those 132, Glial/immune appears in all 6 pairs and each appearance shares
the same 8-sample pain axes of the shared contrasts; the 22 modules within one
pair (e.g. Bone_D7 vs PTX) all share the two reversal vectors. I recomputed T1's
MWU p = 0.0019 on the 102 target + 30 control |z| values; that statistic
treats 132 as independent draws. The manuscript does say "the 132 effects are
grouped within contrasts rather than fully independent, so the p is
approximate" (§3.6), but it never quantifies *how* approximate — i.e. what the
effective number of independent tests is.

**【Why it matters】** If the effective independent units are closer to 22
(modules) or 6 (contrasts) than 132, the MWU p = 0.0019 may be optimistic by
one or two orders of magnitude, and the binomial p in F1 is similarly
uncalibrated (and, per F1, mis-specified regardless). For a study whose point
is the *validity* of a statistical layer, leaving the dependence unquantified is
precisely the gap a methods reviewer must close.

I note in fairness that the BH q-values are not at risk: BH tolerates
(positive) dependence and still controls FDR, and the manuscript uses BH only
for q-values, not for the descriptive MWU/binomial. The exposure is confined to
the two descriptive tests.

**【Specific fix】** Add an effective-independence quantification. Two options,
either acceptable:

(a) Report the MWU p under a permutation null that shuffles module labels across
the 132 rows while preserving the contrast/module block structure — this gives a
dependence-aware p with no independence assumption.

(b) State the effective number of independent tests via the eigenvalues of the
132×132 |z| correlation matrix (or, more simply, bound it at "≤22 modules and
≤6 contrasts").

Then change the wording from "132 independent tests" to "132 deduped hypotheses
(KO/OE mirror removed; residual correlation across contrasts and modules not
assumed independent)."

---

### F3 — T3's reported gene-level MDE (0.488 log₂) is computed on a 5-gene high-SE tail and is mislabeled "median SE of expressed genes"

**【Problem】** The negative-result script computes `median_se_expressed_genes =
0.2048` and `detectable_abs_log2FC_genes = 0.488` from a filter that selects
the five *least*-powered genes, not typical expressed genes; the label and the
resulting MDE contradict the genes the manuscript itself reports as detected.

**【Evidence】** `scripts/33_negative_result_quantification.py:127-129`:
`stable = sc[(sc.q > 0.5) & (sc.logFC.abs() >= 0.2)]`, then
`se = stable.logFC.abs() / stable.t.abs()` and `med_se = median(se) = 0.2048`.
I recomputed this and found exactly **5** genes satisfy the filter — genes with
a large effect size but non-significant call, i.e. the high-variance tail. The
variable is named `median_se_expressed_genes`, but the *true* median SE across
all 34,603 pseudobulk genes is **0.0006**, and the median SE among all
non-significant genes (q > 0.5) is also far below 0.2048.

Consequently `detectable_abs_log2FC_genes = 2.38 × 0.2048 = 0.488` is the
detectable effect for the *noisiest* genes, not for typical genes. The
inconsistency is direct and checkable: `sc_GSE292522_pseudobulk_DE.csv` shows
Gfap detected at |logFC| = 0.0036 (q = 0.006, SE 0.0006), Ifit2 at 0.048 (q =
0.042, SE 0.0123), Oas1c at 0.029 (q = 0.014, SE 0.0059) — all far below
0.488. So the deposited "detectable |log₂FC| = 0.488" cannot be the typical
detectable effect; it is an artefact of a mislabeled, mis-filtered statistic.

The manuscript prose (§3.8) correctly says gene-level DE "detected Ifit2 at
−0.048," so the *text* is fine; the *artefact number* is wrong and mislabeled.
I also verified, in the author's favor, that the myeloid-marker null is genuine:
Cd68 logFC = −0.0005 (SE 0.0010, q = 0.78) and Aif1 logFC = +0.0009 (SE 0.0036,
q = 0.92) — tiny SE, negligible effect, so a large effect *would* have been
detected and "genuine null for large effects" is well-supported by the per-gene
SEs, independent of the broken summary statistic.

**【Why it matters】** T3's role is to calibrate power so that the myeloid-fraction
null is read as "underpowered, not evidence of equality." A correct gene-level
power statement strengthens that calibration. The mislabeled 0.488 instead
implies genes need |logFC| > 0.48 to be called — which the same table refutes
(Gfap at 0.0036). A careful reader cross-checking the deposited artefact against
the prose will find them inconsistent, which erodes trust in the quantification
layer exactly where the paper argues about quantification. The error is silent
because the prose never cites 0.488; but a reviewer or replicator opening
`negative_result_quantification.json` will see a number that contradicts the
detection examples.

**【Specific fix】** In `33_negative_result_quantification.py:127-129`, replace the
filter with the median SE over all expressed genes (or over all genes with a
reliable t, i.e. `sc.logFC.abs() / sc.t.abs()` with t ≠ 0), and rename the
variable to reflect what it actually summarises (e.g. `median_se_all_genes`).
Recompute `detectable_abs_log2FC_genes` from that correct median (expected
≈ 2.38 × 0.0006 ≈ 0.0014). If the intent was to report power for a
*representative* gene, state that explicitly and use the all-gene median, not the
5-gene tail. Leave the myeloid-fraction MDE (0.00123) untouched — that one is
correctly derived from the per-animal fractions.

---

### F4 — T5 per-test MDE is the right frame but rests entirely on bootstrap-SE symmetry and uses a mis-specified t-based geometry

**【Problem】** The one-sample MDE = (z₀.₀₂₅ + z₀.₈₀) × median bootstrap SE is the
correct frame for a single-estimate Δ, but the entire "median target lies below
MDE" conclusion hangs on (a) the bootstrap SEs being reliable and symmetric and
(b) the t-based variant, which uses a two-sample n=2/n=2 geometry that does not
correspond to the actual design.

**【Evidence】** `scripts/33_negative_result_quantification.py:165-171`:
`z_norm = 1.959964 + 0.841621 = 2.801585`;
`se_delta = (pc.ci95_hi - pc.ci95_lo) / (2 * 1.96)`;
`med_se = 0.0185`; `mde_units = 0.0518`;
`d_mde_unpaired = TTestIndPower.solve_power(nobs1=2, ratio=1, ...) = 5.653`;
`mde_units_t = 0.1046`.
I recomputed all five from `etiology_bootstrap_differences.csv` and they match
exactly (med_se 0.0185, normal-approx MDE 0.0518, t-based 0.1046).

The SE-recovery formula assumes the 2,000-resample bootstrap CI is symmetric and
normal — defensible here because the Δ distribution over resamples is unimodal,
but it is never checked (no skewness/kurtosis is reported). The t-based MDE uses
`nobs1=2` (a two-sample n=2 vs n=2 design), which models the *two reversal
arms* as independent observations; but KO/OE are exact mirrors, so there is no
second arm — a one-sample geometry (df = n_boot − 1) is the honest counterpart.
The "conservative" label is correct in direction, but the geometry is
mis-specified rather than merely conservative.

There is also an internal tension the manuscript handles honestly: the exact
permutation already flags 31/132 tests at p < 0.05 (I recomputed 31), so the
layer *does* detect many target effects; the MDE lens says the *median* target
(|Δ| = 0.031) is below the 0.052 threshold. These are consistent (large effects
detected, median not), but the MDE framing can read as "the layer is
underpowered" when the exact test is in fact more powerful at small n than a
z-test on the same SE.

**【Why it matters】** T5 is the paper's "even when we give the layer the benefit
of the doubt on power, it still fails on validity" argument. If the SE recovery
or the conservative t-based number is shaky, an advocate could dismiss the power
caveat. More importantly, the per-test MDE is the right conceptual frame (each Δ
is one estimate) and should be stated as such without leaning on a mis-specified
two-sample correction. The bootstrap SEs themselves are the fragile link: they
are estimated from 2,000 resamples of 4-case/4-control groups, where the
effective resample size is tiny and the SE of the SE is large; the median is
used to stabilize this, which is reasonable, but the stability should be stated.

**【Specific fix】** Keep the normal-approx one-sample MDE (0.0518) as the primary
number — it is correctly framed. Replace the t-based variant's `nobs1=2,
ratio=1` (two-sample n=2/n=2) with a one-sample t reference:
`TTestIndPower` is the wrong tool for a one-sample problem; use the bootstrap df
directly (e.g. the 0.05/0.80 t critical with df ≈ 1,999, which is ≈1.96 and
gives essentially the same 0.0518), or simply drop the t-based number and keep
only the normal-approx. Add a one-line check that the bootstrap Δ distributions
are approximately symmetric (e.g. report the median absolute skewness across the
132 rows), so the `(ci95_hi − ci95_lo)/(2·1.96)` recovery is justified. Keep the
sentence "31/132 tests reach exact permutation p < 0.05" alongside the MDE
statement to avoid over-claiming power failure.

---

### F5 — The "two logFC sources" design is honestly disclosed; residual dual-claim risk is minor and should be stated once up front

**【Problem】** Using muscle-adjusted limma-trend logFC for §3.3 (within-study)
and unadjusted per-sample mean logFC for §3.5 (cross-aetiology Δ) creates two
estimators that a reader must track; the risk of dual claims is low because the
two agree on the headline markers, but the separation deserves one consolidated
statement.

**【Evidence】** `manuscript_final.md §2.3` (line 59) explains the split and
asserts "where the two disagree on a headline marker we quote the
muscle-adjusted value." I checked the agreement on the headline canonical
markers by comparing `etiology_marker_cross.csv` (unadjusted per-sample mean,
Bone_D7) with the muscle-adjusted moderated logFC in Table 1 / §3.3: Cd68
+2.839 vs +2.86 (Δ 0.02), Aif1 +2.192 vs +2.18, Galanin −0.672 vs −0.67. The
two sources agree to within rounding on every headline marker, so there is no
cherry-picked divergence. The split is also justified mechanically: the exact
permutation of §3.5 enumerates case/control splits of raw per-sample means, and
"a per-animal covariate cannot be carried inside that enumeration"
(`manuscript_final.md` line 59). The §3.3 estimator instead uses the moderated t
from the shared 13,601-gene fit (`20_multistage_pipeline.py` STAGE 6,
`module_level_signal.csv`).

**【Why it matters】** This is a real design decision that could have been a
dual-claim trap, but it is handled transparently and the estimators corroborate
each other. The only residual risk is reader confusion and the possibility that
a non-headline marker where the two diverge could be quoted from the more
favorable source without the reader noticing which estimator underlies it.

**【Specific fix】** Add a single up-front sentence in §3.5 (or a footnote to
Table 3/4) stating that the cross-aetiology per-sample-mean Δ and the
within-study muscle-adjusted moderated logFC use different estimators by
necessity, that they agree to <0.03 log₂ on all reported markers, and that any
future divergence will be resolved in favor of the adjusted value. This converts
a latent design choice into an explicit, auditable rule and removes the last
trace of dual-claim risk.

---

## 3. Stands up (verified and sound)

**S1 — Exact-permutation floors and the paired null that omits the control term
are both correct.** I recomputed from `etiology_marker_cross.csv`,
`etiology_module_logfc_cross.csv`, and `etiology_bootstrap_differences.csv`:

- Minimum unpaired p = 0.0006121199755152 = **3/4,901** (manuscript §2.7, §5). ✓
- Minimum paired (marker/module) p = 0.04225352112676056 = **3/71**
  (manuscript §2.7). ✓
- Minimum paired reversal p = 0.028169014084507 = **2/71** (manuscript §2.7). ✓

Worked check of the floors. For the unpaired null the 4,900 differences include
the observed split and its mirror (swap both case/control assignments), so the
minimum attainable `ge` is 2, giving p = (2+1)/(4900+1) = 3/4901. For the paired
marker null the shared control cancels, leaving 70 eight-case splits; the
observed split and its negating complement are both in the null, again `ge` ≥ 2,
so 3/71. For the paired *reversal* null the complement split is *not* an exact
negation because the two pain axes swap asymmetrically, so only the observed
split sits at the extreme (`ge` ≥ 1), giving 2/71. This asymmetry is correctly
reasoned in `22_etiology_comparison.py:283-328` and confirmed by the data.

The paired null legitimately omits the control term for the within-bone pair
(Bone_D14 vs Bone_D7, identical sham) and the PTX/OXA pair (identical vehicle):
the control is literally the same samples (`22_etiology_comparison.py:471-496`,
`534-559`), so `Δ = mean(A_case) − mean(B_case)` is exact. This is a defensible
design choice, not a shortcut, and I would defend it against any reviewer who
assumed a shared-control correction was required.

**S2 — T1 anti-selectivity is a genuine, well-executed finding.** From
`etiology_bootstrap_differences.csv` (per_contrast, knockout): n_target = 102,
n_control = 30; median |z| target = 1.0729, control = 2.4702; two-sided MWU p =
0.001914 (manuscript reports 0.0019). The REMOVE-IFN variant (drop the
interferon control module) gives median |z| control = 2.3327 and MWU p =
0.018697 (manuscript reports 0.019) — anti-selectivity persists without the
single strongest control. This is the paper's strongest methodological result
and it holds. I carry one caveat from F2: the 132 |z| values are correlated, so
the *exact* calibration of 0.0019 is approximate, but even after quantifying
dependence (F2) the effect is large enough that the conclusion is secure.

**S3 — Multiplicity correction is correct.** KO/OE dedup (264→132) in
`22_etiology_comparison.py:687-695` is principled and verified (KO and OE are
exact mirrors). BH applied per scope; `etiology_summary.json` reports
`min_q_permutation_per_contrast = 0.165447` (manuscript "min q = 0.165") and
`n_sig_randset_q05_percontrast_indep = 32` with
`min_q_randset_per_contrast = 0.021956`. I reproduced 32/132 randset-significant
and 14 of those being controls. The correction is not where the paper is weak,
and the "132 independent" language issue (F2) does not touch the q-values
because BH tolerates dependence.

**S4 — The myeloid-fraction "underpowered, not evidence of equality" claim is
defensible, and the Cd68/Aif1 null is genuine.** From
`sc_GSE292522_composition.csv` and `sc_GSE292522_pseudobulk_DE.csv`: Con myeloid
fraction mean 0.00117, B16F10 0.00052; pooled SD (denominator n₁+n₂−2 = 6) =
0.00051; observed Cohen's d = 1.267 vs MDE d = 2.38. MWU on the 4-vs-4
per-animal fractions: U = 14/2, two-sided p = 0.114286, one-sided in the
observed (B16F10 < Con) direction p = 0.057143 — matching the manuscript. With n
= 4 per group, "non-significant ⇒ inconclusive" is the only honest reading;
claiming equality would require a TOST/equivalence design that the study
correctly does not mount. And the myeloid markers are a *bona fide* null: Cd68
logFC = −0.0005 (SE 0.0010, q = 0.78), Aif1 = +0.0009 (SE 0.0036, q = 0.92) —
a large effect would have been detected, so "genuine null for large effects" is
supported by the per-gene SE, not merely asserted. (The only bug in this
neighborhood is the unrelated T3 median-SE mislabel, F3.)

---

## 4. Questions for the authors

1. **On T2 (F1):** Do you agree the binomial is a restatement of T1 rather than
   independent evidence? If you keep it, will you replace the p-value with a
   label-permutation test (or a descriptive cross-tab) so the validity claim
   does not rest on a mis-specified test?

2. **On independence (F2):** Can you report the effective number of independent
   tests among the 132 (e.g. via the correlation structure of the |z| values, or
   a permutation null)? Without it, how should a reader calibrate the MWU
   p = 0.0019?

3. **On T3 (F3):** The `median_se_expressed_genes = 0.2048` figure comes from a
   5-gene filter, not expressed genes. Will you correct the deposited
   `negative_result_quantification.json` and the variable name? The correct
   typical detectable |log₂FC| is ≈0.0014 — do you agree the prose (which uses
   the Ifit2 example) is unaffected?

4. **On T5 (F4):** Will you replace the two-sample n=2/n=2 t-based MDE with a
   one-sample t reference, and add a symmetry check on the bootstrap Δ
   distributions?

5. **On the paired reversal null (S1):** For the within-bone and PTX/OXA pairs
   the shared control is the same animals, so omitting it is correct — but for
   the *reversal* endpoint the two pain axes are equal (c0a = c0b) yet swap
   asymmetrically under the 8-case permutation, which is why the floor is 2/71
   rather than 3/71. Is this asymmetry intended, and should it be stated more
   explicitly in §2.7 so a reader does not expect 3/71 for reversal?

6. **On composition vs biology (§3.5/§3.6):** Since the erythrocyte
   cross-aetiology Δ (−2.48) is offered as a batch/compositional proxy, and
   neuronal modules move in the same direction as erythrocyte, do you plan to
   report the *composition-adjusted* cross-aetiology neuronal Δ, or is the
   current "cannot be separated from composition" caveat the intended final
   word?

---

## 5. What I actually checked

**Scripts read (full):** `scripts/_limma.py` (notably `exact_two_sided_p` lines
420–426, `emp_p_two_sided` 278–300, `cos_batch` 303–337,
`degree_matched_sets` 360–383, `all_splits` 386–402),
`scripts/22_etiology_comparison.py` (full; paired null branches at 283–328,
471–496, 534–559; dedup/BH at 687–695), `scripts/33_negative_result_quantification.py`
(full), `scripts/20_multistage_pipeline.py` (full, for the §3.3
muscle-adjusted estimator and Table 1 provenance).

**Artefacts recomputed against (all numbers below are mine, recomputed from the
CSVs/JSONs with pandas/numpy/scipy):**

- `etiology_bootstrap_differences.csv` — T1 (MWU p = 0.001914; medians
  1.0729/2.4702; REMOVE-IFN p = 0.018697), T2 (32 sig / 14 control / base
  30/132 / binomial p = 0.00671), T5 (med_se 0.0185, z_sum 2.801585, MDE 0.0518,
  t-based 0.1046, max Δ 0.2743, 31/132 exact-p<0.05, 41/102 exceed normal MDE,
  23/102 exceed t-based MDE).

- `etiology_marker_cross.csv` and `etiology_module_logfc_cross.csv` — exact-
  permutation floors: unpaired min p = 3/4,901 (0.0006121199755152); paired min
  p = 3/71 (0.04225352112676056).

- `etiology_summary.json` — `min_q_permutation_per_contrast = 0.165447`,
  `n_sig_randset_q05_percontrast_indep = 32`, `n_independent_tests_per_scope =
  132`; all reproduced.

- `negative_result_quantification.json` — every T1/T2/T3/T5 value reproduced to
  the printed precision; **discrepancy found** in T3 `median_se_expressed_genes
  = 0.2048` (computed on 5 high-SE genes, not expressed genes; true median SE
  across all genes = 0.0006).

- `sc_GSE292522_composition.csv` + `sc_GSE292522_pseudobulk_DE.csv` — T3: Con
  0.00117 / B16F10 0.00052, pooled SD 0.00051, d = 1.267, MDE d = 2.38, MWU U =
  14/2 two-sided p = 0.114286, one-sided(observed dir) p = 0.057143; Cd68 SE
  0.0010 (q 0.78), Aif1 SE 0.0036 (q 0.92) → genuine null.

- `sc_GSE292522_summary.json` — myeloid MWU p two-sided 0.1143 / less 0.0571;
  reproduces the manuscript.

**Discrepancies stated:**

1. T3 `median_se_expressed_genes = 0.2048` (and downstream
   `detectable_abs_log2FC_genes = 0.488`) is computed from a 5-gene (q>0.5 &
   |logFC|≥0.2) tail, not from expressed genes; the label is wrong and the
   number contradicts detected genes (Gfap 0.0036, Ifit2 0.048). True typical
   detectable |log₂FC| ≈ 0.0014. *(F3)*

2. No arithmetic discrepancy in T1, T2, T5, the permutation floors, the
   multiplicity correction, or the myeloid-fraction test — all reproduced
   exactly. The concerns there are about *interpretation/assumption* (constant-p
   binomial, independence of the 132, bootstrap-SE reliance, two-sample t
   geometry), not about the numbers themselves.

**Items explicitly NOT re-run:** I did not re-execute the full `22_*`/`20_*`
pipelines (not required — every reported number is recomputable from the
deposited artefacts, which I did). I did not re-parse the raw FPKM
(`data/raw/GSE149648`, `GSE160543`); the parsing assertions in the scripts were
reviewed by reading, not by re-running the parser on raw files. I did not open
any file in the forbidden set (`docs/review/REVIEW_round*`, `round1..5_*`,
`R3_editorial_review.md`, `author_verification_statement.md`,
`GITHUB_DEPOSIT_SOP.md`, `manuscript_v0.*.md`, or other reviewers' files).

---

## 6. Minor / line-level notes (design-relevant)

- **§2.3 / Table 1 provenance.** Table 1's day-7 module t-values (e.g.
  interferon +24.92) come from the muscle-adjusted limma-trend fit on the
  13,601 shared genes (`20_multistage_pipeline.py` STAGE 6 +
  `module_level_signal.csv`), while Table 3/4 use unadjusted per-sample means.
  The manuscript is explicit about this; I only ask for the one consolidated
  statement suggested in F5.

- **Bootstrap degeneracy.** `etiology_summary.json` reports
  `max_degenerate_boot = 0`, so the 2,000-resample bootstrap never tripped the
  degenerate-resample guard (`22_etiology_comparison.py:353-356`). Good — the
  SEs are not artificially inflated by dropped resamples, which supports the T5
  SE-based MDE.

- **Randset averaging across studies.** For unpaired contrasts the degree-
  matched null uses the mean of the two adjacency matrices
  (`22_etiology_comparison.py:362-364`). This is a reasonable pooled degree
  reference, but it is worth a sentence noting that the two studies have
  separately built networks and the average is a compromise, not a shared graph;
  the null's validity for an unpaired cross-study comparison could be made more
  explicit.

- **BH and dependence.** The q-values are BH within scope; BH is robust to the
  (positive) dependence present here, so the q-values (min 0.165 for exact-
  perm, 0.022 for randset) are safe. The independence caveat applies only to the
  MWU/binomial descriptive tests, as noted in F2.

- **One-sided MWU direction.** The manuscript reports the myeloid one-sided p =
  0.057 in the observed direction (B16F10 fraction lower). My first pass
  computed it in the wrong direction (0.97) because of argument order; the
  manuscript's 0.057 is correct (P(U ≤ 2) = 4/70). No issue — noted only for
  transparency.

- **The reversal null is genuinely exact.** For the reversal endpoint the seed
  and operator P are fixed and only the pain axis varies under permutation, so
  `null = cos(delta, permuted axes)` is exhaustive
  (`22_etiology_comparison.py:270-271`); I confirmed obs always lies inside the
  null (minimum p never reaches 0 across all 528 rows). This is a strength that
  should be highlighted in the methods, because it is the property that makes the
  floors hard rather than arbitrary.

---

## 7. Recommended priority for revision

- **Must fix before publication:** F3 (correct the mislabeled T3 gene-level MDE
  in the deposited JSON and variable name) and F1 (replace or relabel the T2
  binomial so the validity claim does not rest on a mis-specified test).

- **Should fix:** F2 (quantify effective independence / dependence-aware MWU)
  and F4 (one-sample t-based MDE + bootstrap-symmetry check).

- **Nice to have:** F5 (one consolidated two-estimator statement) and the §2.7
  clarification on the 2/71 vs 3/71 reversal floor asymmetry (Q5).

## 8. Cumulative effect of the fixes on the paper's claims

None of the four findings (F1–F4) requires retracting a biological or
methodological conclusion. The cross-aetiology divergence (§3.5), the day-7
interferon event (§3.3), the additive-baseline redundancy (§3.4), and the
anti-selective failure of in-silico perturbation (§3.6) all survive. What the
fixes change is the *rhetoric of validity*: after F1–F2, T1 stands as the
single, defensible statistical demonstration of anti-selectivity, and T2 is
either removed or relabeled as descriptive. After F3, the gene-level power
statement becomes internally consistent with the detected genes. After F4, the
per-test MDE is correctly framed as one-sample. The paper is stronger, not
weaker, for having these tightened — which is the outcome a validity study
should want.

I would recommend acceptance after the must-fix items are addressed; the
should-fix items can be handled in the same revision or as minor revisions.
