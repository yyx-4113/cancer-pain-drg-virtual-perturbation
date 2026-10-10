# R3-A3 Review — single-cell / pseudobulk methodology

**Reviewer role:** single-cell / pseudobulk methodology expert
**Manuscript:** `docs/manuscript/manuscript_final.md` (first-submission posture)
**Scope of this review:** scripts 31/32/34, §2.8, §3.8, Table 5, the single-cell
parts of the Abstract/Limitations, and every §3.8 claim recomputed from
`results/tables/sc_GSE292522_*` and — where the artefacts were insufficient —
directly from the raw 10x matrices in `data/processed/GSE292522/`.
**Independence:** I read no prior-round review files, no other reviewer's
output in this directory, no earlier manuscript version, no `.workbuddy/`
files, and none of the pre-existing `A3_*.py` scripts in `scripts/` (their
provenance is ambiguous under the independence rules, so I deliberately did
not open them). All recomputations below come from my own script
(`_A3_panel_verify_sc.py`, deposited alongside this report) run against the
raw data and the deposited artefacts. The reproducibility gate
(`scripts/98_check_manuscript_numbers.py`) was run by me and ends with
`276 checks, 0 failed`.

---

## 1. Verdict

**Minor Revision.**

The single-cell layer is the most honestly disclosed part of this manuscript.
The animal (not the cell) is genuinely the unit of inference in every test; the
myeloid blind spot is quantified and disclosed rather than hidden; the marker-
bin identity problem is flagged by the authors themselves; and every headline
§3.8 number I checked reproduces — most from the deposited artefacts, and the
composition and annotation numbers to within 1e-16 from an independent
replication of the pipeline on the raw matrices. I found no pseudoreplication
anywhere in the single-cell inference.

The revision needed is real but bounded: one statistical formula error
(observed Cohen's d), one broken methodological promise (leave-one-animal-out
refits promised in §2.8 but absent for §3.8), an off-by-one compartment count,
an unverifiable QC number (11.9% purity) described with an inaccurate
criterion, and an over-generalised Abstract/Discussion statement about
"model-specificity" that §3.8 itself does not license. None of these
undermines the section's core conclusion; all are fixable without new
analysis beyond what the raw data already support.

---

## 2. Summary of the manuscript (as it bears on the single-cell layer)

The paper reanalyses two public bulk rat DRG RNA-seq series (Walker-256 tibia
bone cancer, GSE149648; paclitaxel/oxaliplatin neuropathy, GSE160543) and
shows, with exact permutation tests, that the two aetiologies drive divergent
DRG programmes (myeloid/interferon in bone cancer at day 7; neuronal injury in
CIPN). A second, methodological claim is that in-silico perturbation of these
bulk transcriptomes fails its own negative controls and is anti-selective, so
it cannot rank analgesic targets. §3.8 adds an independent single-cell
validation on GSE292522 (mouse DRG, 10x, 4 B16F10 cancer-nerve-injury animals
vs 4 littermate controls, 55,967 cells): whole-DRG pseudobulk DE (limma-trend,
animal as unit) finds a neuronal-injury programme (Sprr1a, Atf3) and no
up-regulated interferon-stimulated genes; a per-compartment extension
(Usoskin-2015-style marker-bin annotation, per-animal compartment pseudobulk,
BH within compartment) extends the null across the non-myeloid compartments;
the myeloid compartment — where a Walker-256-type programme would reside — is
essentially absent (46 cells total) and therefore untestable; cross-model
Spearman concordance with Walker-256 is weak (ρ = 0.112 / 0.098 / 0.030 for
D1 / D7 / D14). The conclusion drawn is that the interferon programme is a
feature of the Walker-256 model, not of cancer-pain DRG generally.

---

## 3. Major issues

**M1. The observed Cohen's d for the myeloid-fraction drop is computed with a
non-standard pooled SD; the correct value is 1.2680, not 1.37.**
- *Citation:* manuscript §3.8 (line 155) and Table 5 (line 315): "the observed
  Cohen's d is 1.37 against a minimum detectable d of 2.38"; derived in
  `scripts/33_negative_result_quantification.py:117-118`.
- *Evidence:* lines 117–118 of script 33 compute
  `pooled_sd = sqrt(SS_Con/7 + SS_B16F10/7)`, i.e. the total sum of squares
  divided by 7. The standard pooled SD for two independent groups is
  `sqrt((SS_Con + SS_B16F10)/(n1+n2−2)) = sqrt(SS_total/6)`. Recomputed from
  `sc_GSE292522_composition.csv`: mean Con fraction 0.00117, mean B16F10
  0.00052, mean difference 0.00065; correct pooled SD gives **d = 1.2680**;
  the script's /7 variant gives exactly 1.3696 (the deposited 1.37). The
  error shrinks the pooled SD and inflates d by the factor √(7/6).
- *Required fix:* correct the formula (denominator 6, not 7) and update
  §3.8, Table 5 and the JSON. Note the correction is *favourable* to the
  paper's argument (observed d is even further below the MDE of 2.3808, which
  I reproduce exactly), so nothing qualitative changes — but a manuscript
  whose central methodological claim is statistical rigour cannot carry a
  wrong degrees-of-freedom in a reported effect size. The 98 gate does not
  catch this because its expectation for this number is hardcoded
  (`98_check_manuscript_numbers.py:369-370`), unlike the T1/T2/T5 checks,
  which are genuinely runtime-derived; the gate should recompute d from the
  composition CSV.

**M2. §2.8 promises leave-one-animal-out refits for every relied-upon §3.8
result; none exist.**
- *Citation:* manuscript §2.8 (line 97): "any result relied on in the text
  carries its leave-one-animal-out refit outcome."
- *Evidence:* scripts 31, 32 and 34 contain no LOO computation, and
  `results/tables/` contains no single-cell LOO artefact (the LOO artefacts
  that exist belong to the perturbation layer). The results most in need of
  it are precisely the fragile ones: Ifit2 down in the neuron bin at
  q = 0.0459 (a hair above 0.05 in any refit is possible), Ifit1 in
  satellite glia (q = 0.0037) and Isg15 in satellite glia (q = 0.0071),
  where the satellite-glia bin is 25–108 cells per animal and dominated by
  one large sample (108 of 361 cells in Sample 7 alone; my recomputation).
  Sample 7 is visibly atypical in every composition column (neuron fraction
  0.559 vs 0.68–0.72 in the other B16F10 samples).
- *Required fix:* either run and report the LOO refits for the §3.8 claims
  relied upon (the four ISG results, Gal/Scn10a neuron-bin results and the
  myeloid-fraction statistics), or delete the promise. As written, a reader
  assumes a robustness check that was never performed.

**M3. The Abstract and Discussion assert "model-specificity" more strongly
than §3.8's own blind-spot logic licenses, and the B16F10 model is
under-specified for the comparison to be judgeable.**
- *Citation:* Abstract (line 21): "The interferon program is not reproduced in
  the B16F10 model … so the Walker-256 interferon signature is model-specific
  rather than a generic cancer-pain signature." Discussion (line 165): "the
  interferon signal is model-specific." §3.8 itself (line 157) is careful:
  the compartment test "cannot exclude a myeloid-restricted interferon
  programme, and we do not claim that it does."
- *Evidence:* the Walker-256 interferon signal was measured in bulk tissue in
  which myeloid infiltration is a dominant carrier; the B16F10 dataset
  contains 46 myeloid barcodes total (my recomputation: exactly 46, range
  0.0147%–0.1994% per sample), so the one compartment that carries the
  discovery signal is unsampled. A "not reproduced" claim is therefore
  demonstrated only for the non-myeloid compartments. Additionally, the
  manuscript never states (a) where B16F10 was inoculated (the repo's own
  script header says femur/tibia, but the manuscript says only
  "cancer-nerve-injury model"), or (b) at what post-inoculation time point
  the DRG were harvested. A transient day-7-type interferon peak cannot be
  excluded from a single, unspecified time point; this materially weakens
  even the hedged version of the claim.
- *Required fix:* align the Abstract and Discussion with §3.8's hedge ("does
  not generalise at the level detectable in the non-myeloid compartments of
  this dataset"), and state inoculation site, harvest time point, and animal
  age/sex for GSE292522 in §2.8. The Limitations paragraph (line 177) already
  has the right framing — use it consistently.

**M4. A large, statistically significant compositional shift in the B16F10
dataset — neuron fraction down, Schwann fraction up — is present in the
deposited composition table but never mentioned, although it bears directly
on how the whole-DRG pseudobulk should be read.**
- *Citation:* `results/tables/sc_GSE292522_composition.csv`; manuscript §3.8
  reports only the myeloid fraction.
- *Evidence (recomputed, per-animal MWU, the same test the manuscript uses
  for myeloid):* Neuron fraction 0.7859 (Con) vs 0.6695 (B16F10), MWU
  two-sided p = 0.0286; Schwann 0.1885 vs 0.2928, p = 0.0286; Fibroblast
  0.0084 vs 0.0215 (2.6×, p = 0.1143). The paper applies exactly this
  compositional logic to the discovery dataset (§3.3: neuronal module
  declines = one observation about cell composition, not five about targets)
  but does not apply it to the validation dataset, where 703 genes are
  called down-regulated in whole-DRG pseudobulk that is CPM-normalised
  across a shifting cell-type mixture. The up-regulated injury genes
  (Sprr1a, Atf3) and the per-compartment results are not threatened by this
  (within-bin CPM is dilution-free), but the whole-DRG "down" side and the
  flat-ISG statement should be interpreted against the mixture shift, and
  the symmetry with the paper's own §3.3 reasoning is too conspicuous to
  leave implicit.
- *Required fix:* add the neuron/Schwann composition result (and its MWU p)
  to §3.8, and one sentence acknowledging that the whole-DRG down-regulated
  set is partly compositional in this dataset too.

---

## 4. Minor issues

**m1. "Seven testable compartments" is wrong; the artefact contains eight.**
§3.8 (line 157) and Table 5 (line 316) say seven. The rule (≥5 cells in ≥6
animals) excludes exactly four of twelve candidates — T cells (1 cell total,
my recomputation), NP and TH neuron subtypes (0–1 and 0–3 cells per animal,
my recomputation), and myeloid (≥5 cells in 5/8 animals) — leaving **eight**
testable compartments, and `sc_GSE292522_celltype_DE.csv` contains all eight
(Endothelial, Erythrocyte, Fibroblast, Neuron, Neuron.NF, Neuron.PE,
SatelliteGlia, Schwann). 12 − 4 = 8. The substantive conclusion is unaffected
(no ISG is up-regulated in the eighth, Erythrocyte, either — all q ≈ 0.88),
but the count and the arithmetic must be corrected in both places.

**m2. The 11.9% neuron-bin purity figure is not reproducible from the deposit
and its described criterion is not "non-overlapping".**
§2.8 (line 97) says the neuron bin was scored "against an independent,
non-overlapping criterion (Snap25/Tubb3/Syp-positive and myelin-negative)".
Snap25, Tubb3 and Syp are *members of the Neuron marker set used for the
assignment* (`scripts/32_percelltype_sc.py:55-56`), so the criterion overlaps
the assignment markers by construction. No script or artefact in the deposit
computes the 11.9%. My best-faith replications give per-sample strict purity
9.3%–19.4% (mean 15.8%, pooled 15.8%–17.4% depending on the myelin gene set),
so 11.9% is plausible in magnitude but unverifiable in value. Fix: deposit
the QC script, state the exact criterion, and replace "non-overlapping" with
an accurate description.

**m3. "55,967 cells as deposited" (§2.8, line 95) vs "as retained by our
parsing and quality-control pipeline" (§3.8, line 155) — these are different
claims, and the second is wrong.** My independent count of the eight raw
matrices is exactly 55,967 barcodes; no filtering occurred, so there is no
"retained by QC" step to describe. State plainly: all deposited barcodes were
used, no cell-level QC was applied. (This also means empty/damaged droplets
are included — harmless for pseudobulk sums, but say it.)

**m4. The myeloid-fraction range "0.02–0.2% of cells per sample" (line 157)
rounds the lower bound up.** Recomputed range: 0.01468%–0.19937%. Write
"0.015–0.2%".

**m5. Silent marker-gene drops in the annotation sets.** "Blbp" (Fabp7) and
"Oct6" (Pou3f1) are absent from the 10x features (old aliases never mapped),
and the Pgp9.5 symbol is absent in both scripts' casings ("Pgp9.5" in script
31, "PgP9.5" in script 32 — an inconsistency between the two scripts that
happens to be moot). `marker_score()` silently skips missing genes
(`32_percelltype_sc.py:97-102`), so the SatelliteGlia set is effectively six
markers and Schwann seven, with no warning. Log the dropped markers, or map
aliases. (I verified the binning results are otherwise driven by Sox10,
Slc1a3, Aqp4, Kcnj10, Glul, Fgfr3, so the impact is small.)

**m6. Unreported compartment results that complicate the narrative.** In
`sc_GSE292522_celltype_DE.csv`: Cx3cr1 is *up* in the Neuron bin (+0.012,
t = 5.45, q = 0.0087) and Neuron.NF bin (q = 0.0060); Gal is also
significantly up in the Schwann bin (+0.453, q = 0.0061) and Endothelial bin
(+0.671, q = 0.0001); Snap25 is down in the Schwann bin (q = 0.0352) and Lyz2
down in Endothelial (q = 0.0009). None of these is required to change a
conclusion, but the Cx3cr1-in-neurons result is directly relevant to the
disagreement with ref [15] about neuron-intrinsic signalling (it is either a
contamination signal or a genuinely Cx3cr1+ neuronal subpopulation — the text
should say which reading it favours), and selective reporting of only the
convenient compartment results will be noticed by any sc-literate reader.

**m7. "log₂(CPM + 1)" is mislabelled.** Both scripts actually compute
counts per 10,000 (`31_process_sc.py:126`, `32_percelltype_sc.py:154`,
variable name `cp10k`). Constant scaling does not affect t or q, but the
reported logFC values (e.g. Sprr1a +1.55) are on a log₂(CP10K+1) scale, not
log₂(CPM+1); label it correctly.

**m8. Stale artefact inconsistency:** `sc_GSE292522_celltype_summary.json`'s
`celltype_ISG_discrimination` lists seven compartments and omits Erythrocyte,
while `sc_GSE292522_celltype_DE.csv` and `sc_GSE292522_celltype_scores.csv`
both contain it (Erythrocyte has 3 Con + 3 B16F10 ISG-score samples — the
MWU should have been computable). Regenerate the JSON so the artefacts agree.

**m9. Concordance presentation.** The per-day ρ values (0.112 / 0.098 /
0.030) reproduce exactly from `sc_GSE292522_pseudobulk_DE.csv` against the
per-day bulk tables (12,158 shared genes; my recomputation matches
`sc_GSE292522_concordance_per_day.json` to all four decimals). But
"highest at day 1" invites a per-day gradient reading of three tiny
correlations whose pairwise differences were never tested; and the p-values
(5.8 × 10⁻³⁵, 1.4 × 10⁻²⁷, 9.6 × 10⁻⁴) are unreported. Report the p-values
and describe all three simply as negligible.

**m10. Per-compartment ISG testing covered only the 10 ISG genes in the
FOCUS list** (Oas1c, down in whole-DRG pseudobulk at q = 0.0136, was never
tested per compartment because it is not in FOCUS). The sentence "no
interferon-stimulated gene was up-regulated anywhere" is true for what was
tested; state the scope.

---

## 5. Numbers I independently recomputed

Recomputation sources: deposited CSV/JSON artefacts, and — for rows marked
† — an independent re-run of the annotation pipeline on the raw 10x matrices
(`_A3_panel_verify_sc.py`, this directory).

| # | Claim (manuscript §3.8 / Table 5) | Manuscript | Recomputed | Match? |
|---|---|---|---|---|
| 1 | Total cells | 55,967 | 55,967 (raw barcodes, †) | ✓ |
| 2 | Myeloid cells total | ~46 | 46 exactly (†) | ✓ |
| 3 | T cells | 1 | 1 (†) | ✓ |
| 4 | Myeloid fraction Con / B16F10 | ~0.12% / ~0.05% | 0.1170% / 0.0517% | ✓ |
| 5 | Drop in myeloid fraction | 55% | 55.8% (ratio 0.442) | ✓ |
| 6 | Observed Cohen's d | 1.37 | **1.2680** (script 33's SS/7 variant gives 1.3696) | ✗ (M1) |
| 7 | MDE Cohen's d (n=4+4) | 2.38 | 2.3808 | ✓ |
| 8 | MWU U, two-sided / one-sided p | U = 2/14; 0.114 / 0.057 | U1 = 2, U2 = 14; 0.1143 / 0.0571 | ✓ |
| 9 | Myeloid fraction per-sample range | 0.02–0.2% | 0.01468%–0.19937% | ~ (m4) |
| 10 | Myeloid ≥5 cells in N animals | 5 of 8 | 5 of 8 (†) | ✓ |
| 11 | Pseudobulk DE q<0.05 | 1,517 up / 703 down | 1,517 / 703 (n = 34,603; BH reproduced to 6e-15) | ✓ |
| 12 | Sprr1a | +1.55 log₂, t = 54.5 | +1.5539, t = 54.49 | ✓ |
| 13 | Atf3 | +0.84 | +0.8355 (q = 0.0001) | ✓ |
| 14 | Ifit2 pseudobulk | −0.048, q = 0.042 | −0.0476, q = 0.0417 | ✓ |
| 15 | Oas1c pseudobulk | q = 0.014 | −0.0291, q = 0.0136 | ✓ |
| 16 | Gfap pseudobulk | +0.004, q = 0.006 | +0.0036, q = 0.0061 | ✓ |
| 17 | Cd68 / Aif1 unchanged | Δ ≈ 0 | −0.0005 (q = 0.783) / +0.0009 (q = 0.920) | ✓ |
| 18 | Gal whole-DRG | +1.42 | +1.4247 | ✓ |
| 19 | Gal in Neuron bin | +1.47, q < 0.001 | +1.472, q = 0.0001 | ✓ |
| 20 | Gal in Neuron.NF bin | +2.54, q < 0.001 | +2.537, q < 5 × 10⁻⁵ | ✓ |
| 21 | Scn10a down in Neuron bin | q = 0.012 | −0.247, q = 0.0117 | ✓ |
| 22 | SatelliteGlia "no detectable Gfap" | logFC 0.000 | artefact: 0.000 (q = 0.9997); raw data †: 0 UMI in 5/8 animals, ≤3 UMI total in the others, ≤0.19 CPM — near-absolute absence confirmed | ✓ |
| 23 | Ifit1 ↓ SatelliteGlia | q = 0.004 | −0.523, q = 0.0037 | ✓ |
| 24 | Ifit2 ↓ Neuron | q = 0.046 | −0.049, q = 0.0459 | ✓ |
| 25 | Isg15 ↓ SatelliteGlia | q = 0.007 | −0.630, q = 0.0071 | ✓ |
| 26 | No ISG up in any testable compartment | 7 compartments | true, but across **8** compartments (m1) | ~ |
| 27 | Excluded compartments (T cell; NP; TH; myeloid) | — | verified †: NP 0–1 cells/animal, TH 0–3, TCell 1 total, myeloid 5/8 with ≥5 | ✓ |
| 28 | Concordance ρ D1 / D7 / D14 | 0.112 / 0.098 / 0.030 | 0.1115 / 0.0985 / 0.0299 (n = 12,158 shared genes) | ✓ |
| 29 | Neuron-bin purity | 11.9% | not reproducible: 9.3–19.4% per sample, pooled 15.8–17.4% under natural criterion variants (†) | ✗ (m2) |
| 30 | Reproducibility gate | 276 checks, 0 failed | reproduced, ends "276 checks, 0 failed" | ✓ (but T3 d expectation hardcoded, see M1) |
| 31 | Composition shift Neuron / Schwann | not reported | 0.7859→0.6695 and 0.1885→0.2928, both MWU p = 0.0286 (M4) | — |

Annotation-replication check (†): my from-scratch per-animal compartment
fractions match `sc_GSE292522_composition.csv` to max |Δ| = 9.7 × 10⁻¹⁷,
confirming the deposited fractions are exactly what the described pipeline
produces from the raw matrices.

---

## 6. What would make this publishable (concrete, ordered)

1. **Fix the Cohen's d formula** in `scripts/33_negative_result_quantification.py:117-118`
   (denominator n1+n2−2 = 6), regenerate the JSON, and update §3.8 / Table 5
   (d becomes 1.27; the argument gets stronger, not weaker). Make the 98 gate
   recompute d from the composition CSV instead of asserting 1.37 (M1).
2. **Deliver or delete the LOO promise.** Run leave-one-animal-out refits for
   the §3.8 results the text relies on (the three compartment ISG results,
   Gal/Scn10a in neurons, the myeloid statistics), report the ranges in §3.8
   or an artefact table; otherwise remove the sentence at line 97 (M2).
3. **Align Abstract and Discussion with §3.8's own hedge** on
   model-specificity, and add the B16F10 model details (inoculation site,
   harvest time point, animal age/sex) to §2.8 (M3).
4. **Report the neuron/Schwann composition shift** with its MWU p and one
   sentence of interpretation, applying the paper's own §3.3 compositional
   logic to the validation dataset (M4).
5. **Correct the compartment count** ("eight", twice: §3.8 and Table 5) and
   the per-sample myeloid range ("0.015–0.2%") (m1, m4).
6. **Deposit the purity QC computation** with its exact criterion, fix the
   "non-overlapping" description, and reconcile "as deposited" vs "as
   retained by QC" for the cell count (m2, m3).
7. **Log silently dropped marker genes** and resolve the Pgp9.5/PgP9.5
   casing divergence between scripts 31 and 32 (m5).
8. **Add one short paragraph on the unreported compartment results**
   (Cx3cr1 up in the neuron bins; Gal in Schwann/Endothelial bins), stating
   how they are read under bin impurity (m6).
9. **Cosmetics:** relabel log₂(CP10K+1); regenerate
   `sc_GSE292522_celltype_summary.json` so it includes the Erythrocyte
   compartment; report concordance p-values and drop the "highest at day 1"
   gradient framing; state the FOCUS-gene scope of the per-compartment ISG
   claim (m7–m10).

After 1–9 this section would be, in my judgement, fully sound: the
pseudobulk design is correct, the inference unit is the animal throughout,
and the blind spot is disclosed with unusual candour for this literature.
