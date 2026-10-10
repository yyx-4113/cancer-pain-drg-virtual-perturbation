# R3-A4 Review — Provenance / Reproducibility Auditor

**Manuscript:** `docs/manuscript/manuscript_final.md` (treated as first submission)
**Gate script:** `scripts/98_check_manuscript_numbers.py` — run by me, exit 0,
**`276 checks, 0 failed`** (reproduced 2026-10-10).
**Independence:** no prior review files, no other reviewers' outputs, no
`.workbuddy`, no earlier manuscript versions were opened. All numbers below
were recomputed by me from `results/tables/*` (and, where noted, from the raw
`data/processed/GSE292522/` matrices) with the venv Python.

---

## 1. Verdict

**Minor Revision.**

The quantitative layer of this manuscript is in unusually good shape. I
independently recomputed 262 distinct claims (Tables 1–5, every inline
statistic in §3.1–§3.8, Abstract, Limitations, Conclusions): **256 match the
deposited artefacts at the stated precision**. The producing-script chain is
traceable for essentially every headline number. However, I found one claim
that is **falsified by the manuscript's own artefacts** (the 3/71 paired
p-floor), one **internal count error** ("seven testable compartments"; the
artefact has eight), one **unreproducible number with no producing script**
(11.9% neuron purity), a misleadingly attached median (§3.2), a wording
contradiction on the 55,967-cell attribution, and a data-availability URL
that currently returns **404**. The reproducibility gate itself is
majority-hardcoded and contains tautological presence-only checks; it is a
good drift alarm but not a validity check. None of these issues undermines
the paper's two main claims; all are fixable without new analysis.

---

## 2. Summary of the manuscript

A single-author computational reanalysis of public bulk rat DRG RNA-seq
(GSE149648, Walker-256 tibial bone cancer, n=4/group × 4 groups; GSE160543,
paclitaxel/oxaliplatin/vehicle, n=4/group; same platform GPL24782), plus an
independent single-cell validation (GSE292522, mouse B16F10
cancer-nerve-injury DRG, 4+4 animals, 55,967 cells). Two claims: (1)
*biological* — bone cancer pain and CIPN drive qualitatively different DRG
programmes (myeloid/interferon in Walker-256 D7, neuronal-injury in CIPN),
tested with exhaustive exact permutation over C(8,4)² = 4,900 null pairs;
(2) *methodological (negative)* — in-silico perturbation of these bulk
transcriptomes cannot rank analgesic targets because it fails its own
negative controls (an erythrocyte gene set carries the largest standardised
effect, |z| = 7.39; controls are anti-selected over targets, median |z| 2.47
vs 1.07, MWU p = 0.0019; binomial over-representation p = 0.0067) while
remaining adequately powered per test (MDE ≈ 0.105 log-units; largest target
Δ = 0.274 = 2.6× MDE; 31/132 exact p < 0.05). The single-cell validation shows
the interferon programme is Walker-256-model-specific (concordance ρ = 0.098
D7, 0.112 D1, 0.030 D14). Sections 2.9 and the deposited gate script
(`98_check_manuscript_numbers.py`) promise that "a separate checker asserts
every number reported below against the deposited CSV and JSON artefacts."

---

## 3. Major issues

**M1. The claimed paired-contrast p-floor (3/71) is falsified by the
deposited artefacts for the perturbation-reversal endpoint.**
- *Citation:* §2.7 ("the floor is therefore 3/71 ≈ 0.042"), §5 Limitations
  ("3/71 ≈ 0.042 for the two paired contrasts that share a control"), and
  implicitly §2.7's second bullet ("Perturbation reversal difference … under
  the same exhaustive permutation scheme").
- *Evidence:* In `etiology_module_logfc_cross.csv` and
  `etiology_marker_cross.csv` the minimum attained paired
  `p_permutation_exact` is exactly 3/71 = 0.042254 — the floor claim holds
  there. But in `etiology_bootstrap_differences.csv` (per_contrast, knockout)
  **six paired rows attain p = 2/71 = 0.028169**, below the claimed floor:
  CTRL erythrocyte (Bone_D14 vs Bone_D7), CTRL erythrocyte (PTX vs OXA),
  CTRL pan-neuronal (PTX vs OXA), Ca2+ channel a2d (PTX vs OXA),
  Glial / immune (Bone_D14 vs Bone_D7), Injury / regeneration (PTX vs OXA).
  Root cause (code-level): `exact_two_sided_p` (`scripts/_limma.py:420-426`)
  counts ties with a 1e-12 window (`np.abs(null) >= abs(obs) - 1e-12`), while
  the 8-case exchange null's true-split element is recomputed through a
  *different code path* than the point estimate, and the script's own
  self-check only enforces agreement to 1e-9
  (`scripts/22_etiology_comparison.py:314-324`). When the true-split null
  element lands more than 1e-12 below |obs|, the observed statistic fails to
  be counted as a member of its own null, and only the mirror survives
  (ge = 1 → (1+1)/71 = 2/71). For the unpaired 4,900-pair null this cannot
  happen (the null is literally built from the same `sa`/`sb` arrays, so the
  observed pair and its mirror tie exactly; min attained = 7/4901 ≥ 3/4901 ✓).
- *Impact:* none of the quoted headline p/q values change (min
  `q_permutation` = 0.165 is untouched; no quoted paired reversal p appears
  in the text), but a Methods/Limitations statement of exactness is
  contradicted by the deposit — exactly the class of error this paper
  criticises elsewhere.
- *Required fix:* either (a) loosen the tie tolerance in
  `exact_two_sided_p` to match the self-check tolerance (e.g. 1e-9), rerun
  script 22, and verify the floor 3/71 is then attained as minimum; or (b)
  reword §2.7/§5 to state the floor for the logFC/marker endpoints and note
  the reversal endpoint's effective floor of 2/71 under the current
  implementation. Option (a) is preferred and cheap.

**M2. "Seven testable compartments" is contradicted by the deposited
compartment artefact: eight compartments were tested.**
- *Citation:* §3.8 ("the ≥5-cells-in-≥6-animals rule excluded four
  compartments … Among the seven testable compartments, no
  interferon-stimulated gene was up-regulated anywhere") and Table 5 ("7
  testable compartments; myeloid untestable").
- *Evidence:* §3.8's own arithmetic fails: 12 candidate compartments − 4
  excluded (T cells; NP and TH neurons; myeloid) = **8**, not 7.
  `sc_GSE292522_celltype_DE.csv` contains **eight** cell types: Endothelial,
  Erythrocyte, Fibroblast, Neuron, Neuron.NF, Neuron.PE, SatelliteGlia,
  Schwann. The erythrocyte compartment passes the stated rule (≥5 cells in
  6/8 animals: 5, 4, 22, 3, 7, 6, 26, 15 per the per-sample counts in
  `results/review_a3/compartment_cellcounts.json`) and **was** tested by
  limma-trend. The count "7" comes from
  `sc_GSE292522_celltype_summary.json`, whose MWU discrimination lists 7
  compartments only because the erythrocyte bin had <3 Con samples for the
  MWU — a different artefact answering a different question. Additionally
  unreported: the erythrocyte compartment test is run on an **unbalanced**
  2 Con vs 4 B16F10 samples, whereas §2.8 states every retained compartment
  was tested "across the four tumour and four control animals".
- *Impact:* the substantive conclusion (no ISG up-regulated in any
  compartment) survives — I verified the "no ISG up" check over all eight
  compartments in the CSV — but the manuscript misdescribes its own
  exploratory analysis, in the very paragraph whose purpose is to demonstrate
  QC candour.
- *Required fix:* state 8 testable compartments; add erythrocyte to the
  discussion (its DE was computed); disclose the 2-vs-4 imbalance for the
  erythrocyte compartment; reconcile §2.8's "four tumour and four control
  animals" phrasing with the actual per-compartment sample counts.

**M3. The 11.9% neuron-purity claim has no producing script, no artefact,
and does not reproduce under any natural reading of its own criterion.**
- *Citation:* §2.8 ("scored against an independent, non-overlapping criterion
  (Snap25/Tubb3/Syp-positive and myelin-negative), the assigned 'neuron'
  compartment is 11.9% pure").
- *Evidence:* `grep` over all of `scripts/` and `results/` finds no script or
  artefact computing neuron purity (the A3 scripts compute related QC
  quantities but not this one). I recomputed it from the raw 10x matrices
  (`data/processed/GSE292522/`, all 8 samples, 40,537 Neuron-assigned cells
  under the same marker-argmax rule as scripts 31/32):
  - Snap25/Tubb3/Syp **any-transcript-positive** AND Mpz/Pmp22/Mag/Plp1(+Sox10)-negative
    → **17.0%**;
  - Snap25/Tubb3/Syp **all-three-positive** AND myelin-negative → **3.0%**
    (2.99% with Sox10 in the myelin set, 3.04% without).
  Neither reading yields 11.9%; the criterion as written ("positive",
  "myelin-negative") is under-specified (thresholds? any/all? which myelin
  genes?), so the figure is not reproducible as stated.
- *Impact:* the qualitative point (bins are coarse, not validated cell types)
  is strongly supported — by my recomputation and by the fact that the
  "neuron" bin absorbs 40,537/55,967 = 72% of all cells — but a specific,
  unreproducible number in a candour passage is a provenance defect.
- *Required fix:* deposit a small script + JSON that computes the purity with
  an exactly specified criterion, or replace "11.9% pure" with the
  reproducible variant(s) and state the criterion explicitly.

**M4. The reproducibility gate is majority-hardcoded and partly
tautological; it would not catch a wrong-but-consistent derivation.**
- *Citation:* §2.9 ("A separate checker asserts every number reported below
  against the deposited CSV and JSON artefacts"); scripts/98_check_manuscript_numbers.py.
- *Evidence (source audit):* of the 276 checks, roughly **175 are
  hardcoded-literal comparisons** (`chk(artefact_value, <literal from the
  manuscript>)` — e.g. `want_adj`, `want_t`, `want_mod`, `want_mk` at lines
  103–138, 184–196, 206–214, 289–294, 301–317). These catch manuscript↔artefact
  drift only if the literals are updated in lockstep with the manuscript; a
  bug inside script 22 that propagates identically into artefact, manuscript
  and literal sails through. **68 checks are genuinely runtime-derived and
  strong**: Table 2/3 cells are regenerated from the manuscript's own parsed
  markdown and compared to the CSVs at run time (lines 403–457); the
  T1/T2/T5 statistics are recomputed from the CSV and compared to the JSON
  (lines 329–368); the largest-|z| row's identity, KO/OE mirror identity,
  bootstrap-CI coverage and per-compartment "no ISG up" are derived checks.
  **33 checks are tautological** — `in_text()` (lines 58–72) only asserts
  that a literal string appears somewhere in the manuscript, i.e. the
  manuscript compared to itself; for the T1/T2/T5 block the preceding
  recompute gives these meaning, but for the rest ("0.991", "0.0019",
  "46 cells", "0.105", …) presence is not correctness.
- *Gate coverage gaps* (manuscript numbers the gate never checks
  numerically): |z| = 7.39 and z = −5.78 (the gate checks only that the
  max-|z| row is the erythrocyte control, not its value); pan-neuronal
  z = −3.74; seven of eight Table 4 rows (only interferon checked); the
  Table 1 logFC column and 21 of 22 q values; median |logFC|
  0.073/0.158/0.119; sample-correlation range 0.913–0.994; PCA 46.8/17.4/7.2/4.8;
  next-highest muscle share 0.16%; concordance ρ = 0.112 (D1) and 0.030 (D14)
  (only D7 checked); Atf3 +0.84; Ifit2 −0.048/q 0.042; Oas1c q 0.014; Gfap
  +0.004/q 0.006; Gal +1.47/+2.54/+1.42; Ifit1/Isg15 compartment q; U = 2/14;
  Cohen's d = 1.37/2.38 (compared JSON-to-hardcode only, JSON not recomputed
  by the gate); the 46/1 cell counts (text-presence only); "seven testable
  compartments"; 11.9% purity; the 2/71-vs-3/71 floor (M1). I verified all of
  these myself (Section 5); most pass — but two of them (M1, M2) fail,
  which is exactly the demonstration that the gate's coverage claim
  ("asserts every number") overstates.
- *Required fix:* (a) make the paired-reversal null reuse the point-estimate
  code path or widen the tie tolerance (M1), then rerun; (b) add runtime
  checks for the numbers listed above (cheap: all are one CSV lookup each);
  (c) replace the pure-presence `in_text` strings for unverified numbers with
  derived comparisons, or say in §2.9 that the checker covers Tables 2/3 and
  the negative-result layer cell-by-cell plus N literals, rather than
  "every number".

---

## 4. Minor issues

**m1. §3.2 median |logFC| (0.073 / 0.158 / 0.119) are medians over all
13,893 genes, but the sentence attaches them to the q < 0.05 gene sets.**
Recomputed: median |logFC| among q<0.05 genes is 0.340 (D7, n=4,600) and
1.199 (D14, n=235); D1 has **zero** genes at q<0.05, so its "median" over
that set does not exist. The quoted values match the whole-gene median
exactly (0.0731 / 0.1576 / 0.1187). Fix the sentence ("median |logFC| across
all tested genes was …") or quote the correct per-set values.

**m2. 55,967-cell attribution is stated two incompatible ways.** §2.8
(Methods): "55,967 cells as deposited"; §3.8: "55,967 cells as retained by
our parsing and quality-control pipeline". Script 31 performs **no cell-level
QC** — `total_cells` is the raw sum of the eight matrix column counts
(6,514+7,322+6,019+5,817+6,175+7,480+9,826+6,814 = 55,967), i.e. as
deposited. The §3.8 phrasing implies a retention filter that does not exist.
Pick one (the true) description: "as deposited".

**m3. Data-availability GitHub URL is not accessible (HTTP 404).**
`https://github.com/yyx-4113/cancer-pain-drg-virtual-perturbation` returns
404 (the account `yyx-4113` exists; the repo is private or not yet created).
The Zenodo archival ("to be archived at Zenodo upon acceptance") is clearly
flagged as pending — good — but the statement currently implies the
repository is the live deposit while it is unreachable by a reviewer. Either
make the repo public at submission or reword ("code available from the
author on request; will be deposited at Zenodo upon acceptance"). Note this
also means the 46-cell/T-cell counts (M-adjacent) currently have their only
artefact (`results/review_a3/compartment_cellcounts.json`) outside
`results/tables/` and outside any accessible deposit.

**m4. Script 33's JSON carries stale internal text that contradicts its own
computed fields.** `negative_result_quantification.json` →
`T5_reversal_score_power.interpretation` says "**26**/132 exact permutation
p<0.05" while the same object's field says `n_exact_permutation_p_lt_0.05`
= **31** (recomputed: 31 ✓), and the T5 docstring asserts "d = 2.381 at 80%
power" for the unpaired n=2+2 design, which actually gives d = **5.653**
(2.381 is the n=4+4 value used by T3; the code itself uses
`solve_power(nobs1=2)` correctly). Not manuscript-facing, but these are
hand-maintained strings drifting inside a provenance artefact — the exact
failure mode the paper preaches against. Regenerate the interpretation
strings from the computed fields.

**m5. Table 3's "—" cells are unexplained.** The omitted cells exist in
`etiology_module_logfc_cross.csv` and do not contradict anything (e.g. CTRL
erythrocyte Bone_D14 = +0.73, PTX = −0.80; GABAergic Bone_D14 = −0.19), but
no selection criterion for the dashes is stated. One sentence ("— indicates
|logFC| < X" or "contrasts not shown") suffices.

**m6. "55% drop" in the myeloid fraction is 55.8% at full precision**
(0.000517/0.00117 = 0.442; recomputed from the per-sample composition CSV).
"55%" is defensible rounding of the JSON's `myeloid_ratio_B16F10_over_Con`
= 0.44, but "56%" would be exact. Cosmetic.

**m7. Two Table 1 values sit exactly on rounding boundaries**: Prostanoid/COX
mean t = 1.0650 → "+1.07" and K⁺/HCN mean t = −5.6350 → "−5.64" (half-up);
Python's default half-even would give 1.06/−5.63. Not an error, but the
rounding convention should be consistent (the same table uses half-up
elsewhere, so keep it and note it in the deposit README).

**m8. Declarations audit (all present, all plausible):** Funding (none) ✓;
Conflicts (none) ✓; CRediT — single author holding all 14 roles; complete
but a venue may query the full list for a sole author; Acknowledgements ✓;
Generative-AI declaration present, specific ("code debugging and manuscript
language editing"), with responsibility statement ✓. One consistency note:
the AI declaration says "every numerical claim was checked against the
deposited artefacts" — accurate for the gate-covered numbers, overstated for
M1–M3 as of this review; the fixes above make it true.

---

## 5. Numbers I independently recomputed

Format: claim | manuscript | recomputed (4+ dp where relevant) | match?
All from `results/tables/*` unless noted. 262 claims checked; 256 match.
Failing/note rows bolded.

**Table 2 — `etiology_marker_cross.csv` (36 logFC cells + Δ + q, all rows):**

| Claim | MS | Recomputed | Match |
|---|---|---|---|
| Cd68 Bone_D7/D14/PTX/OXA | +2.84/+2.49/−0.58/+0.40 | +2.8448/+2.4870/−0.5847/+0.3989 | ✓ |
| Cd68 Δ (D7 vs PTX), q | +3.42, 0.012 | +3.4183, 0.012487 | ✓ |
| Aif1 4 aetiologies | +2.19/+1.95/−0.88/−0.46 | +2.1859/+1.9547/−0.8766/−0.4616 | ✓ |
| Aif1 Δ, q | +3.08, 0.012 | +3.0757, 0.012487 | ✓ |
| Gal 4 aetiologies | −0.67/+0.01/+1.92/−0.35 | −0.6704/+0.0050/+1.9217/−0.3539 | ✓ |
| Gal Δ, q | −2.59, 0.012 | −2.5921, 0.012487 | ✓ |
| Ccl2 Δ (D7 vs OXA), q | +1.39, 0.017 | +1.3907, 0.017322 | ✓ |
| Il1b Δ (D7 vs PTX), q | +1.23, 0.018 | +1.2305, 0.018219 | ✓ |
| Calca logFCs, Δ, q | −0.27/+0.50/−0.77, 0.027 | −0.2693/+0.5005/−0.7698, 0.027186 | ✓ |
| Scn10a logFCs, Δ, q | −0.46/+0.24/−0.70, 0.026 | −0.4624/+0.2398/−0.7022, 0.026046 | ✓ |
| Atf3 PTX/D7 logFC (§3.5) | +1.13 / +0.22 | +1.128 / +0.217 | ✓ |
| Cd68/Aif1 D14 Δ (§3.5) | +3.07 / +2.84 | +3.0710 / +2.8376 | ✓ |

**Table 3 — `etiology_module_logfc_cross.csv` (all shown cells + Δ + q):**

| Claim | MS | Recomputed | Match |
|---|---|---|---|
| Glial/immune Δ (D7 vs PTX), q | +2.70, 0.0073 | +2.6951, 0.007345 | ✓ |
| Cytokine Δ vs OXA (§3.5) | +1.75 | +1.7487 (Bone_D7 vs OXA) | ✓ |
| Cytokine Δ (D7 vs PTX) | +1.70 | +1.7024 | ✓ |
| TLR Δ | +1.28 | +1.2755 | ✓ |
| Neurotrophin Δ | −0.97 | −0.9717 | ✓ |
| Ca²⁺ α2δ Δ | −0.60 | −0.6001 | ✓ |
| Neuropeptides Δ, q | −1.35, 0.015 | −1.3450, 0.014528 | ✓ |
| Pan-neuronal Δ | −0.69 | −0.6907 | ✓ |
| Glutamate/NMDA Δ | −0.64 | −0.6354 | ✓ |
| Na⁺ Δ | −0.65 | −0.6537 | ✓ |
| GABAergic Δ | −0.49 | −0.4950 | ✓ |
| **Erythrocyte Δ (D7 vs OXA), q** | **−2.48, 0.0073** | **−2.4831, 0.007345 (unpaired)** | **✓** |

**p-floor conventions:**

| Claim | MS | Recomputed | Match |
|---|---|---|---|
| Unpaired null size | 4,900 | n_null = 4,900 in all unpaired rows | ✓ |
| Paired null size | 70 | n_null = 70 in all paired rows | ✓ |
| Unpaired floor | 3/4,901 | attained min = 7/4,901 ≥ 3/4,901; mirror tie exact | ✓ (floor correct) |
| **Paired floor (logFC/marker)** | **3/71** | **min attained = exactly 3/71 = 0.042254** | **✓** |
| **Paired floor (reversal)** | **3/71** | **6 rows attain 2/71 = 0.028169** | **✗ → M1** |

**Negative-result layer — `etiology_bootstrap_differences.csv` +
`negative_result_quantification.json` (all recomputed from the CSV):**

| Claim | MS | Recomputed | Match |
|---|---|---|---|
| Raw tests per scope | 264 | 264 | ✓ |
| Dedup KO universe | 132 | 132 (102 target + 30 control) | ✓ |
| KO/OE mirror r | −1.000 | −1.000000 (corr of Δ_KO, Δ_OE) | ✓ |
| max\|p_KO − p_OE\| | 0 | 0.0 | ✓ |
| No Δ at q<0.05 | 0 | 0 | ✓ |
| Min q | 0.165 | 0.1651 | ✓ |
| Glial/immune Δ, CI, p | 0.274, 0.203–0.287, 0.0014 | 0.2743, 0.2035–0.2866, 0.001429 | ✓ |
| Randset sig raw / dedup | 64 / 32 | 64 / 32 | ✓ |
| Controls / targets among sig | 14 / 18 | 14 / 18 | ✓ |
| All sig CIs exclude 0 | all | 32/32 | ✓ |
| Median \|z\| controls / targets | 2.47 / 1.07 | 2.4702 / 1.0729 | ✓ |
| MWU two-sided p | 0.0019 | 0.001915 | ✓ |
| Targets' percentile | 55.7th | 55.7 (pct of \|z\| rank) | ✓ |
| Binomial p | 0.0067 | 0.00671 (14/32 vs base 30/132 = 22.7%) | ✓ |
| Median bootstrap SE | 0.019 | 0.018505 | ✓ (2 dp) |
| MDE d (n=2+2, 2-sided, 80%) | (implicit) | 5.653 (solve_power) | ✓ chain |
| MDE units | ≈0.105 | 5.653 × 0.018505 = 0.10460 | ✓ |
| Largest target Δ / ratio | 0.274 / 2.6× | 0.2743 / 2.62 | ✓ |
| Exact p<0.05 count | 31/132 | 31 | ✓ |
| Erythrocyte \|z\| D14 vs D7 | 7.39 | 7.3878 (largest \|z\| in table, kind=control) | ✓ |
| Erythrocyte z PTX vs OXA | −5.78 | −5.7821 | ✓ |
| Pan-neuronal z (randset-sig) | −3.74 | −3.7362 (Bone_D7 vs PTX) | ✓ |
| CI covers observed | implied | 264/264 | ✓ |

**Discovery (§3.1–3.4, Table 1, Table 4):**

| Claim | MS | Recomputed | Match |
|---|---|---|---|
| Genes GSE149648 / shared | 13,893 / 13,601 | 13,893 / 13,601 | ✓ |
| DE adj q<0.05 D1/D7/D14 | 0/4,600/235 | 0/4,600/235 | ✓ |
| D7 adj q<0.10 | 5,867 | 5,867 | ✓ |
| **Median \|logFC\| D1/D7/D14** | **0.073/0.158/0.119** | **all-gene median 0.0731/0.1576/0.1187; q<0.05-set median 0.340 (D7), 1.199 (D14), n/a (D1)** | **✓ value / ✗ attachment → m1** |
| Sample-corr range | 0.913–0.994 | 0.9133–0.9935 (off-diagonal) | ✓ |
| PC1–PC4 | 46.8/17.4/7.2/4.8 | 46.85/17.45/7.23/4.76 | ✓ |
| Worst muscle sample/share | GSM4506604, 6.26% | GSM4506604, 6.26% | ✓ |
| Next-highest share | 0.16% | 0.16% | ✓ |
| Samples above 1% | 1 | 1 | ✓ |
| Muscle-drop D7 count, Jaccard, ρ | 3,041, 0.584, 0.912 | 3,041, 0.584, 0.912 | ✓ |
| Table 1 (22 rows: mean t, mean logFC, q) | e.g. IFN +24.92/+2.61/6.8e-5 | +24.923/+2.614/6.8e-5; all 22 rows match | ✓ (m7 boundary notes) |
| §3.3 marker logFC D7 | Cd68 +2.86, Aif1 +2.18, Gfap +1.33, Ccl2 +0.81, Il1b +0.80 | 2.862/2.183/1.331/0.812/0.800 | ✓ |
| Table 4 (8 rows: reversal, z, %oracle, p, q) | e.g. IFN +0.228/+4.61/75.4/0.0020/0.011 | +0.2282/+4.6111/75.427/0.001998/0.010989; all 8 rows match | ✓ |
| ρ(network, additive) | 0.991 | 0.991 | ✓ |
| Mean \|reversal\| net/add | 0.057/0.033 | 0.057/0.033 (ratio 1.727 ≈ "70%") | ✓ |
| Mean % of oracle | 20.5% | 20.5 | ✓ |

**Sensitivity (§3.7):**

| Claim | MS | Recomputed | Match |
|---|---|---|---|
| Grid points | 60 + 15 = 75 | n_grid_points = 75 | ✓ |
| Discovery ρ min/median | 0.959/0.988 | 0.959/0.988 | ✓ |
| Δ ρ min | 0.989 | 0.989 | ✓ |
| Top-5 identical | all settings | min top5 overlap = 5 | ✓ |
| DE q<0.05 across spans | 316–325 | 316–325 (incl. trend off = 323) | ✓ |
| ρ(t vs reference) | ≥0.9999 | 0.999996–1.000000 (`sensitivity_limma_trend.csv`) | ✓ |
| Cross-study ρ net vs add | 0.929 | 0.929 | ✓ |

**Single-cell (§3.8, Table 5) — plus independent recompute from raw matrices
where marked ⟂:**

| Claim | MS | Recomputed | Match |
|---|---|---|---|
| Total cells | 55,967 | 55,967 = Σ matrix columns (⟂ raw mtx) | ✓ (wording → m2) |
| Pseudobulk DE up/down | 1,517/703 | 1,517/703 | ✓ |
| Sprr1a logFC, t | +1.55, 54.5 | +1.5542, 54.5 | ✓ |
| Atf3 logFC | +0.84 | +0.8402 | ✓ |
| Ifit2 logFC, q | −0.048, 0.042 | −0.0478, 0.0423 | ✓ |
| Oas1c q | 0.014 | 0.0138 | ✓ |
| Gfap logFC, q | +0.004, 0.006 | +0.0043, 0.0062 | ✓ |
| Cd68 / Aif1 "unchanged" | Δ≈0 | −0.0005 (q 0.78) / +0.0009 (q 0.92) | ✓ |
| Gal whole-DRG | +1.42 | +1.4207 | ✓ |
| ρ D7 / D1 / D14 | 0.098/0.112/0.030 | 0.0985/0.1115/0.0299; D7 re-derived ⟂ from pseudobulk DE vs DE_D7_adj over 12,158 genes = 0.09846 | ✓ |
| Myeloid fraction Con/B16 | 0.12%/0.05% | 0.117%/0.0517% | ✓ |
| Drop | 55% | 55.8% | ✓ (m6) |
| Cohen's d obs / MDE | 1.37 / 2.38 | 1.365 / 2.381 (both recomputed) | ✓ |
| MWU U, p (2-s / 1-s) | 2/14, 0.114, 0.057 | U = 2, p = 0.114286 / 0.057143 | ✓ |
| Myeloid cells total | 46 | 46 (3+10+12+5+5+3+7+1; per-sample 0.015–0.199% ≈ "0.02–0.2%") | ✓ |
| T cells total | 1 | 1 | ✓ |
| Animals with ≥5 myeloid | 5 of 8 | 5 of 8 | ✓ |
| Compartments tested | **7** | **8** in `sc_GSE292522_celltype_DE.csv` (incl. Erythrocyte) | **✗ → M2** |
| Per-compartment ISG: none up | 0 up | 0 up across all 8 compartments | ✓ |
| Ifit1 SG / Ifit2 Neuron / Isg15 SG q | 0.004/0.046/0.007 | 0.0041/0.0459/0.0071 | ✓ |
| Gal Neuron logFC/q; Gal NF | +1.47, q<0.001; +2.54 | +1.470, 0.0001; +2.540 | ✓ |
| Scn10a Neuron q | 0.012 | 0.0117 | ✓ |
| Gfap SG logFC | 0.000 | 0.000 | ✓ |
| **Neuron-bin purity** | **11.9%** | **17.0% (any-positive) / 3.0% (all-three) ⟂; no script/artefact** | **✗ → M3** |
| Excluded compartments | 4 (T; NP; TH; myeloid) | 4, correct per rule | ✓ (but arithmetic 12−4≠7) |

**Software/tests/declarations:** requirements.txt versions match §2.9
exactly ✓; "test suite of 14 groups" = t1–t14 in `99_unit_tests.py` ✓;
declarations all present (m8); GitHub URL 404 (m3).

**Chain-of-custody summary:** Table 2/3 ← script 22; Table 4 ← script 20/21
pipeline (`virtual_perturbation_scores.csv`) + 22; negative-result layer ←
script 33 (recomputes from script 22's CSV); single-cell ← scripts 30/31/32;
per-day ρ ← script 34; 46-cell counts ← A3_stage2/A3_stage3 (artefact in
`results/review_a3/`, outside `results/tables/`); **11.9% purity ← no
script** (M3); sensitivity ← script 23. Every quoted headline number has a
producing script except the purity figure.

---

## 6. What would make this publishable

Ordered by cost/benefit; nothing here requires new biology.

1. **Fix the paired-permutation tie handling (M1).** Widen the tie tolerance
   in `exact_two_sided_p` to match script 22's 1e-9 self-check (one line),
   rerun script 22, confirm min paired p = 3/71 is then attained, and rerun
   the gate. Alternatively reword §2.7/§5. This removes the only statement in
   the paper that its own artefacts falsify.
2. **Correct the compartment count (M2):** "eight testable compartments",
   add the erythrocyte compartment to the §3.8 narrative, and disclose the
   2-Con-vs-4-B16F10 imbalance for it. Reconcile §2.8's "four tumour and four
   control animals" phrasing.
3. **Make the 11.9% purity reproducible (M3):** deposit a ~40-line script +
   JSON specifying the criterion exactly (which genes, any/all, threshold),
   or replace the figure with the reproducible variants (17.0% / 3.0%).
4. **Upgrade the gate (M4):** add the ~25 uncovered numbers listed in M4 as
   runtime-derived `chk`s (all are single CSV lookups); keep the Table 2/3
   regeneration; reword §2.9 from "asserts every number" to an accurate
   coverage statement. After fix 1, the gate's coverage claim becomes true.
5. **Wording repairs:** m1 (median |logFC| attachment), m2 (55,967 "as
   deposited"), m3 (GitHub 404 — make the repo public or reword),
   m4 (regenerate script 33's stale interpretation strings), m5 (explain
   Table 3 dashes).
6. **Move the compartment cell-count JSON into `results/tables/`** so the
   46/1-cell claims have an artefact inside the deposited table set (they are
   currently only text-presence-checked).
7. Optional polish: state the rounding convention (m7); "56% drop" (m6);
   consider trimming the single-author CRediT list to the roles a venue
   expects (m8).

With items 1–5 done, the manuscript's §2.9 claim that "every number reported
below" is checked would be true rather than aspirational, and I would have no
remaining provenance objection to publication.
