# R4-A3 Review — single-cell / pseudobulk methodology and transcriptomic annotation

**Manuscript:** `docs/manuscript/manuscript_final.md` (treated as a first submission)
**Reviewer role:** independent referee, single-cell / pseudobulk methodology & transcriptomic annotation
**Independence:** I read only `_PANEL_BRIEF.md` under `docs/review/`, the manuscript, the analysis scripts 22/31/32/34/37/98, and the deposited CSV/JSON artefacts. I read no other reviewer output and no earlier manuscript version.

---

## 1. Verdict

**Minor Revision. Publishable.**

The single-cell section (§2.8 / §3.8) is methodologically sound where it matters most: every group
comparison is made at the level of the **animal** (per-animal pseudobulk, per-animal fractions,
per-animal composite scores), the one place where compartment and gene set are defined by the same
markers is disclosed as carrying no compartment information, and the CD45-depletion blind spot is
stated rather than hidden. Across the whole manuscript I recomputed every single-cell number I could
touch and every one matched the artefacts. The gate (`scripts/98_check_manuscript_numbers.py`) runs to
`287 checks, 0 failed`. What remains is a short list of cosmetic-to-minor defects, none of which
affects a conclusion. There are **no major issues**.

---

## 2. Summary of the manuscript (shows you actually read it)

The paper is a single-author reanalysis of two public bulk rat DRG RNA-seq series on the same platform
(GSE149648, Walker-256 tibial bone-cancer pain, sham/d1/d7/d14; GSE160543, paclitaxel/oxaliplatin
neuropathy; n = 4 per group), plus an independent mouse single-cell validation (GSE292522, B16F10
cancer-nerve-injury, 4 + 4 animals, 55,967 cells).

It makes two claims of unequal weight. (1) *Biological:* bone-cancer pain and CIPN drive qualitatively
different DRG programmes — a type-I interferon / myeloid programme peaking at day 7 in the Walker-256
tibia model versus a neuronal-injury programme in CIPN — tested on canonical markers and 22 modules by
an **exhaustive label permutation** (C(8,4)² = 4,900 pairs), with the muscle-contamination covariate
carried through the within-study DE. (2) *Methodological (negative):* the same data do **not** support
ranking analgesic targets by in-silico perturbation; once negative-control modules, degree-matched
nulls, an oracle anchor and an additive/ridge baseline are applied, the layer is reproducible but
**anti-selective** (median |z| 2.47 for controls vs 1.07 for targets, MWU p = 0.0019), correlates 0.991
with the additive model, and ranks an erythrocyte control above every target.

My remit is the single-cell layer (§2.8/§3.8 and the underlying artefacts) plus the bulk limma-trend
claims of §3.2/§3.3. In §3.8 the argument is that the whole-DRG pseudobulk shows a neuronal-injury,
not interferon, response (1,517 up / 703 down at q < 0.05; Sprr1a +1.55, t = 54.5), that the
interferon signal is therefore model-specific, and that a per-compartment per-animal pseudobulk
(8 testable compartments, myeloid untestable at ~46 cells total) cannot find any ISG up-regulated
anywhere — while flagging that CD45 depletion makes a myeloid-restricted interferon programme
undetectable and that the satellite-glia label is not GFAP-validated (7 Gfap UMI across eight animals,
zero in five).

---

## 3. Major issues

**None.** No defect in the single-cell/pseudobulk layer, and none in the bulk limma-trend claims of
§3.2/§3.3 that I could verify, would stop me recommending acceptance at a soundness-oriented SCIE
journal. Specifically:

- **Inference unit.** Every group comparison in the single-cell layer uses the animal as the unit.
  Whole-DRG pseudobulk (script 31) fits limma-trend over the 8 animals; per-compartment pseudobulk
  (script 32, `design = np.column_stack([ones, B16F10_indicator])`, `fit_limma` across the present
  samples) fits over animals; the composition and composite-ISG tests are Mann–Whitney / Welch on
  **per-animal** values (`sc_GSE292522_composition.csv`, `sc_GSE292522_celltype_scores.csv`,
  `sc_GSE292522_compartment_loo.json`). No cell-level test is used as a group comparison. I found **no
  single-cell pseudoreplication**.
- **Per-compartment construction.** `scripts/32_percelltype_sc.py:131-138` builds each compartment's
  pseudobulk as the **sum of counts of that compartment's cells within each animal** (`agg[ctype][prefix] =
  mat[:, mask].sum(axis=1)`), then log₂(CP10k+1) and limma-trend across animals. Cells are assigned by
  argmax to exactly one broad compartment, so no cell is double-counted; neuron subtype bins are nested
  inside Neuron by design (see Minor 9).
- **Double-dipping.** The only place where the bin definition and the tested genes coincide is the
  neuronal-injury localisation (galanin/Scn10a in the neuron bin), and the manuscript explicitly states
  that this "carries no compartment information" (§3.8/line 159, `scripts/32_percelltype_sc.py:55-56`).
  The ISG panel used for the central compartment claim shares no gene with any annotation marker panel
  (`BROAD`/`NEURON_SUB`), so the ISG endpoint is not circular.
- **Bulk limma-trend §3.2/§3.3.** Reconfirmed from `DE_D1/D7/D14_vs_Sham_adj.csv` and
  `module_level_signal.csv` (see §5).

---

## 4. Minor issues

**1. The erythrocyte compartment is described with the tumour/control counts reversed.**
*Citation:* §2.8 (line 97) and §3.8 (line 157): "the erythrocyte bin clears the cell-count rule but is
tested with **two tumour and four control** animals".
*Evidence:* from `sc_GSE292522_composition.csv` × `sc_GSE292522_gfap_purity_check.csv`, the erythrocyte
bin has ≥5 cells in S1 (5 cells), S3 (22) — both **control** — and S5 (7), S6 (5), S7 (26), S8 (15) —
all four **tumour**. So it is **two control and four tumour**; `sc_GSE292522_celltype_scores.csv` rows
122-139 confirm the present samples are {S1,S3} Con + {S5,S6,S7,S8} B16F10, i.e. the DE was fitted on
2 + 4, not 4 + 2.
*Required fix:* swap the labels ("two control and four tumour").

**2. "log₂(CPM + 1)" mislabels the scaling actually used.**
*Citation:* §2.8 (line 95) and §3.8 (line 157).
*Evidence:* both `scripts/31_process_sc.py:126` and `scripts/32_percelltype_sc.py:154` use
`M.div(M.sum(0))*1e4 + 1.0`, i.e. counts **per 10,000**, not per million. The pipeline is internally
consistent, but the label is wrong.
*Required fix:* call it "log₂(counts-per-10,000 + 1)" (or switch the constant to 1e6).

**3. The per-compartment pseudobulk does not collapse duplicate gene symbols, unlike the whole-DRG
pseudobulk.**
*Citation:* `scripts/31_process_sc.py:100` (`s.groupby(s.index).sum()`) versus
`scripts/32_percelltype_sc.py:152` (no `groupby`; raw `canon_sym` index retained).
*Evidence:* each of the 8 feature files has 55,335 rows with **105 duplicated symbols / 107 extra rows**
(mostly `Gm####`, `Ncrna`, `Pakap`, `Rmrp`, `Mir1839`, etc.). In script 32 these rows survive into
`M`, are summed again in the library-size denominator (`M.sum(axis=0)`), and `fit.loc[g]` would return
multiple rows for a duplicated focus gene. None of the 105 duplicates is in `FOCUS`, `MARKERS`, `ISG`,
`INJURY` or `GLIAL`, so **none of the reported numbers is affected** — but the two code paths normalise
non-identically.
*Required fix:* add the same `.groupby(index).sum()` dedup (or explicitly assert the duplicate set is
disjoint from every tested panel).

**4. The composition displacement rests on a coarse argmax proxy and is reported without a
marker-expression caveat.**
*Citation:* §3.8 (line 155): "the neuronal share fell from 79% to 67% … the Schwann share rose from
19% to 29%".
*Evidence:* I reproduce 78.59% → 66.95% and 18.85% → 29.28% with per-group Mann–Whitney two-sided
p = 0.0286 for both (see §5). The manuscript already says the assignment is "coarse" and panels
"overlap", but it does not state that a compositional shift inferred from argmax-over-markers can in
principle be generated by differential expression of the markers themselves. A quick check against
`sc_GSE292522_pseudobulk_DE.csv` argues *against* a pure artefact here (neuron panel mean logFC =
+0.02, i.e. flat; Schwann panel mean = +0.18, up), which is worth saying explicitly. Note also that with
n = 4 vs 4 the MWU floor is 0.0286, and 8 fractions are tested without multiplicity control.
*Required fix:* add one sentence stating the marker-expression confound, note the whole-DRG panel means
are flat (neuron) / up (Schwann), and mark the composition tests exploratory.

**5. "no interferon-stimulated gene was up-regulated anywhere" (per-compartment) is broader than the
genes tested.**
*Citation:* §3.8 (line 157).
*Evidence:* the per-compartment gene-level test in `sc_GSE292522_celltype_DE.csv` reports only the 10
ISGs that appear in `FOCUS` (`Ifit1,Ifit2,Ifit3,Isg15,Oas1a,Oas1b,Bst2,Mx1,Stat1,Irf7`), whereas the
composite ISG score (`sc_GSE292522_celltype_scores.csv`) covers the full ~38-gene `ISG` set. The
whole-DRG claim (all 38, via `sc_GSE292522_summary.json:ISG_up_in_B16F10 = []`) is fine.
*Required fix:* scope the compartment sentence to "among the ISGs reported" or anchor it on the
composite score, which covers the full set.

**6. Myeloid fraction lower bound and drop percentage are slightly off.**
*Citation:* §3.8 (line 155): "~0.12% in controls and ~0.05% in B16F10, a **55%** drop"; line 157:
"**0.02**–0.2% of cells per sample".
*Evidence:* recomputed per-animal fractions span **0.0147%** (S8) to 0.1994% (S3); the drop is
0.052%/0.117% = ratio 0.442 → **55.8%**, i.e. ≈56%.
*Required fix:* "0.01–0.2%" and "a 56% drop".

**7. "U = 2/14" is a non-standard notation for a 4 × 4 Mann–Whitney.**
*Citation:* §3.8 (line 155) and Table 5.
*Evidence:* for n = 4 vs 4 the U statistic ranges 0–16; scipy returns **U = 2.0**, two-sided
p = 0.1143, one-sided ("less") p = 0.0571. "2/14" is not the U/max ratio.
*Required fix:* write "U = 2, n = 4 + 4" (or give the rank-sum statistic).

**8. §3.3 marker values are quoted from the adjusted within-study fit while Table 2 shows the
unadjusted per-sample means — flagged in §2.3 but easy to misread.**
*Citation:* §3.3 (line 115): "Cd68 +2.86, Aif1 +2.18, Gfap +1.33, Ccl2 +0.81, Il1b +0.80"; Table 2:
Cd68 +2.84, Aif1 +2.19, Gfap +1.50, Ccl2 +0.90, Il1b +0.80.
*Evidence:* `etiology_marker_de.csv` (adjusted, gate-checked) gives Cd68 +2.86, Aif1 +2.18, Gfap +1.33,
Ccl2 +0.81; Table 2's per-sample-mean path gives +2.84 / +2.19 / +1.50 / +0.90. The two paths are
deliberately different (§2.3) but the reader hits both without a pointer.
*Required fix:* a one-phrase cross-reference at §3.3 to the two logFC sources.

**9. Nested compartments are counted among the "8 testable compartments".**
*Citation:* §2.8 (line 97): "twelve candidate compartments — neurons (with NF … and TH subtypes) …".
*Evidence:* `Neuron.NF`/`Neuron.PE` cells are a subset of the `Neuron` bucket
(`scripts/32_percelltype_sc.py:136`), so `Neuron`, `Neuron.NF` and `Neuron.PE` are not disjoint and
are not independent for BH purposes.
*Required fix:* note that subtype bins are nested within Neuron when describing the compartment set.

---

## 5. Numbers I independently recomputed

All values recomputed from the deposited artefacts with
`C:/Users/Administrator/.workbuddy/binaries/python/envs/default/Scripts/python.exe`.

| Claim (manuscript) | Manuscript | Recomputed | Match? |
|---|---|---|---|
| Neuron share, control → B16F10 | 79% → 67% | 78.59% → 66.95% | Yes |
| Schwann share, control → B16F10 | 19% → 29% | 18.85% → 29.28% | Yes |
| Neuron / Schwann Mann–Whitney p | 0.029 (both) | 0.0286 (both) | Yes |
| Myeloid fraction, control → B16F10 | 0.12% → 0.05% | 0.1170% → 0.0517% | Yes |
| Myeloid Mann–Whitney | U = 2/14, p = 0.114 two-sided / 0.057 one-sided | U = 2.0, 0.1143 / 0.0571 | Yes (notation aside) |
| Myeloid observed Cohen's d | 1.27 | 1.267 (n−2 pooled) | Yes |
| Myeloid MDE Cohen's d (80% power) | 2.38 | 2.3808 (exact noncentral-t; 1.98 normal-approx) | Yes |
| Myeloid cells in the whole dataset | 46 | 46 (3+10+12+5+5+3+7+1) | Yes |
| T cells in the whole dataset | 1 | 1 (S7 only) | Yes |
| Total cells | 55,967 | 55,967 | Yes |
| Erythrocyte bin group counts | "two tumour, four control" | **2 control (S1,S3), 4 tumour (S5-8)** | **No** |
| "8 testable compartments" | 8 | 8 (T, NP, TH, myeloid excluded) | Yes |
| Satellite-glia Gfap UMI | 7, zero in five animals | 7; zeros in S1,S2,S3,S5,S8 = 5 | Yes |
| Pseudobulk DE, B16F10 vs Con (q < 0.05) | 1,517 up / 703 down | 1,517 / 703 | Yes |
| Sprr1a | +1.55 log₂, t = 54.5 | +1.5539, t = 54.49 | Yes |
| Atf3 | +0.84 | +0.8355 | Yes |
| Ifit2 | −0.048 log₂, q = 0.042 | −0.0476, q = 0.0417 | Yes |
| Oas1c | q = 0.014 | −0.0291, q = 0.0136 | Yes |
| Gfap (pseudobulk) | +0.004, q = 0.006 | +0.0036, q = 0.0061 | Yes |
| Cd68 / Aif1 (pseudobulk) | unchanged | −0.0005 / +0.0009 (q = 0.78 / 0.92) | Yes |
| Concordance ρ (D7 / D1 / D14) | 0.098 / 0.112 / 0.030 | 0.0985 / 0.1115 / 0.0299 | Yes |
| Neuron-bin ISG: Ifit2 q | 0.046 | 0.0459 | Yes |
| Satellite-glia ISG: Ifit1 q / Isg15 q | 0.004 / 0.007 | 0.0037 / 0.0071 | Yes |
| Composite ISG Welch p (neuron / satellite glia) | 0.93 / 0.14 | 0.9255 / 0.1418 | Yes |
| LOO sign-consistency, SG & Schwann | 8/8 drops; significant ≤2 | SG 8/8, p<0.05 in 1; Schwann 8/8, p<0.05 in 2 | Yes |
| Galanin (whole-DRG / neuron / NF-neuron) | +1.42 / +1.47 (q<0.001) / +2.54 | +1.425 / +1.472 (q=1e-4) / +2.537 | Yes |
| Scn10a, neuron bin | q = 0.012 | −0.247, q = 0.0117 | Yes |
| §3.2 DE counts d1/d7/d14 (q<0.05 adj) | 0 / 4,600 / 235 | 0 / 4,600 / 235 | Yes |
| §3.2 median \|logFC\| d1/d7/d14 | 0.073 / 0.158 / 0.119 | 0.0731 / 0.1576 / 0.1187 | Yes |
| §3.2 significant-set median \|logFC\| d7/d14 | 0.34 / 1.20 | 0.3398 / 1.1994 | Yes |
| §3.3 interferon module mean t (D7) | +24.9 (q = 6.8e-5) | +24.923 (q = 6.76e-5) | Yes |
| §3.3 housekeeping control (D7) | +0.65 (q = 0.52) | +0.653 (q = 0.524) | Yes |
| §3.1 sample correlation range | 0.913–0.994 | 0.9133–0.9935 | Yes |
| §3.1 PCA PC1–PC4 | 46.8 / 17.4 / 7.2 / 4.8% | 46.85 / 17.45 / 7.23 / 4.76% | Yes |
| §3.1 muscle share, worst / next | 6.26% / 0.16% | 6.257% / 0.160% | Yes |
| Permutation floors (unpaired / paired logFC / paired reversal) | 3/4901, 3/71, 2/71 | all p are exact k/4901 or k/71; paired minima 0.04225 and 0.02817 | Yes |
| T1 median \|z\| control vs target; MWU p | 2.47 vs 1.07; p = 0.0019 | 2.4702 vs 1.0729; p = 0.001914 | Yes |
| T2 controls among randset-significant; binomial p | 14/32 vs base 30/132; p = 0.0067 | 14/32; base 0.2273; p = 0.00671 | Yes |
| T5 MDE, targets > MDE, exact sig | 0.052; 41/102; 31/132 | 0.0518; 41; 31 | Yes |
| Gate run | 287 checks, 0 failed | 287 checks, 0 failed (exit 0) | Yes |

Reproducibility gate audit: `scripts/98_check_manuscript_numbers.py` derives its statistical
expectations at run time from the CSVs/JSON (T1 MWU, T2 binomial, T5 MDE recomputed from
`etiology_bootstrap_differences.csv`; Tables 2/3 regenerated cell-by-cell from the CSVs). The remaining
literals are manuscript-description constants (sample sizes, gene counts) that are separately
artefact-checked. This is the right design and I did not find a stale literal protecting an error.

---

## 6. Publishability assessment

**The manuscript is publishable in its current scientific substance.** The single-cell/pseudobulk
layer answers the question it poses, at the correct unit of inference, and states its own limits
(CD45-depletion blind spot, non-GFAP-validated satellite-glia label, nested/coarse bins,
model-specificity rather than generalisation, and the unresolved disagreement with ref [15]). The
negative-result layer and the bulk DE counts are reproducible and gate-checked. The one
factual error I found (Minor 1) is a text transposition inside a caveat and changes no conclusion.

**Ordered minimum to reach "ready":**

1. **Minor 1** — correct the erythrocyte bin to "two control and four tumour" (§2.8 line 97, §3.8
   line 157).
2. **Minor 2** — relabel "CPM" as counts-per-10,000 (§2.8/§3.8).
3. **Minor 5** — scope the per-compartment "no ISG up-regulated" sentence to the genes tested, or
   anchor it on the composite score.
4. **Minor 4** — add the marker-expression confound sentence for the composition displacement and mark
   those tests exploratory.
5. **Minor 3** — add the symbol-dedup (or its disjointness assertion) to `scripts/32_percelltype_sc.py`.
6. **Minors 6–9** — wording/notation fixes (56% not 55%; U = 2, n = 4+4; cross-reference the two
   logFC paths; note nested subtype bins).

Nothing on that list requires re-running the analysis; items 1–2, 5–9 are text/prose, and item 5 is a
defensive assertion. On completion I would accept.
