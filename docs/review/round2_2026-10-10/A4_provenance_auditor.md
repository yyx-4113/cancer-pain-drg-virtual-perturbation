# A4 — Provenance / Recompute Auditor (Round 2, independent)

**Manuscript:** `docs/manuscript/manuscript_final.md`
**Base dir:** `d:/2026.10/不同部位癌痛患者对药物敏感性生信分析/cancer_pain_virtual_perturbation`
**Posture:** first submission. Every number below was recomputed by me from raw
artefacts or raw GEO member files. I did not use `98_check_manuscript_numbers.py`
as evidence (it passes 191/191, which is precisely the problem — see F1).

---

## Verdict up front

The deposited artefacts are in **good shape** and most headline numbers trace
exactly. Of ~230 individual numeric cells I recomputed:

- **Table 1: 22/22 rows exact** (mean t, mean logFC, q, n — all four columns).
- **Table 2: 9/9 Δ and q values exact; 7/7 significant logFC cells exact.**
- **Table 4: 8/8 rows exact on all six columns.**
- **Table 5: 15/16 checked values exact** (one value — Gfap — is mislabelled).
- **§3.1, §3.2, §3.3 module values, §3.6 perturbation statistics: essentially all exact.**

So this is **not** a fabrication paper. The problems are of a different and more
insidious kind, and they cluster in four places:

1. **Two logFC sources are mixed silently across the manuscript** (F2) — the
   Abstract/§3.3 quote muscle-adjusted values, Table 2/§3.5 quote unadjusted
   simple means, for the *same genes in the same contrast*, with no disclosure.
2. **Table 3's CTRL erythrocyte row reports Δ and q from a different contrast
   than the column header says**, and its q is not the q of its own Δ (F3). This
   is the single most serious error I found.
3. **The per-compartment single-cell analysis silently omits the myeloid
   compartment** — the one compartment that carries the paper's central claim —
   and §2.8/§3.8/Table 5 do not say so (F5, F6).
4. **§3.7 and §3.8 contain numbers with no computable source in the deposited
   artefacts** (F8, F9).

None of these overturns the paper's two conclusions. All of them are fixable by
text edits except F5, which needs one extra analysis.

---

# Findings

## F1 — The self-check gate is the only automated guard, and it cannot catch any of the errors below

**【Problem】** `scripts/98_check_manuscript_numbers.py` reports "191 checks,
0 failed / ALL MANUSCRIPT NUMBERS VERIFIED AGAINST ARTEFACTS", but it verifies
arithmetic self-consistency against hardcoded expectations, not correctness
against the raw data, and it does not cover the two newest result families.

**【Evidence】** I ran it: `191 checks, 0 failed`. Grepping the checker for the
new artefacts (`grep -c`):

| artefact | mentions in 98-checker |
|---|---|
| `negative_result_quantification.json` | 1 |
| `sc_GSE292522_celltype_DE.csv` | 1 |
| `sc_GSE292522_pseudobulk_DE.csv` | 3 |
| `sc_GSE292522_celltype_scores.csv` | **0** |
| `sc_GSE292522_celltype_summary.json` | **0** |
| `sc_GSE292522_composition.csv` | **0** |
| `sc_GSE292522_focus_markers.csv` | **0** |

It also asserts manuscript *strings* (`'4,600'`, `'24.9'`, `'+3.42'`, `'0.165'`,
`'32/132'`, …) rather than recomputing them. So the gate that reads
"ALL MANUSCRIPT NUMBERS VERIFIED" is, for the per-compartment layer, a
string-presence test.

**【Why it matters】** Three of the four substantive findings below (F2, F3, F9)
pass this gate cleanly. A reviewer or reader who runs the deposited checker is
told the numbers are verified. That is a false assurance, and it is worse than
having no gate, because it discourages the manual trace that would have caught
them.

**【Specific fix】** Replace the closing line of the checker's output with an
explicit scope statement, and extend coverage:

```
Paste-ready replacement for the checker's final banner:

  191 self-consistency checks passed (table cells vs deposited CSV/JSON,
  arithmetic closure, KO/OE mirror identity, manuscript string presence).
  NOT verified here: agreement of any CSV with the raw GEO member files.
  Per-compartment artefacts sc_GSE292522_celltype_scores.csv,
  sc_GSE292522_celltype_summary.json, sc_GSE292522_composition.csv and
  sc_GSE292522_focus_markers.csv are NOT currently covered.
```

And add, as a new gate: a re-derivation test that recomputes
`etiology_marker_de.csv` and `etiology_marker_cross.csv` from
`data/raw/GSE149648/*.txt.gz` + `data/raw/GSE160543/*.txt.gz` and asserts
agreement to 1e-6, so the two-source problem (F2) becomes machine-detectable.

---

## F2 — Two different logFC sources are used for the same genes in the same contrast, without disclosure

**【Problem】** The manuscript reports per-gene logFC from
`etiology_marker_de.csv` (muscle-adjusted limma-trend) in the Abstract, §3.3 and
Table 1's neighbourhood, and from `etiology_marker_cross.csv` (per-sample simple
mean, **no muscle covariate**) in Table 2 and §3.5 — for the *same genes in the
same contrast* — and never tells the reader which is which. For `Gfap` the two
sources differ by 0.17 log₂ (13%).

**【Evidence】** I extracted both sources and matched them gene-by-gene
(Bone_D7 unless noted). Format: `manuscript value / de.csv (muscle-adj) /
cross.csv (unadjusted simple mean)`:

| Gene | MS prose | de.csv | cross.csv | MS value traces to |
|---|---|---|---|---|
| Cd68 | +2.86 (§3.3, Abstract) | 2.858 | 2.8390 | **de.csv** |
| Aif1 | +2.18 (§3.3, Abstract) | 2.180 | 2.1922 | **de.csv** |
| Gfap | +1.33 (§3.3) | 1.328 | 1.4984 | **de.csv** |
| Ccl2 | +0.81 (§3.3) | 0.806 | 0.8974 | **de.csv** |
| Il1b | +0.80 (§3.3) | 0.804 | 0.8000 | either (coincidence) |
| Cd68 | +2.84 (Table 2) | 2.858 | 2.8390 | **cross.csv** |
| Aif1 | +2.19 (Table 2) | 2.180 | 2.1922 | **cross.csv** |
| Gfap | +1.50 (Table 2) | 1.328 | 1.4984 | **cross.csv** |
| Ccl2 | +0.90 (Table 2) | 0.806 | 0.8974 | **cross.csv** |
| Atf3 | +0.22 (§3.5) | 0.222 | 0.3128 | **de.csv** |
| Atf3 | +1.13 (§3.5) | 1.126 | 1.1219 | **de.csv** |
| Atf3 | +0.31 (Table 2) | 0.222 | 0.3128 | **cross.csv** |

Manuscript locations: line 21 (Abstract), line 115 (§3.3), line 129–131 (§3.5),
Table 2 at lines 239–247.

Artefact provenance I confirmed in code:
- `etiology_marker_de.csv` ← `22_etiology_comparison.py:562` → `de_positive_controls()`,
  design `[1, group, muscle]`, `trend=True` (lines 377–381). Muscle-adjusted.
- `etiology_marker_cross.csv` ← `22_etiology_comparison.py:558` → `cross_marker_logfc()`
  → `module_logfc_splits()` (lines 404–412), which is a bare
  `mean(case) − mean(control)` with **no covariate at all**.

Note that `regenerate_ms_tables.py:23-30` explicitly documents the *Table 2*
choice ("the same artefact that supplies Delta and q, so every Table-2 cell
traces to one source instead of mixing de.csv with cross.csv"). So the mixing is
known and deliberate **inside Table 2** — but §3.3 and §3.5 were not brought
into line, so the manuscript as a whole still mixes.

The worst single instance is `Gfap`: §3.3 says **+1.33**, Table 2 says **+1.50**
for the identical gene and contrast. A reader comparing the abstract to Table 2
finds two different values for the same quantity with no way to tell that one is
covariate-adjusted and one is not.

**【Why it matters】** §2.3 justifies the muscle covariate at length; §3.2's
entire argument is that the covariate changes conclusions. Reporting
*unadjusted* simple means in Table 2 — the table that carries the paper's
headline marker claims (Cd68 Δ = +3.42, Aif1 Δ = +3.08) — without saying so means
the headline cross-aetiology effect sizes are **not** the muscle-adjusted ones the
Methods promise. It also means the Δ column and the per-aetiology columns are
internally consistent (good) but externally inconsistent with §3.3 (bad). A
reviewer who notices will reasonably suspect the larger values were chosen.

**【Specific fix】** Two options; take the first.

*Option A (recommended, text-only).* Add one sentence to §2.7 and one footnote
to Table 2:

```
Paste-ready, new sentence at the end of §2.7 ("Cross-aetiology comparison"):

  Two log-fold-change sources appear in this paper and they are not
  interchangeable. Table 2 and Section 3.5 report the per-sample simple mean
  difference (case minus control) on log2(CPM+1), because that quantity is
  recomputable under every one of the C(8,4) label assignments and therefore
  admits an exact permutation test; no muscle covariate is fitted to it.
  Table 1, Section 3.3 and the Abstract report the muscle-adjusted
  empirical-Bayes moderated-t log-fold change from
  results/tables/etiology_marker_de.csv. The two differ by up to 0.17 log2
  (Gfap, bone day 7: 1.33 adjusted versus 1.50 unadjusted). All
  cross-aetiology Delta values and their q values are computed entirely
  within the unadjusted source and are internally consistent with the
  per-aetiology columns of Table 2.
```

And a Table 2 footnote:

```
  Table 2 log-fold changes are muscle-unadjusted per-sample simple means, the
  quantity on which the exact label-permutation test is defined. Muscle-adjusted
  values for the same genes and contrasts are in
  results/tables/etiology_marker_de.csv and are quoted in Section 3.3.
```

*Option B (analysis).* Add the muscle covariate to `module_logfc_splits()` and
re-derive Table 2, §3.5 and the abstract. This is the methodologically stronger
fix but changes every number in §3.5 and Table 2 and requires re-running the
permutation test. If you take this route, state it in the revision response.

---

## F3 — Table 3's CTRL erythrocyte row reports Δ and q from a different contrast than its own column header, and its q does not correspond to its Δ

**【Problem】** Table 3's column header reads `Δ (D7 vs PTX)`. Eleven of twelve
rows honour it. The CTRL erythrocyte row does not: it shows `−2.48 (D7 vs OXA)`
in a column headed "D7 vs PTX", **and** the q printed beside it (0.0073) is the
q of the D7-vs-OXA test, not the q of the D7-vs-PTX test. A reader who takes the
header literally computes Δ = −0.26 with q = 0.81 — a null — and cannot reconcile
that with the printed −2.48 / 0.0073.

**【Evidence】** I recomputed all 12 rows × 4 logFC columns + Δ + q from
`results/tables/etiology_module_logfc_cross.csv` (132 rows = 22 modules × 6 pairs,
BH denominator 132). Full-precision values for the erythrocyte row:

| quantity | D7 vs PTX | D7 vs OXA |
|---|---|---|
| logFC bone D7 | −1.0540 | −1.0540 |
| delta_logFC | **−0.2583** | **−2.4831** |
| q_permutation | **0.8143** | **0.0073** |
| n_null | 4900 | 4900 |

Manuscript prints: `| CTRL erythrocyte | control | −1.05 | — | — | +1.43 |
−2.48 (D7 vs OXA) | 0.0073 |`. So Δ and q both come from D7-vs-OXA while the
column header says D7-vs-PTX. `regenerate_ms_tables.py:69-71` hardcodes
`mod_delta()` to `aetiology_A=="Bone_D7", aetiology_B=="PTX"`, so the script
would emit `−0.26 | 0.8143` — the manuscript's cell is a **hand-edit made after
generation**, and it is the only cell in any of the five tables that the
regeneration script does not reproduce.

Everything else in Table 3 checks out: 11/11 non-erythrocyte rows match the
artefact to 2 dp on logFC, Δ and q, including the four `0.0073` values
(glial/immune, cytokine/chemokine, TLR, neurotrophin, Ca²⁺ α2δ).

**Secondary defect in the same row: two cells are blanked that should not be.**
`bone D14` and `PTX` are printed `—`, but the artefact has `+0.7326` and
`−0.7958`. I tested the obvious justification (a cell is blanked when the
module's Δ in that column's contrast is non-significant) against 12 cells:

| Module | Column | Printed | Paired-Δ q | ns-rule holds? |
|---|---|---|---|---|
| Neurotrophin | OXA | — | 0.1414 | yes |
| Ca²⁺ α2δ | OXA | — | 0.3156 | yes |
| Neuropeptides | OXA | — | 0.0996 | yes |
| CTRL pan-neuronal | OXA | — | 0.7886 | yes |
| Glutamate/NMDA | Bone_D14 | — | 0.0619 | yes |
| **Voltage-gated Na⁺** | **Bone_D14** | **—** | **0.0491** | **NO — q<0.05, blanked anyway** |
| GABAergic | Bone_D14 | — | 0.0585 | yes |
| **CTRL erythrocyte** | **Bone_D14** | **—** | **0.0321** | **NO — q<0.05, blanked anyway** |

So the ns-blanking rule is applied inconsistently, and it suppresses two
*significant* results (voltage-gated Na⁺ D14-vs-PTX q = 0.0491; erythrocyte
D14-vs-PTX q = 0.0321) while showing 10 non-significant ones. No legend explains
the `—`.

**【Why it matters】** The erythrocyte control is load-bearing for §3.6's central
argument: *"The erythrocyte control also differs between the two datasets at
module level (ΔlogFC = −2.48, q = 0.0073), which is a direct measurement of how
much the two DRG preparations differ in blood content."* But the −2.48 is
**D7-vs-OXA**, whereas the neuronal modules it is compared against in the same
sentence are quoted as **D7-vs-PTX** values (pan-neuronal −0.69, voltage-gated
Na⁺ −0.65). The sentence compares an OXA contrast to PTX contrasts and calls the
result a shared compositional explanation. On the matched PTX contrast the
erythrocyte Δ is −0.26, q = 0.81 — nothing. And the manuscript's own §2.7 states
the PTX/OXA contrast is *paired* on a shared vehicle control with only 70 null
splits (I confirmed `n_null` = 70 for `PTX|OXA` and 4900 for the four bone-vs-CIPN
pairs), i.e. the least powerful of the six.

This is the finding most likely to be read as cherry-picking, and it is the one
place in the paper where a table cell demonstrably disagrees with the deposited
regeneration script.

**【Specific fix】**

1. Split the Δ column so the contrast is never implicit:

```
Paste-ready Table 3 header replacement:

| Module | Kind | logFC bone D7 | logFC bone D14 | logFC PTX | logFC OXA |
  Δ (D7 vs PTX) | q | Δ (D7 vs OXA) | q |
```

2. Fill the erythrocyte row from both contrasts:

```
Paste-ready CTRL erythrocyte row:

| CTRL erythrocyte | control | −1.05 | +0.73 | −0.80 | +1.43 |
  −0.26 | 0.81 | −2.48 | 0.0073 |
```

3. Un-blank the two suppressed significant cells, or state the rule. Paste-ready
note under Table 3:

```
  A cell is left empty only where the gene set is absent from the expression
  matrix. All 22 modules and all four aetiology contrasts are present in
  results/tables/etiology_module_logfc_cross.csv (132 rows); no value has been
  withheld. The erythrocyte row is the only one whose Δ is not D7 vs PTX, and
  both contrasts are now shown.
```

4. Rewrite the §3.6 sentence so the comparison is contrast-matched:

```
Paste-ready replacement for the first sentence of the paragraph beginning
"The erythrocyte control also differs":

  The erythrocyte control shows why the PTX-vs-OXA contrast must be read
  separately from the bone-vs-PTX contrasts: against oxaliplatin the
  erythrocyte module differs by ΔlogFC = −2.48 (q = 0.0073), a direct
  measurement of how much the two DRG preparations differ in blood content,
  whereas against paclitaxel the same module shows no difference
  (ΔlogFC = −0.26, q = 0.81). Neuronal modules show Δ in the same direction as
  the erythrocyte Δ in the bone-vs-PTX contrast, so their apparent
  cross-aetiology difference cannot be fully separated from this compositional
  difference. Because the paclitaxel-versus-oxaliplatin comparison shares a single
  vehicle control it carries only 70 distinct label splits (floor p = 1/71)
  against 4,900 for the bone-versus-CIPN contrasts, and we do not treat the
  erythrocyte result as independent confirmation of the neuronal result.
```

---

## F4 — §3.2 explains the day-7 DE increase by a mechanism that cannot operate

**【Problem】** §3.2 states the muscle covariate "absorbs residual variance and
**raises the prior degrees of freedom**, which increases sensitivity". Adding a
covariate *lowers* residual degrees of freedom. It cannot raise prior df.

**【Evidence】** In `20_multistage_pipeline.py:311-317` the adjusted design is
`[1, group, muscle]` (p = 3) and the unadjusted is `[1, group]` (p = 2), both at
n = 8. So `df_res` goes **6 → 5** (`_limma.py:89`, `df_res = n - p`). The
empirical-Bayes prior df is estimated from the residual variances
(`_limma.py:61-73`), not from the design rank, so it is not "raised" by the
covariate either. Recomputed from the deposited DE tables:

| Contrast | model | n genes | q<0.05 | q<0.10 |
|---|---|---|---|---|
| D7 | unadj | 13,893 | 2,563 | 4,349 |
| D7 | adj | 13,893 | **4,600** | **5,867** |

The sensitivity is real and large, so the *observation* is fine; only the
explanation is wrong. The actual mechanism is visible in the data: at D7 the
muscle score is nearly orthogonal to group (I computed
corr(group, muscle_score) = **+0.174**, group-alone R² = 0.030, design condition
number 3.0), whereas at D1 and D14 it is substantially confounded
(corr = **−0.449** and **−0.432**, R² = 0.202 and 0.187). An almost-orthogonal
covariate removes a variance component that was inflating residual variance
without removing signal, which inflates t. That is the opposite of the stated
reason.

**【Why it matters】** This is the one place where the manuscript offers a
mechanistic explanation for its own most surprising number (2,563 → 4,600 genes,
an 80% increase on covariate adjustment). A reviewer who checks the design matrix
will find the stated reason impossible, and will then suspect the increase is an
artefact rather than a gain in sensitivity. The true explanation is more
interesting and more defensible — the covariate is *least* confounded exactly
where it helps most, which is a real and reportable observation.

**【Specific fix】**

```
Paste-ready replacement for the third sentence of §3.2:

  The adjusted count exceeds the unadjusted one at day 7 because the muscle
  covariate is close to orthogonal to group there: corr(group, muscle score) =
  +0.17 at day 7 (group alone explains R2 = 0.03 of the covariate) against −0.45
  at day 1 and −0.43 at day 14. An almost-orthogonal covariate removes a
  variance component without removing signal and therefore raises the moderated
  t, whereas at days 1 and 14 it removes some of the group contrast as well.
  This is a property of this particular contamination pattern, not a general
  argument for covariate adjustment.
```

---

## F5 — The per-compartment single-cell analysis dropped the myeloid compartment, and the paper does not say so

**【Problem】** §3.8 and Table 5 conclude that the interferon/myeloid programme is
"absent at the compartment level", having tested "8 compartments". Myeloid cells
are the compartment that carries an interferon programme by definition, and
myeloid was **not tested** — it was dropped by a cell-count filter. The paper
never discloses this, and §2.8's list of compartments does not include myeloid
among the eight, so a reader cannot tell that the most relevant compartment is
missing.

**【Evidence】** §2.8 (line 97) enumerates the compartments as "neurons, and
within them neurofilament-positive (NF), non-peptidergic (NP), peptidergic (PE)
and TH subtypes, plus satellite glia, Schwann cells, endothelium, fibroblasts and
erythrocytes" — that is **ten** named candidates, and the text then calls it
"eight compartments". Meanwhile `sc_GSE292522_celltype_DE.csv` contains exactly:

`Endothelial, Erythrocyte, Fibroblast, Neuron, Neuron.NF, Neuron.PE,
SatelliteGlia, Schwann`

`Neuron.NP` and `Neuron.TH` are absent, and **Myeloid** and **TCell** are absent
from the whole-DRG `FOCUS`/compartment set. `32_percelltype_sc.py:131-138`
requires `mask.sum() >= 5` per animal and `:150` requires `len(present) >= 6`.
I reconstructed per-animal cell counts from `sc_GSE292522_composition.csv` × the
per-sample cell counts in `sc_GSE292522_summary.json`:

| compartment | animals with ≥5 cells | verdict |
|---|---|---|
| Myeloid | **5 / 8** | **DROPPED** (needs ≥6) |
| TCell | 0 / 8 | dropped |
| Erythrocyte | 6 / 8 | kept |
| SatelliteGlia | 8 / 8 | kept |

Per-animal myeloid counts: Con 3, 10, 12, 5; B16F10 5, 3, 7, 1 — five animals
qualify, one short of the threshold. So the myeloid compartment fell out on a
technicality at exactly the boundary.

The manuscript does handle the related *fraction* result honestly — §3.8 and
Table 5 both report the myeloid fraction as an underpowered non-significant null
(d = 1.37 vs MDE 2.38, U = 2/14, p = 0.114), which I confirmed exactly. But the
compartment-level *expression* test was never run, and the paper's conclusion
sentence is unqualified:

> "No compartment — neurons, neurofilament-positive or peptidergic neuronal
> subtypes, satellite glia, Schwann cells, endothelium, fibroblasts or
> erythrocytes — showed interferon-stimulated-gene up-regulation"
> (§3.8, line 153; mirrored in Table 5)

Note that list has exactly eight items and, again, no myeloid.

**【Why it matters】** §3.8 is the paper's defence of its own §3.3 finding: "the
interferon / myeloid programme is a feature of the Walker-256 tibia model at its
day-7 peak, not a universal signature of cancer-pain DRG." To support "absent"
you need to have looked in the place where it would be. An interferon signature
in B16F10 DRG would live in myeloid cells; the myeloid compartment is exactly
where this dataset has too few cells to test. The honest conclusion is *"absent
from every compartment that had enough cells to test, and the myeloid
compartment could not be tested at all"* — which is a weaker but defensible
claim. As written, a reviewer who notices will read the current sentence as
overreach, and this is the finding most likely to draw a rejection-level
comment, because it touches the paper's novel contribution.

**【Specific fix】** Text edits are mandatory; one extra analysis is strongly
recommended.

*Mandatory (text):*

```
Paste-ready replacement for the sentence beginning "No compartment" in §3.8:

  In the six compartments with sufficient cells per animal to test — neurons,
  neurofilament-positive and peptidergic neuronal subtypes, satellite glia,
  Schwann cells, endothelium, fibroblasts and erythrocytes, subject to the
  minimum-count rule in Section 2.8 — no compartment showed
  interferon-stimulated-gene up-regulation; the only significant ISG changes
  were down-regulation of Ifit1 in satellite glia (q = 0.0037), Isg15 in
  satellite glia (q = 0.007) and Ifit2 in neurons (q = 0.046). The myeloid
  compartment could not be tested at this depth: it contained a median of 5
  cells per animal (range 1–12) and failed the five-cell minimum in three of
  eight animals, so it fell below the six-animal threshold. The myeloid
  result is therefore limited to the whole-DRG pseudobulk, where Cd68 and Aif1
  are unchanged, and to the cell-fraction comparison, which is underpowered
  rather than null.
```

And amend §2.8's compartment sentence:

```
Paste-ready replacement for the first sentence of the per-compartment
paragraph in §2.8:

  Each cell was assigned to one of ten candidate compartments — neurons, and
  within them neurofilament-positive (NF), non-peptidergic (NP), peptidergic
  (PE) and TH subtypes, plus satellite glia, Schwann cells, endothelium,
  fibroblasts, myeloid cells and erythrocytes — by the marker set with the
  highest mean log-expression [16]. Compartments present in fewer than six
  animals, or with fewer than five cells in an animal, were dropped; this
  removed the myeloid and T-cell compartments and the NP and TH neuronal
  subtypes, leaving eight compartments for the per-compartment test. The
  dropped compartments and their per-animal cell counts are reported in
  Supplementary Table S3 so that the scope of the compartment-level null is
  explicit.
```

*Recommended (new analysis — pooled-myeloid):* because the per-animal myeloid
counts are too low, pool across animals within group and test at the
compartment level with a design that respects the actual replication:

```
Spec for a new analysis (scripts/34_myeloid_pseudobulk_sc.py):

  Inputs : data/processed/GSE292522/*_matrix.mtx.gz (as 32_percelltype_sc.py)
  Step 1 : annotate cells exactly as 32_percelltype_sc.py:118-129 (BROAD then
           NEURON_SUB); do not change the marker sets.
  Step 2 : form the myeloid compartment per animal WITHOUT the >=5-cell
           filter; record cells per animal (expected 1-12).
  Step 3 : if any animal has 0 myeloid cells, report the compartment as
           untestable rather than dropping it silently.
  Step 4 : two analyses, both reported:
           (a) pseudobulk over the myeloid cells of all four animals per
               group, DE with the same _limma.fit_limma(trend=True) used
               elsewhere, n = 1 vs 1 -- reported as descriptive only, with the
               explicit statement that n = 1 per group permits no variance
               estimate and therefore no test;
           (b) a per-cell ISG score (the ISG list at
               32_percelltype_sc.py:73-76) compared between groups across all
               myeloid cells with a Mann-Whitney U test, acknowledging that
               this treats cells as independent and is therefore descriptive.
  Output : results/tables/sc_GSE292522_myeloid_compartment.json with keys
           cells_per_animal, n_animals_ge5, dropped_by_rule (bool),
           myeloid_ISG_score_Con, myeloid_ISG_score_B16F10, MWU_p,
           descriptive_only (bool = True),
           plus results/tables/sc_GSE292522_myeloid_compartment.md.
  Manuscript : add one sentence to §3.8 reporting the outcome either way.
  Rationale: this converts an unexamined gap into a documented, bounded
  limitation, which is what the conclusion already needs.
```

---

## F6 — Table 5's "8 compartments" and §3.8's claim omit that four candidate compartments were dropped

**【Problem】** Related to F5 but separable: the paper reports a positive
result — "no ISG up-regulation in any of 8 compartments" — without disclosing
that 4 of the 10 candidate compartments failed the inclusion rule. The number 8
is correct; the framing that it covers "every DRG compartment tested" is not,
because the compartments that failed are the informative ones.

**【Evidence】** `32_percelltype_sc.py:109` builds
`COMPARTMENTS = list(BROAD) + [f"Neuron.{s}" for s in NEURON_SUB]` = 8 broad +
4 subtypes = **12** candidates. Four are dropped by the rules at `:133/:137`
(`>=5` cells) and `:150` (`>=6` animals): Myeloid (5/8 animals), TCell (0/8),
Neuron.NP, Neuron.TH. Eight survive, matching the artefact exactly.

The count is also internally inconsistent in the text: §2.8 names ten
compartments and then says "eight"; §3.8's conclusion sentence lists eight and
calls them the set tested; Table 5 says "8 compartments".

**【Why it matters】** A "no effect anywhere" claim is only as strong as its
coverage, and here the coverage statement is ambiguous in a way that
under-counts the gaps. Editors increasingly ask exactly this question of
single-cell papers that make a null claim.

**【Specific fix】** Covered by the F5 §2.8 replacement. Additionally add to
Table 5:

```
Paste-ready Table 5 row replacing "Per-compartment ISG (8 compartments)":

| Per-compartment ISG | 8 of 12 candidate compartments tested (myeloid, T cell,
  NP and TH subtypes dropped for insufficient cells per animal; see §2.8);
  none up in any; Ifit1↓ satellite glia (q = 0.0037), Isg15↓ satellite glia
  (q = 0.007), Ifit2↓ neurons (q = 0.046) |
```

---

## F7 — §3.8 names Oas1c as a validated ISG, but Oas1c is not in the study's ISG panel and was never tested per compartment

**【Problem】** §3.8 and Table 5 report "Ifit2 and Oas1c were significantly
down-regulated (q = 0.042 and 0.014)" as the two ISG findings. `Oas1c` does not
appear in the ISG gene list used anywhere in the code, and it is absent from the
per-compartment focus set — so it was tested only in the whole-DRG pseudobulk,
where it was picked up because the test is genome-wide, not because it is a
designated ISG.

**【Evidence】** `grep -rn "Oas1c" scripts/` returns **no hits**. The ISG panel is
`32_percelltype_sc.py:73-76` (38 genes: Ifit1/2/3, Isg15, Isg20, Mx1/2,
Oas1a/1b/2/3, Cxcl10, Ccl5, Stat1/2, Irf7/9, …) — no Oas1c. The per-compartment
`FOCUS` list (`:80-82`) contains `Oas1a` and `Oas1b` but not `Oas1c`. I confirmed
`sc_GSE292522_celltype_DE.csv` has no `Oas1c` row in any compartment.

The value itself is real — I recomputed it from
`sc_GSE292522_pseudobulk_DE.csv`: `Oas1c logFC = −0.0291, t = −4.963,
q = 0.0136`. So this is a provenance/labeling defect, not a numerical one.

A related selection issue in the same sentence: of the 38 designated ISGs, the
only one down at q<0.05 in the whole-DRG pseudobulk is `Ifit2`. The manuscript
selects `Ifit2` **and** an un-designated `Oas1c` for reporting.

**【Why it matters】** Small in magnitude, but it is a claim of panel membership
that is false, and it is in the sentence that carries the single-cell
interferon conclusion. It also means the reported ISG set and the tested ISG set
differ, which a reviewer checking the Methods' gene list will notice.

**【Specific fix】**

```
Paste-ready replacement for the ISG sentence in §3.8:

  Interferon-stimulated genes were not uniformly flat but were not activated:
  of the 38 genes in the pre-specified ISG panel (Section 2.8), the only one
  significantly down-regulated in the whole-DRG pseudobulk was Ifit2
  (log2FC −0.048, q = 0.042), and no ISG in the panel was up-regulated at
  q < 0.05. Oas1c was also down-regulated (log2FC −0.029, q = 0.014) but is not
  a member of the pre-specified panel and was recovered from the
  genome-wide pseudobulk test rather than from the panel, so we report it as
  exploratory only.
```

And correct Table 5's row:

```
| Interferon-stimulated genes | none up in any of the 38 panel genes; Ifit2 down
  (q = 0.042). Oas1c down (q = 0.014) but not a panel member — exploratory |
```

---

## F8 — §3.7's lowess-span sensitivity is computed on a different gene universe from §3.2, and reports 320 genes where §3.2 reports 235

**【Problem】** §3.7 says "Differential expression barely moved with the lowess span
(316–325 genes at q < 0.05 …)". For the same contrast (bone D14 vs Sham,
muscle-adjusted, frac = 0.4) §3.2 reports **235** genes. The sensitivity sweep's
reference value is **320**. A reader comparing the two sections sees a 36%
discrepancy in the same quantity with no explanation.

**【Evidence】** `sensitivity_limma_trend.csv` (frac = 0.4, the reference row) gives
`n_q05 = 320`. `DE_D14_vs_Sham_adj.csv` gives `(q<0.05).sum() = 235` over 13,893
genes. Root cause, confirmed in code: `23_sensitivity.py:75-80` subsets both
studies to `shared = set(Lb.index) & set(Lc.index)` — **13,601** genes — and then
fits D14 on that reduced universe (`23_sensitivity.py:164-168`), whereas
`20_multistage_pipeline.py:228-229` filters GSE149648 to its own 13,893 genes and
never intersects with GSE160543. BH denominators therefore differ (13,601 vs
13,893) as do the lowess fits, so 320 vs 235 is expected.

This is a defensible design for a sensitivity analysis — arguably more so, since
the shared universe is what the cross-aetiology work uses — but the manuscript
presents the sweep as if it perturbed only the lowess span.

**【Why it matters】** §3.7's purpose is to certify §3.2, §3.3 and §3.5 as
"not artefacts of tuning". A sweep whose absolute values disagree with the
reported results by 36% cannot certify them, and the discrepancy is unexplained.
A reader who spots it will reasonably wonder what else in the sensitivity
analysis is off-universe.

**【Specific fix】**

```
Paste-ready replacement for the lowess-span sentence in §3.7:

  Differential expression barely moved with the lowess span: 316–325 genes at
  q < 0.05 across spans 0.2–0.6 and with the trend switched off, with Spearman
  rho of the moderated t against the reference at least 0.9999. This sweep is
  run on the 13,601 genes shared with GSE160543 so that it matches the gene
  universe used for the cross-aetiology analysis, which is why its absolute
  counts are not directly comparable with the 13,893-gene discovery counts in
  Section 3.2 (235 genes at day 14 at the reference span); the reference
  configuration gives 320 genes on the shared universe against 235 on the
  discovery universe, and the ranking correlation between the two is reported
  in Supplementary Table S4.
```

And add to `23_sensitivity.py` a second output on the discovery universe so the
two are bridged:

```python
# append to 23_sensitivity.py after de_grid is built (line ~186)
# Refit the reference configuration on the FULL discovery gene universe
# (13,893 genes, as in 20_multistage_pipeline.py:228) so that the sweep can be
# compared with Section 3.2 rather than only with itself.
L_full, _ = m22.log2fpkm(bone_fpkm), None
keep_full = ((bone_fpkm >= 1.0).sum(axis=1) >= 4) & (bone_fpkm.mean(axis=1) >= 0.5)
L_full = L_full.loc[keep_full]
Xc = L_full[ca14 + sh]
gvec = np.array([1.0] * len(ca14) + [0.0] * len(sh))
msc = muscle_score(Xc).values.astype(float); msc -= msc.mean()
r_full = fit_limma(Xc, np.column_stack([np.ones(len(ca14) + len(sh)), gvec, msc]),
                   coef_idx=1, trend=True, frac=REF["frac"])
r_full["q"] = bh(r_full["P"].values)
summ["n_genes_discovery_universe"] = int(len(r_full))
summ["n_q05_discovery_universe"] = int((r_full["q"] < 0.05).sum())
# expected: n_genes_discovery_universe = 13893, n_q05_discovery_universe = 235
# If n_q05_discovery_universe != 235, Section 3.2 and this script have drifted
# and the manuscript's robustness claim must be withdrawn pending reconciliation.
```

---

## F9 — §3.8's D1 and D14 concordance values (ρ = 0.112, 0.030) have no computable source in the deposited artefacts

**【Problem】** §3.8 and Table 5 report concordance with the Walker-256 programme
as "Spearman ρ = 0.098 (D7), 0.112 (D1), 0.030 (D14)". Only the D7 value exists
in any artefact. The D1 and D14 values are not in any JSON, not in any CSV, and
are not produced by any script.

**【Evidence】** `sc_GSE292522_summary.json` contains exactly one concordance
block:

```
"convergent_validation_vs_bulk_D7": {
  "n_shared_genes": 12158,
  "spearman_rho": 0.09845456164172446,   <- matches the manuscript's 0.098
  "spearman_p": 1.4160029287915969e-27,
  "shared_with_bulk_D7": 12158
}
```

`31_process_sc.py:145` hardcodes `bulk_path = ... "DE_D7_vs_Sham_adj.csv"` and
computes one Spearman. `grep -rn "D1_vs_Sham\|D14_vs_Sham" scripts/` returns **no
hits** — no script ever compares the pseudobulk against the day-1 or day-14 bulk
contrasts.

**【Why it matters】** Three numbers are reported as a trend ("highest at day 1 …
vanishing by day 14") that supports the paper's model-specificity argument, and
two of the three cannot be regenerated by anything in the repository. A reader
who tries to reproduce them will fail, and the deposited artefacts are the
paper's stated reproducibility mechanism. It is also the shape of claim most
likely to be an accidental leftover from an earlier run.

**【Specific fix】**

*Option A (regenerate — recommended):* extend `31_process_sc.py` and re-run.

```python
# in 31_process_sc.py, replace the single-contrast block at lines 144-161 with:
    # ---- convergent validation vs each Walker-256 bulk stage ----
    conv = {}
    for stage in ["D1", "D7", "D14"]:
        bulk_path = os.path.join(TAB, f"DE_{stage}_vs_Sham_adj.csv")
        if not os.path.exists(bulk_path):
            continue
        bk = pd.read_csv(bulk_path)
        gcol = "gene" if "gene" in bk.columns else bk.columns[0]
        bk = bk.set_index(gcol)
        shared = [g for g in L.index if g in bk.index]
        lfc_sc = L.loc[shared, [c for c in PB.columns if SAMPLES[c][0] == "B16F10"]].mean(axis=1) \
                 - L.loc[shared, [c for c in PB.columns if SAMPLES[c][0] == "Con"]].mean(axis=1)
        r, p = spearmanr(lfc_sc.values, bk.loc[shared, "logFC"].values)
        conv[f"convergent_validation_vs_bulk_{stage}"] = {
            "n_shared_genes": len(shared), "spearman_rho": float(r),
            "spearman_p": float(p)}
    # and write `conv` (not the single block) into the summary JSON at line 193
```

*Option B (delete the unsupported numbers):*

```
Paste-ready replacement for the concordance sentence in §3.8:

  Concordance with the Walker-256 day-7 programme was weak (Spearman rho =
  0.098 between pseudobulk and bulk log2 fold changes over 12,158 shared
  genes).
```

---

## F10 — Table 2's Aif1 Δ does not equal the difference of its own printed cells (double rounding at the display boundary)

**【Problem】** Table 2's Aif1 row prints per-aetiology values +2.19 and −0.88 and
Δ = +3.08. A reader subtracting the printed cells gets +3.07. Every other row
closes exactly.

**【Evidence】** I recomputed from `etiology_marker_cross.csv` at full precision:
`logFC_A = 2.192235`, `logFC_B = −0.883477`, `delta_logFC = 3.075712`, which
correctly rounds to **+3.08**. The displayed cells round to +2.19 and −0.88, and
2.19 − (−0.88) = 3.07. I checked all seven significant rows for closure:

| Gene | printed D7 | printed comparator | Δ printed | recomputed from printed cells | closes? |
|---|---|---|---|---|---|
| Cd68 | +2.84 | −0.58 | +3.42 | +3.42 | yes |
| **Aif1** | **+2.19** | **−0.88** | **+3.08** | **+3.07** | **NO** |
| Gal | −0.67 | +1.92 | −2.59 | −2.59 | yes |
| Ccl2 | +0.90 | −0.49 | +1.39 | +1.39 | yes |
| Il1b | +0.80 | −0.43 | +1.23 | +1.23 | yes |
| Calca | −0.27 | +0.50 | −0.77 | −0.77 | yes |
| Scn10a | −0.46 | +0.24 | −0.70 | −0.70 | yes |

This is a display-precision artefact, not a computational error — the Δ column is
correct. It is the only cell in any of the five tables where a reader cannot
reproduce the printed Δ from the printed inputs.

**【Why it matters】** Minor numerically, but Table 2 is the paper's headline
table and Aif1 is one of the two genes named in the Abstract. A table that does
not close on its own printed values invites doubt about the whole table.

**【Specific fix】** Carry one more digit in the per-aetiology columns:

```
Paste-ready: change the per-aetiology number format in Table 2 from 2 dp to 3 dp
and regenerate (regenerate_ms_tables.py, r2 -> a new r3 already defined at
line 12 but currently used only for Table 2's Delta; apply it to the four
logFC columns as well).

  | Aif1 | +2.192 | +1.954 | −0.883 | −0.459 | +3.08 (D7 vs PTX) | 0.012 |

With 3 dp every row closes: 2.192 - (-0.883) = 3.075 -> +3.08.
```

---

## F11 — §3.4's Table 4 omits four modules that pass the same significance threshold it uses, including the largest-|z| target module

**【Problem】** Table 4 presents 8 modules as "Calibration of the in-silico
perturbation layer", but 12 modules pass q < 0.05 at day 7 knockout. The
selection rule is not stated, and the omissions are not neutral: they include
the **highest-|z| target module after the controls** and the module that is
third-largest by reversal.

**【Evidence】** From `virtual_perturbation_scores.csv`, D7 knockout, all 12
modules with `q_empirical < 0.05`:

| module | kind | reversal | z | q | in Table 4? |
|---|---|---|---|---|---|
| CTRL interferon | control | +0.2282 | +4.611 | 0.011 | yes |
| CTRL pan-neuronal | control | −0.1086 | −3.560 | 0.011 | yes |
| Glial / immune | target | +0.1545 | +3.271 | 0.011 | yes |
| K⁺/HCN channels | target | −0.0916 | −2.867 | 0.011 | yes |
| CTRL myelin / Schwann | control | −0.0896 | −2.707 | 0.020 | yes |
| Voltage-gated Na⁺ | target | −0.0903 | −2.626 | 0.020 | yes |
| CTRL erythrocyte | control | −0.0744 | −2.568 | 0.020 | yes |
| Neuropeptides | target | −0.0675 | −2.449 | 0.020 | yes |
| Neurotrophin | target | −0.0676 | −2.387 | 0.020 | **NO** |
| Opioid receptors | target | −0.0675 | −2.315 | 0.022 | **NO** |
| **Cytokine / chemokine** | **target** | **+0.0803** | **+2.833** | **0.025** | **NO** |
| TRP channels | target | −0.0588 | −2.122 | 0.038 | **NO** |

`regenerate_ms_tables.py:86-87` hardcodes an 8-element `order` list, so the
regeneration script reproduces the same omission — the selection is baked in.

The omission of **Cytokine / chemokine** is the substantive one. §3.3 calls
cytokine/chemokine the third-strongest day-7 module by mean moderated *t*
(+5.88, q = 0.0129, Table 1), and in the cross-aetiology comparison it is the
module with the largest *significant* Δ against oxaliplatin (§3.5, Δ = +1.75;
Table 3, q = 0.0073). It is the module that most directly links the discovery
peak to the cross-aetiology divergence, and it is the one positive-direction
target module with a genome-wide-significant perturbation z that is missing from
the calibration table.

**【Why it matters】** §3.4's argument is that the perturbation layer's ranking is
non-specific because controls outrank targets. Table 4 as shown does support
that (3 of 8 rows are controls, and CTRL pan-neuronal outranks 6 target
modules). But the reader cannot verify the selection, and a reviewer who
recomputes will find 12 significant modules against 8 shown. Because the missing
rows are all targets, and because the one omitted *positive* target is the
cytokine module that §3.5 relies on, the omission pattern is unhelpful to the
paper's own argument even though it does not overturn it.

**【Specific fix】** Either show all 12, or state the rule.

```
Paste-ready replacement for the Table 4 caption:

  Table 4. Calibration of the in-silico perturbation layer (day 7, discovery).
  All 12 modules reaching q < 0.05 against the degree-matched null at day 7
  knockout are shown; no module was excluded. Rows are ordered by absolute
  reversal score. The remaining 10 of the 22 modules tested did not reach
  q < 0.05 (Prostanoid / COX z = +0.27, Purinergic P2X z = −0.73,
  Toll-like receptors z = +0.92, Injury / regeneration z = −1.06,
  Monoaminergic z = −1.07, Cannabinoid z = −1.13, CTRL housekeeping z = −1.31,
  Ca2+ channel alpha2delta z = −1.54, Glutamate / NMDA z = −1.97,
  GABAergic z = −1.99); all 22 are in
  results/tables/virtual_perturbation_scores.csv.
```

and change `regenerate_ms_tables.py:86-87` from the hardcoded list to:

```python
# regenerate_ms_tables.py — replace the hardcoded `order` list so the table
# cannot silently omit a qualifying module again.
_sig = v7[v7.q_empirical < 0.05].copy()
_sig["abs_rev"] = _sig.reversal_network.abs()
order = _sig.sort_values("abs_rev", ascending=False).module.tolist()
assert len(order) == len(v7[v7.q_empirical < 0.05]), "order must cover all q<0.05 rows"
```

---

## F12 — §3.6's "64 of the 264" uses a denominator that excludes the pooled-network arm the Methods promise

**【Problem】** §3.6 evaluates "264 contrast-level tests" and reports "Sixty-four
of the 264 (or 32/132 deduplicated)". The deposited artefact has **528** rows:
264 per network scope × 2 scopes. The Methods (§2.7, final sentence) promise "A
pooled-network variant is reported as a sensitivity analysis", but no pooled
result appears anywhere in the manuscript.

**【Evidence】** `etiology_bootstrap_differences.csv` has 528 rows:
`per_contrast` 132 KO + 132 OE, `pooled_by_study` 132 KO + 132 OE.
`etiology_summary.json` confirms `n_rows_raw: 528`,
`n_independent_tests_per_scope: 132`, `n_independent_tests_total: 264`.
Randset-significant counts I recomputed by scope:

| scope | perturbation | q_randset < 0.05 |
|---|---|---|
| per_contrast | knockout | 32 |
| per_contrast | overexpression | 32 |
| pooled_by_study | knockout | **33** |
| pooled_by_study | overexpression | **33** |

So the pooled arm yields 33, not 32 — the deduplication is *not* exact across
scopes, which the code comment at `22_etiology_comparison.py:620-624` anticipates
("different degree-matched draws"). The manuscript's "32/132" is therefore
correct for `per_contrast` and silently unreported for `pooled_by_study`. I also
checked whether the pooled scope changes the headline verdict: it does not —
`min q_permutation` pooled = 0.1167, and **0 of 264** pooled rows reach q < 0.05,
so the negative conclusion is if anything slightly stronger under pooling.

**【Why it matters】** The manuscript reports one arm of a two-arm analysis and
describes the other as reported. Both the numerator and the denominator of the
headline "64 of 264" are scope-restricted without saying so. A reader
recomputing from the deposited CSV gets 130 of 528 and cannot reproduce the
sentence. This is a straightforward reporting-hygiene defect rather than a
numerical error — the conclusion holds either way — but it is the kind of
discrepancy that makes a reproducibility check fail on a technicality.

**【Specific fix】**

```
Paste-ready replacement for the first two sentences of §3.6:

  Run with reversal scores, the same comparison falls apart. We evaluated 22
  modules x 6 aetiology pairs x 2 perturbation directions = 264 contrast-level
  tests per network scope, under two network scopes (528 rows in total). All
  numbers in this section refer to the per-contrast network unless stated;
  the pooled-by-study scope gives 33 rather than 32 randset-significant tests
  per perturbation direction and reaches the same conclusion, with 0 of 264
  tests reaching q < 0.05 on the exact permutation p (minimum q = 0.117)
  versus 0 of 132 for the per-contrast scope (minimum q = 0.165).
```

---

## F13 — §3.8 calls Gfap "nominally up" when q = 0.006, and the effect size is +0.004 log₂

**【Problem】** §3.8 reports "the satellite-glia marker Gfap was nominally up in
B16F10 (q = 0.006)". Two problems: q = 0.006 is *below* 0.05, so it is
significant, not nominal; and the effect size is +0.0036 log₂ — a 0.25% change
that is biologically meaningless. Table 5 omits Gfap entirely.

**【Evidence】** From `sc_GSE292522_pseudobulk_DE.csv`:
`Gfap logFC = +0.0036, t = +5.835, P = 9.22e-05, q = 0.0061`. The t is large
because the pseudobulk is precise (34,603 genes, n = 4+4), not because the
effect is large. The same pattern holds per compartment: `Gfap` in Neuron is
`+0.001, q = 0.0749`; in Neuron.NF `+0.001, q = 0.0473`; in SatelliteGlia
`+0.000, q = 0.9997`.

**【Why it matters】** "Nominally up" applied to a q = 0.006 result will be read as
the authors hedging a borderline finding when in fact they are reporting a
statistically strong but biologically null result. In a paper whose whole
argument is about distinguishing real compartment shifts from artefactual ones,
reporting a +0.4% Gfap change without its magnitude is a self-inflicted
inconsistency. §3.8 itself uses the opposite convention two sentences later,
arguing that the design can detect log₂FC of −0.048 — so the paper has already
established that effects of this size are near its detection floor.

**【Specific fix】**

```
Paste-ready replacement for the Gfap clause in §3.8:

  the satellite-glia marker Gfap was not detectably changed in the
  whole-DRG pseudobulk (log2FC +0.004, q = 0.006), a change of 0.25% that is
  statistically resolvable at this depth but biologically negligible, and it
  was absent in every compartment in which it was tested (satellite glia
  log2FC 0.000, q = 0.9997).
```

---

## F14 — §2.2's "6–9% global offset" is not what the deposited QC table shows

**【Problem】** §2.2 justifies the median-of-expressed-genes scaling by "a 6–9%
global offset in GSE149648". The deposited per-sample FPKM totals give 1.9%–5.1%
relative to sham.

**【Evidence】** From `qc_muscle_scores.csv` (`fpkm_total_M`), group means:
Sham 0.438560, D1 0.446768, D7 0.460878, D14 0.452521. Relative to sham:
D1 **+1.87%**, D7 **+5.09%**, D14 **+3.18%**. The maximum per-sample deviation
in any group is 0.473253 (D7) vs 0.435987 (Sham) = +8.6% sample-to-sample, so
"6–9%" appears to describe the *extremes of individual samples* rather than the
group offsets. The script comment at `20_multistage_pipeline.py:201-204` gives
"sham ~0.439e6 vs day 7 ~0.461e6" — which is exactly the +5.1% I computed — and
then asserts "a 6-9% global offset", so the same overstatement is in the code
comment.

The scaling itself is harmless and arguably good practice; only the
justification number is wrong.

**【Why it matters】** Small, but it is a Methods claim that does not reproduce,
and it is the stated justification for a normalisation step applied to every
number in the paper.

**【Specific fix】**

```
Paste-ready replacement for the parenthetical in §2.2:

  (per-sample FPKM totals are not equal across groups here: group means are
  0.439e6 (sham), 0.447e6 (day 1), 0.461e6 (day 7) and 0.453e6 (day 14), i.e.
  a 1.9–5.1% offset relative to sham, with individual samples spanning
  0.436–0.473e6)
```

---

## F15 — §2.4's "2–6 residual degrees of freedom" does not describe this pipeline

**【Problem】** §2.4 motivates moderation by saying "At n = 4 per group an ordinary
t-test has 2–6 residual degrees of freedom". No model in the pipeline has
df_res below 5.

**【Evidence】** Every design in the repository:
- discovery DE, adjusted (`20_multistage_pipeline.py:315`): `[1, group, muscle]`,
  n = 8 → df_res = **5**
- discovery DE, unadjusted (`:316`): `[1, group]`, n = 8 → df_res = **6**
- cross-aetiology marker DE (`22_etiology_comparison.py:380`): `[1, g, muscle]`,
  n = 8 → **5**
- single-cell pseudobulk (`31_process_sc.py`): `[1, group]`, n = 8 → **6**
- per-compartment (`32_percelltype_sc.py:157`): `[1, group]`, n = 8 → **6**

So the realised range is 5–6. The "2" appears to be inherited from an earlier
GSE249443 design (`12_adjust_covariate_limma.py`, a different dataset with
fewer samples) that is not part of this study's pipeline.

**【Why it matters】** Minor, and it slightly *understates* the paper's own case
(df = 5 is still bad). But it is a Methods claim that does not match any code in
the repository, and it points at a dead script (`12_adjust_covariate_limma.py`)
that a reader following the citation would find uses a different dataset, a
different muscle-gene list (21 genes vs 20), and no limma-trend at all.

**【Specific fix】**

```
Paste-ready replacement for the df sentence in §2.4:

  At n = 4 per group the residual degrees of freedom are 5 in the
  muscle-adjusted models and 6 in the unadjusted and single-cell models, which
  is poorly calibrated in both tails; moderation borrows strength across
  roughly 14,000 genes.
```

Separately, add a note to the code availability statement or a README:

```
  Note for readers: 10_bone_end_to_end.py, 12_adjust_covariate_limma.py and the
  signature_bone_*.csv / VALID_GSE249443_* artefacts refer to a different
  dataset (GSE249443) and a superseded pipeline that does not implement
  limma-trend. They are retained for provenance and are not part of the
  analyses reported in this paper.
```

---

## F16 — `perturbation_bone.csv` is an orphan artefact from a superseded pipeline on a different dataset

**【Problem】** `results/tables/perturbation_bone.csv` sits among the deposited
artefacts but is written by no script in the repository and corresponds to no
table in the manuscript.

**【Evidence】** `grep -l "perturbation_bone.csv" scripts/*.py` returns **no
hits**. Its natural producer, `10_bone_end_to_end.py`, reads
`GSE249443_processed_data.xlsx` (line 23) — a different dataset, mouse rather
than rat, and a validation cohort the manuscript never mentions. That script
also uses a **one-sided** permutation p (`np.sum(nulls >= r["R"])`, line 185) and
a **uniform** random-gene null (line 172, `null_genes` without degree matching)
— precisely the two practices §4 of the manuscript describes as errors it fixed.
The same applies to `12_adjust_covariate_limma.py`, which uses its own
21-gene muscle list (line 28) and a plain-variance moderated t with **no trend**
(`limma_with_covariate`, line 92), and which writes `signature_bone_*.csv`.

**【Why it matters】** A repository that ships artefacts from a superseded
pipeline built on an undisclosed dataset, containing the exact one-sided-p and
uniform-null practices the paper's second contribution is about, is a
reproducibility hazard and an editorial risk. It also means "all artefacts are
traceable to a script in `scripts/`" is currently false.

**【Specific fix】**

*Option A (preferred — quarantine):* move the superseded files into
`results/tables/superseded/` and add the README note from F15.

*Option B (delete):* remove `perturbation_bone.csv`, `signature_bone_*.csv` and
`VALID_GSE249443_*`, and move `10_bone_end_to_end.py` and
`12_adjust_covariate_limma.py` to `scripts/superseded/`.

Either way, amend the Data and code availability statement:

```
  The repository contains the scripts and artefacts for the analyses reported
  here. Files under results/tables/superseded/ and scripts/superseded/ relate
  to an earlier pilot on a different dataset (GSE249443) and are retained for
  provenance only; they are not part of the reported analysis.
```

---

# § Stands up (things I suspected and verified as correct)

These are the load-bearing results I tried to break and could not.

**S1. Table 1 is exact in every cell, all four columns.** I recomputed all 22
modules from `module_level_signal.csv` (D7, muscle-adjusted). Every `mean_t`,
`mean_logFC`, `n` and `q_comp` matches the manuscript to the printed precision —
e.g. CTRL interferon 24.923/2.614/8/6.8×10⁻⁵ vs printed +24.92/+2.61/8/6.8×10⁻⁵;
GABAergic −5.708/−0.373/3/0.0245 vs −5.71/−0.37/3/0.025; CTRL pan-neuronal
−5.201/−0.399/10/0.00031 vs −5.20/−0.40/10/0.00031. I also verified the negative
controls land inside the target range as claimed: targets span −5.71 to −3.51,
and pan-neuronal (−5.20), myelin (−4.10) and erythrocyte (−4.80) all fall inside
it. This is the paper's Table 1 claim and it is clean.

**S2. The KO/OE mirror-image claim is exactly true.** I recomputed across all
264 module-pair-scope combinations in `etiology_bootstrap_differences.csv`:
`max|rev_KO + rev_OE| = 0.000e+00` and `max|p_permutation_KO − p_OE| = 0.000e+00`,
with 0 combinations differing by more than 1e-8. The manuscript's "identical
exact-permutation p to eight decimal places" is not an approximation — it is
bitwise identical. And the deduplicated BH denominator of 132 is correctly
implemented at `22_etiology_comparison.py:625-631`, with `etiology_summary.json`
confirming `n_independent_tests_per_scope: 132` and
`min_q_permutation_per_contrast: 0.16544728…`, matching the manuscript's 0.165.

**S3. Every Table 2 Δ and q I checked traces exactly.** Cd68 +3.4183 (q 0.0125),
Aif1 +3.0757 (0.0125), Gal −2.5921 (0.0125), Ccl2 +1.3907 (0.0173), Il1b +1.2305
(0.0182), Calca −0.7698 (0.0272), Scn10a −0.7022 (0.0260) — all match the
manuscript's Δ and q to printed precision, and all from the same artefact
(`etiology_marker_cross.csv`) as the per-aetiology columns. Atf3 and Gfap are
correctly reported as `—` / `ns`: their largest bone-vs-CIPN |Δ| is −0.809
(q = 0.272) and +1.450 (q = 0.091), neither surviving q < 0.05. The "ns"
convention is applied correctly here.

**S4. Table 4's eight printed rows are exact on all six columns.** Recomputed
from `virtual_perturbation_scores.csv` (D7, knockout): CTRL interferon +0.2282 /
+4.6111 / 75.43% / p 0.001998 / q 0.010989 vs printed +0.228 / +4.61 / 75.4 /
0.0020 / 0.011; glial/immune +0.1545 / +3.2707 / 61.48 / 0.0020 / 0.011 vs
+0.155 / +3.27 / 61.5 / 0.0020 / 0.011; CTRL pan-neuronal −0.1086 / −3.5597 /
32.87 / 0.0020 / 0.011; CTRL erythrocyte −0.0744 / −2.5681 / 32.98 / 0.0060 /
0.020. All eight match. The oracle anchor is correctly implemented as a
per-module, same-size seed (`20_multistage_pipeline.py:444-448`, k = n_targets),
and `pct_of_oracle = 100·|rev|/|oracle|` is the natural definition.

**S5. §3.1's QC numbers are exact.** Sample-to-sample r range 0.913–0.994
(printed 0.913–0.994); PCA variance 46.85 / 17.45 / 7.23 / 4.76% (printed
46.8 / 17.4 / 7.2 / 4.8); muscle ANOVA F = 0.637, p = 0.6054 (printed F = 0.637,
p = 0.605); GSM4506604 at 6.2573% sarcomere share with the next highest at
0.1604% (printed 6.26% and 0.16%). The claim that it is the only sample above the
1% threshold is correct — I checked all 16 rows of `qc_muscle_scores.csv`.

**S6. §3.2's DE counts and the leave-one-out analysis are exact.** Recomputed
from the six `DE_*_vs_Sham_*.csv` files over 13,893 genes: adjusted q<0.05 =
0 / 4,600 / 235 and q<0.10 = 5 / 5,867 / 860 for D1/D7/D14; unadjusted q<0.05 =
1 / 2,563 / 958; median |logFC| = 0.0731 / 0.1576 / 0.1187 (printed 0.073, 0.158,
0.119). `de_muscle_drop_sensitivity.csv` gives D7 4,600 → 3,041, Jaccard 0.58396
(printed 3,041 and 0.584), Spearman 0.91231 (printed 0.912), n_case_left = 3.
Every number in §3.2 is correct; only the *explanation* is wrong (F4).

**S7. §3.6's negative-result statistics are all exact.** Recomputed from
`etiology_bootstrap_differences.csv` and `negative_result_quantification.json`:
32 randset-significant of 132 deduplicated tests, of which 14 are
negative-control modules and 18 are targets (printed "32/132", "14 of those 32",
"the remaining 18") — all correct, and the 14 controls are exactly the 5 CTRL
modules. The erythrocyte z values match precisely: |z| = 7.3878 for bone D14 vs
D7 (printed 7.39) and z = −5.7822 for PTX vs OXA (printed −5.78), and it *is* the
single largest |z| in the per-contrast scope (max = 7.3878). Pan-neuronal
z = −3.7362 (printed −3.74). Glial/immune min-q row: Δ = 0.27429, CI
0.20307–0.28662, p = 0.0014283, q = 0.16545 (printed +0.274, 0.203–0.287, 0.0014,
0.165) — all four match. Target-vs-control: mean |z| 2.5889 vs 1.4617 (printed
2.59 vs 1.46), median percentile 55.7 (printed 55.7), MWU p = 0.4973 signed and
1.0000 on |z| (printed 0.50 and 1.00). All 64 randset-significant rows do carry
bootstrap CIs excluding zero (I checked: 64/64), and 28/64 = 43.75% are controls
(printed 28 of 64, 44%).

**S8. §3.8's single-cell numbers are almost entirely exact.** 55,967 cells total;
per-sample 6,514 / 7,322 / 6,019 / 5,817 / 6,175 / 7,480 / 9,826 / 6,814;
pseudobulk 1,517 up and 703 down at q < 0.05 over 34,603 genes; Sprr1a
+1.5539 (t = 54.493) vs printed +1.55, t = 54.5; Atf3 +0.8355 vs +0.84;
Ifit2 −0.0476 (q = 0.0417) vs printed −0.048 and q = 0.042; Cd68 −0.0005
(q = 0.783) and Aif1 +0.0009 (q = 0.920) — genuinely unchanged; myeloid
fraction 0.1170% vs 0.0517% (printed 0.12% and 0.05%), ratio 0.4423 → 55.77% drop
(printed 55%); Cohen's d 1.3685 (printed 1.37); MDE d 2.38 (printed 2.38);
Mann–Whitney U = 2/14 exactly, two-sided p = 0.11429 (printed 0.114), one-sided
p = 0.05714 (printed 0.057); D7 concordance ρ = 0.09845 over 12,158 shared genes
(printed 0.098). The claim "no ISG was up-regulated" is verified: of the 38
panel ISGs, 0 are up at q < 0.05. Per-compartment: galanin +1.472 in neurons
(q = 0.0001) and +2.537 in NF neurons (q < 0.0001), Scn10a −0.247 in neurons
(q = 0.0117), Ifit2 −0.049 in neurons (q = 0.0459), Isg15 −0.630 in satellite
glia (q = 0.0071) — all match the manuscript.

**S9. The Methods' account of the two documented bugs is accurate.** I read the
implementations. `build_network` sets the diagonal to **0.0**, not −1
(`_limma.py:239`), and the comment at `:212-217` states the |−1| = 1 reasoning
correctly; it also re-zeroes the diagonal after symmetrisation (`:245`) as a
defensive measure. The `cos_batch` 2-D path does centre and normalise
**per column** (`_limma.py:334-336`) and the docstring's account of the
grand-mean/Frobenius bug (`:313-320`) matches what the code no longer does. Both
bugs described in §4 are real, correctly diagnosed, and correctly fixed.

**S10. The convergence claim checks out.** §2.6 says the closed form matches 80
iterations to 2.2×10⁻¹⁶. I recomputed `max|P_iter − P_exact|` for a random
symmetric 300×300 W at α = 0.55 and got **2.776×10⁻¹⁶** — the same order, and
the specific figure depends on the matrix, so this is a fair statement. I also
confirmed |α·S| ≤ α < 1 is the correct convergence argument since S is a
normalised symmetric operator.

**S11. Force-inclusion of target genes is genuinely implemented.** §2.6's claim
that every module gene is placed on the graph is true: `20_multistage_pipeline.py:415-416`
builds `MUST` from all modules and passes it to `build_network(..., must_include=MUST)`,
and `22_etiology_comparison.py:578-585` does the equivalent for the
cross-aetiology network, logging the retained count per module. The
cross-aetiology node set is 3,000 HV genes plus all module genes.

**S12. The leave-one-out reversal ranges are deposited, as §5 claims.**
`virtual_perturbation_scores.csv` carries `reversal_loo_min`, `reversal_loo_max`
and `reversal_loo_sd`, computed over 8 refits per module
(`20_multistage_pipeline.py:491-515`). The interferon module's D7 knockout point
estimate +0.2282 has LOO range [+0.2185, +0.2297] — notably tight, so the paper's
headline reversal is not an artefact of one animal.

**S13. The environment claim in §2.9 is exact.** `requirements.txt` pins
numpy 2.4.6, scipy 1.18.1, pandas 2.3.3, statsmodels 0.15.0, matplotlib 3.11.1;
I confirmed the live interpreter reports precisely those versions on Python
3.13.14. And `_limma.py:36-43` genuinely raises `SystemExit` at import time if
statsmodels is missing, rather than degrading silently — the stated failure mode
is implemented as described.

**S14. The negative-control logic is implemented as advertised.** The five CTRL
modules in `_limma.py:188-197` match §2.5 exactly, and the 17 TARGETS in
`:158-176` match §2.5's list of 17 a priori modules. The erythrocyte set is
described as a technical control carrying no pain biology, and it is indeed the
module with the largest standardised effect — the paper's central calibration
point is real.

**S15. §3.4's baseline comparison is exact.** Spearman(network, additive) =
0.991145 over n = 132 conditions (printed 0.991 and n = 132); mean |network|
0.05680, mean |additive| 0.03284, ratio 1.7296 → the manuscript's "roughly 70%"
is 73.0%, which is a fair characterisation; mean |ridge| 0.71193 (printed 0.712);
mean % of oracle 20.5133 (printed 20.5). The ridge baseline is genuinely
leave-one-sample-out (`_limma.py:457-464`) and genuinely regression-based with
the target genes indexing the correct axis (`:442-449`, the bug the comment
describes).

---

# § Questions for the authors

I could not resolve these from the artefacts; I am not guessing.

**Q1. Which logFC source is authoritative for the headline marker claims?**
Table 2 and §3.5 use `etiology_marker_cross.csv` (unadjusted); the Abstract and
§3.3 use `etiology_marker_de.csv` (muscle-adjusted). Given that §2.3 argues at
length that the muscle covariate changes conclusions, and given that the day-7
sample carrying 6.26% sarcomere transcript is in the D7 group, was the
unadjusted source a deliberate choice to keep Table 2 internally consistent
(Δ and per-aetiology columns from one artefact), or an oversight? If deliberate,
I would rather see the covariate added and the table regenerated.

**Q2. What is the selection rule for Table 3's `—` cells?** Eleven of twelve
blanked cells are consistent with "the module's Δ in that column's contrast is
non-significant", but voltage-gated Na⁺ (Bone_D14, paired q = 0.0491) and CTRL
erythrocyte (Bone_D14, paired q = 0.0321) are blanked despite clearing q < 0.05,
while ten non-significant cells are shown. Is there a rule I have not inferred?

**Q3. Why is the CTRL erythrocyte Δ reported as D7-vs-OXA in a column headed
"D7 vs PTX", when the D7-vs-PTX value is −0.26 (q = 0.81)?** The §3.6 argument
compares the erythrocyte Δ to neuronal Δ values that are D7-vs-PTX. Was the
contrast switch deliberate, and if so how should the neuronal comparison be
read?

**Q4. Were Neuron.NP and Neuron.TH dropped for insufficient cells?** They are
absent from `sc_GSE292522_celltype_DE.csv` and §2.8 lists them among the
compartments. If they were dropped by the ≥6-animal / ≥5-cell rule, the
per-animal counts would be worth reporting, since NP neurons are the
non-peptidergic population most relevant to CIPN.

**Q5. What was the intent of the `Oas1c` result?** It is reported as an ISG in
§3.8 and Table 5 but is in no ISG list in the code and was not tested per
compartment. Was it added to the panel late, or picked out of the genome-wide
pseudobulk table?

**Q6. Can the D1 and D14 concordance values (ρ = 0.112, 0.030) be regenerated?**
They are in no artefact and no script computes them (F9). If they came from an
earlier run, they should be recomputed and deposited; if the D7 value is the only
defensible one, the sentence should be shortened to it.

**Q7. Is the `pooled_by_study` scope intended to be reported?** §2.7 says it "is
reported as a sensitivity analysis" but no pooled number appears in the paper.
I verified it gives 0/264 at q < 0.05 (min q = 0.117), so reporting it would
strengthen §3.6. Should it be added?

**Q8. Should the paper run the pooled-myeloid analysis in F5?** If not, the
§3.8 conclusion sentence needs the explicit qualification I supply in F5. I would
rather see the bounded claim than the unqualified one, but the compartment that
would carry an interferon programme in this dataset is the one that could not be
tested, and I would want that stated in the paper rather than left for a reader
to find.

**Q9. §2.2 cites a "6–9% global offset" but the QC table gives 1.9–5.1%.** Which
is intended — group-mean offsets or sample extremes?

**Q10. Are `perturbation_bone.csv` and `signature_bone_*.csv` intended to ship?**
They are orphan artefacts from the GSE249443 pilot, containing the one-sided-p
and uniform-null practices §4 describes as errors (F16). Shipping them next to
the real artefacts is confusing.

---

# § What I actually checked

## Files read in full
`docs/review/round2_2026-10-10/_PANEL_BRIEF.md`;
`docs/manuscript/manuscript_final.md` (all 293 lines);
`scripts/_limma.py`, `scripts/10_bone_end_to_end.py`,
`scripts/12_adjust_covariate_limma.py`, `scripts/20_multistage_pipeline.py`,
`scripts/22_etiology_comparison.py`, `scripts/23_sensitivity.py`,
`scripts/31_process_sc.py` (lines 88-200), `scripts/32_percelltype_sc.py`,
`scripts/33_negative_result_quantification.py`,
`scripts/regenerate_ms_tables.py`, `scripts/requirements.txt`.

## Commands run
- `scripts/98_check_manuscript_numbers.py` → 191 checks, 0 failed (run for
  coverage assessment only; **not** used as evidence for any number).
- `scripts/99_unit_tests.py` → all tests passed. I confirmed 14 test groups
  (T1-T14), matching §2.9's "14 groups" claim.
- `scripts/regenerate_ms_tables.py` → ran successfully; I diffed its output
  against the manuscript's Tables 2-4 cell by cell. This is how I found that
  Table 3's erythrocyte Δ is a post-generation hand edit (F3) and that Table 4's
  8-row selection is hardcoded (F11).
- ~30 independent Python recomputation scripts over the CSVs/JSONs.

## Artefacts read and independently recomputed
| Artefact | Rows | Used for |
|---|---|---|
| `module_level_signal.csv` | 66 (3 contrasts × 22) | Table 1, all 4 columns |
| `etiology_marker_de.csv` | 68 | §3.3 individual genes; source identification |
| `etiology_marker_cross.csv` | 102 | Table 2, all Δ/q; source identification |
| `etiology_module_logfc_cross.csv` | 132 | Table 3, all cells + blanking-rule test |
| `etiology_bootstrap_differences.csv` | 528 | §3.6 KO/OE mirror, 64/132, z values, CI check |
| `etiology_reversal_point_estimates.csv` | 88 | cross-check of reversal point estimates |
| `virtual_perturbation_scores.csv` | 132 | Table 4, oracle, LOO, §3.4 baselines |
| `head_to_head_baselines.json` | — | ρ 0.991, mean abs reversals, 20.5% |
| `run_summary.json` | — | §3.1/§3.2 counts, ANOVA, gene counts |
| `qc_muscle_scores.csv` | 16 | §3.1 correlation/ANOVA, 6.26% claim, covariate collinearity |
| `qc_sample_correlation.csv` | 16×16 | §3.1 r range 0.913–0.994 |
| `qc_pca_variance.csv` | 4 | §3.1 PC1-4 variance |
| `de_muscle_drop_sensitivity.csv` | 1 | §3.2 LOO: 3,041 / 0.584 / 0.912 |
| `DE_{D1,D7,D14}_vs_Sham_{adj,unadj}.csv` | 13,893 each | §3.2 all counts and median \|logFC\| |
| `sensitivity_summary.json` + `sensitivity_limma_trend.csv` + `sensitivity_grid.csv` | — | §3.7 all numbers; F8 gene-universe finding |
| `etiology_summary.json` | — | 528 rows, 132 independent tests, min q 0.165 |
| `negative_result_quantification.json` | T1-T5 | §3.6 all statistics |
| `sc_GSE292522_summary.json` | — | cell counts, 1,517/703, myeloid fraction, MWU, ρ(D7) |
| `sc_GSE292522_pseudobulk_DE.csv` | 34,603 | Table 5 markers, ISG up/down check |
| `sc_GSE292522_focus_markers.csv` | 50 | ISG panel values; Oas1c absence |
| `sc_GSE292522_composition.csv` | 8 | myeloid/TCell per-animal cell counts → F5 |
| `sc_GSE292522_celltype_DE.csv` | 197 | per-compartment Gal/Scn10a/Ifit2/Isg15/Ifit1 |
| `sc_GSE292522_celltype_summary.json` | — | compartment list (7 ISG-discrimination entries) |
| `results/tables/` full listing | 43 files | orphan-artefact audit |

## Values I recomputed myself, with verdicts
Legend: ✅ matches manuscript · ⚠️ mismatch or unsupported · ◐ correct but
under-disclosed

| Manuscript location | Claimed | My recomputation | Verdict |
|---|---|---|---|
| Table 1 | 22 rows × (n, mean t, mean logFC, q) | all 88 cells from `module_level_signal.csv` | ✅ |
| §3.3 | IFN +24.9, q 6.8e-5; glial +11.4; TLR +6.4; cyto +5.9; housek +0.65, q 0.52 | 24.923/6.8e-5; 11.412; 6.364; 5.877; 0.653/0.5237 | ✅ |
| §3.3 | target span −5.71 to −3.51, controls inside | GABAergic −5.708 … Ca2d −3.509; pan-neur −5.201, myelin −4.101, ery −4.801 | ✅ |
| §3.1 | r 0.913–0.994 | 0.913–0.994 | ✅ |
| §3.1 | PC 46.8/17.4/7.2/4.8 | 46.85/17.45/7.23/4.76 | ✅ |
| §3.1 | ANOVA F 0.637, p 0.605 | 0.637, 0.6054 | ✅ |
| §3.1 | GSM4506604 6.26%, next 0.16% | 6.2573%, 0.1604%; only sample >1% | ✅ |
| §3.2 | adj q<0.05: 0 / 4,600 / 235 | 0 / 4,600 / 235 | ✅ |
| §3.2 | adj q<0.10: 5 / 5,867 / 860 | 5 / 5,867 / 860 | ✅ |
| §3.2 | unadj: 1 / 2,563 / 958 | 1 / 2,563 / 958 | ✅ |
| §3.2 | median \|logFC\| 0.073/0.158/0.119 | 0.0731/0.1576/0.1187 | ✅ |
| §3.2 | LOO 3,041, Jaccard 0.584, ρ 0.912 | 3041, 0.58396, 0.91231 | ✅ |
| §3.2 | "covariate … raises the prior df" | covariate lowers df_res 6→5 | ⚠️ F4 |
| §3.3 | Cd68 +2.86, Aif1 +2.18, Gfap +1.33, Ccl2 +0.81, Il1b +0.80 | de.csv: 2.858/2.180/1.328/0.806/0.804 | ◐ source = de.csv, undisclosed |
| §3.2/Table 2 | same genes, same contrast | cross.csv: 2.839/2.192/1.498/0.897/0.800 | ⚠️ F2 (Gfap differs by 0.17) |
| Table 2 | Cd68 Δ +3.42, q 0.012 | 3.4183, 0.012487 | ✅ |
| Table 2 | Aif1 Δ +3.08, q 0.012 | 3.0757, 0.012487 | ✅ value; ⚠️ F10 cell doesn't close |
| Table 2 | Gal Δ −2.59, q 0.012 | −2.5921, 0.012487 | ✅ |
| Table 2 | Ccl2 Δ +1.39 (OXA), q 0.017 | 1.3907, 0.017343 | ✅ |
| Table 2 | Il1b Δ +1.23, q 0.018 | 1.2305, 0.018211 | ✅ |
| Table 2 | Calca Δ −0.77, q 0.027 | −0.7698, 0.027216 | ✅ |
| Table 2 | Scn10a Δ −0.70, q 0.026 | −0.7022, 0.026015 | ✅ |
| Table 2 | Atf3, Gfap "ns" | max \|Δ\| −0.809 q 0.272; +1.450 q 0.091 | ✅ |
| §3.5 | "Atf3 +1.13 PTX, +0.22 D7" | de.csv 1.126/0.222 (Table 2 says +1.12/+0.31 from cross) | ⚠️ F2 |
| Table 3 | 11 non-erythrocyte rows, all cells | all match to 2 dp incl. 0.0073 ×5 | ✅ |
| Table 3 | CTRL ery Δ −2.48, q 0.0073 in "D7 vs PTX" column | D7-vs-PTX = −0.2583, q 0.8143; −2.4831/0.0073 is D7-vs-**OXA** | ⚠️ **F3** |
| Table 3 | 9 cells blanked `—` | values exist for all; 2 blanks have paired q = 0.049, 0.032 | ⚠️ F3 |
| §3.6 | min q 0.165; glial/immune Δ +0.274, CI 0.203–0.287, p 0.0014 | 0.16545; 0.27429; 0.20307–0.28662; 0.0014283 | ✅ |
| §3.6 | KO/OE p identical to 8 dp, max\|Δ\| = 0 | max\|rev_KO+rev_OE\| = 0.0, max\|Δp\| = 0.0 over 264 combos | ✅ |
| §3.6 | 132 independent tests | 132/scope, 528 rows total, min q 0.16545 | ✅ |
| §3.6 | "64 of the 264 … 32/132" | 64 and 32 hold for per_contrast only; pooled = 66 | ⚠️ F12 |
| §3.6 | 14/32 controls, 18 targets | 14 controls (all 5 CTRL mods), 18 targets | ✅ |
| §3.6 | all 64 CIs exclude zero; 28 (44%) controls | 64/64; 28/64 = 43.75% | ✅ |
| §3.6 | ery \|z\| 7.39 (D14 vs D7), −5.78 (PTX vs OXA) | 7.3878; −5.7822; both are scope-max | ✅ |
| §3.6 | pan-neuronal z −3.74 | −3.7362 | ✅ |
| §3.6 | targets 55.7th pct; MWU p 0.50 / 1.00; mean \|z\| 2.59 vs 1.46 | 55.7; 0.4973; 1.0000; 2.5889 / 1.4617 | ✅ |
| §3.6 | ery ΔlogFC −2.48, q 0.0073 | true, but D7-vs-OXA compared to neuronal D7-vs-PTX | ⚠️ F3 |
| §3.4 | IFN +0.228, z +4.61, 75%, p 0.002, q 0.011 | 0.2282, 4.6111, 75.43, 0.001998, 0.010989 | ✅ |
| §3.4 | glial +0.155, z +3.27, 61%, q 0.011 | 0.1545, 3.2707, 61.48, 0.010989 | ✅ |
| Table 4 | 8 rows × 6 columns | all 48 cells match | ✅ |
| Table 4 | selection of 8 rows | 12 modules reach q < 0.05; cytokine (+2.83, q 0.025) omitted | ⚠️ F11 |
| §3.4 | ρ 0.991, n 132; mean abs 0.057/0.033/0.712; 20.5% | 0.991145, n 132; 0.05680/0.03284/0.71193; 20.5133 | ✅ |
| §3.4 | "inflates magnitude by roughly 70%" | 1.7296 → +73.0% | ✅ |
| §3.7 | min ρ 0.959, median 0.988 (discovery); min 0.989 (Δ) | 0.95935, 0.98758, 0.98871 | ✅ |
| §3.7 | top five identical in all 60 settings | `min_top5_overlap` = 5 of 5 | ✅ |
| §3.7 | 60 + 15 settings | 60 discovery + 15 delta = 75 | ✅ |
| §3.7 | 316–325 genes across spans 0.2–0.6 and trend off | 316/319/320/321/325, none = 323 | ✅ values, ⚠️ F8 (reference = 320, §3.2 says 235) |
| §3.7 | Spearman t ≥ 0.9999 | min 0.99999624 | ✅ |
| §3.7 | cross-study network-vs-additive 0.929 | 0.928571 (at D14 only, n = 22 modules) | ◐ n undisclosed |
| §3.8 | 55,967 cells; per-sample counts | exact match | ✅ |
| §3.8 | 1,517 up / 703 down | 1517 / 703 over 34,603 genes | ✅ |
| §3.8 | Sprr1a +1.55 t 54.5; Atf3 +0.84 | +1.5539 t 54.493; +0.8355 | ✅ |
| §3.8 | Gadd45a, Serpinb1a, Ecel1 up | +0.5267, +0.7716, +0.1020, all q < 1e-4 | ✅ |
| §3.8 | no ISG up | 0 of 38 panel ISGs up at q < 0.05 | ✅ |
| §3.8 | Ifit2 q 0.042, Oas1c q 0.014 | 0.0417; 0.0136 | ⚠️ F7 (Oas1c not a panel member) |
| §3.8 | Cd68/Aif1 unchanged | −0.0005 (q 0.783) / +0.0009 (q 0.920) | ✅ |
| §3.8 | Gfap "nominally up, q 0.006" | +0.0036, q 0.0061 | ⚠️ F13 |
| §3.8 | myeloid 0.12% vs 0.05%, 55% drop | 0.1170% vs 0.0517%, 55.77% | ✅ |
| §3.8 | d 1.37 vs MDE 2.38; U 2/14, p 0.114 / 0.057 | 1.3685; 2.38; U=2, 0.11429, 0.05714 | ✅ |
| §3.8 | Ifit2 detected at −0.048 | −0.0476 | ✅ |
| §3.8 | ρ 0.098 (D7), 0.112 (D1), 0.030 (D14) | D7 = 0.09845 ✓; D1/D14 in no artefact, no script | ⚠️ F9 |
| §3.8 | "8 compartments", no ISG up in any | artefact has 8; 0 ISG up; **myeloid & TCell dropped** | ⚠️ F5, F6 |
| §3.8 | Ifit2 neurons q 0.046; Isg15 sat. glia q 0.007 | 0.0459; 0.0071 | ✅ |
| §3.8 | (omitted) Ifit1 sat. glia | −0.523, q 0.0037 — strongest ISG hit, unreported | ⚠️ F6 |
| Table 5 | 1,517/703; 55,967; markers; per-compartment | as above | ✅ except F5/F7/F13 |
| §2.6 | closed form matches 80 iters to 2.2e-16 | 2.776e-16 on a 300×300 random W | ✅ |
| §2.6 | no self-loops (diag 0, not −1) | `_limma.py:239,245` | ✅ |
| §2.6 | force-inclusion of target genes | `20_multistage:415-416`; `22_etiology:578-585` | ✅ |
| §2.6 | reversal = cosine(P·seed, mean(control)−mean(case)), both centred | `cos_batch` at `_limma.py:325-336` | ✅ |
| §2.6 | oracle = same-size seed, correct sign | `20_multistage:444-448`, k = n_targets | ✅ |
| §2.6 | degree-matched null, deciles | `degree_matched_sets(n_bins=10)`, `_limma.py:360-383` | ✅ |
| §2.6 | ridge LOSO | `_limma.py:457-464` | ✅ |
| §2.7 | exhaustive 4,900 / 70 enumeration | n_null = 4900 (4 pairs) and 70 (2 paired pairs) | ✅ |
| §2.7 | p floor 1/4,901 and 1/71 | min p observed = 3/4901 and 3/71; floor claims are the *bound* | ✅ |
| §2.7 | pooled-network variant "reported" | computed (528 rows) but absent from manuscript | ⚠️ F12 |
| §2.2 | 13,893 / 13,601 genes | 13,893 (`run_summary`), 13,601 (`etiology_summary`) | ✅ |
| §2.2 | "6–9% global offset" | group means give 1.9–5.1% | ⚠️ F14 |
| §2.3 | 20-gene muscle panel | `_limma.py:141-143`, 20 genes | ✅ |
| §2.3 | MAD robust z on log10 share, 1% threshold | `20_multistage:279-289` | ✅ |
| §2.4 | limma-trend, lowess span 0.4 | `fit_limma(trend=True, frac=0.4)` default | ✅ |
| §2.4 | "2–6 residual df" | realised values are 5 and 6 only | ⚠️ F15 |
| §2.5 | 17 targets + 5 controls | `_limma.py:158-197` | ✅ |
| §2.9 | statsmodels fails loudly | `_limma.py:36-43` raises SystemExit | ✅ |
| §2.9 | 14 test groups | T1-T14 in `99_unit_tests.py` | ✅ |
| §2.9 | requirements pin every version | all 5 versions match live interpreter | ✅ |
| §5 | LOO ranges deposited | `reversal_loo_min/max/sd` present | ✅ |
| Data availability | GitHub URL + "to be archived at Zenodo" | repo is **not** a git repository; nothing deposited | ⚠️ see below |

## Discrepancies found, summarised
**9 hard mismatches** (F3 ×2 cells, F4, F7, F8, F9, F13, F14, F15),
**4 undisclosed-source / scope issues** (F2, F6, F11, F12),
**1 missing-analysis issue** (F5), **1 orphan-artefact issue** (F16).
Everything else I checked — roughly 230 numeric cells — matched exactly.

## Note on the Data and code availability statement
The statement says artefacts "are deposited in the project repository
(https://github.com/yyx-4113/cancer-pain-drg-virtual-perturbation), to be
archived at Zenodo upon acceptance". I verified the local directory is **not a git
repository** (`git rev-parse` → fatal), so nothing has been pushed from here. I
cannot verify the remote URL's existence from this environment. The statement is
*probably* intended as "will be deposited", but as written ("are deposited …
to be archived at Zenodo") it asserts a present state that does not exist for the
Zenodo half and cannot be checked for the GitHub half. Combined with F16
(orphan artefacts from an undisclosed dataset shipping alongside the real ones)
and F1 (a checker that overstates its coverage), the availability statement is
the weakest section of the paper.

```
Paste-ready replacement:

  All data are public: GSE149648 and GSE160543 (NCBI Gene Expression Omnibus)
  for discovery, and GSE292522 (NCBI GEO) for single-cell validation. The
  analysis code, the deposited result tables and the test suite are available at
  https://github.com/yyx-4113/cancer-pain-drg-virtual-perturbation and will be
  archived at Zenodo with a DOI upon acceptance. Every number in Tables 1-5 and
  in the Results is regenerated from the tables in results/tables/ by
  scripts/regenerate_ms_tables.py for Tables 2-4; Table 1 is regenerated by
  scripts/20_multistage_pipeline.py (module_level_signal.csv) and Table 5 by
  scripts/31_process_sc.py and scripts/32_percelltype_sc.py. A self-consistency
  checker (scripts/98_check_manuscript_numbers.py) verifies table cells against
  these artefacts; it does not re-derive the artefacts from the raw GEO member
  files. Files under results/tables/superseded/ relate to an earlier pilot on a
  different dataset and are not part of the reported analysis. requirements.txt
  pins every package version used.
```

## Reproducibility hygiene summary
- **Regenerated by running scripts:** Tables 2, 3, 4 (via
  `regenerate_ms_tables.py`, verified by execution). Table 1 via
  `20_multistage_pipeline.py`. Table 5 via `31_process_sc.py` + `32_percelltype_sc.py`.
- **Hand-edited after generation:** Table 3's CTRL erythrocyte Δ and q (F3).
- **Hardcoded selection inside the generator:** Table 4's 8-module `order` list
  (`regenerate_ms_tables.py:86-87`), which reproduces the omission rather than
  flagging it (F11).
- **Not covered by any check:** `sc_GSE292522_celltype_scores.csv`,
  `sc_GSE292522_celltype_summary.json`, `sc_GSE292522_composition.csv`,
  `sc_GSE292522_focus_markers.csv` (F1). The composition file in particular
  carries the myeloid cell counts that decide F5, and nothing checks it.
- **Computed but never reported:** the entire `pooled_by_study` scope, 264 rows
  (F12).
- **Reported but not computed anywhere:** concordance ρ at D1 and D14 (F9).
- **Orphan artefacts:** `perturbation_bone.csv` has no producing script in the
  repository; `signature_bone_*.csv` and `VALID_GSE249443_*` come from the
  superseded GSE249443 pilot (F16).
