# Round 4 Independent Review Panel — Shared Brief

**Manuscript under review:** `docs/manuscript/manuscript_final.md`
(single-author computational reanalysis; ~6,150 words to References, 295-word abstract)

**Round purpose:** This is a **post-revision confirmation round**. The author has
completed an internal revision cycle and now needs an *independent* panel to
decide whether the manuscript has reached publishable state, or whether any
substantive problem remains. Render a decisive verdict.

---

## 1. What this manuscript is

A single-author computational reanalysis of public bulk rat DRG (dorsal root
ganglion) RNA-seq, comparing two cancer-pain aetiologies, plus an independent
single-cell validation.

**Discovery data (bulk, rat DRG, same platform GPL24782):**
- GSE149648 — Walker-256 bone cancer pain (tibial implantation); sham, d1, d7, d14; n=4/group
- GSE160543 — paclitaxel / oxaliplatin / vehicle neuropathy; n=4/group

**Validation data (single-cell, mouse DRG, 10x Genomics):**
- GSE292522 — B16F10 cancer-nerve-injury, 4 animals + 4 controls, 55,967 cells

**The two claims the paper makes:**
1. *Biological (positive):* bone cancer pain and CIPN drive qualitatively
   different DRG programmes — myeloid/interferon-driven in the Walker-256
   tibia model, neuronal-injury-driven in the other — demonstrable on canonical
   markers with exact permutation tests.
2. *Methodological (negative):* the same data do **not** support ranking
   analgesic targets by in-silico perturbation. The paper's position is that
   the layer fails its own negative controls and is anti-selective (negative
   controls score higher than analgesic targets), while remaining adequately
   powered for large effects.

**Target venue class:** a SCIE journal tolerant of exploratory / neutral /
computational evidence. Reporting standard in play: ARRIVE-style animal
reporting plus reproducibility discipline.

---

## 2. Independence discipline (MANDATORY)

**Forbidden to read:**
- Any file under `docs/review/` other than THIS brief (`round4_2026-10-10/_PANEL_BRIEF.md`).
  Do **not** open `round1_*`, `round2_*`, `round3_*`, `REVIEW_round*.md`,
  `R3_editorial_review.md`, or any other reviewer's output.
- Any earlier manuscript version (`manuscript_v0.*`). Only `manuscript_final.md`.
- Any task-status / project-overview / memory file under `.workbuddy/`
- `SUBMISSION_MANIFEST.md`, deposit SOP, author verification statement

**Required posture:**
- Treat the manuscript as a **first submission**. No credit for effort.
- Every judgement must come from text or source data you read yourself.
- Any claim you CAN verify, you **MUST** verify against scripts/artefacts.
- Do not soften findings; do not invent problems either. Each issue needs a
  file/line citation and, where possible, a recomputed number.
- **Be decisive about publishability.** If a problem is cosmetic, say Minor. If
  it would stop you recommending acceptance at a soundness-oriented SCIE
  journal, say Major and specify exactly what must change.

---

## 3. Source data you must recompute from (mandatory)

All paths relative to the project root
`cancer_pain_virtual_perturbation/`. Python:
`C:/Users/Administrator/.workbuddy/binaries/python/envs/default/Scripts/python.exe`
(pandas / numpy / scipy / statsmodels installed).

**Artefact tables (`results/tables/`):**
- `etiology_marker_cross.csv` — per-marker cross-aetiology ΔlogFC with exact
  permutation p (`p_permutation_exact`), BH q (`q_permutation`), `paired_control`
- `etiology_module_logfc_cross.csv` — same for module mean logFC
- `etiology_bootstrap_differences.csv` — per-contrast per-module Δ with
  bootstrap CI, `z_vs_randset`, `p_randset`, `q_randset`, exact permutation p,
  `kind` (target/control), `perturbation` (knockout/overexpression)
  — knockout and overexpression are exact mirrors (r = −1.000); the
  deduplicated hypothesis universe is `network_scope=='per_contrast'` &
  `perturbation=='knockout'` (n=132)
- `negative_result_quantification.json` — T1/T2/T3/T5 summaries
- `etiology_marker_loo.json`, `etiology_marker_de.csv` — marker leave-one-animal-out
- `sc_GSE292522_celltype_DE.csv`, `sc_GSE292522_pseudobulk_DE.csv`,
  `sc_GSE292522_composition.csv`, `sc_GSE292522_compartment_loo.json`,
  `sc_GSE292522_gfap_purity_check.csv`, `sc_GSE292522_concordance_per_day.json`
- `DE_D1_vs_Sham_adj.csv`, `DE_D7_vs_Sham_adj.csv`, `DE_D14_vs_Sham_adj.csv`

**Figures (`docs/manuscript/figures/`):** fig1/fig2/fig3 in .png and .pdf.
Check every plotted number is consistent with the CSV it cites.

**Scripts (`scripts/`):** 01–37 analysis; **98_check_manuscript_numbers.py** —
the reproducibility gate (run it yourself: it should end with
`287 checks, 0 failed`). Audit whether the gate's assertions actually derive
expectations from artefacts (vs hardcoded literals).

**Key statistical conventions to audit (re-derive, don't trust):**
- Exact permutation p: null over C(8,4)=70 split pairs × 70 partner splits =
  4,900; the two-sided convention counts both the observed statistic and its
  exact mirror, giving floor 3/4,901 for unpaired contrasts. For the two
  shared-control (paired) contrasts the enumeration is 70 split pairs; the
  manuscript claims a floor of 3/71 for the module/marker logFC endpoint and
  2/71 for the reversal endpoint. **Verify this endpoint dependence is real**
  (implementation in `scripts/22_etiology_comparison.py`), not a tolerance
  artefact.
- Cohen's d denominator, MDE derivation, dedup logic (r = −1.000), binomial
  base rate — recompute each.

---

## 4. What each reviewer must produce

Write your full report to `docs/review/round4_2026-10-10/<your-file>.md`
using this heading structure:

```
# R4-X Review — <role>
## 1. Verdict (one of: Accept / Minor Revision / Major Revision / Reject)
## 2. Summary of the manuscript (shows you actually read it)
## 3. Major issues (numbered, each with: citation, evidence, required fix)
## 4. Minor issues (numbered, same format)
## 5. Numbers you independently recomputed (table: claim | manuscript | recomputed | match?)
## 6. Publishability assessment (is it ready? if not, the ordered minimum to get there)
```

Be exact: manuscript line numbers or section numbers, script file/line for
code issues, recomputed numbers with 4+ decimals where relevant.

---

## 5. Round-4-specific audit targets

Beyond a normal first-submission review, this round must explicitly audit:

1. **Internal consistency of the negative-result layer** — do §3.6, the
   Abstract, the Limitations and the Conclusions tell one coherent story about
   anti-selectivity and power? Recompute T1 (Mann–Whitney on |z|, dedup KO
   universe), T2 (binomial on randset-significant controls vs base rate),
   T5 (MDE = 2.802 × median bootstrap SE; how many of 132/102 tests clear it).
2. **The single-cell section (§3.8)** — check every compartment claim against
   `sc_GSE292522_celltype_DE.csv`, the composition displacement, the
   myeloid-blind-spot disclosure (CD45 depletion), the satellite-glia / Gfap
   caveat, and the Baruch et al. disagreement, against artefacts.
3. **Table 2/3/5 provenance** — verify every cell against the CSVs; verify the
   gate does this correctly (runtime-derived expectations).
4. **Journal-submission readiness** — abstract length/keywords/running head vs
   a real SCIE venue; declarations (funding/COI/CRediT/AI-use/data
   availability); the data-availability statement's placeholder URL.
5. **Survivorship of old errors** — grep for suspicious literals and recompute
   any number you can; report any that contradict the current artefacts.
6. **Figures** — do they plot what the captions say, from the cited files?
