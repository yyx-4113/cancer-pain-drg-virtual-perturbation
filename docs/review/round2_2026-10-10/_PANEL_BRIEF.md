# Round 2 Independent Review Panel — Shared Brief

**Manuscript under review:** `docs/manuscript/manuscript_final.md`
(~7,440 words; single-author computational reanalysis)

---

## 1. What this manuscript is

A single-author computational reanalysis of public bulk rat DRG (dorsal root
ganglion) RNA-seq, comparing two cancer-pain aetiologies, plus an independent
single-cell validation.

**Discovery data (bulk, rat DRG, same platform GPL24782):**
- GSE149648 — Walker-256 bone cancer pain (tibial implantation); sham, d1, d7, d14; n=4/group
- GSE160543 — paclitaxel / oxaliplatin / vehicle neuropathy; n=4/group

**Validation data (single-cell, mouse DRG, 10x Genomics):**
- GSE292522 — B16F10 cancer-nerve-injury, 4 animals + 4 littermate controls, 55,967 cells

**The two claims the paper makes:**
1. *Biological (positive):* bone cancer pain and CIPN drive qualitatively
   different DRG programmes — myeloid/interferon-driven in the Walker-256
   tibia model, neuronal-injury-driven in the other — demonstrable on canonical
   markers with exact permutation tests.
2. *Methodological (negative):* the same data do **not** support ranking
   analgesic targets by in-silico perturbation. The paper's position is that this
   layer fails its own negative controls, degree-matched nulls and additive
   baseline, and is additionally underpowered at n=4.

**Target venue class:** a SCIE journal that is tolerant of exploratory / neutral /
computational evidence (the author's constraint). Reporting standard in play:
ARRIVE-style animal reporting plus reproducibility discipline.

---

## 2. Independence discipline (MANDATORY)

**Forbidden to read:**
- `docs/review/REVIEW_round1_2026-10-10.md`
- `docs/review/R3_editorial_review.md`
- `docs/review/round1_2026-10-10/*` (all Round-1 expert files, incl. `_PANEL_BRIEF.md`)
- `docs/manuscript/manuscript_v0.1.md`, `manuscript_v0.2.md`, `manuscript_v0.3.md`
- Any task-status / project-overview / memory file under `.workbuddy/`
- `SUBMISSION_MANIFEST.md`, deposit SOP, author verification statement

**Also forbidden:** reading any other reviewer's output file in this directory.

**Required posture:**
- Do **not** assume the manuscript is mature, polished, or has passed prior review.
- Treat it as a **first submission**. No credit for effort, no "they already
  handled this".
- Every judgement must come from text or source data you read yourself.
- Any claim in the manuscript that you CAN verify, you **MUST** verify.

---

## 3. Source data you must recompute from (mandatory)

Do not trust the manuscript's numbers. Recompute from these artefacts:

| Artefact | Contents |
|---|---|
| `results/tables/etiology_marker_de.csv` | per-gene, per-contrast logFC/t (muscle-adjusted limma-trend) |
| `results/tables/etiology_marker_cross.csv` | per-sample simple-mean logFC — source for Δ columns |
| `results/tables/etiology_bootstrap_differences.csv` | 264 contrast-level perturbation tests: z_vs_randset, q_randset, ci95_lo/hi, p_permutation_exact |
| `results/tables/etiology_reversal_point_estimates.csv` | reversal/add/ridge point estimates per module × aetiology |
| `results/tables/sc_GSE292522_pseudobulk_DE.csv` | whole-DRG pseudobulk DE |
| `results/tables/sc_GSE292522_composition.csv` | per-sample myeloid fraction |
| `results/tables/sc_GSE292522_summary.json` | single-cell summary |
| `results/tables/sc_GSE292522_celltype_DE.csv` | **NEW** per-compartment DE (8 compartments) |
| `results/tables/sc_GSE292522_celltype_scores.csv` | **NEW** per-compartment marker/ISG scores |
| `results/tables/negative_result_quantification.json` | **NEW** T1/T2/T3/T4/T5 |
| `scripts/` | all analysis code — read the implementations, do not assume |

**A review that does not recompute numbers is incomplete.** State which values you
recomputed and whether they matched.

---

## 4. Environment / interpreter notes (do not rediscover these)

- Python: `C:/Users/Administrator/.workbuddy/binaries/python/envs/default/Scripts/python.exe`
  (has numpy, scipy, pandas, statsmodels, matplotlib).
- Run: `cd "d:/2026.10/不同部位癌痛患者对药物敏感性生信分析/cancer_pain_virtual_perturbation" && <python> your_script.py`
- **Heredoc works**: `"<python>" - <<'PY' ... PY`
- The 10x matrices are **genes × cells** (6514 × ~55k), gzipped `.mtx`.
- `scipy` sparse `.sum(axis=1)` returns a `np.matrix` — it has **no** `.toarray()`;
  use `np.asarray(...).ravel()`.
- Reading all 8 single-cell matrices takes ~8 min. Budget for it; run long jobs in background.
- If a dispatch fails on a tool error, retry with `run_in_background: true`.
- Write files with the Write/Edit tool; do not try to end the whole task.
- Do not mention what tools you use; write review comments only.

---

## 5. Output contract (MANDATORY — four parts per item)

Every finding must have all four:

- **【Problem】** one sentence
- **【Evidence】** pinned to file:line, or table/section + exact numbers;
  numbers you cite must be ones **you recomputed yourself**
- **【Why it matters】** concrete effect on conclusions / credibility / acceptance
- **【Specific fix】** a paste-ready English replacement sentence, **or** an explicit
  spec for a new analysis (variables, strata, output columns)

> "Consider strengthening the discussion" is **banned**.

### Also required in every report

- **§ Stands up (≥3, with evidence)** — explicitly mark things you suspected but
  found to be correct. This is a deliverable, not filler.
- **§ Questions for the authors** — state what you need to know; do not guess answers.
- **§ What I actually checked** — files read, commands run, values recomputed vs
  the manuscript's, with any discrepancy stated explicitly.

Write your report to `docs/review/round2_2026-10-10/<your_codename>_<role>.md`.
Expect 300–650 lines.
