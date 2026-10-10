# Round 5 Independent Review Panel — Shared Brief (confirmation round)

**Manuscript under review:** `docs/manuscript/manuscript_final.md`
(single-author computational reanalysis; ~6,150 words to References; 295-word abstract;
16 references; 3 figures)

**Round purpose:** This is a **final confirmation round**. The author has run several
internal review-and-revise cycles and needs a fresh, independent panel to decide
whether the manuscript is now **publishable at a soundness-oriented SCIE venue**, or
whether any substantive problem remains. Be decisive. If you find no blocking issue,
say so plainly and give your acceptance-level verdict.

---

## 1. What this manuscript is

A single-author computational reanalysis of public bulk rat DRG RNA-seq comparing two
cancer-pain aetiologies, plus an independent single-cell validation.

- **Discovery (bulk rat DRG, same platform GPL24782):** GSE149648 (Walker-256 bone cancer
  pain; sham, d1, d7, d14; n=4/group) and GSE160543 (paclitaxel/oxaliplatin; n=4/group).
- **Validation (single-cell mouse DRG, 10x):** GSE292522 (B16F10 cancer-nerve-injury,
  4 animals + 4 controls, 55,967 cells).

**Claims:** (1) *Biological:* bone cancer pain and CIPN drive qualitatively different DRG
programmes (myeloid/interferon-driven in Walker-256; neuronal-injury-driven in the other),
demonstrable on canonical markers with exact permutation tests. (2) *Methodological
negative:* in-silico perturbation cannot rank analgesic targets — it fails its own negative
controls and is anti-selective, while remaining powered for large effects.

**Venue class:** SCIE journal tolerant of exploratory/neutral/computational evidence;
ARRIVE-style animal reporting plus reproducibility discipline.

---

## 2. Independence discipline (MANDATORY)

**Forbidden to read:** any file under `docs/review/` other than THIS brief; any earlier
manuscript version (`manuscript_v0.*`); any `.workbuddy/` file; `SUBMISSION_MANIFEST.md`
and deposit documents.

**Posture:** treat as a **first submission**. Every judgement from text/source you read.
Any claim you can verify, you **must** verify against the artefacts. No softening, no
invented problems. Each issue needs a file/line citation and, where possible, a
recomputed number. Distinguish genuine blocking issues (Major) from cosmetic ones (Minor).

---

## 3. Source data (mandatory)

Project root `cancer_pain_virtual_perturbation/`. Python:
`C:/Users/Administrator/.workbuddy/binaries/python/envs/default/Scripts/python.exe`

`results/tables/`: `etiology_marker_cross.csv`, `etiology_module_logfc_cross.csv`,
`etiology_bootstrap_differences.csv` (dedup universe = `network_scope=='per_contrast'` &
`perturbation=='knockout'`, n=132), `negative_result_quantification.json` (T1–T5),
`etiology_marker_loo.json`, `sc_GSE292522_*` (pseudobulk DE, celltype DE, composition,
compartment LOO, gfap purity, concordance), `DE_D*_vs_Sham_adj.csv`.
Figures: `docs/manuscript/figures/fig1..3.{png,pdf}`.
Gate: `scripts/98_check_manuscript_numbers.py` — run it; it should end
`291 checks, 0 failed`. Audit whether its expectations are runtime-derived or hardcoded.

**Re-derive, don't trust:** exact permutation p and its endpoint-dependent floor
(3/4,901 unpaired, 3/71 paired logFC, 2/71 paired reversal); Cohen's d denominator;
T5 MDE = 2.802 × median bootstrap SE; T4 target percentile; dedup mirror r = −1.000;
the erythrocyte bin's 2-control/4-tumour imbalance; the log₂(counts-per-10,000+1) scale.

---

## 4. Output

Write your report to `docs/review/round5_2026-10-10/<your-file>.md`:

```
# R5-X Review — <role>
## 1. Verdict (Accept / Minor Revision / Major Revision / Reject)
## 2. Summary of the manuscript
## 3. Major issues (numbered: citation, evidence, required fix)
## 4. Minor issues (numbered)
## 5. Numbers you independently recomputed (claim | manuscript | recomputed | match?)
## 6. Publishability assessment (ready? if not, the ordered minimum)
```

---

## 5. Round-5 focus

Confirm there is **no remaining blocking issue**. Specifically re-check: (a) the negative
layer's internal consistency (Abstract / §3.6 / Limitations / Conclusions) and its numbers;
(b) §2.8/§3.8 single-cell claims against artefacts (compartment counts, composition,
Gfap, CD45 blind spot); (c) every Table 2/3/5 cell vs the CSVs; (d) reference numbering
now runs in first-appearance order and every DOI resolves; (e) figures match their captions
and cited files; (f) declarations (Funding/COI/CRediT/Acknowledgements/Ethics/AI-use/data).
Report any stale number, misstatement, or new problem. If none, say so.
