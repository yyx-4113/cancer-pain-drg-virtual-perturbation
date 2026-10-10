# Divergent DRG programs in bone cancer pain and chemotherapy-induced peripheral neuropathy

**Analysis code and deposited artefacts for the manuscript:**

> *Divergent dorsal root ganglion programs in bone cancer pain and chemotherapy-induced peripheral neuropathy: an exact cross-study reanalysis and a negative-control calibration of in-silico perturbation*
> Yongxin Yang (single author). Manuscript in preparation, 2026.

This repository holds everything needed to reproduce every number in that manuscript:
the analysis pipeline, the CSV/JSON artefacts each number is traced to, the figures,
a runtime number-checker, and a unit-test suite. It also documents, as a **negative
result**, why in-silico knockout/overexpression perturbation of bulk DRG transcriptomes
does **not** rank analgesic targets — and the negative controls that make that failure
visible.

---

## What the study does

1. **Cross-aetiology reanalysis.** Two rat DRG RNA-seq studies on the same platform —
   GSE149648 (tibial Walker-256 bone cancer pain; sham and days 1/7/14) and GSE160543
   (paclitaxel/oxaliplatin neuropathy) — are compared with empirical-Bayes moderated *t*
   (limma-trend) plus a muscle-contamination covariate. Differences are tested by an
   **exhaustive label permutation** over all C(8,4)² = 4,900 case/control split pairs,
   giving an exact p rather than an asymptotic one.
2. **Negative-control calibration of in-silico perturbation.** Seventeen analgesic-target
   modules are perturbed in silico and calibrated against five negative-control modules,
   degree-matched nulls, an oracle anchor, and additive/ridge baselines. The layer is
   shown to add nothing over an additive model and to fail its own controls.
3. **Independent single-cell check.** GSE292522 (mouse DRG, B16F10 melanoma vs
   keratinocyte control, n = 4 + 4) is analysed as per-animal pseudobulk; the bone-cancer
   interferon signature does not reproduce, and the model's response is neuronal-injury-driven.

---

## Repository layout

```
cancer_pain_virtual_perturbation/
├── README.md                     ← this file
├── CITATION.cff                  ← machine-readable citation metadata
├── LICENSE                       ← MIT
├── requirements.txt              ← pinned environment (exact versions used)
├── MANIFEST.sha256               ← SHA-256 of every deposited artefact
├── GITHUB_DEPOSIT_SOP.md         ← 中文：建仓 / 打 tag / Zenodo 归档操作手册
├── author_verification_statement.md
│
├── scripts/                      ← analysis pipeline + audit helpers
│   ├── 00–02   exploratory proof-of-concept (earliest virtual-perturbation layer)
│   ├── 10–12   bone-pain end-to-end, limma-trend, muscle-covariate adjustment
│   ├── 20–23   final manuscript pipeline: discovery, validation, cross-aetiology,
│   │           sensitivity analysis
│   ├── 30–37   single-cell validation: fetch, process, per-cell-type DE,
│   │           negative-result quantification, per-day concordance, figures,
│   │           marker and compartment leave-one-animal-out
│   ├── 98      runtime number-checker (the "gate")
│   ├── 99      unit-test suite
│   ├── _limma.py  shared limma-trend implementation (used by 11/12/20/22/23)
│   └── _r3_* _r4_* _r5_*  / A3_*  ← session audit helpers, kept for provenance only
│
├── results/
│   ├── tables/                   ← DEPOSITED ARTEFACTS (CSV/JSON) — every number traces here
│   ├── figures/                  ← figure raster/vector outputs
│   └── logs/ , review_a3/        ← run logs and reviewer recompute outputs
│
├── docs/
│   ├── manuscript/
│   │   ├── manuscript_final.md   ← the manuscript
│   │   ├── manuscript_v0.{1,2,3}.md  ← earlier drafts (provenance)
│   │   └── figures/              ← fig1–fig3 (.png + .pdf)
│   └── review/
│       ├── round1…round5_2026-10-10/  ← independent reviewer reports
│       └── REVIEW_round{3,4,5}_2026-10-10.md  ← consolidated verdicts
│
└── review_A2/                    ← round-2 statistician independent recompute scripts
```

> `scripts/_*.py` and `scripts/A3_*.py` are **not** part of the pipeline. They are the
> ad-hoc checks and renumbering helpers used during review rounds 3–5, retained so the
> audit trail is complete. The pipeline proper is the numbered `00–99` set.

---

## Data sources

All data are public; nothing is redistributed here (see `.gitignore`).

| Accession | Organism / model | Role | Fetch |
|---|---|---|---|
| **GSE149648** | rat DRG, tibial Walker-256 bone cancer pain (sham, D1/D7/D14) | discovery | `data/raw/` |
| **GSE160543** | rat DRG, paclitaxel / oxaliplatin neuropathy | cross-aetiology comparator | `data/raw/` |
| **GSE292522** | mouse L3–L5 DRG, B16F10 melanoma vs MPEK-BL6 control, 14 d | single-cell validation | `data/processed/GSE292522/` |
| GSE249443 | rat DRG, bone cancer pain (earlier cohort) | proof-of-concept / pipeline validation | `data/raw/` |

GEO is the canonical source. Download the series matrices / RAW tarballs into `data/raw/`
(or run `scripts/30_fetch_sc.py` for the GSE292522 single-cell matrices). The two bulk
expression matrices that are small enough to redistribute (`data/processed/*.csv`,
≈12 MB total) are committed so the bulk pipeline can be re-run without re-downloading.

---

## Reproducing the analysis

### Environment

```bash
python -m venv .venv
# Windows:  .venv/Scripts/pip install -r requirements.txt
# POSIX:    .venv/bin/pip    install -r requirements.txt
```

Pinned versions (identical to those used for the manuscript):

```
numpy==2.4.6   scipy==1.18.1   pandas==2.3.3   statsmodels==0.15.0   matplotlib==3.11.1
```

`statsmodels` is **not optional**: it supplies the lowess fit used by limma-trend. If it
is missing the run raises `SystemExit` rather than silently degrading to a constant prior
variance — an earlier revision did the latter, which would have changed every moderated
*t* in the study with no trace in the output.

### Run order

```bash
# 1. Cross-aetiology core result  (produces results/tables/etiology_*, DE_*)
python scripts/22_etiology_comparison.py

# 2. Negative-result quantification  (produces negative_result_quantification.json)
python scripts/33_negative_result_quantification.py

# 3. Single-cell validation layers
python scripts/31_process_sc.py            # per-animal pseudobulk
python scripts/32_percelltype_sc.py        # per-cell-type DE + composition
python scripts/34_concordance_per_day.py   # per-day cross-model concordance
python scripts/36_marker_loo.py            # marker leave-one-animal-out
python scripts/37_compartment_loo.py       # compartment leave-one-animal-out

# 4. Figures
python scripts/35_figures.py               # reads results/tables/ at run time

# 5. Checks — MUST both pass
python scripts/98_check_manuscript_numbers.py   # expect: "N checks, 0 failed"
python scripts/99_unit_tests.py                 # expect: ALL PASSED
```

### The number-checker ("the gate")

`scripts/98_check_manuscript_numbers.py` asserts every number reported in the manuscript
against the deposited artefacts and exits non-zero on any mismatch. Its expectations are
**derived from the artefacts at run time** (the perturbed quantities T1/T2/T5 and the
Table 2/3 cells are recomputed from the CSV/JSON), not hardcoded — so a stale artefact
fails the check instead of being "protected" by a stale literal. It is the single
reproducible pass/fail gate for the study.

---

## Caveats and known limitations

These are stated in full in the manuscript; the short version:

- **Negative result, by design.** The in-silico perturbation layer provides no target
  ranking. This is the finding, not a defect. Its failure is visible only because
  negative-control modules, degree-matched nulls and an additive baseline are reported.
- **Two logFC sources.** `etiology_marker_de.csv` (muscle-covariate limma-trend) is used
  for the §3.3 discovery statements; `etiology_marker_cross.csv` (per-sample mean logFC)
  is used for the Δ columns, q values and the exhaustive permutation, because a per-animal
  covariate cannot be carried inside the split enumeration. The two differ in the third
  decimal; §2.3 discloses this.
- **Exact-p floor is endpoint-dependent.** The paired permutation floor is 2/71 for the
  reversal endpoint and 3/71 for the marker-logFC endpoint (the reversal statistic's
  complement is not an exact mirror because the control mean is held fixed).
- **Single-cell caveats.** Cell-type assignment is argmax over overlapping marker sets
  (coarse, no independent validation); the myeloid compartment is unmeasurable after
  CD45⁺ depletion before sequencing; the erythrocyte bin is retained with an unequal
  2 control / 4 tumour split; Gfap is essentially absent from the satellite-glia bin.
- **Data availability.** Code and artefacts are deposited here. The manuscript's
  Data-and-code-availability statement points to this repository; the archived Zenodo DOI
  is added on acceptance.

---

## Software and hardware

- Python 3.13 (see `requirements.txt` for the exact numeric stack).
- All analyses are CPU-only and single-machine; no GPU or cluster is required.
- The single-cell matrices are large (~679 MB compressed); peak RAM for the pseudobulk
  step is a few GB.

---

## Citation

If you use this code or the derived artefacts, please cite the manuscript (see
`CITATION.cff`). The repository itself is citable via its version tag / Zenodo DOI.

## License

MIT — see `LICENSE`. The **data** remain under the terms of their GEO submitters.

## Author

Yongxin Yang — Department of Anesthesiology, The Second Affiliated Hospital of Fujian
University of Traditional Chinese Medicine, Fuzhou, Fujian 350003, China.
ORCID [0009-0004-9698-6552](https://orcid.org/0009-0004-9698-6552) · 960856791@qq.com
