# Supporting Information

## S1 Text. Data and code availability

This manuscript is a reanalysis of publicly deposited transcriptome data; no new primary data were generated.

- **Discovery datasets.** GSE149648 (rat DRG, Walker-256 tibial cancer pain model, sham and days 1/7/14; n = 4 per group) and GSE160543 (rat DRG, paclitaxel/oxaliplatin/vehicle neuropathy; n = 4 per group), both on platform GPL24782, retrieved from NCBI GEO.
- **Validation dataset.** GSE292522 (mouse DRG, 10x Genomics single-cell RNA-seq; B16F10 melanoma versus keratinocyte control, n = 4 + 4 animals), retrieved from NCBI GEO.
- **Analysis code and artefacts.** The complete analysis pipeline (scripts 00-99), all intermediate CSV/JSON artefacts, and the test suite are deposited at https://github.com/yyx-4113/cancer-pain-drg-virtual-perturbation (release v1.0.0). `requirements.txt` pins every software version used (Python 3.13.12, numpy 2.4.6, scipy 1.18.1, pandas 2.3.3, statsmodels 0.15.0, matplotlib 3.11.1).
- **Reproducibility gates.** `scripts/98_check_manuscript_numbers.py` re-derives every reported number from the deposited artefacts at run time (291 checks, all passing); `scripts/99_unit_tests.py` asserts the properties that past bugs violated (14 test groups, all passing). `MANIFEST.sha256` records the SHA-256 of every deposited artefact.
- **Version control.** The repository carries a single tagged release, v1.0.0, which is the version cited throughout the manuscript.

All materials required to reproduce the results are publicly available at the URL above; no login or request is required.
