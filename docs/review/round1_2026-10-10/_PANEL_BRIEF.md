# Independent Review Panel — Round 1 Brief

**Manuscript under review:** `docs/manuscript/manuscript_final.md`
"Divergent dorsal root ganglion programs in bone cancer pain and chemotherapy-induced peripheral neuropathy: a calibrated in-silico perturbation analysis"
Single author: Yongxin Yang. Target: a methods-aware computational-biology / pain-neuroscience journal (e.g. Journal of Pain Research / BMC Medical Genomics tier). Article type: original research.

## Independence discipline (mandatory)
Forbidden to read, in whole or in part:
- `docs/review/R3_editorial_review.md` and any `REVIEW_*.md`, `RESPONSE_*.md`, `REVISION_*.md` in this repo
- any task-status file, project-overview file, `SUBMISSION_MANIFEST.md`, deposit SOP, `author_verification_statement.md`
- **other reviewers' output files in `docs/review/round1_2026-10-10/`** — do NOT read them
- any prior review/response written by the author or editor

Do not assume the manuscript is mature or has passed previous review. Treat it as a first submission.
Every judgement must come from text or source data you read yourself.
**Any claim in the manuscript that you CAN verify, you MUST verify** with your own recomputation.

## Output contract (every item needs all four parts)
- 【Problem】 one sentence
- 【Evidence】 pinned to `file:line`, or table/section + exact numbers; numbers you cite must be ones you recomputed yourself
- 【Why it matters】 concrete effect on conclusions / credibility / acceptance
- 【Specific fix】 a paste-ready English replacement sentence, or an explicit spec for a new analysis (variables, strata, output columns)

"Consider strengthening the discussion" is banned.

## Also required
- § Stands up (≥3, with evidence) — explicitly mark things you suspected but found to be correct.
- § Questions for the authors — state what you need to know; do not guess answers.
- § What I actually checked — files read, commands run, values recomputed vs the manuscript, with the discrepancy stated.

## Where the data live (recompute against these)
- `results/tables/run_summary.json`, `qc_muscle_scores.csv`, `de_muscle_drop_sensitivity.csv`, `module_level_signal.csv`, `virtual_perturbation_scores.csv`, `head_to_head_baselines.json`
- `results/tables/etiology_summary.json`, `etiology_module_logfc_cross.csv`, `etiology_marker_cross.csv`, `etiology_marker_de.csv`, `etiology_bootstrap_differences.csv`
- `results/tables/sc_GSE292522_pseudobulk_DE.csv`, `sc_GSE292522_focus_markers.csv`, `sc_GSE292522_composition.csv`, `sc_GSE292522_summary.json`
- `results/tables/sensitivity_summary.json`
- Source GEO series: GSE149648 (rat DRG, Walker-256 tibia TCI), GSE160543 (rat DRG, paclitaxel/oxaliplatin), GSE292522 (mouse DRG, B16F10 cancer-nerve-injury, 10x).

## Environment traps for this dataset
- GSE160543 ships inconsistent per-sample column layouts (six columns for most members, thirteen for some); a fixed column index silently parses the Entrez-id column and yields zero fold change. Confirm the parser asserts non-integer-valued matrices.
- ychr The in-silico perturbation uses a k-NN co-expression graph; target genes must be force-included or most ion-channel/opioid modules collapse to a single node.
- KO and OE are exact mirror images (reversal(KO) = −reversal(OE)); the 264 contrast-level tests therefore collapse to 132 independent hypotheses per network scope. Verify the BH denominator was corrected accordingly and that the permutation endpoint is still null after correction.
- Single-cell validation (GSE292522) is a different species (mouse), model (B16F10 nerve-injury), and quantitation (10x) — its role is to show the interferon programme does NOT generalise, not to confirm it.

## Manuscript's own gate
The repo ships `scripts/98_check_manuscript_numbers.py` (178 checks vs the CSV/JSON). It passing means arithmetic/provenance are internally consistent — it does NOT mean the design is sound. Review the design layer independently.
