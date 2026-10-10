# Panel brief — Round 6 independent review (cancer-pain virtual perturbation)

## Independence discipline (MANDATORY)

Treat this as a **first submission**. You have not seen earlier review rounds.

**Forbidden to read (do not open, do not grep, do not reference):**
- `docs/review/REVIEW_round1_2026-10-10.md` … `REVIEW_round5_2026-10-10.md`
- `docs/review/R3_editorial_review.md`
- `docs/review/round1_2026-10-10/`, `round2/`, `round3/`, `round4/`, `round5/`
- `author_verification_statement.md`
- `GITHUB_DEPOSIT_SOP.md`
- `docs/manuscript/manuscript_v0.1.md`, `manuscript_v0.2.md`, `manuscript_v0.3.md`
- Any other reviewer's output file inside `docs/review/round6_2026-10-10/` (you write one file; you do not read the others).

Every judgement must come from text or source data YOU read. Any number the manuscript cites, you MUST attempt to recompute or trace to a deposited artefact yourself before trusting it. Do NOT assume the manuscript is mature or "passed review".

## What the manuscript is

Single-author computational reanalysis: *Divergent dorsal root ganglion programs in bone cancer pain and chemotherapy-induced peripheral neuropathy: an exact cross-study reanalysis and a negative-control calibration of in-silico perturbation* (author Yongxin Yang).

- **Study 1 (biological):** two same-platform rat DRG bulk RNA-seq studies — GSE149648 (Walker-256 tibial bone cancer pain, sham + D1/D7/D14, n=4/group) and GSE160543 (paclitaxel/oxaliplatin CIPN, n=4/group) — compared by exhaustive label permutation (C(8,4)²=4,900 split pairs) for an *exact* p. Claim: the two aetiologies drive qualitatively different DRG programs (myeloid/IFN up in bone-cancer, down in paclitaxel; galanin opposite), readable from individual markers.
- **Study 2 (methodological, the headline negative result):** 17 analgesic-target modules + 5 negative-control modules perturbed in silico; calibrated against degree-matched nulls, an oracle anchor, an additive baseline. Claim: the perturbation layer adds nothing, is anti-selective (negative controls carry LARGER standardised effects than targets; median |z| 2.47 vs 1.07, MWU p=0.0019), and fails its own controls.
- **Study 3 (single-cell validation):** GSE292522 mouse B16F10 DRG (n=4+4), per-animal pseudobulk. Claim: the Walker-256 day-7 interferon programme does NOT generalise (neuronal-injury programme instead); interferon is model-specific.

## Source data & environment

- Manuscript: `docs/manuscript/manuscript_final.md` (read in full).
- Pipeline scripts: `scripts/00-99` (key: `20_multistage_pipeline.py`, `22_etiology_comparison.py` [core], `33_negative_result_quantification.py` [T1-T5], `_limma.py` [shared stats/perturbation primitives], `31/32/34/36/37_*` [single-cell], `98_check_manuscript_numbers.py` [the gate], `99_unit_tests.py`).
- Deposited artefacts (every number traces here): `results/tables/*.csv` and `*.json` (esp. `etiology_module_logfc_cross.csv`, `etiology_marker_cross.csv`, `etiology_bootstrap_differences.csv`, `negative_result_quantification.json`, `etiology_summary.json`, `sc_GSE292522_*.csv`, `sc_GSE292522_*summary.json`).
- Raw data ARE present: `data/raw/GSE149648`, `data/raw/GSE160543`.
- Python (managed, required): `C:/Users/Administrator/.workbuddy/binaries/python/envs/default/Scripts/python.exe`. Packages pinned in `requirements.txt` (numpy 2.4.6, scipy 1.18.1, pandas 2.3.3, statsmodels 0.15.0, matplotlib 3.11.1). Use this venv; do NOT use a system python.
- Two logFC sources (disclosed in §2.3): `etiology_marker_de.csv` = muscle-adjusted limma-trend fit on the 13,601 shared genes (used for §3.3 single-gene statements); `etiology_marker_cross.csv` = per-sample mean logFC (used for Δ columns, q, the exact permutation). The §3.3 per-gene values (e.g. Cd68 +2.86) therefore differ in the third decimal from the §3.5 Table-3 per-sample values (Cd68 +2.84) — this is by design, not an error, but check it is coherent.

## Mandatory checks (not recomputing = review incomplete)

- Recompute or trace at least THREE key numbers to the artefacts/scripts yourself (pick the ones central to your layer). The gate `98_check_manuscript_numbers.py` is green — but a green gate only proves arithmetic/provenance, NOT design validity. State explicitly which numbers you personally verified and against what.
- Critically scrutinise the **exact-permutation** machinery in `22_etiology_comparison.py` + `_limma.py`: the unpaired null is 4,900 (70×70); the *paired* (shared-control) null is a 70-split "8-case swap" built in `run()` (reversal) and in `cross_module_logfc`/`cross_marker_logfc` (module/marker Δ). Confirm whether `obs` actually lies inside each null and whether the stated floors (3/4,901 unpaired; 3/71 paired marker; 2/71 paired reversal) are correct. (Editor note: the reversal paired null carries the control inside the axis and is consistent; the module/marker paired null omits the control term but is applied only to the within-bone pair whose controls are identical samples, so the control cancels — judge this for yourself rather than trusting this hint.)
- Check the multiplicity correction: KO/OE are exact mirrors, so 264 raw tests → 132 independent; BH applied per scope. Confirm the "no reversal Δ survives correction (min q=0.165)" and "32/132 randset-significant, 14 controls" claims.
- Check the negative-result quantification (T1-T5 in `negative_result_quantification.json`): anti-selectivity MWU p=0.0019; remove-IFN p=0.019; binomial p=0.0067; Cohen's d 1.27 vs MDE 2.38 for myeloid fraction; MDE≈0.052 (normal) / ≈0.10 (t-based) for reversal Δ; 41/102 targets exceed MDE.
- Single-cell: verify Table 5 numbers (1,517 up/703 down; Sprr1a +1.55 t=54.5; myeloid 0.12% vs 0.05%, U=2/14 p=0.114; concordance ρ=0.098) against `sc_GSE292522_*.csv`/`*summary.json`. Assess the argmax-over-marker-set cell assignment and the CD45-depletion / myeloid-blind-spot caveats.
- Data/code availability: manuscript states code at `https://github.com/yyx-4113/cancer-pain-drg-virtual-perturbation` (release v1.0.0). Verify the claim is accurate: is the repo actually public and does tag v1.0.0 exist with the deposited artefacts? Check `git -C . tag`, `git -C . log --oneline -1`, and whether the local tree matches `MANIFEST.sha256`. Report if you cannot confirm the remote is live.

## Output contract (per review item)

Four mandatory parts:
- 【Problem】 one sentence
- 【Evidence】 pinned to file:line, or table/section + exact numbers; numbers you cite must be ones YOU recomputed or traced
- 【Why it matters】 concrete effect on conclusions / credibility / acceptance
- 【Specific fix】 a paste-ready English replacement sentence, or an explicit spec for a new analysis

"Consider strengthening the discussion" is banned.

## Also required

- **§ Stands up** (≥3, with evidence) — things you suspected but found correct. This is a deliverable.
- **§ Questions for the authors** — what you need to know; do not guess answers.
- **§ What I actually checked** — files read, commands/scripts run, values recomputed vs the manuscript, with any discrepancy stated.

Write your report to `docs/review/round6_2026-10-10/<codename>_<role>.md`. Do NOT mention what tools you used; write only review comments.
