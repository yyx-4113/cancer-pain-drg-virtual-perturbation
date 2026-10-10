# A3 — Implementation / Provenance Audit (independent, fresh review)

**Manuscript:** *Divergent dorsal root ganglion programs in bone cancer pain and chemotherapy-induced peripheral neuropathy: an exact cross-study reanalysis and a negative-control calibration of in-silico perturbation*
**Reviewer:** A3 (implementation / provenance)
**Date:** 2026-10-10
**Scope:** I audited only the implementation layer — reproducibility gate, table traceability to deposited artefacts, two-logFC coherence, cross-file drift, git/repo availability, and single-cell QC. I did not re-read other reviewers' files.

---

## 1. Reproducibility gate (run exactly as specified)

I ran the three scripts from the project root (`d:/.../cancer_pain_virtual_perturbation`) with the managed interpreter.

### 1.1 `scripts/98_check_manuscript_numbers.py`
- **Result:** `291 checks, 0 failed` — `ALL MANUSCRIPT NUMBERS VERIFIED AGAINST ARTEFACTS`, exit 0.
- The checker regenerates 68 table cells at run time and asserts ≥30 key manuscript strings (4,600; 3,041; 0.991; 24.9; +3.42; −2.59; 6.26; 3/4,901; 3/71; 2/71; 0.165; Jaccard 0.584; 0.959; 0.989; 264; 132; 20.5; 0.929; 55,967; 32/132; anti-selectiv; 0.0019; binomial p = 0.0067; 46 cells; eight testable compartments; 1.27; ≈0.052; two control and four tumour; etc.).
- **Verdict:** PASS.

### 1.2 `scripts/99_unit_tests.py`
- **Result:** `ALL TESTS PASSED` (14 groups: T1–T14). Covers no-self-loops, KO/OE exact-mirror, identifier-matrix rejection, D14≠D1, degree-matched bins, exhaustive splits contain the true assignment, batched-cosine per-column normalisation, muscle-score direction, ridge LOSO, etc.).
- **Verdict:** PASS.

### 1.3 `scripts/33_negative_result_quantification.py`
- **Result:** ran to completion (exit 0), wrote `results/tables/negative_result_quantification.json`.
- **Diff vs deposited reference `deposited_neg_ref.json`:** I flattened both to 54 leaf keys and compared. **0 differences** (all 54 keys identical, including every nested value).
- **Verdict:** PASS — the freshly regenerated quantification is byte-for-byte equivalent in content to the deposited reference.

**Gate summary:** all three pass; the deposited `negative_result_quantification.json` is confirmed identical in content to what script 33 produces on this tree.

---

## 2. Table traceability — recomputed from the deposited artefacts

I independently recomputed at least five (in practice every) table cell from the CSV/JSON artefacts and compared to the manuscript prose. Source files are listed.

### 2.1 Table 1 — Day-7 module mean moderated *t* (muscle-adjusted)
- **Source:** `results/tables/module_level_signal.csv`, written by `scripts/20_multistage_pipeline.py:401` (the script-20/10/11 pipeline output). Extracted the `contrast == "D7"` rows.
- **Recomputed → manuscript:**
  - CTRL interferon: t = 24.923 → **+24.92** (ms: +24.92 ✓); q = 6.7591e-5 → **6.8 × 10⁻⁵** (ms: 6.8e-5 ✓).
  - Glial/immune: +11.412 → +11.41 ✓; q 0.002667 → 0.0027 ✓.
  - Toll-like receptors: +6.364 → +6.36 ✓; q 0.0245 ✓.
  - Cytokine/chemokine: +5.877 → +5.88 ✓; q 0.0129 ✓.
  - CTRL housekeeping: +0.653 → +0.65 ✓; q 0.524 ✓.
  - CTRL erythrocyte: −4.801 → −4.80 ✓; q 0.0186 → 0.019 ✓.
  - CTRL pan-neuronal: −5.201 → −5.20 ✓; q 0.000314 → 0.00031 ✓.
  - Neurotrophin: −5.292 → −5.29 ✓; q 0.0335 → 0.033 ✓.
  - K⁺/HCN: −5.635 → −5.64 ✓; q 0.00534 → 0.0053 ✓.
  - GABAergic: −5.708 → −5.71 ✓; q 0.0245 ✓.
  - All 22 rows match the manuscript Table 1 within rounding. **Verified.**

### 2.2 Table 2 — Calibration of the in-silico perturbation layer (day-7, discovery)
- **Source:** `results/tables/virtual_perturbation_scores.csv` (within-study reversal scores, D7 knockout rows).
- **Recomputed → manuscript (top-8 rows):**
  - CTRL interferon: rev 0.2282 → **+0.228** ✓; z 4.611 → **+4.61** ✓; %oracle 75.43 → **75.4** ✓; p 0.0020 ✓; q 0.011 ✓.
  - Glial/immune: 0.1545 → +0.155 ✓; z 3.271 → +3.27 ✓; 61.48 → 61.5 ✓; q 0.011 ✓.
  - CTRL pan-neuronal: −0.1086 → −0.109 ✓; z −3.560 → −3.56 ✓; q 0.011 ✓.
  - K⁺/HCN: −0.0916 → −0.092 ✓; z −2.867 → −2.87 ✓; q 0.011 ✓.
  - Voltage-gated Na⁺: −0.0903 → −0.090 ✓; z −2.626 → −2.63 ✓; q 0.020 ✓.
  - CTRL myelin/Schwann: −0.0896 → −0.090 ✓; z −2.707 → −2.71 ✓; q 0.020 ✓.
  - CTRL erythrocyte: −0.0744 → −0.074 ✓; z −2.568 → −2.57 ✓; q 0.020 ✓.
  - Neuropeptides: −0.0675 → −0.067 ✓; z −2.449 → −2.45 ✓; q 0.020 ✓.
- **Footer:** Spearman ρ between network and additive reversal over all 132 rows = **0.9911 → 0.991** ✓; mean |reversal| network **0.0568 → 0.057** ✓, additive **0.0328 → 0.033** ✓, ridge **0.7119 → 0.712** ✓. (Computed directly: `spearmanr(reversal_network, reversal_additive)` over the 132-row file = 0.9911; mean-abs as above.)
- **Verified.**

### 2.3 Table 3 — Cross-aetiology marker difference (exact permutation)
- **Source:** `results/tables/etiology_marker_cross.csv` (per-sample mean logFC; the Δ and q columns).
- **Recomputed → manuscript (key cells):**
  - Cd68 Bone_D7 2.839 → **+2.84** ✓; Δ(D7 vs PTX) 3.418 → **+3.42** ✓; q 0.012487 → **0.012** ✓.
  - Aif1 Bone_D7 2.192 → +2.19 ✓; Δ 3.076 → +3.08 ✓; q 0.012 ✓.
  - Gal Bone_D7 −0.672 → −0.67 ✓; PTX 1.920 → +1.92 ✓; Δ −2.592 → −2.59 ✓; q 0.012 ✓.
  - Ccl2 D7 vs OXA Δ 1.391 → +1.39 ✓; q 0.017 ✓.
  - Il1b D7 vs PTX Δ 1.230 → +1.23 ✓; q 0.018 ✓.
  - Calca D7 vs PTX Δ −0.770 → −0.77 ✓; q 0.027 ✓.
  - Scn10a D7 vs PTX Δ −0.702 → −0.70 ✓; q 0.026 ✓.
  - Atf3 q = 0.1197 (>0.05) → **ns** ✓; Gfap q = 0.0978 (>0.05) → **ns** ✓.
- All 9 rows reconcile. **Verified.**

### 2.4 Table 4 — Cross-aetiology module mean logFC difference
- **Source:** `results/tables/etiology_module_logfc_cross.csv`.
- **Recomputed → manuscript:**
  - Glial/immune: bone D7 2.053 → +2.05 ✓; PTX −0.643 → −0.64 ✓; Δ(D7 vs PTX) 2.695 → **+2.70** ✓; q 0.007345 → **0.0073** ✓.
  - Cytokine/chemokine: Δ 1.702 → +1.70 ✓; q 0.0073 ✓.
  - Toll-like receptors: Δ 1.275 → +1.28 ✓; q 0.0073 ✓.
  - Neurotrophin: Δ −0.972 → −0.97 ✓; q 0.0073 ✓.
  - Ca²⁺ α2δ: Δ −0.600 → −0.60 ✓; q 0.0073 ✓.
  - Neuropeptides: Δ −1.345 → −1.35 ✓; q 0.0145 ✓.
  - CTRL pan-neuronal: Δ −0.691 → −0.69 ✓; q 0.0145 ✓.
  - Glutamate/NMDA: Δ −0.635 → −0.64 ✓; q 0.0162 ✓.
  - Voltage-gated Na⁺: Δ −0.654 → −0.65 ✓; q 0.0310 ✓.
  - GABAergic: Δ −0.495 → −0.49 ✓; q 0.0310 ✓.
  - CTRL erythrocyte: Δ(D7 vs OXA) −2.483 → **−2.48** ✓; q 0.0073 ✓.
- All 11 rows reconcile. **Verified.**

### 2.5 Table 5 — Independent single-cell validation
- **Source:** `results/tables/sc_GSE292522_summary.json`, `sc_GSE292522_pseudobulk_DE.csv`, `sc_GSE292522_composition.csv`.
- **Recomputed → manuscript:**
  - Pseudobulk DE: n_up 1517, n_down 703 → **1,517 up / 703 down** ✓.
  - Sprr1a: logFC 1.554 → **+1.55** ✓; t 54.49 → **54.5** ✓.
  - Atf3: logFC 0.835 → +0.84 ✓ (ms prose).
  - Ifit2 q 0.04165 → 0.042 ✓; Oas1c q 0.01356 → 0.014 ✓; Gfap +0.004, q 0.006 ✓.
  - Myeloid fraction Con 0.001170 → **0.12%** ✓; B16F10 0.000517 → **0.05%** ✓.
  - Myeloid MWU two-sided p 0.114286 → **0.114** ✓; one-sided 0.057 → 0.057 ✓.
  - Convergent ρ vs D7: 0.098455 → **0.098** ✓ (ms also lists 0.112 D1, 0.030 D14 in the same json).
  - Total cells 55,967 → **55,967** ✓.
- **Verified.**

### 2.6 Pivotal anti-selectivity claim (§3.6) — independently recomputed
From `etiology_bootstrap_differences.csv` (per_contrast, knockout, 132 dedup rows: 30 controls, 102 targets):
- median |z| controls = **2.470** → ms 2.47 ✓; targets = **1.073** → ms 1.07 ✓.
- Mann–Whitney two-sided p = **0.0019** ✓.
- n significant vs randset null (q<0.05) = **32** ✓; of which **14 controls / 18 targets** ✓.
- binomial p (14/32 vs base rate 30/132) = **0.0067** ✓.
- min q on exact permutation = 0.1654 → ms 0.165 ✓ (from `etiology_summary.json`).
- **Verified** — the paper's central methodological conclusion is fully reproducible from the artefacts.

### 2.7 §3.1 / §3.2 QC and DE-count verification
Beyond the five tables, I also recomputed the headline QC and DE-count claims from `qc_muscle_scores.csv`, `qc_pca_variance.csv`, `qc_sample_correlation.csv` and the `DE_*_vs_Sham_adj.csv` / `_unadj.csv` files.

- **Muscle score ANOVA (§3.1):** one-way ANOVA on the muscle *score* (mean z-score of the 20 sarcomere genes) across Sham/D1/D7/D14 → **F = 0.637, p = 0.605**, exactly matching the manuscript. (Note: the absolute FPKM-share measure gives a different F=0.984/p=0.433; the manuscript correctly uses the z-score "muscle score," not the share — verified by column.)
- **GSM4506604 contamination (§3.1):** `muscle_fpkm_pct` = 6.257% → **6.26%** ✓, flagged `MUSCLE_CONTAMINATED`; next-highest share is 0.16%, an order of magnitude lower ✓.
- **PCA variance (§3.1):** PC1–PC4 = 46.85 / 17.45 / 7.23 / 4.76% → **46.8%, 17.4%, 7.2%, 4.8%** ✓.
- **Sample correlation (§3.1):** pairwise correlation range = **0.913–0.994** ✓.
- **DE counts at q<0.05, muscle-adjusted 13,893 genes (§3.2):** D1 = 0, D7 = 4,600, D14 = 235 → manuscript "0, 4,600 and 235" ✓. At q<0.10: 5 / 5,867 / 860 ✓.
- **Median |logFC| across tested genes (§3.2):** 0.073 / 0.158 / 0.119 → manuscript "0.073, 0.158 and 0.119" ✓; within the significant sets 0.34 (D7) and 1.20 (D14) ✓.
- **Unadjusted DE counts (§3.2):** D1/D7/D14 = 1 / 2,563 / 958 → manuscript "1, 2,563 and 958" ✓.
- **Leave-one-out drop of GSM4506604 (§3.2):** the manuscript states this reduces the day-7 count to 3,041 (Jaccard 0.584). This figure is produced by `scripts/22b_patch_artifact.py` / `de_muscle_drop_sensitivity.csv`; the LOO sensitivity artefact is deposited and the gate's "3,041" / "Jaccard 0.584" string checks pass, confirming the number is reproduced by the checker. I did not re-run the full LOO here (it requires re-fitting with the flagged sample removed), but the deposited sensitivity artefact and the passing gate are consistent with the stated value.

All §3.1/§3.2 claims reconcile with the artefacts. **No discrepancy.**

---

## 3. Two-logFC coherence (§3.3 +2.86 vs Table 3 +2.84)

- §3.3 Cd68 **+2.86** traces to `etiology_marker_de.csv` (Bone_D7, logFC 2.858). This file is the **muscle-adjusted single-study limma-trend fit on 13,893 genes** — it matches `DE_D7_vs_Sham_adj.csv` (13,893 genes, muscle-adjusted) Cd68 = 2.86199 exactly.
- Table 3 Cd68 **+2.84** traces to `etiology_marker_cross.csv` (Bone_D7, 2.83902), which is the **unadjusted per-sample mean logFC**. It matches `DE_D7_vs_Sham_unadj.csv` (13,893 genes, unadjusted) Cd68 = 2.83902 exactly.
- The 0.02 gap is therefore the **muscle-covariate adjustment** (adjusted 2.862 vs unadjusted per-sample mean 2.839), recomputed on the *same* 13,893 discovery gene set — **not** a stale or contradictory number, and **not** a gene-set nesting effect.
- Cross-gene check: Aif1 de 2.180 (→+2.18 §3.3) vs cross 2.192 (→+2.19 Table 3); Gal de −0.688 (→−0.67) vs cross −0.672 (→−0.67). Both consistent (de slightly lower due to covariate). **No contradiction found.**
- **However**, see Finding F1: the manuscript's own §2.3 prose mis-attributes the §3.3 source.

---

## 4. Stale-number / cross-file drift

- **Repository URL + version:** manuscript §6, `README.md` (lines 176–177), and `CITATION.cff` (lines 17–20) all cite `https://github.com/yyx-4113/cancer-pain-drg-virtual-perturbation` at **release v1.0.0**. `CITATION.cff` version = `1.0.0`, date-released `2026-10-10`. **Agree everywhere.**
- **Accessions:** GSE149648 / GSE160543 / GSE292522 identical across manuscript, README, CITATION.cff. n = 4 per group consistent.
- **Key statistics:** 0.991, 0.929, 4,600/3,041, 55,967, 46 cells, 132/264 — all identical between manuscript and the gate/artefacts (gate passed).
- **One minor naming mismatch (F3):** README §"Data sources" calls the GSE292522 control "MPEK-BL6 control", whereas the manuscript §2.8 and §3.8 call it "non-tumorigenic keratinocytes" / "keratinocyte-injected". MPEK-BL6 is a mouse epidermal keratinocyte line, so the entities are the same, but the label differs and could confuse a reader cross-referencing.
- No other version-string, n-value, or statistic drift detected.

---

## 5. Git / repository availability

- **Local tree at v1.0.0:** `git rev-parse HEAD` = `8882852afc11da21f156136da000c9278f78d30b`; `git rev-list -n 1 v1.0.0` = same commit. **The working tree checked out is exactly the tagged v1.0.0 commit.**
- **Tag:** `v1.0.0` present, message "docs: drop Zenodo from data-availability statement…". Branch `main`.
- **Remote:** `origin` = `https://github.com/yyx-4113/cancer-pain-drg-virtual-perturbation.git`.
- **Remote liveness — UNVERIFIED:** I attempted `git ls-remote origin` (15 s timeout); the remote helper aborted the session (`fatal: remote helper 'https' aborted session`), i.e. the sandbox has no outbound network. **I could not confirm the GitHub repository is actually public/live.** The Data-availability statement therefore remains UNVERIFIED for liveness. (All other deposit integrity checks below pass locally.)
- **Untracked working-tree files (`git status -s`):** `?? .github/`, `?? _review_recompute.py`, `?? deposited_neg_ref.json`, `?? docs/review/round6_2026-10-10/`, `?? err.txt`. Of these, `.github/workflows/release.yml` is a CI release workflow that is **not part of the v1.0.0 tag** (it is untracked, so the tagged tree is clean of it). `_review_recompute.py`, `deposited_neg_ref.json`, `docs/review/round6_2026-10-10/`, `err.txt` are review/working artifacts, also absent from the tag. **Assessment:** none of these affect the reproducibility or integrity of the deposited v1.0.0 tree; the untracked `.github/` is harmless but means the working tree is not byte-identical to the release. Not a credibility blocker, but note it.
- **MANIFEST.sha256 internal consistency:** 190 entries. I recomputed SHA-256 for 8 files (`CITATION.cff`, `README.md`, `requirements.txt`, `results/tables/negative_result_quantification.json`, `etiology_marker_cross.csv`, `virtual_perturbation_scores.csv`, `sc_GSE292522_summary.json`, `docs/manuscript/manuscript_final.md`) and **all 8 matched** the recorded hashes — including the JSON that script 33 just overwrote, confirming the script is deterministic and the manifest is current. **Manifest is internally consistent.**

---

## 6. Single-cell QC (§3.8)

- **Cell counts:** from `sc_GSE292522_composition.csv` × per-sample cell counts in `sc_GSE292522_summary.json` (cells_per_sample): myeloid = **46** total (ms: "46 myeloid cells" ✓); T cells = **1** total (ms: "1 T cell" ✓); total = **55,967** (ms ✓).
- **Myeloid ≥5-cells-in-≥6-animals rule:** myeloid counts per sample = 3, 10, 12, 5, 5, 3, 7, 1 → **≥5 cells in exactly 5 of 8 animals**, matching the manuscript's "with ≥5 cells in only 5 of 8 animals" and the exclusion of myeloid from the testable set. ✓
- **8 testable compartments:** of 12 candidate compartments, the four excluded are T cells (1 cell total), NP neurons, TH neurons, and myeloid (5/8 rule) → **8 retained**. ✓
- **Erythrocyte 2-control / 4-tumour split:** per-sample erythrocyte cell counts = Con: 5, 4, 22, 3; B16F10: 7, 5, 26, 15 → samples with ≥5 erythrocytes = **2 control / 4 tumour**, exactly as the manuscript states ("two control and four tumour animals"). ✓
- **Concordance ρ = 0.098** (D7) recomputed from `sc_GSE292522_summary.json` (0.098455). ✓
- **Verified** — the single-cell QC claims are all backed by the deposited composition/summary artefacts.

---

## 7. Findings (four-part structure)

### F1 — §2.3 mis-states the gene set underlying the §3.3 per-gene values
【Problem】 The Methods (§2.3) claims the §3.3 per-gene values "come from … the muscle-adjusted limma-trend fit on the 13,601 genes shared by both studies," but the actual §3.3 numbers are taken from the 13,893-gene single-study fit.
【Evidence】 §2.3 prose vs `etiology_marker_de.csv` (Bone_D7 Cd68 2.858) which equals `DE_D7_vs_Sham_adj.csv` (13,893 genes, muscle-adjusted) Cd68 = 2.86199. The 2.86-vs-2.84 gap is muscle-covariate (adjusted 2.862 vs unadjusted per-sample mean 2.839), not gene-set nesting. The `README.md` §"Two logFC sources" correctly attributes §3.3 to `etiology_marker_de.csv`; only the manuscript §2.3 is wrong.
【Why it matters】 A reader following §2.3 would believe the 2.86/2.84 difference is a gene-set (13,601 vs 13,893) effect and that §3.3 uses the *shared* fit; neither is true. It is a provenance/clarity defect, not a numerical error, but it undermines the paper's own coherence argument ("same model on nested gene sets, not conflicting measurements") because the two numbers are actually adjusted-vs-unadjusted on the same gene set.
【Specific fix】 Rewrite the §2.3 sentence to: "The per-gene values of Section 3.3 come from the muscle-adjusted single-study limma-trend fit on the 13,893-gene discovery set (`etiology_marker_de.csv`; the same fit as Section 3.2). The cross-aetiology differences of Section 3.5 use the unadjusted per-sample mean logFC (`etiology_marker_cross.csv`), because the exhaustive permutation of Section 2.7 enumerates case/control splits of raw per-sample means and a per-animal covariate cannot be carried inside that enumeration. The muscle score did not differ between groups, so the 2.86 vs 2.84 gap for Cd68 is the covariate adjustment, not a gene-set effect."

### F2 — Remote GitHub repository liveness is unverified
【Problem】 The Data-availability statement points to a GitHub release v1.0.0, but I could not confirm the repository is actually public/live from this environment.
【Evidence】 `git ls-remote origin` aborted (`fatal: remote helper 'https' aborted session`); no outbound network in the sandbox. Local tag `v1.0.0` exists and the tree is at it, but liveness of the *remote* was not established.
【Why it matters】 Reviewers and readers cannot act on an availability statement whose target may be private, renamed, or empty. For a reproducibility-focused submission this is a gap in the verifiable evidence, even though the local deposit is internally consistent (gate passes, MANIFEST consistent).
【Specific fix】 Before acceptance, the author should confirm the repository is public and the v1.0.0 release tag resolves to the deposited tree (e.g., a public clone + `git checkout v1.0.0` + run order in README). State explicitly that the release tag, not just the default branch, is the cited artefact.

### F3 — GSE292522 control label differs between README and manuscript
【Problem】 The README labels the GSE292522 control "MPEK-BL6 control" while the manuscript describes it as "non-tumorigenic keratinocytes" / "keratinocyte-injected."
【Evidence】 `README.md` line 91 ("B16F10 melanoma vs MPEK-BL6 control") vs manuscript §2.8 ("non-tumorigenic keratinocytes, the control") and §3.8 ("keratinocyte-injected animals").
【Why it matters】 MPEK-BL6 is indeed a mouse keratinocyte line, so the entities are the same, but the divergent naming invites a reader to think two different controls are in play and undercuts the cross-document consistency this paper otherwise maintains.
【Specific fix】 Use one consistent descriptor throughout — e.g., "non-tumorigenic keratinocytes (MPEK-BL6)" in both README and manuscript.

### F4 — Working tree contains untracked files not in the v1.0.0 tag
【Problem】 `git status` shows untracked `.github/`, `_review_recompute.py`, `deposited_neg_ref.json`, `docs/review/round6_2026-10-10/`, `err.txt`, so the working tree is not identical to the tagged release.
【Evidence】 `git status -s` output; `.github/workflows/release.yml` present but untracked (not in `v1.0.0`).
【Why it matters】 Not a credibility blocker — these are CI/review/working artifacts and the tagged tree is clean. But it means "the repository at this commit" ≠ "the v1.0.0 release contents," and a reviewer cloning the default branch would also pull in review scaffolding. Worth a one-line note in the deposit SOP.
【Specific fix】 Either commit `.github/workflows/release.yml` into the release (if the CI is intended to be part of the deposit) or add it to `.gitignore` and note in `GITHUB_DEPOSIT_SOP.md` that review/CI scaffolding is intentionally excluded from the tagged artefact set.

---

## 7b. Appendix — full row-by-row reconciliation

For transparency I list every manuscript table row against the recomputed artefact value (rounded to the manuscript's displayed precision). "match" means the artefact value rounds to the manuscript figure.

### Table 1 (all 22 D7 rows, source `module_level_signal.csv`)
| module | artefact t | ms t | artefact q | ms q | verdict |
|---|---|---|---|---|---|
| CTRL interferon | 24.923 | +24.92 | 6.76e-5 | 6.8e-5 | match |
| Glial/immune | 11.412 | +11.41 | 0.00267 | 0.0027 | match |
| Toll-like receptors | 6.364 | +6.36 | 0.0245 | 0.0245 | match |
| Cytokine/chemokine | 5.877 | +5.88 | 0.0129 | 0.0129 | match |
| Prostanoid/COX | 1.065 | +1.07 | 0.315 | 0.315 | match |
| Purinergic P2X | 0.692 | +0.69 | 0.272 | 0.272 | match |
| CTRL housekeeping | 0.653 | +0.65 | 0.524 | 0.524 | match |
| Cannabinoid | −0.646 | −0.65 | 0.934 | 0.934 | match |
| Monoaminergic | −1.079 | −1.08 | 0.961 | 0.961 | match |
| Injury/regeneration | −1.399 | −1.40 | 0.791 | 0.791 | match |
| TRP channels | −2.467 | −2.47 | 0.221 | 0.221 | match |
| Neuropeptides | −2.630 | −2.63 | 0.136 | 0.136 | match |
| Ca²⁺ α2δ | −3.509 | −3.51 | 0.136 | 0.136 | match |
| Opioid receptors | −3.761 | −3.76 | 0.0772 | 0.077 | match |
| CTRL myelin/Schwann | −4.101 | −4.10 | 0.0174 | 0.017 | match |
| Voltage-gated Na⁺ | −4.347 | −4.35 | 0.0211 | 0.021 | match |
| Glutamate/NMDA | −4.353 | −4.35 | 0.0601 | 0.060 | match |
| CTRL erythrocyte | −4.801 | −4.80 | 0.0186 | 0.019 | match |
| CTRL pan-neuronal | −5.201 | −5.20 | 3.14e-4 | 0.00031 | match |
| Neurotrophin | −5.292 | −5.29 | 0.0335 | 0.033 | match |
| K⁺/HCN | −5.635 | −5.64 | 0.00534 | 0.0053 | match |
| GABAergic | −5.708 | −5.71 | 0.0245 | 0.0245 | match |

### Table 2 (all 8 displayed D7-KO rows, source `virtual_perturbation_scores.csv`)
| module | rev | ms rev | z | ms z | %or | ms %or | q | ms q |
|---|---|---|---|---|---|---|---|---|
| CTRL interferon | 0.2282 | +0.228 | 4.611 | +4.61 | 75.43 | 75.4 | 0.0110 | 0.011 |
| Glial/immune | 0.1545 | +0.155 | 3.271 | +3.27 | 61.48 | 61.5 | 0.0110 | 0.011 |
| CTRL pan-neuronal | −0.1086 | −0.109 | −3.560 | −3.56 | 32.87 | 32.9 | 0.0110 | 0.011 |
| K⁺/HCN | −0.0916 | −0.092 | −2.867 | −2.87 | 36.46 | 36.5 | 0.0110 | 0.011 |
| Voltage-gated Na⁺ | −0.0903 | −0.090 | −2.626 | −2.63 | 35.93 | 35.9 | 0.0195 | 0.020 |
| CTRL myelin/Schwann | −0.0896 | −0.090 | −2.707 | −2.71 | 30.91 | 30.9 | 0.0195 | 0.020 |
| CTRL erythrocyte | −0.0744 | −0.074 | −2.568 | −2.57 | 32.98 | 33.0 | 0.0195 | 0.020 |
| Neuropeptides | −0.0675 | −0.067 | −2.449 | −2.45 | 26.85 | 26.9 | 0.0195 | 0.020 |

### Table 3 (all 9 rows, source `etiology_marker_cross.csv`)
| gene | Bone_D7 | Bone_D14 | PTX | OXA | Δ | ms Δ | q | ms q |
|---|---|---|---|---|---|---|---|---|
| Cd68 | 2.839 | 2.492 | −0.579 | 0.404 | 3.418 | +3.42 | 0.0125 | 0.012 |
| Aif1 | 2.192 | 1.954 | −0.883 | −0.459 | 3.076 | +3.08 | 0.0125 | 0.012 |
| Gal | −0.672 | 0.013 | 1.920 | −0.347 | −2.592 | −2.59 | 0.0125 | 0.012 |
| Ccl2 | 0.897 | 0.843 | −0.230 | −0.493 | 1.391 | +1.39 | 0.0173 | 0.017 |
| Il1b | 0.800 | 0.740 | −0.431 | −0.305 | 1.230 | +1.23 | 0.0182 | 0.018 |
| Calca | −0.272 | −0.182 | 0.498 | 0.041 | −0.770 | −0.77 | 0.0272 | 0.027 |
| Scn10a | −0.461 | −0.207 | 0.241 | −0.044 | −0.702 | −0.70 | 0.0260 | 0.026 |
| Atf3 | 0.313 | 0.317 | 1.122 | 0.109 | (ns) | — | 0.1197 | ns |
| Gfap | 1.498 | 0.594 | 0.086 | 0.049 | (ns) | — | 0.0978 | ns |

### Table 4 (all 11 rows, source `etiology_module_logfc_cross.csv`)
| module | D7 | D14 | PTX | OXA | Δ | ms Δ | q | ms q |
|---|---|---|---|---|---|---|---|---|
| Glial/immune | 2.053 | 1.302 | −0.643 | −0.132 | 2.695 | +2.70 | 0.00735 | 0.0073 |
| Cytokine/chemokine | 1.396 | 1.456 | −0.306 | −0.353 | 1.702 | +1.70 | 0.00735 | 0.0073 |
| Toll-like receptors | 0.854 | 0.685 | −0.421 | −0.136 | 1.275 | +1.28 | 0.00735 | 0.0073 |
| Neurotrophin | −0.566 | −0.292 | 0.405 | — | −0.972 | −0.97 | 0.00735 | 0.0073 |
| Ca²⁺ α2δ | −0.275 | −0.019 | 0.325 | — | −0.600 | −0.60 | 0.00735 | 0.0073 |
| Neuropeptides | −0.515 | −0.241 | 0.830 | — | −1.345 | −1.35 | 0.0145 | 0.0145 |
| CTRL pan-neuronal | −0.444 | −0.210 | 0.247 | — | −0.691 | −0.69 | 0.0145 | 0.0145 |
| Glutamate/NMDA | −0.437 | — | 0.198 | — | −0.635 | −0.64 | 0.0162 | 0.0162 |
| Voltage-gated Na⁺ | −0.492 | — | 0.162 | — | −0.654 | −0.65 | 0.0310 | 0.0310 |
| GABAergic | −0.430 | — | 0.065 | — | −0.495 | −0.49 | 0.0310 | 0.0310 |
| CTRL erythrocyte | −1.054 | — | — | 1.429 | −2.483 | −2.48 | 0.00735 | 0.0073 |

### Table 5 (single-cell, source json/csv as noted)
- Samples/cells: 4 tumour + 4 control; 55,967 — `sc_GSE292522_summary.json` `total_cells` = 55967 ✓.
- DE q<0.05: 1,517 up / 703 down — `pseudobulk_DE` = 1517/703 ✓.
- Sprr1a +1.55 (t 54.5) — logFC 1.554, t 54.49 ✓.
- Atf3 +0.84 — logFC 0.835 ✓.
- Ifit2 down q 0.042 / Oas1c q 0.014 — q 0.04165 / 0.01356 ✓.
- Gfap +0.004 q 0.006 — logFC 0.004, q 0.00615 ✓.
- Cd68/Aif1 unchanged — logFC −0.001/−0.001, q 0.783/0.920 ✓.
- Myeloid fraction 0.12%/0.05% — Con 0.001170, B16F10 0.000517 ✓.
- MWU U=2/14, p 0.114 (one-sided 0.057) — `myeloid_MWU_p_twosided` 0.11429, `_less` 0.05714 ✓.
- 8 testable compartments; erythrocyte 2 control / 4 tumour — composition counts ✓.
- ρ 0.098 (D7), 0.112 (D1), 0.030 (D14) — `convergent_validation_vs_bulk_D7.spearman_rho` 0.09845 ✓.

---

## 8. Stands up (what is solid)

1. **The reproducibility gate is real and effective.** `98_check_manuscript_numbers.py` (291 checks, 0 failed) and `99_unit_tests.py` (14 groups, all passed) both pass on this tree, and the gate recomputes table cells at run time rather than comparing to hardcoded literals — so a stale artefact would fail. This is stronger provenance practice than most computational submissions.
2. **Every table number I recomputed (T1, T2, T3, T4, T5 and the §3.6 anti-selectivity statistics) reconciles with the deposited artefacts to the rounding shown in the manuscript.** No stale or contradictory table value was found.
3. **The negative-result quantification is deterministic and self-consistent:** script 33's output is byte-equivalent in content to the deposited `deposited_neg_ref.json` (0/54 key diffs), and the MANIFEST hashes for 8 spot-checked files (including that JSON, README, CITATION, and three table CSVs) all match. The deposit is internally coherent.

---

## 9. Questions for the authors

1. Can you confirm the GitHub repository `yyx-4113/cancer-pain-drg-virtual-perturbation` is public and that the `v1.0.0` tag (not just `main`) resolves to the exact deposited tree? (Network was blocked here, so I could not verify liveness — F2.)
2. In §2.3 you attribute the §3.3 per-gene values to the "13,601 genes shared" muscle-adjusted fit, but the actual §3.3 numbers equal the 13,893-gene single-study adjusted fit (`etiology_marker_de.csv`). Was a separate 13,601-shared adjusted fit ever computed for §3.3, or should §2.3 be corrected to say §3.3 uses the single-study 13,893 fit (F1)?
3. The README calls the GSE292522 control "MPEK-BL6" while the manuscript says "keratinocyte." Are these the same cell preparation, and can the wording be unified (F3)?
4. Minor: the manuscript §3.6 states "32 of the 264 (or 32/132 deduplicated)" significant against the randset null, and also "64 of the 264 … 14 of those 32 were negative-control modules." The de-duplication arithmetic (264 → 132, and 32 controls+targets) is consistent in the artefacts; please confirm the "64 of the 264" figure is meant as the pre-deduplication KO+OE count (32 × 2) so a reader is not confused by the two denominators.

---

## 10. What I actually checked (commands, outputs, discrepancies)

**Commands run (from project root, managed interpreter):**
- `python scripts/98_check_manuscript_numbers.py` → `291 checks, 0 failed; ALL MANUSCRIPT NUMBERS VERIFIED AGAINST ARTEFACTS` (exit 0).
- `python scripts/99_unit_tests.py` → `ALL TESTS PASSED` (T1–T14, exit 0).
- `python scripts/33_negative_result_quantification.py` → exit 0; wrote `results/tables/negative_result_quantification.json`.
- Python one-liners to diff the regenerated JSON vs `deposited_neg_ref.json` → **0/54 key differences**.
- Python recompute of T1 from `module_level_signal.csv` (script 20:401) → all 22 D7 rows match manuscript Table 1.
- Python recompute of T2 from `virtual_perturbation_scores.csv` (D7 knockout) → top-8 rows + footer ρ=0.9911, mean|net|=0.057, mean|add|=0.033, mean|ridge|=0.712 all match.
- Python recompute of T3/T4 from `etiology_marker_cross.csv` and `etiology_module_logfc_cross.csv` → all rows match.
- Python recompute of T5 from `sc_GSE292522_summary.json`, `sc_GSE292522_pseudobulk_DE.csv`, `sc_GSE292522_composition.csv` → 1,517/703, Sprr1a +1.55/t54.5, myeloid 0.12%/0.05%, MWU p=0.114, ρ=0.098, 46 myeloid/1 T cell/55,967 total, erythrocyte 2/4 split all match.
- Python recompute of §3.6 anti-selectivity from `etiology_bootstrap_differences.csv` (132 dedup KO rows) → median |z| 2.47/1.07, MWU p=0.0019, 32 sig (14 ctrl/18 tgt), binomial p=0.0067 — all match.
- Two-logFC check: `etiology_marker_de.csv` Cd68 2.858 = `DE_D7_vs_Sham_adj.csv` (13,893, adj) 2.86199; `etiology_marker_cross.csv` Cd68 2.839 = `DE_D7_vs_Sham_unadj.csv` 2.83902 — confirms 2.86/2.84 is covariate, not gene-set.
- `git log --oneline -1`, `git tag`, `git status -s`, `git remote -v`, `git rev-parse HEAD`/`git rev-list -n 1 v1.0.0` → HEAD == v1.0.0 commit 8882852; untracked `.github/`, review files; remote configured but `git ls-remote` blocked by network.
- MANIFEST.sha256: recomputed SHA-256 of 8 files → all 8 matched.

**Discrepancies found:**
- F1: §2.3 gene-set attribution for §3.3 is inaccurate (says 13,601 shared; actual 13,893 single-study fit). Numbers themselves correct.
- F2: remote repo liveness unverified (network blocked).
- F3: GSE292522 control named "MPEK-BL6" (README) vs "keratinocyte" (manuscript).
- F4: untracked working-tree files not in the v1.0.0 tag (harmless).

**No numerical contradictions, no stale table values, no failed gate, no MANIFEST hash mismatches were found.** The implementation layer is, with the four minor exceptions above, reproducible and internally consistent.

### 10b. Verbatim gate / unit-test excerpts (evidence)

`98_check_manuscript_numbers.py` tail (selected lines):
```
  [PASS] Table 3 Ca2+ channel a2d PTX logFC: 0.3246
  [PASS] text: '4,600'
  [PASS] text: '0.991'
  [PASS] text: '+3.42'
  [PASS] text: '46 cells'
  [PASS] text: 'eight testable compartments'
========================================================================
291 checks, 0 failed
ALL MANUSCRIPT NUMBERS VERIFIED AGAINST ARTEFACTS
```

`99_unit_tests.py` tail (selected lines):
```
  [PASS] T10c the true assignment is one of the splits
  [PASS] T12b LOSO differs from in-sample (overfitting check)
  [PASS] T14e differs from the grand-mean / Frobenius implementation
======================================================================
ALL TESTS PASSED
```

`33_negative_result_quantification.py` diff vs `deposited_neg_ref.json` (my own comparison script, flattened to 54 leaf keys):
```
KEYS total fresh=54 deposited=54
DIFFS: 0
```

`git` state:
```
HEAD  = 8882852afc11da21f156136da000c9278f78d30b
v1.0.0= 8882852afc11da21f156136da000c9278f78d30b   (HEAD == tag commit)
git ls-remote origin -> fatal: remote helper 'https' aborted session  (network blocked)
git status -s -> ?? .github/  ?? _review_recompute.py  ?? deposited_neg_ref.json
                 ?? docs/review/round6_2026-10-10/  ?? err.txt
```

MANIFEST spot-check (8 of 190 entries):
```
OK   CITATION.cff
OK   README.md
OK   requirements.txt
OK   results/tables/negative_result_quantification.json
OK   results/tables/etiology_marker_cross.csv
OK   results/tables/virtual_perturbation_scores.csv
OK   results/tables/sc_GSE292522_summary.json
OK   docs/manuscript/manuscript_final.md
sum: ok=8 bad=0
```

These excerpts are the raw evidence underlying the verdicts above; they were produced on the project tree as delivered and are not edited.
