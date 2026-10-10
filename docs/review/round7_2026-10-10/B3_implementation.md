# Reviewer B3 — Implementation / Provenance Audit (fresh, independent)

**Manuscript:** `docs/manuscript/manuscript_final.md` (title: *Divergent dorsal root ganglion programs in bone cancer pain and chemotherapy-induced peripheral neuropathy: an exact cross-study reanalysis and a negative-control calibration of in-silico perturbation*)
**Role:** Implementation / provenance auditor. I treat this as a first submission and assume flaws are possible. I ran the gate and unit-test scripts myself, recomputed table values from the deposited CSV/JSON artefacts, independently resolved the three new DOIs on the live web, and inspected the git state.
**Bottom line:** The computational results are exceptionally well traced: the gate reports **291/0** and the unit suite reports **ALL TESTS PASSED**, and my independent recomputation of ≥15 cells across Tables 1–5 agrees with the manuscript to the printed precision. The three new references [17][18][19] are all real and correctly detailed. The one genuine provenance defect is that the **`v1.0.0` release tag cited in the Data-availability statement does not contain the submitted manuscript** (the working tree is ahead of the tag, and the *remote* `v1.0.0` tag points to a different, older commit). This must be fixed before the release is cited as the reproducible snapshot.

---

## 1. Reproducibility gate and unit tests (mandatory focus 1)

**Commands run (managed python, `PYTHONIOENCODING=utf-8`):**
- `python scripts/98_check_manuscript_numbers.py manuscript_final.md` → **291 checks, 0 failed**
- `python scripts/99_unit_tests.py` → **ALL TESTS PASSED**

**Gate output (excerpt, verbatim):**
```
[PASS] Table 3 Ca2+ channel a2d PTX logFC: 0.3246
[PASS] Table 3 CTRL erythrocyte Bone_D7 logFC: -1.054
[PASS] Table 3 CTRL erythrocyte OXA logFC: 1.429
... 68 table cells regenerated
[PASS] text: '4,600'  ... '0.991' ... '24.9' ... '6.26' ... '3/4,901' ...
[PASS] text: '55,967' ... '32/132' ... '0.0019' ... '46 cells' ... 'eight testable compartments'
[PASS] all references carry a DOI: 19
[PASS] figure file present (fig1.png / fig2.png / fig3.png): True
========================================================================
291 checks, 0 failed
ALL MANUSCRIPT NUMBERS VERIFIED AGAINST ARTEFACTS
```

**Unit-test output (excerpt, verbatim):** T1 (moderated t/logFC recovery), T2 (BH against hand-computed), T3/T4 (diffusion closed-form, max|P80−Pex|=2.22e-16), T5 (no self-loops/symmetry), T6 (KO≡−OE), T7 (identifier matrix rejected), T8 (D14≠D1 parsing), T9 (degree-matched bins), T10 (C(8,4)=70 splits), T11 (two-sided p), T12 (ridge LOSO), T14 (batched cosine 2-D), T13 (muscle score) — **all groups PASS, final line "ALL TESTS PASSED".**

【Problem】 The author reports 291/0 and all unit tests passing after revision; this needed independent confirmation.
【Evidence】 I executed both scripts with the managed interpreter; both exited 0. Gate: `291 checks, 0 failed`. Unit suite: every group T1–T14 reports `[PASS]` and ends `ALL TESTS PASSED`. The gate also regenerated 68 table cells at run time rather than comparing to hardcoded literals (its own design), so the pass is meaningful, not circular.
【Why it matters】 If the gate or unit suite had failed, every downstream number in the manuscript would be unverified and the central "negative-control-calibrated" claim would be unsubstantiated. A green gate is the single most important credibility signal for a computation-heavy paper.
【Specific fix】 None required. Recommend the author keep the gate as the acceptance criterion and re-run it in CI on any future change (see §6 — the current gate is *not* wired into the untracked `.github/workflows/release.yml`, and the released tag predates the current manuscript, so the published release would not re-run it on the exact submitted text).

---

## 2. Table traceability — independent spot-checks (mandatory focus 2)

I recomputed values directly from the deposited artefacts, not from the gate. The two logFC sources the manuscript distinguishes are `etiology_marker_de.csv` (muscle-adjusted limma-trend) and `etiology_marker_cross.csv` (unadjusted per-sample mean), plus `module_level_signal.csv`, `etiology_module_logfc_cross.csv`, and the `sc_GSE292522_*` JSON/CSV.

**Table 1** (`module_level_signal.csv`, D7 rows):
- CTRL interferon: `mean_t=24.923`→+24.9; `mean_logFC=2.614`→+2.61; `q=6.76e-05`→6.8×10⁻⁵. ✓
- Glial/immune: `mean_t=11.412`→+11.4; `mean_logFC=2.012`→+2.01; `q=0.002667`→0.0027. ✓
- CTRL erythrocyte: `mean_t=-4.801`→−4.80; `mean_logFC=-1.147`→−1.15; `q=0.0186`→0.019. ✓
- K⁺/HCN: `mean_t=-5.635`→−5.64; `q=0.005339`→0.0053. ✓
- **Table 1 footnote (new):** glial/immune "six genes (Gfap, Aif1, Itgam, Cd68, S100a8, S100a9); one gene did not pass the expression filter … leaving five effective genes (n=5)." The D7 `genes_used` string is `Gfap;Aif1;Itgam;Cd68;S100a9` = exactly 5 genes; the omitted sixth is S100a8. The ISG module `n=8` lists `Irf1;Ifitm3;Ifi47;Bst2;Rsad2;Gbp2;Ubd;Ly6e` = 8 genes. ✓ **The new footnote is accurate.**

**Table 3** (`etiology_marker_cross.csv`): Cd68 Bone_D7=2.839→+2.84, PTX=−0.579→−0.58, Δ=3.418→+3.42, q=0.0125→0.012 ✓; Aif1 Bone_D7=2.192→+2.19, PTX=−0.883→−0.88, Δ=3.076→+3.08 ✓; Gal Bone_D7=−0.672→−0.67, PTX=+1.920→+1.92, Δ=−2.592→−2.59 ✓; Ccl2 D7/OXA Δ=1.391→+1.39, q=0.0173→0.017 ✓; Il1b Δ=1.230→+1.23, q=0.0182→0.018 ✓; Calca Δ=−0.770→−0.77, q=0.0272→0.027 ✓; Scn10a Δ=−0.702→−0.70, q=0.0260→0.026 ✓; Atf3 Bone_D7=0.313→+0.31, PTX=1.122→+1.12, OXA=0.109→+0.11 ✓; Gfap Bone_D7=1.498→+1.50, Bone_D14=0.594→+0.59, PTX=0.086→+0.09, OXA=0.049→+0.05 ✓. All 9 rows recomputed and matched.

**Table 4** (`etiology_module_logfc_cross.csv`): Glial/immune D7/PTX `logFC_A=2.053`→+2.05, `logFC_B=−0.643`→−0.64, Δ=2.695→+2.70, q=0.00735→0.0073 ✓; Cytokine/chemokine D7/PTX +1.40/−0.31/Δ+1.70 ✓; Toll-like D7/PTX +0.85/Δ+1.28 ✓; Neurotrophin D7/PTX −0.57/+0.41/Δ−0.97, q=0.0073 ✓; Ca²⁺ α2δ D7/PTX −0.28/+0.32/Δ−0.60 ✓; Neuropeptides D7/PTX −0.52/+0.83/Δ−1.35, q=0.0145 ✓; CTRL pan-neuronal D7/PTX −0.44/+0.25/Δ−0.69 ✓; Glutamate/NMDA D7/PTX −0.44/+0.20/Δ−0.64, q=0.0162 ✓; Voltage-gated Na⁺ D7/PTX −0.49/+0.16/Δ−0.65, q=0.0310 ✓; GABAergic D7/PTX −0.43/+0.06/Δ−0.49, q=0.0310 ✓; CTRL erythrocyte D7/OXA −1.054/+1.429/Δ−2.483→−2.48, q=0.0073 ✓. All 11 rows matched.

**Table 5** (JSON/CSV below, focus 7): 55,967 cells; 1,517 up / 703 down; myeloid fraction 0.12 % / 0.05 %; MWU p=0.114; concordance ρ 0.098/0.112/0.030; Sprr1a +1.55 (t=54.5); Atf3 +0.84; Gfap +0.004 (q=0.006); Ifit2 q=0.042, Oas1c q=0.014; galanin neurons +1.47 (q<0.001) and NF-neuron +2.54 (q<0.001). All matched.

【Problem】 Whether Tables 1–5 numbers actually trace to the deposited CSV/JSON, including the new Table 1 glial/immune n=5 footnote.
【Evidence】 I independently recomputed ≥15 distinct cells across Tables 1–5 from `module_level_signal.csv`, `etiology_marker_cross.csv`, `etiology_module_logfc_cross.csv`, `sc_GSE292522_summary.json`, `sc_GSE292522_pseudobulk_DE.csv`, `sc_GSE292522_celltype_DE.csv`, and `sc_GSE292522_concordance_per_day.json`. Every value matched the printed figure to the stated rounding. The Table 1 footnote n=5 is confirmed by the 5-gene `genes_used` string; the ISG n=8 by the 8-gene list.
【Why it matters】 Table traceability is the core of provenance. A single mismatched cell would mean the manuscript was hand-edited away from the artefacts. No mismatch was found — this is a strong positive signal.
【Specific fix】 None required. The footnote is accurate; keep it.

---

## 3. New references [17][18][19] (mandatory focus 3)

I resolved each DOI on the live web (publisher pages / PubMed / CrossRef mirrors) and checked authors, journal, year, volume/issue/pages.

- **[17]** `doi:10.1007/s00011-022-01663-2` → *Inflammation Research* 2023;72(1):117–132. Title *"Activation of the STING pathway induces peripheral sensitization via neuroinflammation in a rat model of bone cancer pain."* Authors (manuscript: Zhang Y, Wang W, Gong Z, Peng Y, Li X, Zhang Z, Zhang X, You X, Wu J) match the publisher/PubMed author list (Yuxin Zhang, Wei Wang, Zhihao Gong, Yuan Peng, Xin Li, Zuojing Zhang, Xiaoxuan Zhang, Xingji You, Jingxiang Wu) one-for-one. ✓ Real, correctly detailed.
- **[18]** `doi:10.1073/pnas.2503779122` → *PNAS* 2025;**122(34):e2503779122**. Title *"Unraveling the neuroimmune mechanisms in cancer-induced bone pain: new horizons for therapeutic intervention of the two-phase paradigm."* Authors (manuscript: Wang C, Ji H, Wang H, Chen Z, Zhou L, Yang Y, Jiang Y, Yu G, Jiang L, Tang Z) match the PNAS/CrossRef author list (Changming Wang, Haiwang Ji, Hanwen Wang, Ziyun Chen, Lan Zhou, Yan Yang, Yucui Jiang, Guang Yu, Ling Jiang, Zongxiang Tang). The article number `e2503779122` is the DOI suffix itself — internally consistent. ✓ Real, correctly detailed.
- **[19]** `doi:10.1016/j.brainres.2007.06.066` → *Brain Research* 2007;**1168:46–59**. Title *"An evolving cellular pathology occurs in dorsal root ganglia, peripheral nerve and spinal cord following intravenous administration of paclitaxel in the rat."* Authors (manuscript: Peters CM, Jimenez-Andrade JM, Kuskowski MA, Ghilardi JR, Mantyh PW) match PubMed (Christopher M. Peters, Juan Miguel Jimenez-Andrade, Michael A. Kuskowski, Joseph R. Ghilardi, Patrick W. Mantyh). The manuscript's use of [19] to support "Atf3 … elevated in paclitaxel-treated DRG neurons" is faithful to the source (the paper's abstract reports ATF3 up-regulation in DRG neurons after paclitaxel). ✓ Real, correctly detailed.

【Problem】 Whether the three added references are real and the DOIs resolve to the papers as described, with correct citation metadata.
【Evidence】 All three DOIs resolved on the live web. [17]→*Inflammation Research* 72(1):117–132 (2023); [18]→*PNAS* 122(34):e2503779122 (2025); [19]→*Brain Research* 1168:46–59 (2007). Author lists, journals, years, volumes/issues/pages, and titles all match the manuscript verbatim. No fabricated or mismatched DOI.
【Why it matters】 Fabricated or mismatched DOIs are a hard credibility/integrity failure that can trigger retraction. Here there is no such risk; the new citations are genuine and support the claims they are attached to (STING/innate-immune programme [17]; macrophage infiltration in CIBP [18]; ATF3 as a paclitaxel nerve-injury marker [19]).
【Specific fix】 None required.

---

## 4. Two-logFC coherence (mandatory focus 4)

The manuscript keeps two logFC sources apart: §3.3 uses the muscle-adjusted limma-trend *moderated* logFC from `etiology_marker_de.csv`; Table 3 uses the *unadjusted per-sample mean* logFC from `etiology_marker_cross.csv`. §2.3 (manuscript lines ~59, ~115) states the expected pairs and attributes residual gaps to moderated shrinkage.

From `etiology_marker_de.csv` (§3.3 source, Bone_D7):
- Cd68 `2.858`→**2.86**; Aif1 `2.180`→**2.18**; Gfap `1.328`→**1.33**; Ccl2 `0.806`→**0.81**.

From `etiology_marker_cross.csv` (Table 3 source, Bone_D7):
- Cd68 `2.839`→**2.84**; Aif1 `2.192`→**2.19**; Gfap `1.498`→**1.50**; Ccl2 `0.897`→**0.90**.

These reproduce the manuscript's stated pairs exactly: **Cd68 2.86 vs 2.84; Aif1 2.18 vs 2.19; Gfap 1.33 vs 1.50; Ccl2 0.81 vs 0.90**. The larger gap for the lower-abundance markers (Gfap, Ccl2) is consistent with the manuscript's explanation that "moderated shrinkage pulls low-abundance genes toward zero." §2.3 explicitly accounts for the divergence and states the muscle-adjusted value is quoted as primary when they disagree.

【Problem】 Whether §3.3 (muscle-adjusted moderated logFC) and Table 3 (unadjusted per-sample mean) are internally coherent and whether §2.3 explains the gap.
【Evidence】 Recomputed both sources from `etiology_marker_de.csv` and `etiology_marker_cross.csv`. Cd68 2.86/2.84, Aif1 2.18/2.19, Gfap 1.33/1.50, Ccl2 0.81/0.90 all match the printed numbers; the direction and magnitude of the residual (larger for lower-abundance genes) is exactly as §2.3 predicts. The §2.3 clarification (manuscript lines 59 and ~115) is present and correct.
【Why it matters】 Conflating the two estimators would look like a contradiction (e.g., "Gfap +1.33" in §3.3 but "+1.50" in Table 3). The manuscript avoids this and explains it; a reader who misses §2.3 could otherwise suspect an error.
【Specific fix】 None required. The clarification is adequate; I would only suggest adding one cross-reference sentence in §3.3 itself (e.g., "see §2.3 for the two-logFC convention") so the two sections are explicitly linked at the point of reading.

---

## 5. Stale-number / cross-file drift (mandatory focus 5)

I grepped `README.md`, `CITATION.cff`, `requirements.txt`, and the manuscript for version strings, accession numbers, n-values, and key statistics.

- **GSE292522 control wording.** README line 91 now reads *"B16F10 melanoma vs non-tumorigenic keratinocytes (MPEK-BL6), 14 d"* (the uncommitted diff changed it from "vs MPEK-BL6 control" to "vs non-tumorigenic keratinocytes (MPEK-BL6)"). The manuscript §2.8 (line 95) reads *"B16F10 melanoma or of non-tumorigenic keratinocytes, the control."* The two are **consistent** — README now spells out the cell line (MPEK-BL6) that the manuscript leaves implicit. No contradiction.
- **Accessions.** GSE149648, GSE160543, GSE292522, GSE249443 appear identically in README and manuscript; n=4 per group, n=4+4 single-cell stated consistently.
- **Key statistics.** "55,967 cells", "0.991", "264 / 132", "0.0019", "6.26 %", "3,041", "Jaccard 0.584" appear consistently across README caveats and manuscript.
- **Release URL.** Manuscript Data-availability: *"https://github.com/yyx-4113/cancer-pain-drg-virtual-perturbation (release v1.0.0)"*. README §Caveats: *"The manuscript's Data-and-code-availability statement points to this repository (release v1.0.0)."* Consistent in wording. (Whether v1.0.0 *contains* this manuscript is the separate git finding in §6.)

【Problem】 Whether version strings, accessions, n-values, and headline statistics agree across README, CITATION.cff, requirements.txt, and the manuscript.
【Evidence】 Compared the four files. The GSE292522 control now reads "non-tumorigenic keratinocytes (MPEK-BL6)" in README and "non-tumorigenic keratinocytes, the control" in the manuscript — consistent. Accessions, group sizes, and the cited statistics match. No stale or conflicting number was found.
【Why it matters】 Cross-file drift (e.g., an old n or an old accession lingering in README) is a classic sign of an under-maintained repository and would undermine the reproducibility claim.
【Specific fix】 None required for content. One housekeeping note: `CITATION.cff` still says `journal: "Manuscript in preparation"` and `year: 2026`; once the manuscript is accepted this should be updated to the real citation, but that is expected post-acceptance, not a defect now.

---

## 6. Git / repository state — the one genuine defect (mandatory focus 6)

**Commands run:**
```
git -C . log --oneline -1   → 8882852 docs: drop Zenodo from data-availability statement and repo docs (Zenodo archiving paused per author)
git -C . tag                → v1.0.0
git -C . status -s         → M README.md | M docs/manuscript/manuscript_final.md | M scripts/98_check_manuscript_numbers.py
                                ?? .github/ | ?? docs/review/round6_2026-10-10/ | ?? docs/review/round7_2026-10-10/ | ... (tmp files, logs)
git -C . ls-remote <remote> → SUCCESS: HEAD/main = 1ff1fd6… ; refs/tags/v1.0.0 = 8d34f89…
git -C . rev-list -n 1 v1.0.0 (local) → 8882852…
git -C . diff --stat v1.0.0      → README.md | 2 +- ; manuscript_final.md | 33 +- ; scripts/98 | 2 +-
```

**Findings:**
1. **The repo is live** — `git ls-remote` succeeded; the remote `v1.0.0` tag exists and `main` exists. Data-availability liveness is therefore **confirmed**, not unverified.
2. **But the remote `v1.0.0` tag (8d34f89…) is a DIFFERENT commit from the local `v1.0.0` tag / HEAD (8882852…).** The published release therefore does **not** contain the current HEAD.
3. **The working tree is ahead of the local `v1.0.0` tag.** Three tracked files are modified but *uncommitted*: `README.md`, `docs/manuscript/manuscript_final.md`, and `scripts/98_check_manuscript_numbers.py`. The `manuscript_final.md` diff is **substantive**, not cosmetic: it adds the §2.3 two-logFC clarification (the exact Cd68 2.86-vs-2.84 etc. sentence), re-titles §3.3, inserts the new [17][18] citations supporting the myeloid/ISG reading, and adds the [19] citation for Atf3 — i.e., precisely the content added "after revisions."
4. **`.github/` is untracked** (contains only `workflows/release.yml`). It is not part of the `v1.0.0` snapshot and does not affect its correctness; it is simply a CI helper that was never committed.

**Consequence:** A reader who opens the cited "release v1.0.0" gets commit `8d34f89` (remote) — which is *older* than even the local tagged `8882852`, which in turn is *older* than the submitted `manuscript_final.md` in the working tree. The released snapshot therefore does **not** contain the manuscript (and the §2.3 clarification, the Table 1 footnote, the new references) that this review is assessing. The gate/unit tests I ran are green only because I ran them on the working tree; the published release would run them on a different, earlier manuscript and code.

【Problem】 The `v1.0.0` release cited in the Data-availability statement does not contain the submitted manuscript; the working tree is ahead of the tag and the remote tag points to a different, older commit.
【Evidence】 `git tag` = v1.0.0 (local = 8882852 = HEAD). `git status -s` shows `M manuscript_final.md` (among others) — i.e., the submitted manuscript is an *uncommitted* modification on top of the tagged commit. `git ls-remote` shows the *remote* `v1.0.0` = `8d34f89…`, a different SHA than the local `8882852…`. `git diff --stat v1.0.0` shows `manuscript_final.md | 33 +-` of substantive edits (§2.3 clarification, §3.3 re-title, [17][18][19] insertions). So the cited release ≠ submitted manuscript.
【Why it matters】 The entire reproducibility claim rests on the released snapshot matching the paper. If a reviewer/reader clones "release v1.0.0" they cannot reproduce the numbers tied to the *current* manuscript, and the gate they run there would be checking a different text. This is a real provenance gap that must be closed before the release URL is trusted as the canonical artefact.
【Specific fix】 Before citing "release v1.0.0", (a) `git add` and `git commit` the three modified files (README.md, docs/manuscript/manuscript_final.md, scripts/98_check_manuscript_numbers.py) so the working tree is clean and equals a commit; (b) move the `v1.0.0` tag (or cut a new tag, e.g. `v1.0.0-final`) onto that final commit **locally and on GitHub** (`git tag -f v1.0.0 <final> && git push -f origin v1.0.0`); (c) verify by `git ls-remote` that the remote `v1.0.0` SHA now equals the submitted-HEAD SHA. Paste-ready sentence for the Data-availability section: *"Analysis code, artefacts and the test suite are available at https://github.com/yyx-4113/cancer-pain-drg-virtual-perturbation (release v1.0.0, tag `v1.0.0` at commit <final-SHA>), with requirements.txt pinning every version used."*

---

## 7. Single-cell QC (mandatory focus 7)

From `sc_GSE292522_composition.csv` and `results/review_a3/compartment_cellcounts.json` (per-sample cell counts):

- **46 myeloid cells total.** Per-sample myeloid counts: 3, 10, 12, 5 (Con) and 5, 3, 7, 1 (B16F10) → sum = **46**. ✓
- **1 T cell total.** Only Sample7_DRG_B16F10_rep3 has 1 T cell; all others 0 → **1**. ✓
- **55,967 total cells** (`sc_GSE292522_summary.json` `total_cells: 55967`). ✓
- **8 testable compartments.** The 12 candidate compartments minus 4 excluded: T cells (1 cell), NP neurons (max 1), TH neurons (max 3), and myeloid (≥5 cells in only 5 of 8 animals: samples 2,3,4,5,7). 12−4 = **8** (Neuron, SatelliteGlia, Schwann, Endothelial, Fibroblast, Neuron.NF, Neuron.PE, Erythrocyte). ✓
- **Erythrocyte 2-control / 4-tumour split.** Erythrocyte per-sample counts: 5, 4, 22, 3 (Con) and 7, 6, 26, 15 (B16F10). The ≥5-cell rule is met by 6 samples: Sample1(5)+Sample3(22) = **2 Con**, and Sample5(7)+Sample6(6)+Sample7(26)+Sample8(15) = **4 B16F10**. This exactly matches the manuscript's "two control and four tumour animals" flag. ✓
- **Myeloid fraction test.** Summary JSON: Con 0.0011699 (~0.12 %), B16F10 0.0005175 (~0.05 %); MWU two-sided p=0.1142857 ("U = 2/14" = the two U statistics 2 and 14, since n₁n₂=16), one-sided p=0.0571429. ✓ Matches Table 5.
- **Per-compartment ISG.** `sc_GSE292522_celltype_DE.csv`: Neuron,Ifit2 q=0.0459→0.046 ✓; SatelliteGlia,Isg15 q=0.0071→0.007 ✓; SatelliteGlia,Ifit1 q=0.0037→0.004 ✓. "none up in any" is supported by `ISG_up_in_B16F10: []` in the summary. ✓
- **Neuronal-injury expression (no compartment information).** Neuron,Gal logFC=1.472, q=0.0001 (<0.001) ✓; Neuron.NF,Gal logFC=2.537→+2.54, q=0.0 (<0.001) ✓; Neuron,Scn10a q=0.0117→0.012 ✓. All match Table 5.
- **Concordance.** `sc_GSE292522_concordance_per_day.json`: D1 ρ=0.1115→0.112, D7 ρ=0.0985→0.098, D14 ρ=0.0299→0.030. ✓ Matches Table 5.
- **Pseudobulk DE.** Summary: 1,517 up / 703 down ✓. `sc_GSE292522_pseudobulk_DE.csv`: Sprr1a logFC=1.5539→+1.55, t=54.49→54.5 ✓; Atf3 logFC=0.8355→+0.84 ✓; Gfap logFC=0.00358→+0.004, q=0.00615→0.006 ✓; Ifit2 q=0.0417→0.042 ✓; Oas1c q=0.0136→0.014 ✓.

【Problem】 Whether the single-cell composition numbers (46 myeloid, 1 T cell, 55,967 total, 8 testable compartments, erythrocyte 2-control/4-tumour split) and the Table 5 DE/ISG/concordance values are internally consistent with the deposited artefacts.
【Evidence】 Recomputed myeloid total (3+10+12+5+5+3+7+1 = 46) and T-cell total (1) from per-sample counts; confirmed 8 testable compartments and the erythrocyte 2-Con/4-B16F10 split from `compartment_cellcounts.json`; verified every Table 5 value (cells, DE counts, fractions, MWU p, per-compartment ISG q's, galanin neurons/NF, Sprr1a/Atf3/Gfap logFC, Ifit2/Oas1c q, concordance ρ) against the JSON/CSV. All matched.
【Why it matters】 The single-cell arm is the paper's "model-specific" negative-control result; if its QC counts were wrong the whole §3.8 conclusion would collapse. They are correct and the blind-spot (CD45 depletion → myeloid untestable) is honestly stated.
【Specific fix】 None required. The erythrocyte imbalance is already flagged in text; keep it.

---

## 8. Version consistency (mandatory focus 8)

§2.9 states: *"Python 3.13 with numpy 2.4.6, scipy 1.18.1, pandas 2.3.3, statsmodels 0.15.0 and matplotlib 3.11.1 (`requirements.txt`)."*

- **requirements.txt:** `numpy==2.4.6`, `scipy==1.18.1`, `pandas==2.3.3`, `statsmodels==0.15.0`, `matplotlib==3.11.1`. ✓ exact match.
- **README.md (line 114):** `numpy==2.4.6   scipy==1.18.1   pandas==2.3.3   statsmodels==0.15.0   matplotlib==3.11.1`. ✓ exact match; README §Software also states "Python 3.13".
- **CITATION.cff:** does not list the numeric stack (it records repo `version: 1.0.0`, `date-released: 2026-10-10`, and a preferred-citation block). The numeric versions are therefore not duplicated there — this is **not a mismatch**, just an absence; only §2.9, requirements.txt, and README carry them, and those three agree exactly.

Empirically, the managed interpreter executed both `scripts/98_*.py` and `scripts/99_*.py` to green, which demonstrates the pinned environment is internally coherent and runnable.

【Problem】 Whether the §2.9 software versions exactly match requirements.txt, README.md and CITATION.cff.
【Evidence】 §2.9 (numpy 2.4.6, scipy 1.18.1, pandas 2.3.3, statsmodels 0.15.0, matplotlib 3.11.1, Python 3.13) is character-for-character identical to requirements.txt and to README line 114; README §Software also states Python 3.13. CITATION.cff omits the numeric stack (records only repo `version: 1.0.0`) — an absence, not a conflict. Both gate and unit scripts ran successfully under the managed interpreter.
【Why it matters】 Version drift between the methods section and the pinned environment is a common reproducibility failure; here there is none among the three sources that carry versions.
【Specific fix】 None required for consistency. Optional: add the same five version lines to `CITATION.cff` (e.g., a `dependencies`/Notes entry) so the machine-readable citation carries the runtime stack too. (Note: the version numbers themselves — e.g. numpy 2.4.6, scipy 1.18.1 — are unusual/future-dated relative to public release trains; this is outside the consistency check and I did not treat it as a defect since the managed environment executed both scripts successfully. If these are aspirational pins, the author should confirm the environment actually resolves them.)

---

## Findings that need action (summary)

1. **[HIGH — provenance]** The cited `release v1.0.0` does not contain the submitted manuscript (§6). Fix by committing the three modified files and re-tagging `v1.0.0` at the final commit both locally and on GitHub; verify via `git ls-remote`.
2. **[LOW — optional]** Add a one-line cross-reference in §3.3 to §2.3 for the two-logFC convention (§4).
3. **[LOW — optional]** Mirror the five pinned versions into `CITATION.cff` (§8).
4. **[LOW — post-acceptance]** Update `CITATION.cff` `journal`/`year` from "Manuscript in preparation"/2026 once published (§5).

Everything else — gate (291/0), unit tests (all pass), Tables 1–5 traceability, the three new DOIs, two-logFC coherence, cross-file drift, and single-cell QC — verified clean.

---

## § Stands up

1. **The number gate is real and meaningful.** `scripts/98_check_manuscript_numbers.py` regenerates 68 table cells and re-derives key perturbed quantities at run time (its own design note says expectations are *not* hardcoded), then asserts them against the artefacts; it reported **291 checks, 0 failed**. This is far stronger than a typical "we re-ran the script" claim because the artefacts, not the manuscript, are the source of truth.
2. **Tables 1–5 trace exactly to deposited CSV/JSON.** I independently recomputed ≥15 cells (Cd68, Aif1, Gal, Ccl2, Il1b, Calca, Scn10a, Atf3, Gfap across Tables 3/4; interferon/glial/erythrocyte/Toll-like across Table 1; every Table 5 QC/DE/concordance value) and every one matched to the printed rounding, including the new Table 1 glial/immune `n=5` footnote.
3. **The three new references [17][18][19] are genuine and correctly cited.** All three DOIs resolve on the live web to the exact papers, with author lists, journals, years, and volume/issue/pages matching the manuscript verbatim — no fabricated or mismatched DOI, which is the single biggest credibility risk for added citations.
4. **The negative-control / single-cell conclusions are internally honest.** The erythrocyte 2-control/4-tumour imbalance, the CD45-depletion myeloid blind spot, the Gfap-absent satellite-glia label, and the "model-specific" interferon conclusion are all stated with their limiting QC facts, and the underlying counts (46 myeloid, 1 T cell, 55,967 total, 8 testable compartments) check out against the composition artefacts.

---

## § Questions for the authors

1. The remote `v1.0.0` tag currently points to commit `8d34f89…`, which is **not** the submitted HEAD `8882852…`, and your working tree has uncommitted edits to `manuscript_final.md` (the §2.3 clarification, §3.3 re-title, and [17][18][19] citations). **When you cut "release v1.0.0", will it point at the exact submitted manuscript, or at an earlier draft?** Please confirm the final tag SHA equals the submitted-HEAD SHA via `git ls-remote`.
2. §2.9 pins `numpy 2.4.6`, `scipy 1.18.1`, `pandas 2.3.3`, `statsmodels 0.15.0`, `matplotlib 3.11.1`, `Python 3.13`. These version numbers are ahead of the public release trains I would expect. Are these the versions the managed environment actually resolved, or are they aspirational pins? (The gate and unit tests did run green, so the environment is at least self-consistent.)
3. In Table 3, Gfap shows +1.50 (unadjusted cross) but §3.3 quotes +1.33 (muscle-adjusted moderated). You explain this in §2.3 as shrinkage, which is correct — but for the *within-study* §3.3 programme statements, is the muscle-adjusted value always the one intended as "primary" when the two disagree (as §2.3 now says), and is that consistently what Tables 1 and 3 report? I found no contradiction, but I want to confirm the convention is applied uniformly.
4. The erythrocyte compartment is tested with a 2-control / 4-tumour split (Sample1+Sample3 Con; Sample5–8 B16F10, per the ≥5-cell rule). Was any within-compartment test or visualisation produced for that imbalanced design, and is the imbalance explicitly carried into the per-compartment ISG statement in §3.8 (which currently reads "none up in any")?
5. `CITATION.cff` still lists `journal: "Manuscript in preparation"` and `year: 2026`. Is this intentional pending acceptance, or should it already reflect a preprint/submitting journal?

---

## § What I actually checked (commands, outputs, discrepancies)

**Environment:** managed interpreter `C:/Users/Administrator/.workbuddy/binaries/python/envs/default/Scripts/python.exe`, `PYTHONIOENCODING=utf-8`, run from project root `D:/2026.10/不同部位癌痛患者对药物敏感性生信分析/cancer_pain_virtual_perturbation`.

**Scripts executed:**
- `python scripts/98_check_manuscript_numbers.py manuscript_final.md` → exit 0; tail: `291 checks, 0 failed  /  ALL MANUSCRIPT NUMBERS VERIFIED AGAINST ARTEFACTS`. Regenerated 68 table cells; all `[PASS]`.
- `python scripts/99_unit_tests.py` → exit 0; all groups T1–T14 `[PASS]`; final `ALL TESTS PASSED`.

**Web verification of DOIs (live):**
- [17] `10.1007/s00011-022-01663-2` → *Inflammation Research* 72(1):117–132 (2023); authors match.
- [18] `10.1073/pnas.2503779122` → *PNAS* 122(34):e2503779122 (2025); authors/title match.
- [19] `10.1016/j.brainres.2007.06.066` → *Brain Research* 1168:46–59 (2007); authors match.

**Table recomputation (independent of the gate):** read `results/tables/module_level_signal.csv`, `etiology_marker_cross.csv`, `etiology_module_logfc_cross.csv`, `etiology_marker_de.csv`, `sc_GSE292522_composition.csv`, `sc_GSE292522_summary.json`, `sc_GSE292522_celltype_summary.json`, `sc_GSE292522_pseudobulk_DE.csv`, `sc_GSE292522_celltype_DE.csv`, `sc_GSE292522_concordance_per_day.json`, and `results/review_a3/compartment_cellcounts.json`. Reconciled ≥15 cells across Tables 1–5 and all Table-5 QC/DE/concordance values — **zero discrepancies**.

**Git state:**
- `git log --oneline -1` → `8882852 docs: drop Zenodo from data-availability statement and repo docs (Zenodo archiving paused per author)`.
- `git tag` → `v1.0.0`.
- `git status -s` → `M README.md`, `M docs/manuscript/manuscript_final.md`, `M scripts/98_check_manuscript_numbers.py`; untracked `.github/`, `docs/review/round6_2026-10-10/`, `docs/review/round7_2026-10-10/`, `tmp_recompute_*.py`, `err.txt`, `_review_recompute.py`, `deposited_neg_ref.json`.
- `git ls-remote` → **succeeded** (repo live); remote `HEAD`/`main` = `1ff1fd6…`; remote `refs/tags/v1.0.0` = `8d34f89…`.
- Local `v1.0.0` = `8882852…` (= HEAD). Remote `v1.0.0` = `8d34f89…` (**different**). `git diff --stat v1.0.0` → `README.md | 2 +-`, `manuscript_final.md | 33 +-`, `scripts/98 | 2 +-`. The `manuscript_final.md` diff is substantive (§2.3 two-logFC clarification, §3.3 re-title, [17][18][19] insertions).
- **Discrepancy found:** the cited `release v1.0.0` does not contain the submitted manuscript (working tree ahead of tag; remote tag is an older, different commit). This is the only substantive defect. `.github/` (untracked, `workflows/release.yml` only) does **not** affect the v1.0.0 snapshot integrity.

**Cross-file drift checks:** README line 91 "non-tumorigenic keratinocytes (MPEK-BL6)" vs manuscript §2.8 "non-tumorigenic keratinocytes, the control" — consistent. Accessions, n-values, and headline statistics agree across README/CITATION/requirements/manuscript. No stale numbers found.

**Version consistency:** §2.9 == requirements.txt == README line 114 exactly (numpy 2.4.6, scipy 1.18.1, pandas 2.3.3, statsmodels 0.15.0, matplotlib 3.11.1, Python 3.13). CITATION.cff omits the numeric stack (records repo version 1.0.0) — absence, not conflict.

**Overall verdict:** computationally sound and exceptionally well traced; the only required fix before the release URL can be trusted as the canonical reproducible snapshot is to commit the three modified files and re-tag `v1.0.0` so it equals the submitted HEAD (§6).
