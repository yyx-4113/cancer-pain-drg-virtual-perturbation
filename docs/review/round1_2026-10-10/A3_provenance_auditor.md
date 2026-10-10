# A3 — Provenance / Reproducibility Audit

**Manuscript:** `docs/manuscript/manuscript_final.md`
**Auditor role:** Independent provenance auditor, round-1 panel
**Artefact root:** `results/tables/`
**Independence statement:** I did **not** read any other reviewer output (`A1_*`, `A2_*`, `A4_*`, `REVIEW_*/RESPONSE_*/REVISION_*`, `_PANEL_BRIEF.md`) or any `docs/review/round1_2026-10-10/` file other than this one. Every judgement below was derived by me from the manuscript and the deposited CSV/JSON artefacts, recomputed with Python 3.13.12. I did read the deposited scripts only as code (not as reviewer commentary) to confirm artefact structure — never another auditor's text.

---

## Summary verdict

The headline numbers in the Abstract, §3.2, §3.3, §3.4, §3.6 and §3.8 **do** trace exactly to the deposited artefacts, and the gate's arithmetic is sound. **However**, the supporting *tables* (Table 2 and Table 3) and several *prose* claims contain numbers that do **not** trace to any deposited artefact, and there is one genuine cross-section contradiction (§6 over-claims the interferon programme in a way §3.8 explicitly retracts). The gate cannot catch these because it can only validate within-table arithmetic and the Δ/q columns; the per-contrast logFC columns in Table 2/3 and the un-deposited D1/D14 concordance are exactly the class of number it does not bind.

Severity ranking: **F1, F2, F8** = high (tables contradict the artefacts / conclusion drift); **F3, F4, F5** = medium; **F6, F7, F9** = low/cosmetic.

---

## Stands up (verified against deposited artefacts)

These were recomputed from the raw artefacts and match the manuscript **exactly** (or to the stated precision):

- **4,600 / 3,041 D7 genes.** `run_summary.json` → `n_sig_q05_adj.D7 = 4600`; `de_muscle_drop_sensitivity.csv` → `n_q05_reduced = 3041`, `jaccard_q05 = 0.583955` (→0.584), `spearman_logFC = 0.912307` (→0.912). ✓
- **0.165 minimum q.** `etiology_summary.json` → `min_q_permutation_per_contrast = 0.165447…`; `etiology_bootstrap_differences.csv` (per_contrast KO) min `q_permutation = 0.165447`. ✓
- **32/132 randset-significant with 14 controls.** `etiology_summary.json` → `n_sig_randset_q05_percontrast_indep = 32`, `randset_sig_breakdown_percontrast = {target:18, control:14}`. ✓
- **55,967 cells; 1,517 up / 703 down.** `sc_GSE292522_summary.json` → `total_cells = 55967`, `pseudobulk_DE.n_up_B16F10_q05 = 1517`, `n_down = 703`. ✓ (also confirmed by recounting `sc_GSE292522_pseudobulk_DE.csv`: 1,517 / 703).
- **Sprr1a +1.55 (t = 54.5).** `sc_GSE292522_pseudobulk_DE.csv` → logFC = 1.5539 (→1.55), t = 54.493 (→54.5). ✓
- **Myeloid fraction Con 0.00117 vs B16F10 0.000517, MWU p = 0.971.** `sc_GSE292522_summary.json` → `myeloid_fraction_mean.Con = 0.0011699` (→0.117%→"0.12%"), `B16F10 = 0.0005175` (→0.0517%→"0.05%"), `myeloid_MWU_p_greater = 0.9714` (→0.97). The "0.12%"/"0.05%" rounding of 0.00117/0.000517 is acceptable double-rounding. ✓
- **Cd68 +2.86, Aif1 +2.18, Gfap +1.33 at D7.** `etiology_marker_de.csv` (Bone_D7) → 2.858/2.180/1.328 → +2.86/+2.18/+1.33. ✓
- **Galanin Δ = −2.59; Cd68 Δ = +3.42.** `etiology_marker_cross.csv` → Gal maxΔ = −2.592, Cd68 maxΔ = 3.418. ✓
- **Table 1 (all 22 rows) matches `module_level_signal.csv`** (D7) to the digit: e.g. interferon t = 24.923 (→+24.92), q = 6.759e−5 (→6.8×10⁻⁵); pan-neuronal t = −5.201 (→−5.20), q = 3.138e−4 (→0.00031); erythrocyte logFC = −1.147 (→−1.15). ✓
- **Head-to-head / calibration.** `head_to_head_baselines.json` → ρ = 0.9911 (→0.991), n = 132, mean |net| = 0.0568 (→0.057), mean |add| = 0.0328 (→0.033), mean %oracle = 20.51 (→20.5). ✓
- **Sensitivity.** `sensitivity_summary.json` → discovery min ρ = 0.9593 (→0.959), median 0.9876 (→0.988), Δ min 0.9887 (→0.989), 75 grid points (60+15). ✓ `sensitivity_limma_trend.csv` → `spearman_t_vs_reference` = 1.000/0.999996/0.999999 (≥0.9999). ✓
- **QC.** `qc_muscle_scores.csv` → max `muscle_fpkm_pct = 0.06257` (→6.26%, sample GSM4506604), next = 0.001604 (→0.16%); `qc_pca_variance.csv` → PC1–4 = 46.85/17.45/7.23/4.76 (→46.8/17.4/7.2/4.8); sample-to-sample correlation off-diagonal range = 0.913–0.994. ✓
- **ISG magnitude.** Of the 37 ISGs in `sc_GSE292522_focus_markers.csv`, max |logFC| = 0.105 (< 0.11) — so "flat in magnitude" holds. (But see **F5** for the q claim.)

Abstract vs §3.6 consistency check (focus task 2, part 1): the Abstract's "minimum q = 0.165" and "32/132 … 14 negative controls" match §3.6 line-for-line and match the artefact. **No drift there.**

---

## F1 — Table 1 and Table 3 report the same quantity (module mean logFC at D7) with different numbers, and Table 3 disagrees with the artefact

【Problem】 Table 3's "logFC bone D7" column does not equal Table 1's "mean logFC" column for the same modules/contrast, although both are defined in §2.7 as "mean logFC (case versus control) within each study"; Table 3 also disagrees with the deposited `module_level_signal.csv` and `etiology_module_logfc_cross.csv`.

【Evidence】
- `module_level_signal.csv` (D7) and `etiology_module_logfc_cross.csv` (logFC_A for Bone_D7) **agree** on, e.g., Glial/immune = 2.012, Cytokine/chemokine = 1.382, CTRL erythrocyte = −1.147.
- Table 1 (manuscript_final.md:202–223) reproduces these (Glial/immune +2.01, Cytokine/chemokine +1.38, CTRL erythrocyte −1.15) — i.e. Table 1 = artefact.
- Table 3 (manuscript_final.md:241–253) reports the **same** Glial/immune = **+2.05**, Cytokine/chemokine = **+1.40**, CTRL erythrocyte = **−1.05** — none of which match the artefact.
- Full per-module comparison (Table 3 "logFC bone D7" vs artefact):
  - Glial/immune 2.05 vs **2.012**; Cytokine/chemokine 1.40 vs **1.382**; TLR 0.85 vs **0.870**; Neurotrophin −0.57 vs **−0.519**; Ca²⁺α2δ −0.28 vs **−0.237**; Neuropeptides −0.52 vs −0.524 (only match); CTRL pan-neuronal −0.44 vs **−0.399**; Glutamate/NMDA −0.44 vs **−0.376**; Voltage-gated Na⁺ −0.49 vs **−0.424**; GABAergic −0.43 vs **−0.373**; CTRL erythrocyte −1.05 vs **−1.147**.
  - 10 of 11 modules disagree; offsets 0.02–0.10.
- The Δ and q columns of Table 3 **do** match the cross-aetiology artefact, so only the per-contrast logFC columns are wrong/mis-sourced.

【Why it matters】 A reader citing a Table 3 per-contrast module logFC (e.g. "CTRL erythrocyte −1.05 at D7") is citing a number that traces to no deposited artefact, and the manuscript contradicts itself between two adjacent tables. This is precisely the class the gate misses (it checks Δ = D7−PTX inside Table 3's own cells and the q, not cross-table or against the module-level signal).

【Specific fix】 Regenerate Table 3's per-contrast logFC columns directly from `module_level_signal.csv` (bone D7/D14) and `etiology_module_logfc_cross.csv` (PTX/OXA as logFC_A/logFC_B), which will make Table 3 consistent with Table 1 and the artefact. At minimum, set Glial/immune D7 = +2.01 (not +2.05), CTRL erythrocyte D7 = −1.15 (not −1.05), etc.

---

## F2 — Table 2 per-contrast marker logFCs do not trace to `etiology_marker_de.csv`

【Problem】 Several per-contrast marker logFCs in Table 2 disagree with the deposited `etiology_marker_de.csv` (the only source for the CIPN PTX/OXA contrasts, and the bone D7/D14 contrasts corroborated by §3.3). The Δ and q columns in Table 2 **do** match `etiology_marker_cross.csv`, so the fault is isolated to the per-contrast logFC cells.

【Evidence】 Recomputed from `etiology_marker_de.csv`:
- **Ccl2 D7: +0.90 (Table 2, line 232) vs 0.806 (artefact); §3.3 line 111 states "Ccl2 +0.81", corroborating the artefact, not Table 2.** Clear internal contradiction.
- **Gal D14: +0.01 (Table 2, line 231) vs 0.191 (artefact).** Off by 0.18.
- **Il1b OXA: −0.31 (Table 2, line 234) vs −0.160 (artefact).** Off by 0.15.
- **Aif1 OXA: −0.46 (Table 2, line 230) vs −0.291 (artefact).** Off by 0.17.
- **Calca D7: −0.27 vs −0.212; Scn10a D7: −0.46 vs −0.399; Scn10a D14: −0.21 vs −0.180** — all off ~0.05–0.06.
- Confirmed NOT explained by the raw `DE_D14_vs_Sham_adj.csv`: e.g. Gal D14 raw = 0.236, still ≠ +0.01. The Table 2 D14/OXA cells match neither the marker-DE artefact nor the raw DE file consistently.

【Why it matters】 Table 2 is the canonical "canonical-marker cross-aetiology" table. Cells that cannot be reproduced from any deposited artefact undermine the exact-permutation claim the whole section rests on, and one of them (Ccl2 D7) is directly contradicted by the manuscript's own §3.3 prose.

【Specific fix】 Rebuild Table 2's four per-contrast columns from `etiology_marker_de.csv` (Bone_D7/Bone_D14/PTX/OXA). Correct at least: Ccl2 D7 → +0.81, Gal D14 → +0.19, Il1b OXA → −0.16, Aif1 OXA → −0.29. Then re-verify the "Largest Δ"/q columns still hold (they should, since they come from `etiology_marker_cross.csv`).

---

## F3 — §3.5 cytokine/chemokine Δ "against oxaliplatin" is arithmetically wrong

【Problem】 §3.5 (line 129) states "cytokine/chemokine Δ = +1.81 against oxaliplatin and +1.70 against paclitaxel", but the deposited cross-aetiology artefact gives Δ(D7 vs OXA) = +1.75, not +1.81.

【Evidence】 `etiology_module_logfc_cross.csv` → Cytokine/chemokine Bone_D7|OXA: `delta_logFC = 1.749` (q = 0.0073). The paclitaxel figure +1.70 matches the artefact (Δ(D7 vs PTX) = 1.702 → 1.70). The oxaliplatin figure +1.81 does not. (The error looks like the D14 logFC 1.462 was used instead of D7 1.382: 1.462 − (−0.353) = 1.815 ≈ 1.81.)

【Why it matters】 It is a single wrong number in a biological conclusion sentence that the Discussion leans on ("differ in direction-specific ways … large effect sizes"). Internally it also conflicts with Table 3, which (correctly per artefact) shows Cytokine/chemokine D7 = +1.38 and OXA = −0.35.

【Specific fix】 Change "cytokine/chemokine Δ = +1.81 against oxaliplatin" to "**+1.75** against oxaliplatin" (matching `etiology_module_logfc_cross.csv` Δ = 1.749).

---

## F4 — D1 and D14 convergent-validation correlations are reported but not deposited in any artefact

【Problem】 §3.8 (line 145) and Table 5 (line 280) report Walker-256↔B16F10 pseudobulk concordance ρ = 0.112 (D1) and 0.030 (D14), but no deposited CSV/JSON contains these values.

【Evidence】 The only convergent-validation object in the artefacts is `sc_GSE292522_summary.json → convergent_validation_vs_bulk_D7 = {spearman_rho: 0.09845, …}` (the D7 value, which the manuscript also reports as 0.098 — ✓). A full-text search of `results/tables/` for "concordance", "bulk_D1", "bulk_D14" and for the literal strings "0.112"/"0.030" (excluding coincidental substrings inside other numbers) finds **no** per-day concordance artefact. The D1/D14 values are therefore computed somewhere (likely the validation script) but never saved as a numbered, checkable artefact.

【Why it matters】 §2.9 explicitly promises "A separate checker asserts every number reported below against the deposited CSV and JSON artefacts." Two reported numbers violate that guarantee. This is the exact gap the gate cannot see — the numbers are internally consistent (they appear in both §3.8 and Table 5) but untethered to a deposited file.

【Specific fix】 Deposit a small artefact, e.g. `sc_GSE292522_concordance.json` = `{D1:0.112, D7:0.098, D14:0.030, n_shared_genes:…, spearman_p:…}`, and point the §2.9 checker at it.

---

## F5 — "every ISG q > 0.15" is false (Ifit2 q = 0.042)

【Problem】 §3.8 (line 145) claims "Every interferon-stimulated gene tested was flat (|logFC| < 0.11, **q > 0.15**)". The |logFC| < 0.11 part holds, but the q > 0.15 part does not.

【Evidence】 `sc_GSE292522_focus_markers.csv` (37 ISGs): max |logFC| = 0.105 (<0.11, OK). But **Ifit2 has q = 0.042** and Sp100 has q = 0.150, i.e. Ifit2 fails the stated q > 0.15 threshold. (Ifit2's magnitude is indeed flat, logFC = −0.048, so the biological point stands; only the universal "q > 0.15" qualifier is wrong.)

【Why it matters】 A categorical "every ISG q > 0.15" claim is falsifiable and false; a careful reviewer or replicator will flag it. It also slightly undercuts the strength of the "interferon programme is absent" conclusion, which should rest on effect size, not on an over-broad q-statement.

【Specific fix】 Soften to "Every interferon-stimulated gene tested was flat in magnitude (|logFC| < 0.11; all q ≥ 0.04, most q > 0.15)" or simply drop the q-threshold clause and keep "|logFC| < 0.11".

---

## F6 — §3.6 "Sixty-five of the 264" does not match the deposited count (68)

【Problem】 §3.6 (line 133) says "Sixty-five of the 264 (or 32/132 deduplicated) tests were 'significant'", but `etiology_summary.json` reports `n_sig_randset_q05_percontrast = 68` (raw) with `n_sig_randset_q05_percontrast_indep = 32` (dedup). The dedup 32 matches; the raw 65 does not match 68.

【Evidence】 `etiology_summary.json`: `n_sig_randset_q05_percontrast = 68`, `n_sig_randset_q05_percontrast_indep = 32`. Direct count of `etiology_bootstrap_differences.csv` (per_contrast, knockout, `q_randset<0.05`) = 32 per direction; the JSON's raw total is 68.

【Why it matters】 Low severity, but it is an arithmetic/provenance mismatch the gate should have caught if it checked the raw count. Either the manuscript rounded/filtered differently or the artefact was patched after the prose was written.

【Specific fix】 Either change "Sixty-five" → "**Sixty-eight**" (68) to match the artefact, or document the exact denominator/filter that yields 65. The dedup figure 32/132 is correct as written.

---

## F7 — "erythrocyte … single largest standardised effect (z = −5.78)" is not the largest |z|

【Problem】 §3.6 (line 133) calls the erythrocyte set "the single largest standardised effect (z = −5.78 for paclitaxel versus oxaliplatin)", but the erythrocyte module attains larger |z| in two other pairs.

【Evidence】 `etiology_bootstrap_differences.csv` (per_contrast, knockout, `z_vs_randset`): erythrocyte Bone_D7|OXA = **−6.169**, erythrocyte Bone_D14|Bone_D7 = **+7.388**. Both exceed |−5.78| of PTX|OXA. (The −5.78 PTX|OXA value itself is correct and matches the artefact.)

【Why it matters】 Cosmetic, but "single largest" is a factual over-statement; the erythrocyte's largest |z| is actually the within-bone D14-vs-D7 pair (7.39) or the Bone_D7|OXA pair (6.17). Reader will notice if they sort the table.

【Specific fix】 Rephrase to "the erythrocyte set — a technical control with no pain biology — produced the largest standardised effect among the negative controls (e.g. z = −5.78 for paclitaxel vs oxaliplatin, and even larger in other pairs)" or simply drop "single largest".

---

## F8 — §6 over-claims the interferon programme as pan-bone-cancer, contradicting §3.8's model-specificity qualifier

【Problem】 §6 (line 169) concludes the two aetiologies drive "qualitatively different programmes … **myeloid and interferon-driven in the first** [bone cancer pain]", with no qualification, whereas §3.8 (line 147) and the Abstract (line 23) explicitly state the interferon programme is **model-specific** to the Walker-256 tibia model and is **not** a universal cancer-pain signature.

【Evidence】
- §3.8 line 147: "the interferon / myeloid programme is a feature of the **Walker-256 tibia model at its day-7 peak**, not a universal signature of cancer-pain DRG."
- Abstract line 23: "An independent single-cell model shows the interferon programme is **model-specific rather than a universal cancer-pain signature**."
- §4 line 153 correctly qualifies: "**In the Walker-256 tibia model**, bone cancer pain at its transcriptional peak is an innate-immune event…".
- §6 line 169 drops the qualifier entirely: "myeloid and **interferon-driven in the first**" (i.e. bone cancer pain as a disease).

【Why it matters】 This is a genuine cross-section contradiction. §3.8 is the paper's own methodological punchline — that a single-model signature must not be generalised. §6's flat "bone cancer pain is interferon-driven" reverts exactly the over-generalisation §3.8 warns against, and a reader who lands on the Conclusions will take away the opposite of the qualified body.

【Specific fix】 Align §6 with §3.8/Abstract: e.g. "…myeloid- and interferon-driven in the **Walker-256 tibia model** (model-specific, not a universal cancer-pain signature), and neuronal-injury-driven in CIPN…".

---

## F9 — D14 q < 0.10 DE count (860) is in the artefact but omitted from §3.2

【Problem】 §3.2 (line 107) reports q < 0.05 as "0, 4,600 and 235" for days 1/7/14 and then "(5 and 5,867 at q < 0.10)" — only two of the three days are given for q < 0.10, so D14 q < 0.10 = 860 (deposited) is silently dropped.

【Evidence】 `run_summary.json` → `n_sig_q10_adj = {D1:5, D7:5867, D14:860}`. §3.2 gives "5 and 5,867 at q < 0.10" (D1 and D7) and never states D14 = 860.

【Why it matters】 Low severity, but it is an artefact number not reflected in the manuscript (the inverse of the usual problem) and makes the q < 0.10 sentence ambiguous (a reader cannot tell whether D14 was 860 or omitted deliberately).

【Specific fix】 Either add "and 860 at day 14" to the q < 0.10 clause, or explicitly state "D14 not shown".

---

## Double-rounding / precision assessment (focus task 4)

- **Myeloid %:** 0.00117 → "0.12%" and 0.000517 → "0.05%" (Table 5) and "~0.1%" (§3.8) are acceptable one-significant-figure renderings of 0.117%/0.0517%. No error.
- **Galanin Δ:** Table 2 −2.59 vs recomputed −2.592 — correct rounding.
- **ISG claim (F5)** is *not* a rounding artefact: Ifit2 q = 0.042 is genuinely below the stated 0.15 threshold; the claim is substantively wrong, not double-rounded.
- **Min q:** 0.165447 → "0.165" is correct rounding (not 0.17).
- No other double-rounding hazards found in the verified set.

## Internal-contradiction hunt (focus task 2) — full result

- **Abstract ↔ §3.6:** MATCH. "minimum q = 0.165" and "32/132 … 14 negative controls" appear identically in both and match `etiology_summary.json`. No drift.
- **§4 ↔ §3.8:** PARTIAL MATCH with one real drift — see **F8** (§6 un-qualifies the interferon programme that §3.8/Abstract qualify as model-specific). §4's own prose (line 153) is correctly qualified, so the drift is specifically between §4/§3.8/Abstract and §6.
- **Table 1 ↔ Table 3:** CONTRADICTION — see **F1** (same module-D7 logFC, different numbers; Table 3 wrong).
- **§3.5 ↔ Table 2:** CONTRADICTION on Ccl2 D7 (+0.81 in §3.3 prose vs +0.90 in Table 2) — see **F2**.
- **§3.5 ↔ artefact:** cytokine/chemokine OXA Δ wrong — see **F3**.

---

## Questions for the authors

1. What is the provenance of Table 2's D14 and OXA per-contrast marker columns and Table 3's per-contrast module columns? They do not match `etiology_marker_de.csv` / `module_level_signal.csv` / `etiology_module_logfc_cross.csv`. Were they hand-entered from an earlier (pre-patch) artefact version? (See `22b_patch_artifact.py` in the script dir — was an artefact regenerated after the tables were filled?)
2. Where are the D1 (ρ = 0.112) and D14 (ρ = 0.030) concordance values saved? They are reported in §3.8 and Table 5 but absent from every deposited CSV/JSON. Will you deposit them (F4)?
3. For F6: should the raw randset-significant count be 65 or 68? Which denominator/filter produces 65?
4. The erythrocyte randset z reaches +7.39 (Bone_D14|Bone_D7) and −6.17 (Bone_D7|OXA); do you want to keep the "single largest" wording (F7)?
5. Does the Ifit2 q = 0.042 ISG alter your "absent interferon programme" conclusion, or is effect-size (|logFC|<0.11) sufficient on its own (F5)?

---

## What I actually checked

**Files read (manuscript + artefacts only):**
- `docs/manuscript/manuscript_final.md` (full, 281 lines).
- `results/tables/run_summary.json`, `head_to_head_baselines.json`, `etiology_summary.json`, `sc_GSE292522_summary.json`, `sensitivity_summary.json`.
- `results/tables/module_level_signal.csv` (66 rows, D7 subset), `virtual_perturbation_scores.csv` (132 rows), `etiology_marker_cross.csv` (102 rows), `etiology_module_logfc_cross.csv` (132 rows), `etiology_bootstrap_differences.csv` (528 rows), `etiology_marker_de.csv` (68 rows), `qc_muscle_scores.csv`, `qc_pca_variance.csv`, `qc_sample_correlation.csv`, `de_muscle_drop_sensitivity.csv`, `sc_GSE292522_pseudobulk_DE.csv` (34,603 rows), `sc_GSE292522_composition.csv`, `sc_GSE292522_focus_markers.csv` (54 rows), `sensitivity_limma_trend.csv`, `sensitivity_grid.csv`.
- Raw DE files `DE_D7_vs_Sham_adj.csv`, `DE_D14_vs_Sham_adj.csv` (cross-check only).
- I did **not** open any `docs/review/round1_2026-10-10/` file other than this one, nor any `A1/A2/A4/*`, `REVIEW/RESPONSE/REVISION`, `_PANEL_BRIEF.md`, manifest, SOP, or author statement.

**Values recomputed from the artefacts (Python 3.13.12):**
- Table 1 (22 module rows) — all match `module_level_signal.csv` to <0.5% (F1 used this to prove Table 1 = artefact).
- Table 2 per-contrast marker logFCs and Δ/q — matched Δ/q to `etiology_marker_cross.csv`; flagged 6 per-contrast cells that don't match `etiology_marker_de.csv` (F2).
- Table 3 per-contrast module logFCs and Δ/q — matched Δ/q to `etiology_module_logfc_cross.csv`; flagged 10/11 per-contrast D7 cells that contradict `module_level_signal.csv` (F1).
- Table 4 — interferon KO rev 0.2282 / z 4.611 / 75.4% / p 0.0020 / q 0.011, glial/immune 0.1545 / 3.271 / 61.5% — all match `virtual_perturbation_scores.csv`.
- §3.6 — min q 0.165447 (matches), erythrocyte PTX|OXA z −5.782 (matches), pan-neuronal Bone_D7|PTX z −3.736 (matches the −3.74 claim — my first extraction had grabbed the wrong pair, PTX|OXA +4.85), erythrocyte max |z| = 7.388 (F7), raw randset count 68 vs ms 65 (F6).
- §3.8 — 55,967 cells, 1,517/703, Sprr1a 1.5539/54.49, Atf3 0.8355, myeloid 0.00117/0.000517/p 0.9714, ISG max |logFC| 0.105 but Ifit2 q 0.042 (F5), Cd68/Aif1 in sc unchanged (logFC ≈ 0, q > 0.78), Gal logFC 1.425 (neuronal-injury, as expected).
- QC — muscle 6.257%/0.160%, PCA 46.85/17.45/7.23/4.76, sample corr 0.913–0.994 (off-diagonal), DE-drop 3,041/Jaccard 0.584/ρ 0.912.
- Concordance — only D7 (0.09845) present in `sc_GSE292522_summary.json`; D1/D14 (0.112/0.030) absent from all artefacts (F4).
- Markdown table syntax — programmatic pipe/column-count check across Tables 1–5: **no unbalanced rows**.

**Reused (already correct) artefacts verified:** `head_to_head_baselines.json` (ρ 0.991, 0.057/0.033, 20.5%), `sensitivity_summary.json` (0.959/0.988/0.989, 75 points), `sensitivity_limma_trend.csv` (ρ ≥ 0.9999).
