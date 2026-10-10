# REVIEW Round 2 — Independent Panel Consolidated Report

**Date:** 2026-10-10 · **Manuscript:** `docs/manuscript/manuscript_final.md` (~7,440 words)
**Panel:** 5 independent experts · **Editor/consolidator:** WorkBuddy agent (adjudicated, not averaged)

---

## 1. Independence statement

Five reviewers (A1 domain neuroimmunology, A2 biostatistics, A3 single-cell methodology, A4 provenance audit, A5 venue/reporting) were briefed under a first-submission posture with an explicit forbidden-file list (`REVIEW_round1_*`, `R3_editorial_review.md`, all `manuscript_v0.*`, each other's outputs, `.workbuddy/`). Each wrote 600–1,500 lines with mandatory four-part findings and recomputed values.

**Evidence independence worked:** three reviewers with no knowledge of each other hit the same defect from three different angles — **the myeloid compartment was dropped by the cell-count filter and the paper's compartment-null claim silently excludes the one compartment where an interferon programme would reside** (A1 §3.1 with per-animal cell counts 3/10/12/5/5/3/7/1; A3 Finding 8 with the same table; A4 F5 by reconstructing counts from deposited artefacts). A second three-way convergence: the omission of `Ifit1` in satellite glia (q = 0.0037, the *most* significant compartment ISG result) from §3.8's "the only significant ISG changes" list (A1 §3.4, A2 F10, A5 Finding 6 evidence). A third: the lowess-sweep gene-universe mismatch (A2 F9 ≡ A4 F8).

## 2. Verdict table

| Expert | Layer | Verdict | One-line rationale |
|---|---|---|---|
| A1 | Domain | **Major revision** | Compartment-null excludes myeloid; galanin not neuron-restricted; conclusion contradicted by the source publication's own abstract (Ifnar1 causal result) |
| A2 | Statistics | **Major revision** | p-floor claims wrong ×3; shared-control null wrong; **the round-2 power analysis tests a comparison that does not exist and its conclusion is backwards**; T1/T2 stated opposite to the data |
| A3 | sc-methods | **Major revision** | Compartments are not cell types (Neuron bin 11.9% pure; SatelliteGlia bin contains no Gfap); compartment logFC ≈ whole-DRG logFC (no localisation); MWU results are at the n=4 resolution floor |
| A4 | Provenance | **Major revision** | 16 findings incl. one hand-edited table cell that the regeneration script cannot reproduce, and §3.8 numbers (ρ=0.112/0.030) with no computable source |
| A5 | Venue | **Major revision (2 desk-reject items as-is)** | Repo URL is a live 404; five mandatory front-matter statements absent; abstract 628 words |

## 3. Cross-verification (editor-replicated rows marked ✅)

| # | Location | Manuscript says | Independently recomputed | Who | Verdict |
|---|---|---|---|---|---|
| 1 | §2.7/§5 | p floor = 1/4,901 and 1/71 | min p in `etiology_marker_cross.csv` = **0.000612 = 3/4,901**; shared-control floor **3/71 = 0.04225** | A2; **✅ editor** | **WRONG ×3** |
| 2 | §3.6 (T1) | targets "statistically indistinguishable" from controls | dedup two-sided \|z\| MWU **p = 0.0019**, controls *larger* (median 2.47 vs 1.07) | A2; **✅ editor** | **CONTRARY TO DATA** |
| 3 | §3.6 (T2) | randset null "flags targets and controls equally" | 28/64 = 43.8% controls vs 22.7% base rate → binomial **p = 0.00024**: anti-selective | A2; **✅ editor** | **UNDERSTATED** |
| 4 | §3.6/§5/§6 (T5) | reversal MDE d = 5.65; "only powered to detect >0.42 units"; largest observed ≈0.27 | design premise false (KO/OE r = −1.000, single scalar per contrast); correct MDE (bootstrap SE 0.0185 × d 2.38) = **0.044 units**; **26/132** exact p < 0.05 (min 0.00143) — impossible if the layer could only see 0.42 | A2; **✅ editor** (0.3073/0.2743; SE 0.0185; 26/132) | **BACKWARDS** |
| 5 | §3.8 | "the only significant ISG changes were Ifit2 in neurons and Isg15 in satellite glia" | **4** ISG rows at q<0.05; omitted **Ifit1 SG q=0.0037** (largest) and **Neuron.NF Ifit2 q=0.041** | A1, A2, A5; **✅ editor** | **INCOMPLETE** |
| 6 | §3.8 | Gal localises to neurons | whole-DRG +1.425 vs Neuron +1.472 (**ratio 1.03** = no enrichment); Gal also up in endothelium (q=1e-4) and Schwann (q=0.006); Gal is itself a Neuron marker set member (circular) | A3, A1; **✅ editor** | **UNSUPPORTED** |
| 7 | §3.8/Table 5 | Sprr1a/Atf3 injury programme "neuron-intrinsic" | Sprr1a, Atf3, Gadd45a **absent from the compartment DE table** (never tested) | A1; **✅ editor** | **NOT TESTED** |
| 8 | §3.8/§5 | "no compartment…shows ISG up-regulation"/"every DRG compartment" | **Myeloid dropped** (46 cells total; 5/8 animals pass ≥5-cell rule, needs 6); TCell/NP/TH also dropped | A1, A3, A4 (3×) | **SILENTLY UNTESTED** |
| 9 | Table 3 | erythrocyte row Δ −2.48/q 0.0073 in "D7 vs PTX" column | that cell's D7-vs-PTX values are Δ −0.258, q **0.814**; printed values are the D7-vs-OXA test; regeneration script emits −0.26/0.814 → cell is a **hand-edit** | A4 | **PROVENANCE BREAK** |
| 10 | §3.8 | ρ = 0.112 (D1), 0.030 (D14) | no computable source in any deposited artefact | A4 | **UNSOURCED** |
| 11 | §2.8 | "55,967 cells **after quality filtering**" | 55,967 = sum of raw barcodes; no filtering code exists | A3 | **FALSE STATEMENT** |
| 12 | Front matter | — | repo URL returns **404**; no funding/COI/contributions/AI-use statements; abstract **628 words**; 8 keywords; running head 67 chars | A5 | **DESK-REJECT ×2** |

## 4. Graded consolidated issue list

### Tier 0 — conclusion-invalidating (must fix before any submission)

- **T0-1 · Compartment layer fails as framed** (A1 §3.1–3.7, A3 F1–F13, A4 F5, three-way convergence). Five independent failures: (a) myeloid compartment (the IFN site) dropped at 46 cells/5-of-8 animals; (b) "SatelliteGlia" bin contains no satellite glia (Gfap raw counts 0,0,0,3,0,2,2,0; top gene Plp1) so its ISG results cannot be attributed; (c) "Neuron" bin is 11.9% pure by an independent marker criterion; (d) compartment logFC ≈ whole-DRG logFC (Gal ratio 1.03; Neuron slope 0.855, ρ 0.904) → **no localisation information**; (e) all compartment MWU results sit at the n=4 rank-test floor (2/70 = 0.0286) and none survive leave-one-out (Isg15 SG: 1/8; Ifit1 SG: 0/8) while neuronal Gal/Scn10a survive 8/8. **Fix:** withdraw the compartment-null and localisation claims as currently written; either re-engineer (cluster-based annotation, disjoint markers, detection floor ≥10 counts, denominator column, LOO disclosure, myeloid reported descriptively as untestable) or reduce §3.8 to the whole-DRG pseudobulk result, which is clean.
- **T0-2 · The round-2 power analysis (T5) is wrong and its conclusion is backwards** (A2 F3; editor-verified). There is no 2-vs-2 comparison — reversal is one scalar per contrast, KO/OE are sign-flip mirror images (r = −1.000). The dispersion used (between-module SD 0.074) is not the SE of Δ (bootstrap median 0.0185). Corrected MDE ≈ 0.044 units; largest target Δ 0.274 = 6× the threshold; 26/132 exact p<0.05. **The layer is adequately powered; its failure is non-separation (controls outrank targets), which is a different — and stronger — criticism.** Delete the 5.65 sentence from §3.6, §5, §6; replace with A2's paste-ready passage; fix `T5` in `33_negative_result_quantification.py` (use bootstrap SE, not between-module SD).
- **T0-3 · "Indistinguishable"/"non-selective" are the opposite of what the data show** (A2 F4/F5; editor-verified). Controls have *larger* perturbation effects than targets (two-sided p=0.0019; module-label permutation p=0.033) and the randset null flags controls at ~2× base rate (p=0.00024) — an **anti-selective** null. This *strengthens* the paper's thesis but the current sentences assert the weaker, false claim. Also fix the one-sided `alternative="greater"` orientation in `33_*.py`.
- **T0-4 · Compartment-null vs source publication** (A1 §3.5). GSE292522's own paper (Baruch et al., ref [15]) reports that injured neurons in this exact model autonomously initiate a type-I IFN response, causally reversed by Ifnar1 knockout. The manuscript cannot present this dataset as a clean IFN refutation without a reconciliation subsection stating what the source concluded and why the reanalysis differs (compartment coverage, panel overlap = 1/8 genes, dynamic range).
- **T0-5 · Desk-reject items** (A5 F1–F3). Live 404 repo link + "archived upon acceptance" promise; missing Funding/COI/Author contributions/AI-use/Acknowledgements; abstract 628 words, 8 keywords, 67-char running head.

### Tier 1 — analyses/numbers to fix

- **T1-1 · p-floor statements** (A2 F1; ✅editor): replace 1/4,901→**3/4,901 = 6.12×10⁻⁴** and 1/71→**3/71 = 0.0423**, in §2.7, §5, code comment, **and the 98-gate enforced literal "1/4,901"** — the gate currently protects the wrong number (the exact failure mode documented 2026-09-27).
- **T1-2 · Shared-control permutation null** (A2 F2): index-paired null has ~2/3 the spread of the exchangeable-8-case null (Cd68: p 0.155 → 0.465). No current conclusion flips, but the Methods text presents it as a refinement. Enumerate the 8 case labels instead, or disclose.
- **T1-3 · Single-cell "well powered" argument is circular** (A2 F7): median-SE yardstick dominated by zero-variance genes; detected 2,220 genes at median |logFC| 0.045 — 40× the claimed MDE. State Cd68/Aif1 intervals directly (−0.0005 q=0.78; +0.0009 q=0.92). Myeloid d should use conventional pooled SD (1.27 not 1.37) or better, report log-scale ratio (−1.18); drop the post-hoc one-sided p.
- **T1-4 · "55,967 after quality filtering"** (A3 F5): false; reword or implement filtering.
- **T1-5 · A4 F3 Table-3 hand-edited cell**; **A4 F9 unsourced ρ=0.112/0.030** (add per-day concordance computation to `31_process_sc.py` or delete the numbers); **A4 F2 logFC source mixing** (Gfap +1.33 vs +1.50 for the same contrast) — add one Methods sentence declaring de.csv vs cross.csv per section.
- **T1-6 · §2.2 sham framing** (A1 §2.2): if sham is day-14 PBS, "D1 = 0 genes" is a D1-vs-D14 time-matched result presented as a baseline — verify against the series matrix and reword; **A1 §2.4**: oxaliplatin shows no neuronal-injury programme — the PTX/OXA lumping needs qualification.

### Tier 2 — wording/reporting

T2-1 A1 §1.2 Gfap three values (§3.3 +1.33 de.csv / Table 2 +1.50 cross.csv / sc +0.004) — label sources; T2-2 A1 §1.5 interferon module is 8 genes, none of Isg15/Mx1/Cxcl10/Ifit3, and it is a *negative control* whose +24.9 contradicts the compartment story — acknowledge; T2-3 A2 F6 Table 1 `q` is a rank test vs transcriptome, not a test of mean t — rename column + footnote; T2-4 A2 F8/A4 F3-secondary Table 2/3 selection rule + blanked cells — caption note; T2-5 A4 F11 Table 4 omits 4 qualifying modules — fix `regenerate_ms_tables.py` order list; T2-6 A4 F13 Gfap "nominally up" at +0.004 log2 — delete or reword; T2-7 A4 F14/F15 (6–9% offset; 2–6 residual df); T2-8 A5 F5 Methods names 10 compartments under "eight"; T2-9 A5 F6 "confirms absence"→"no up-regulation detected; bounded below detectable size"; T2-10 A5 F7/A5 F8 retitle + scope the universal negative; T2-11 A5 F10 add ≥6 refs (Ernst 2004; Cowen 2017; Gloor 2017; Robinson & Oshlack 2010; Ritchie 2015; Lotfollahi 2023); T2-12 A5 F9 ethics sentence completeness; T2-13 A2 F11 `U = 2/14` → `U = 2 of 16 possible`; T2-14 A5 flagged orphan `VALID_GSE249443` artefacts — deposit-note or mention.

### Tier 3 — format

T3-1 abstract ≤300 words (cut KO/OE bookkeeping + 32/132 detail); T3-2 6 keywords; T3-3 running head ≤50 chars; T3-4 A5 F12 rebuild Table 5 as compartment × score matrix or delete; T3-5 A4 F16 delete orphan `perturbation_bone.csv` from deposit or explain.

## 5. Consensus / complementarity / disagreement

**Consensus (≥3 experts):** myeloid drop (A1/A3/A4); ISG-list incompleteness (A1/A2/A5); "confirms absence" over-claim (A1/A3/A5); lowess universe (A2/A4); desk-reject front matter (A5 alone, no dissent).

**Complementary:** A3's regression/purity/LOO machinery explained *why* A1's biological suspicions were right; A4's regeneration-script test converted A2's arithmetic into a provenance fact (hand-edited cell); A5 located the strategic cost of the framing that A1/A3 showed was empirically over-claimed.

**Disagreements adjudicated:**
1. **A2 F3 vs A5 F8 (power).** A5 accepted MDE 5.65 and asked to scope the universal negative; A2 showed the MDE premise is false and the corrected MDE is 6× below the observed effect. **Ruling: adopt A2.** A5's fix is then wrong in direction — the honest statement is *not* "underpowered, cannot rank" but "adequately powered, and the failure is non-separation (controls outrank targets)". A5's strategic point (foreground the biology) survives independently.
2. **Severity: A4 F4/F13 Minor-ish vs A1 Major framing.** Adopted the stricter reading where the sentence carries a conclusion (§3.2 mechanism claim; satellite-gliosis attribution), lenient where purely numeric display.
3. **A3 "delete SatelliteGlia rows" vs A1 "report as unconfirmed".** Ruled with A3's stricter option for the *named* findings (delete from narrative), A1's option for the *table* (retain with LOO columns).

## 6. Priority must-fix list

**Must reword only:** T0-3, T1-1 (incl. gate literal), T2-6/8/9/13, T3-1..3.
**Must add analysis/recompute:** T0-1 (re-engineer or withdraw compartment layer), T0-2 (re-derive T5), T0-4 (reconciliation subsection), T1-2 (8-case enumeration), T1-3, T1-5 (source ρ or delete), A4 F3 (regenerate the erythrocyte cell), A2 F6 footnote.
**Must do before submission (author actions):** T0-5 — create the repo, reserve Zenodo DOI; add 5 front-matter statements; cut abstract.
**Gate work:** update 98-gate expectations derived from corrected artefacts (3/4,901; new T5; T1 two-sided; ISG list of 4) — never hardcode; add a gate that the regeneration script reproduces **every** table cell (would have caught #9).

## 7. What stands up (do not change)

A2 S1–S11 + A4 §Stands-up + A5 §Stands-up: exhaustive enumeration is exhaustive; **Cd68 D7-vs-PTX reproduces exactly from raw counts**; KO/OE mirror exact (max|Δ|=0); 264→132 collapse correct; min q 0.165/p 0.0014/glial-immune row exact; erythrocyte is the best negative control and holds; 32/132 and 14-control split reproduce; whole-DRG pseudobulk layer is technically clean and reproduces to the last digit (A3 Finding 6); myeloid-fraction honesty (d 1.37 vs MDE) praised by A5 as "the paper's best editorial asset"; Table 1/2/3 headline numbers match artefacts; network ρ 0.991 stands (A5 recomputed 0.992 on a different row aggregation).

## 8. Recommended handling paths

- **Path A (recommended): restructure, same article type.** Foreground the aetiology divergence (positive, exact, survives) with the model-specificity claim properly scoped (T0-4); demote the perturbation calibration to an explicitly secondary "cautionary calibration" section carrying the design-independent failures (erythrocyte |z|=7.39; ρ=0.991) **and the corrected power statement (T0-2)**; withdraw or re-engineer the compartment layer (T0-1). Fix all Tier-1 numbers. This is honest and requires no new data.
- **Path B: downgrade compartment layer to a supplementary note** with the A3 fixes, keep main text bulk-only. Lower risk, less novel.
- **Path C: wording-only.** **Not viable** — Tier 0 items are conclusions, not phrasing.

## 9. Process lessons (what gates could not catch)

1. The 98-gate verified arithmetic and provenance but the wrong p-floor literal **was enforced as truth** — the 2026-09-27 hazard recurred in a new form. New gates needed: (a) every table cell must be regenerated by script (catches hand-edits); (b) every ρ/derived statistic in prose must have a named artefact source; (c) gate expectations must be derived at run time, never stored as literals.
2. An argmax-annotation can only be falsified by an *independent* assignment criterion (A3's) — add a mandatory "compartment purity + canonical-marker-inside-compartment" QC artefact.
3. LOO robustness must be a deposited artefact, not an option — every named result should carry its 8/8 survival count.
4. Per-compartment BH families silently changed the meaning of a cross-compartment claim; family definition must be stated wherever a claim spans families.

## 10. Note to the author

The strongest claim that does not hold — "the interferon programme is absent in every DRG compartment, and the injury programme is neuron-intrinsic" — fails for structural reasons that are not carelessness: the dataset contains **46 myeloid cells in total**, so the compartment that carries interferon biology was unanalysable, and an argmax annotator built on overlapping marker sets cannot localise a programme whose whole-tissue and compartment effect sizes are identical to within 3%. Meanwhile the findings that *do* survive — the myeloid/interferon vs neuronal-injury divergence on canonical markers with exact permutation tests (q ≤ 0.027), the erythrocyte control outranking every analgesic target (|z| = 7.39), and the network model redundant with an additive baseline (ρ = 0.991) — are the paper. The corrected power analysis makes the methods story *stronger*, not weaker: the perturbation layer is adequately powered and still fails, because controls outrank targets. Swapping the headline to the biology and presenting the calibration as an actionable cautionary layer converts this from a fragile negative into a stable structural finding.
