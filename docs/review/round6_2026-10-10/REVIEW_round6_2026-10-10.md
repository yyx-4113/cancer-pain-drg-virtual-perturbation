# Editorial Synthesis — Round 6 Independent Review (2026-10-10)

**Manuscript:** "Divergent dorsal root ganglion programs in bone cancer pain and chemotherapy-induced peripheral neuropathy: an exact cross-study reanalysis and a negative-control calibration of in-silico perturbation"
**Mode:** Fresh 4-reviewer panel. Independence enforced — each reviewer was forbidden from opening `docs/review/round1..5_*`, `REVIEW_round*`, `R3_editorial_review.md`, `author_verification_statement.md`, and each other's files. (The prior rounds' "gate all green" status is treated as *arithmetic* verification only; the panel's remit is the design/biological/implementation layers gates cannot reach.)
**Panel:** A1 Domain (pain neuroscience / anaesthesiology), A2 Design statistician, A3 Implementation/provenance auditor, A4 Venue/reporting editor.
**Outcome of checks:** A3 re-ran `98_check_manuscript_numbers.py` (291/0), `99_unit_tests.py` (all pass), `33_negative_result_quantification.py` (byte-equivalent JSON), and independently traced every cell of Tables 1–5 to the deposited CSV/JSON. A2/A4 recomputed T1–T5 from artefacts. All four agree the **permutation machinery and T1–T5 statistics are correct and internally consistent.**

---

## Verdict

**Revise (moderate) — not a reject.** The statistical core is sound and was independently re-verified. The blocking issues are **interpretive over-statements in the biological narrative (A1)** and a short list of **cross-reference numeric drifts (A4)**; the statistical-design critiques (A2) are real but secondary framing issues. None require new experiments. Estimated effort: a focused rewrite of §3.3/§3.8/§4/§5 plus a §2.3 clarification and three number-consistency edits.

---

## Consolidated findings by severity

### P0 — Must fix (affect conclusions or would block at a careful journal)

**P0-1. "Day 7 is an interferon event, *not* an ion-channel event" is over-stated (A1 Item 1).**
- Evidence: §3.3 (manuscript line 113). The data show neuronal/ion-channel genes **down-regulated**, not absent — Nav1.8 (Scn10a) is −0.46 in bone cancer (Table 3, line ≈291), and §3.3 line 117 itself attributes this to compositional dilution. The source study Zhai 2021 (ref [10]) titled its own paper "ion channels and cytokines" and is cited only for provenance.
- Fix: re-frame §3.3 heading and first sentence as "a prominent interferon-stimulated-gene programme and myeloid marker up-regulation, on a background of neuronal-transcript dilution." Full paste-ready wording in A1 §Item 1.

**P0-2. The ref [14] (Baruch 2025) "conflict our design cannot adjudicate" sentence is factually wrong (A1 Item 2).**
- Evidence: the paper's own per-compartment test found Ifit2 *down* in neurons (q=0.046, §3.8 line 157; Table 5), which is *negative* evidence against Baruch's neuronal-intrinsic IFN claim — so the neuronal arm IS partly adjudicated. Baruch's IFN-I is explicitly an ATF3-driven nerve-healing programme (congruent with the paper's own Sprr1a/Atf3 neuronal-injury reading), which the paper omits and which would *strengthen* its thesis. A1 verified the Baruch abstract (Nature 2025;646:462–473, doi:10.1038/s41586-025-09370-8): the attribution "neuron-intrinsic type-I-interferon response in B16F10" is accurate (not a strawman), so only the "cannot adjudicate" framing is wrong.
- Fix: replace lines 159–161 with the narrowed, honest wording in A1 §Item 2 (neuronal arm reads negative; myeloid arm genuinely untestable because CD45⁺ depleted; note Baruch's IFN is an injury programme).

**P0-3. "Glial activation (Gfap up)" is listed among the *secure* conclusions, but the paper's own companion data undermine it (A1 Item 3).**
- Evidence: §4 line 167 lists "glial activation (Gfap up)" as secure; yet §2.8/§3.8 state the satellite-glia bin "carries almost no Gfap transcript" (7 UMIs across 8 animals). GFAP is a weak DRG glial marker. Gfap is also the weakest of the strongly-up glial/immune genes (Cd68 +2.84, Aif1 +2.19, Gfap +1.33/+1.50).
- Fix: drop "(Gfap up)" from the §4 secure list and add the reconciling caveat in A1 §Item 3.

**P0-4. The §4 clinical-pharmacology sentence is asserted, not supported (A1 Item 4).**
- Evidence: §4 line 167 claims bone pain "responds to anti-resorptive and anti-inflammatory strategies" and CIPN "responds partly to agents acting on neuronal excitability," citing [1][2][16]. Refs [1][2] are generic analgesia ladders (ESMO/NCCN) that do not stratify by DRG mechanism; ref [16] (Luger 2005) emphasizes osteoclast/NGF/endothelin, not IFN; the evidence-based CIPN agent is duloxetine (an SNRI), not primarily a "neuronal excitability" drug. This contradicts §5 line 181 ("nothing here is validated in humans").
- Fix: convert to hypothesis-generating language (A1 §Item 4) and tie the §5 limitation to it.

### P1 — Should fix (consistency / honesty / framing)

**P1-1. §3.6 says MDE ≈0.052 but §5 says ≈0.05 — numeric drift (A4 F1, verified).**
- Evidence: line 143 "≈0.052 log-units" vs line 179 "≈0.05 log-units". The precise value is 2.802×0.0185 = 0.0518. Fix: use "≈0.052" in both (or "≈0.05" in both) — editor recommends "≈0.052" everywhere since the gate asserts 0.052.

**P1-2. §2.3 dual-logFC disclosure must clarify the estimator difference, not imply only "few hundredths" (A4 F3, root cause corrected by editor).**
- Editor's verification: §3.3 values come from `etiology_marker_de.csv`, which `script 22` builds on the **13,601 shared** genes (line 580 restricts `Lb/Lc` to `shared` before `de_positive_controls`); §3.2 single-study counts use 13,893. So A3 F1's claim that "§2.3 wrongly attributes §3.3 to the 13,601 fit" is a **FALSE POSITIVE** — §2.3 is correct. The genuine gap A4 flagged is **§3.3 (muscle-adjusted limma-trend moderated logFC) vs Table 3 (unadjusted per-sample mean logFC)**, both on the same 13,601 set: Cd68 2.86 vs 2.84 (0.02 ✓), Aif1 2.18 vs 2.19 (0.01 ✓), but **Gfap 1.33 vs 1.50 (0.17)** and **Ccl2 0.81 vs 0.90 (0.09)** — far beyond "few hundredths." This is method/expression-driven (moderation shrinks low-expression genes more), not gene-set nesting.
- Fix: in §2.3, state explicitly that §3.3 reports the muscle-adjusted limma-trend *moderated* logFC while Table 3 reports the *unadjusted per-sample mean* logFC on the same gene set; the two agree within a few hundredths for well-expressed markers (Cd68, Aif1) but can differ by ~0.1–0.2 log₂ for lower-expression markers (Gfap, Ccl2), and the muscle-adjusted value is quoted as primary when they disagree.

**P1-3. T2 binomial framed as a formal test though its success probability is not constant across trials (A2 F1).**
- Evidence: "14 of 32 randset-significant are controls, vs base rate 30/132, binomial p=0.0067" (line 139). The per-trial success probability depends on module kind, so the binomial is better read as a descriptive over-representation than an independent test. The recomputed p=0.0067 is correct; only the interpretation/framing needs softening ("descriptive over-representation … binomial p≈0.0067 as an approximate gauge").

**P1-4. "132 independent tests" overstates independence; effective n not quantified (A2 F2).**
- Evidence: the 132 (module×contrast) stats are correlated within contrast; T1/T2 already admit "approximate p" but the headline "132 independent tests" is stated absolutely (line 137). Fix: add "≈132 tests (correlated within contrast; effective n smaller)" or report a conservative effective-n estimate.

**P1-5. Strengthen the AI-disclosure statement (A4 A1).**
- Evidence: § "Declaration of generative AI in the writing process" says an LLM was used for "code debugging and manuscript language editing." Many journals (PLOS, Nature, Elsevier, Wiley) require the tool named, an attestation that AI was absent from analysis design/interpretation, and human responsibility. Fix: name the tool, explicitly attest no AI in statistical design/interpretation, confirm human authorship of all conclusions.

**P1-6. Add the missing mechanistic citations (A1 Item 5).**
- (i) STING/TBK1/NF-κB + IL-1β/IL-6/TNF-α in Walker-256 BCP DRG neurons (direct mechanistic support for the interferon claim); (ii) an explicit sentence in §3.3 engaging Zhai 2021's "ion channels and cytokines" conclusion; (iii) a DRG-resident-macrophage-in-pain citation (justifies Cd68/Aif1 as DRG myeloid markers); (iv) a CIPN DRG neuropeptide/ATF3 citation.

### P2 — Minor / optional

- **P2-1. Table 1 Glial/immune n=5 but module defined with 6 genes (A1 R1).** Footnote that one gene (Itgam or S100a9) dropped at the expression filter; effective n=5.
- **P2-2. §3.6 vs §5 rounding: 3/4,901 = 6.1×10⁻⁴ (line 87) vs 6×10⁻⁴ (line 179) (A4 F2).** Use 6.1×10⁻⁴ in both.
- **P2-3. T5 t-based MDE "≈0.10" may mis-specify a two-sample geometry (A2 F4).** The normal-approx 0.052 is the correct headline; either drop the t-based aside or derive it correctly. Low priority since it is flavor text.
- **P2-4. T3 `median_se_expressed_genes` is computed on a biased subset (A2 F3).** `script 33` line ≈6755 selects `q>0.5 & |logFC|≥0.2` — i.e. non-significant genes with sizeable effect, which have *inflated* SE — yet labels the result `median_se_expressed_genes`. The bias is conservative (over-states the gene-level MDE), so it does not weaken the negative conclusion, but the label is misleading. Fix: relabel (e.g. "median SE of stable, non-DE genes with |logFC|≥0.2") or compute on a representative expressed-gene set.
- **P2-5. README calls the GSE292522 control "MPEK-BL6" while the manuscript says "keratinocyte" (A3 F3).** Same entity, mismatched label — reconcile.
- **P2-6. Untracked working-tree files incl. `.github/` are not in the v1.0.0 tag (A3 F4).** Harmless; note that the release snapshot excludes them.
- **P2-7. Remote GitHub liveness UNVERIFIED (A3 F2).** Network was blocked in the sandbox (`git ls-remote` aborted), so the Data-availability statement's "release v1.0.0" URL could not be confirmed live. Editor recommends the author confirm the repo is public before submission. Note: per the project's standing decision, no Zenodo DOI is expected; GitHub release tag satisfies PLOS ONE / BMC / Scientific Reports, but a few journals require an archived DOI — check the target journal.

---

## Findings the editor REJECTED as false positives (verified)

1. **A3 F1 — "§2.3 wrongly attributes §3.3 to the 13,601 shared fit."** FALSE. `script 22` line 580 restricts `Lb/Lc` to `shared` *before* `de_positive_controls` runs, so `etiology_marker_de.csv` (§3.3 source) is correctly the 13,601 shared fit. §2.3 is right; do not "fix" it on this basis.
2. **A4 F4 — "3/4,901 denominator is internally contradictory."** FALSE. Line 87 states the observed split and its mirror are added to the 4,900 enumerated null pairs, giving 4,901 total — standard for an exact permutation test that includes the observed statistic. The denominator is correct.
3. **A4's reading that "few hundredths" applies to §3.3-vs-Table-3.** PARTIALLY misread: §2.3's "few hundredths" refers to §3.3 (13,601 shared) vs §3.2 (13,893 full) — a gene-set-only difference that is genuinely small. The larger §3.3-vs-Table-3 gap is a *different* comparison (moderated vs unadjusted estimator) and is the real, separately-fixed issue under P1-2.

---

## What the editor independently verified (trust-but-verify)

- Ran `grep`/`python` to confirm Gfap §3.3 +1.33 (de.csv 1.328) vs Table 3 +1.50 (cross.csv 1.498) → gap 0.17; Ccl2 +0.81 vs +0.90 → gap 0.09. Confirms P1-2.
- Read `script 22` lines 560–625 and confirmed `etiology_marker_de.csv` is built on the shared 13,601 genes → A3 F1 rejected.
- Read manuscript lines 87 and 179 and confirmed the 4,901 denominator rationale and the 0.052 vs 0.05 drift → P1-1 and A4 F4 rejected.
- Cross-checked A1 R1 (Table 1 Glial/immune n=5 vs 6-gene definition) against Table 1 line 244 → confirmed.

---

## Prioritized revision plan

1. **P0-1 → P0-4** (biological narrative): rewrite §3.3 heading/first sentence, §3.8 ref [14] paragraph, §4 secure-list + clinical sentence, and add the §5 tie-in. All prose; no number changes.
2. **P1-1, P1-2, P2-2**: three numeric/consistency edits in §2.3, §3.6, §5.
3. **P1-3, P1-4**: soften T2 binomial + "132 independent" framing in §3.6.
4. **P1-5, P1-6**: AI-disclosure strengthening + 4 citations.
5. **P2-1, P2-3…P2-7**: minor footnotes/labels/README sync.
6. After edits: re-run `98_check_manuscript_numbers.py` (gate will need its literal expectations updated for the §3.6/§5 0.052 consistency and any reworded numbers) and re-run A3's trace before re-submission.

**Recommendation:** Apply P0 and P1 as a single revision pass (these are the submission blockers), fold P2 in opportunistically, then re-gate and re-run the panel or a targeted P0-check before submitting to a pain/computational-biology journal (recommended targets given the negative-result emphasis: *Journal of Pain Research*, *BMC Medical Genomics*, or *PLOS ONE* with the strengthened AI disclosure).
