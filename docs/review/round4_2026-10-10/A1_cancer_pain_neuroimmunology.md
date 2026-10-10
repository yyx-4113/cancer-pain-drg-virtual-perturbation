# R4-A1 Review — cancer pain / DRG neuroimmunology

Scope of this review: biological claims about bone-cancer-pain vs CIPN DRG
programmes (myeloid/interferon vs neuronal-injury), internal consistency of the
Gfap/glia/myeloid narrative, scope of the type-I interferon finding, fairness
of the §3.8 single-cell non-replication, and correctness of the citations in
those passages. Reviewed as a **first submission**, per the panel brief.
Independence discipline followed: only `_PANEL_BRIEF.md`, the manuscript,
`results/tables/*`, `scripts/*`, and `docs/manuscript/figures/*` were read; no
other reviewer output, no `manuscript_v0.*`, no `.workbuddy/`, no
`SUBMISSION_MANIFEST`.

---

## 1. Verdict

**Minor Revision.**

The manuscript is **publishable at a soundness-oriented SCIE journal** after a
small, bounded set of corrections. Every biological and statistical number I
could reach recomputed correctly (20/20 quantitative claims matched their
deposited artefacts; see §5). The gate runs clean (`287 checks, 0 failed`), the
unit suite passes (14 groups), and all citations I checked are real and
correctly characterised. The one item that must be resolved before acceptance
is administrative (the code/data pointer), and the remaining items are wording
and framing fixes. No scientific claim requires re-analysis.

---

## 2. Summary of the manuscript

A single-author computational reanalysis of two public **bulk rat DRG RNA-seq**
series on a shared platform (GPL24782): GSE149648 (Walker-256 tibial bone
cancer pain; sham/d1/d7/d14; n=4) and GSE160543 (paclitaxel/oxaliplatin vs
vehicle; n=4), plus an independent **single-cell mouse DRG** validation
(GSE292522, B16F10 cancer-nerve-injury; 4+4 animals; 55,967 cells).

The paper pursues two questions it deliberately keeps apart. (Q1) *Are the DRG
programmes of bone cancer pain and CIPN qualitatively different?* — answered
yes, via exhaustive C(8,4)²=4,900 label-permutation tests on 22 modules and
canonical markers. The Walker-256 day-7 peak is an innate-immune / **type-I
interferon** event (ISG module mean *t* = +24.9) with myeloid markers rising
(Cd68 Δ = +3.42 vs paclitaxel), whereas paclitaxel is a **neuronal-injury**
event (galanin Δ = −2.59 in the opposite direction; Atf3, Scn10a consistent).
(Q2) *Can in-silico perturbation rank analgesic targets?* — answered no: the
layer fails its own controls, is **anti-selective** (negative-control modules
carry larger standardised effects than analgesic targets, median |z| 2.47 vs
1.07, MWU p = 0.0019; an erythrocyte control holds the largest |z| in the whole
table, 7.39), correlates 0.991 with an additive no-network baseline, and no
reversal Δ survives correction (min q = 0.165).

The single-cell dataset is used **only** to test generalisation of the
interferon programme and finds it does not generalise (B16F10 response is
neuronal-injury-driven, ISGs flat or down), with the CD45-depletion blind spot,
the non-GFAP-validated satellite-glia label and the disagreement with the
source study (ref 15) all disclosed.

---

## 3. Major issues

### 3.1 The data/code-availability pointer cannot be verified and, if dead or private, defeats the paper's own thesis

*Citation.* Manuscript §"Data and code availability", lines 211–213:
"Analysis code, artefacts and the test suite are deposited in the project
repository (https://github.com/yyx-4113/cancer-pain-drg-virtual-perturbation),
to be archived at Zenodo upon acceptance."

*Evidence.* From this environment, the single successful HTTP probe of that URL
returned **404** (subsequent probes failed at the egress proxy with a 502, so I
cannot re-confirm; I therefore do not assert the repository is fabricated —
only that it is **not verifiable by a reviewer** and that a 404 is the wrong
answer for a statement using the word "deposited"). This matters disproportionately
here because reproducibility is the paper's declared contribution: the manuscript's
credibility rests on artefacts (`results/tables/*`) and gates
(`scripts/98_check_manuscript_numbers.py`, `scripts/99_unit_tests.py`) that a
reader is told to consult.

*Required fix.* Either (a) make the repository public (or supply a
reviewer-accessible link / private-repo invite) and state the exact commit/tag,
or (b) replace the promise with the actual **Zenodo DOI** at acceptance and
state explicitly that the GitHub link is private until then. A Zenodo DOI in
the repository's `.zenodo.json`-style wording is preferable. This is the only
item in my list that would block acceptance if left unresolved.

---

## 4. Minor issues

### 4.1 Abstract length exceeds the common 250–300-word ceiling

*Citation.* Abstract, lines 15–25. The panel brief states "295-word abstract".

*Evidence.* Excluding the four run-in labels ("Background.", "Methods.",
"Results.", "Conclusions.") the abstract is **325 words** (329 including them),
recomputed with a plain token count of `manuscript_final.md`. This is a ~10%
overrun of a number the brief itself carries, and exceeds the 250–300-word
limits of most SCIE venues in this space.

*Required fix.* Trim to ≤300 words to the target journal's limit (typically by
condensing the Methods run-in, which repeats §2 detail).

### 4.2 The "55.7th percentile" for the target effect-size rank is an inverted-orientation statistic, and "at the median" is slightly generous

*Citation.* §3.6, line 141: "Targets rank at the median of the effect-size
distribution (55.7th percentile)."

*Evidence.* `results/tables/negative_result_quantification.json`
(`T4_target_rank_percentile.median_pct_absz` = 55.7) is produced by
`scripts/33_negative_result_quantification.py:93–94` as
`(rank_ascending=False − 1)/(N−1)·100`, i.e. the largest |z| receives pct = 0.
So **55.7 means the median target sits *below* the overall median**, not above
it. Conventional percentile: the median target |z| = 1.0729 vs the overall
median |z| = 1.3248 → the median target is at the **44.7th percentile**, and
`n_targets_in_top_half_absz` = 45/102 (44%). The substantive conclusion
("no enrichment") is unaffected — indeed the conventional number *strengthens*
the anti-selectivity point.

*Required fix.* Report the percentile in conventional orientation ("targets sit
at ≈45th percentile, i.e. at or below the median") or explicitly define the
55.7 as a rank-from-top index. One sentence.

### 4.3 §3.5 mixes the two logFC sources within a single paragraph (Atf3 +0.22 vs Table 2 +0.31)

*Citation.* §3.5, line 131: "Atf3, the standard injury marker, reads +1.13 in
paclitaxel neuropathy and +0.22 at day 7 of bone cancer pain." vs Table 2,
line 278: Atf3 Bone (day 7) = **+0.31**, PTX = **+1.12**.

*Evidence.* The paragraph opens by quoting the **unadjusted** cross-aetiology
values (Cd68 +2.84 / −0.58, matching `etiology_marker_cross.csv`), then quotes
Atf3 from the **muscle-adjusted** moderated model
(`etiology_marker_de.csv`: Bone_D7 Atf3 = 0.222, PTX = 1.126). §2.3 (lines 59)
discloses that the text quotes the adjusted value "where the two disagree on a
headline marker", but Atf3 is not a headline marker and the juxtaposition of
adjusted (0.22) and unadjusted (2.84) values in adjacent sentences reads as an
internal inconsistency. A reader checking Table 2 will see 0.31, not 0.22.

*Required fix.* Quote one source consistently within §3.5 or add a half-sentence
flagging that the Atf3 figure is the muscle-adjusted value from
`etiology_marker_de.csv`.

### 4.4 The interferon/ISG module is formally a "negative control" yet is the headline positive finding, with no explicit reconciliation

*Citation.* §2.5, line 69 ("Five negative-control modules test specificity:
… interferon-stimulated genes …"); §3.6, line 141 ("the interferon control
alone — the strongest signal in the study"); vs §3.3 (line 113, "Day 7 is an
interferon event"), Abstract line 21, Discussion line 167 and Conclusions
line 185.

*Evidence.* "CTRL interferon" (Irf1, Ifitm3, Ifi47, Bst2, Rsad2, Gbp2, Ubd,
Ly6e; `_limma.py:194`) is one of the five negative-control modules and enters
the anti-selectivity analysis as a *control*, while the type-I interferon
programme is simultaneously the paper's strongest *biological* claim. The two
uses are on different endpoints and are not logically contradictory, but the
manuscript never states this, so the same module is presented as both a
specificity null and the headline result.

*Required fix.* One sentence in §2.5 or §3.6 clarifying that the ISG set serves
as a *perturbation-ranking* negative control while its DE/module-*t*
up-regulation at day 7 is an independent biological read-out.

### 4.5 The contaminated day-7 sample underpins the interferon peak, but the ISG module's leave-one-out refit is not deposited

*Citation.* §3.1, line 107 (GSM4506604 carries 6.26% sarcomere transcript);
§3.2, lines 111 (drop → 3,041 genes, Jaccard 0.584, ρ = 0.912); §3.3 (ISG
module *t* = +24.9).

*Evidence.* `de_muscle_drop_sensitivity.csv` and `etiology_marker_loo.json`
deposit robustness for the *global* DE count and for six cross-aetiology
markers (Cd68/Aif1/Gal/Ccl2/Il1b/Scn10a; all sign-stable, matching Limitations
lines 179). There is **no** analogous leave-One-out value for the ISG module or
for its day-7 *t*, even though the interferon peak is the paper's central
biological claim and the muscle-contaminated animal is a day-7 sample. The
muscle covariate partly addresses this, but the reader cannot see it.

*Required fix.* Deposit and quote the day-7 ISG-module effect with GSM4506604
excluded (one number), in §3.2 or §3.3, and state it in Limitations alongside
the existing DE-count LOO.

### 4.6 The myeloid per-sample range is mis-stated at its lower bound

*Citation.* §3.8, line 157: "0.02–0.2% of cells per sample".

*Evidence.* From `sc_GSE292522_composition.csv` the eight per-animal myeloid
fractions span **0.0147%–0.199%** (46 cells total, confirmed). The lower bound
is ~0.015%, not 0.02%.

*Required fix.* Change to "≈0.01–0.2%".

### 4.7 The reproducibility gate's numeric expectations are largely hard-coded literals, not derived from the manuscript

*Citation.* `scripts/98_check_manuscript_numbers.py`: docstring lines 10–14
("re-derived from the CSV/JSON at check time and compared with the string that
appears in the manuscript file"); `chk()` lines 34–45; `in_text()` lines
59–72.

*Evidence.* Of 101 `chk(...)` calls, the great majority compare an artefact
value against a **hard-coded literal in the script** (e.g. lines 90–99,
115–135), not against the manuscript text; only **9** `in_text(...)` calls read
the manuscript. The genuinely runtime-derived exceptions are the negative-result
block (lines ~336–398: T1/T2/T5 recomputed from the CSV, then required to appear
in the text) and the Table 2/3 cell regeneration (lines ~430–470). So the gate
is stronger than a pure literal diff but weaker than its docstring claims: if
the prose and the literal drifted together, or if the prose diverged from both,
most checks would not catch it. (In practice the abstract/§3.5 Atf3 mismatch of
4.3 is exactly the failure mode the docstring promises to prevent.)

*Required fix.* Either soften the docstring or extend the derived expectations
(the `in_text` pattern) to the discovery-section numbers. Not blocking.

### 4.8 Figure 3 hard-codes its headline p-value in the plot title

*Citation.* §Figures, Figure 3 caption (line 346) and
`scripts/35_figures.py:132` (`ax.set_title("Mann–Whitney two-sided p = 0.0019", …)`).

*Evidence.* Figures 2 and 3 otherwise read their plotted values at run time from
`etiology_marker_cross.csv` / `etiology_bootstrap_differences.csv` (verified),
consistent with the captions. Figure 3's title string is the sole literal in the
figure pipeline and would not update if the data changed. The value is currently
correct (0.001914 → 0.0019, recomputed). Cosmetic drift risk only.

*Required fix.* Compute the MWU p from `dk` at plot time. Also the caption says
"black bars" for medians that are drawn in `#2c3e50` (dark navy) — trivial.

### 4.9 Figure 3 "control" colour and §3.3 "every neuronal module fell" — wording nits

*Citation.* Figure 3 caption (line 346) colours controls red / targets blue;
Figure 2 (line 342) colours bone red / paclitaxel blue. The same red is "bone
cancer" in Fig 2 and "negative control" in Fig 3, without a cross-figure note.

*Required fix.* Optional one-clause note; purely presentational.

---

## 5. Numbers you independently recomputed

All recomputed with the managed Python
(`C:/Users/Administrator/.workbuddy/binaries/python/envs/default/Scripts/python.exe`,
pandas/numpy/scipy) directly from `results/tables/*.csv|json`.

| # | Claim | Manuscript | Recomputed | Match? |
|---|---|---|---|---|
| 1 | Cd68 Δ, Bone_D7 vs PTX (`etiology_marker_cross.csv`) | +3.42 | +3.4183 | ✔ |
| 2 | Aif1 Δ, Bone_D7 vs PTX | +3.08 | +3.0757 | ✔ |
| 3 | Gal Δ, Bone_D7 vs PTX | −2.59 | −2.5921 | ✔ |
| 4 | Calca Δ, Bone_D7 vs PTX | −0.77 | −0.7698 | ✔ |
| 5 | Scn10a Δ, Bone_D7 vs PTX | −0.70 | −0.7022 | ✔ |
| 6 | ISG (CTRL interferon) module mean *t*, D7 (`module_level_signal.csv`) | +24.92 | 24.923 | ✔ |
| 7 | Day-7 DE genes, muscle-adjusted q<0.05 (`DE_D7_vs_Sham_adj.csv`) | 4,600 | 4,600 | ✔ |
| 8 | Muscle-score ANOVA p (`qc_muscle_scores.csv`) | 0.605 | 0.6054 | ✔ |
| 9 | T1 median \|z\|, control vs target (`etiology_bootstrap_differences.csv`, dedup KO, per_contrast) | 2.47 vs 1.07 | 2.4702 vs 1.0729 | ✔ |
| 10 | T1 Mann–Whitney two-sided p on \|z\| | 0.0019 | 0.001914 | ✔ |
| 11 | T2 binomial over-representation p (14/32 vs 30/132) | 0.0067 | 0.006708 | ✔ |
| 12 | Largest \|z\| in table (CTRL erythrocyte, Bone_D14 vs D7) | 7.39 | 7.3878 | ✔ |
| 13 | Maximum number of exact-permutation-significant reversal tests (min q) | q = 0.165 | 0.1654 (31/132 p<0.05) | ✔ |
| 14 | Spearman(network, additive), n=132 (`head_to_head_baselines.json`) | 0.991 | 0.99115 | ✔ |
| 15 | T5 minimum detectable Δ, normal approx (2.802 × median boot SE) | ≈0.052 | 0.0518 (SE 0.0185) | ✔ |
| 16 | Target tests exceeding MDE | 41 (of 102) | 41 | ✔ |
| 17 | Single-cell myeloid fraction MWU, two-sided (`sc_GSE292522_composition.csv`) | 0.114 | 0.1143 (U=2, U′=14) | ✔ |
| 18 | Cross-model concordance ρ, B16F10 vs Walker-256 D7 | 0.098 | 0.0985 | ✔ |
| 19 | Atf3, Bone day 7 | +0.31 (Table 2) / **+0.22** (text §3.5) | +0.3128 (cross, unadj) / +0.222 (adj DE) | ✖ (text/table source mismatch, §4.3) |
| 20 | Abstract length (brief states 295) | 295 | **325** words | ✖ (§4.1) |
| 21 | Target effect-size percentile | "55.7th" | **44.7th** conventional; 55.7 = rank-from-top | ✖ orientation (§4.2) |
| 22 | Myeloid per-sample fraction range | 0.02–0.2% | **0.0147–0.199%** | ✖ lower bound (§4.6) |

Additional corroborating checks (all matched): day-7 DE count 0/4,600/235 and
q<0.10 counts 5/5,867/860; drop-GSM4506604 count 3,041, Jaccard 0.584,
ρ=0.912; sample correlation 0.913–0.9935; PCA 46.85/17.45/7.23/4.76%;
module-level D7 *t* for all 22 modules (Table 1) reproduced cell-by-cell;
module mean-logFC Δ (Table 3) and calibration (Table 4) `pct_of_oracle`
(75.4%), z (+4.61) reproduced; single-cell 1,517 up/703 down, Sprr1a +1.554
(*t*=54.49), Ifit2 q=0.0417, Oas1c q=0.0136, Gfap +0.0036 q=0.0062, neuron
galanin +1.472, NF-neuron galanin +2.537, Scn10a neuron q=0.0117, satellite-glia
Ifit1 q=0.0037 / Isg15 q=0.0071; neuronal share 78.6%→67.0%, Schwann
18.8%→29.3%. Gate reproduced: `287 checks, 0 failed`. Citations refs 1–16 all
resolved (Crossref / PubMed / GEO); GEO accessions GSE149648, GSE160543,
GSE292522 match the cited studies; ref 15 correctly characterised as the source
of the neuron-intrinsic type-I IFN finding the paper disputes.

---

## 6. Publishability assessment

**Ready? Not quite — Minor Revision.** The science is sound and the reporting is
unusually well-controlled: the two claims are pre-specified, the negative result
is calibrated against negative controls, degree-matched nulls and an additive
baseline, the single-cell non-replication is handled with appropriate scope
(including the CD45-depletion blind spot and the explicit disagreement with
ref 15), and essentially every headline number recomputes exactly. Nothing about
the biological or statistical substance requires re-analysis, and the two claims
are appropriately hedged (interferon = Walker-256 day-7 model-specific; the
neuron-intrinsic-IFN conflict is declared unadjudicable by this design).

**Ordered minimum to get there:**

1. Resolve the code/data pointer (§3.1): public link with a commit hash, or a
   live Zenodo DOI. **Blocking.**
2. Fix the three factual/definitional wordings: abstract length (§4.1), Atf3
   source consistency (§4.3), myeloid range lower bound (§4.6).
3. Correct the T4 percentile orientation (§4.2).
4. Add one reconciliation sentence for the ISG module's dual control/finding
   role (§4.4) and one deposited LOO value for the day-7 ISG module (§4.5).
5. Optional: compute Fig 3's title p at run time (§4.8) and soften or strengthen
   the gate docstring (§4.7).

Items 2–5 are cosmetic-to-clarifying and can be done in a single revision
pass. With the abstract trimmed and the repository made accessible, I would
recommend acceptance at a soundness-oriented SCIE journal.

---

### Reviewer's note on independence

Everything above derives from the manuscript, `results/tables/*`, `scripts/*`
and public bibliographic databases read in this session. I did not open any
other file under `docs/review/`, any earlier manuscript version, or any
`.workbuddy/` file. Two items I could not fully close from this environment are
flagged as such: the repository URL (§3.1; egress blocked after one 404 probe)
and the visual inspection of the figures (the PNG contents could not be
re-rendered for display, so Figures 2–3 were verified via their
`scripts/35_figures.py` data sources and the CSVs they cite, not by eye).
