# R4-A4 Review — Data provenance, reproducibility and computational audit

Reviewer role: independent peer reviewer (data provenance / traceability / computational audit).
Manuscript: `docs/manuscript/manuscript_final.md` (single-author computational reanalysis).
Review round: Round 4. Posture: first submission, no credit for effort.
Independence: obeyed the brief. Read only `_PANEL_BRIEF.md` outside the manuscript; no other
`docs/review/` file, no `manuscript_v0.*`, no `.workbuddy` file was opened.
Python: `C:/Users/Administrator/.workbuddy/binaries/python/envs/default/Scripts/python.exe` (3.13.14).

---

## 1. Verdict

**Major Revision.**

The provenance and traceability of this manuscript are, on the whole, *exceptionally* strong — the
numbers I could check all trace to deposited artefacts, the certificate gate regenerates Table 2/3
cells cell-by-cell rather than pattern-matching them, and the software/test-suite claims are exact.
I found **no contradicted or untraceable number** among the ~45 I recomputed.

The verdict is nonetheless **Major Revision**, not Accept, on a single blocking item: the
**Data and code availability statement is not true as written and is not verifiable**. The URL it
gives (`https://github.com/yyx-4113/cancer-pain-drg-virtual-perturbation`) returns HTTP **404**
against both the GitHub HTML endpoint and the GitHub REST API, from a network that *does* reach
github.com. The sentence uses the present tense — "Analysis code, artefacts and the test suite
**are deposited in the project repository** (URL)" — which asserts an accessible deposit that does
not resolve. For a paper whose central methodological selling point is "every number was checked
against the deposited artefacts," an unreachable deposit is a material gap that a
reproducibility-oriented SCIE venue will (correctly) treat as blocking. The fix is one line and
changes no result, so this is a *scoped* Major, not a rejection.

---

## 2. Summary of the manuscript

A single-author reanalysis of public bulk rat DRG RNA-seq (GSE149648, Walker-256 tibial bone cancer
pain, sham/d1/d7/d14; GSE160543, paclitaxel/oxaliplatin CIPN; n=4/group, same platform GPL24782),
with an independent mouse single-cell validation (GSE292522, B16F10 cancer-nerve-injury, 4+4
animals, 55,967 cells). Two claims: (1) the two cancer-pain aetiologies drive qualitatively
different DRG programmes — myeloid/type-I-interferon at the Walker-256 day-7 peak, neuronal-injury
in CIPN — demonstrable on canonical markers with exhaustive label permutation (C(8,4)²=4,900 for
unpaired contrasts, 70 splits for shared-control contrasts); and (2) the same data do **not**
support ranking analgesic targets by in-silico perturbation, the layer being *anti-selective*
(negative controls carry larger standardised effects than analgesic targets) while still powered
for large effects. The perturbation failure is triangulated against negative-control modules,
degree-matched nulls, an oracle anchor, and additive/ridge baselines. Robustness to four free
parameters (α, k, gene-count, lowess span; 60+15 settings) is reported, and the single-cell
section discloses hard limits (CD45 depletion blinding the myeloid compartment; an unvalidated
satellite-glia label; a disagreement with the source study on neuron-intrinsic interferon).

---

## 3. Major issues

### M1. Data-availability statement asserts a deposit that does not resolve (blocking)

- **Citation:** `manuscript_final.md`, "Data and code availability" (line 213):
  "Analysis code, artefacts and the test suite are deposited in the project repository
  (https://github.com/yyx-4113/cancer-pain-drg-virtual-perturbation), to be archived at Zenodo
  upon acceptance".
- **Evidence:** I probed the URL directly. Network access to GitHub is confirmed working (the
  request returned a real GitHub TLS/session handshake and server headers, not a proxy failure):
  - `curl -sI https://github.com/yyx-4113/cancer-pain-drg-virtual-perturbation` →
    `HTTP/1.1 404 Not Found` (Server: github.com).
  - `curl -s https://api.github.com/repos/yyx-4113/cancer-pain-drg-virtual-perturbation` →
    `{"message": "Not Found", "status": "404"}`.
  I **cannot** rule out that the repository exists but is *private*; I state only what I observed:
  anonymous access returns 404. A private repository is not a valid data-availability deposit for
  a reader or reviewer.
- **Why it blocks:** the paper's stated value is reproducibility "against the deposited artefacts".
  The claim "are deposited" (present tense) is falsified/unverifiable as accessed. The second clause
  ("to be archived at Zenodo upon acceptance") is a future action and is fine — the problem is
  purely the present-tense assertion.
- **Required fix (exactly one of):** (a) make the repository public now and confirm the URL resolves
  to the deposited code + `results/tables/` artefacts + test suite; **or** (b) reword to the honest
  state, e.g. "Analysis code, artefacts and the test suite will be released publicly in the project
  repository (URL) and archived at Zenodo upon acceptance; they are available from the author on
  request during review." If the venue mandates open code/data at submission, (a) is compulsory.

---

## 4. Minor issues

### m2. The verification gate is partly tautological — 13% of checks are pure string-presence

- **Citation:** `scripts/98_check_manuscript_numbers.py:59-73` (`in_text`), and the 30-item
  key-statements loop at lines 481-488.
- **Evidence:** I ran the gate: it ends `287 checks, 0 failed` (`ALL MANUSCRIPT NUMBERS VERIFIED
  AGAINST ARTEFACTS`). Of the **287 runtime checks**, 250 are `chk(...)` calls and **37 are
  `in_text(...)`**. Every `in_text` call only asserts that a substring occurs *somewhere* in ~1,200
  lines; it never validates the number's context or relation to the artefact. 30 of the 37 use
  fully hardcoded needles (`"4,600"`, `"0.991"`, `"24.9"`, `"+3.42"`, `"32/132"`, …). If a value
  were moved to the wrong sentence but left present anywhere, these would still pass.
- **Notably under-verified:** the §2.7 floors `"3/4,901"`, `"3/71"`, `"2/71"`, and the §3.8
  `"46 cells"` and `"eight testable compartments"`, are **only presence-checked** — never derived.
  I therefore verified the floors by hand (see §5/§6): the shared-control logFC floor is genuinely
  `3/71` (I recomputed the minimum paired `p_permutation_exact` in
  `etiology_marker_cross.csv` and `etiology_module_logfc_cross.csv` as 0.0423 = 3/71), the unpaired
  floor is `3/4,901` given `exact_two_sided_p` (`scripts/_limma.py:420-426`, `p=(ge+1)/(n+1)` with
  the mirror in the null), and the `2/71` reversal floor follows from the reversal projection not
  being antisymmetric at the complement split. No error — but the gate does not prove this itself.
- **Fix:** derive these five values at check time (e.g. assert `min(paired.p_permutation_exact)
  == 3/71`) instead of asserting the literal is present.

### m3. 59% of gate assertions freeze the expected value as a literal rather than recomputing it

- **Citation:** the `chk(label, got, want)` calls where `want` is a hardcoded number (e.g. lines
  92-170, 178-266).
- **Evidence:** counting the 250 `chk` calls: **68** are the Table 2/3 cell regeneration (line
  447-477; both sides read at run time) and **13** are the T1/T2/T5 recomputations (lines 337-384,
  where `got`/`want` are computed from CSV/JSON). The remaining **~169 (≈59%)** compare an
  artefact-derived `got` against a **hardcoded literal** `want`. These catch *artefact drift*
  (good) but the "truth" is a frozen constant, so a coordinated change of both artefact and
  manuscript to a new value would fail rather than pass — conservative, but not self-recomputing.
- **Fix:** for the headline quantities, recompute the expectation from the raw analysis inputs
  where feasible; otherwise document that the gate is a *drift detector* anchored to a signed-off
  snapshot, which is a legitimate but weaker guarantee.

### m4. Figure 3 hard-codes its significance annotation

- **Citation:** `scripts/35_figures.py:132` — `ax.set_title("Mann–Whitney two-sided p = 0.0019", ...)`.
- **Evidence:** the plotted points and median bars are read at run time from
  `etiology_bootstrap_differences.csv` (lines 117-127), but the p-value in the title is a string
  literal, unlike every other plotted quantity. It currently matches (0.0019, see §5), but it would
  silently go stale if the CSV changed.
- **Fix:** compute the p in the script (it already loads the CSV) and format the title from it.

### m5. Two logFC sources give the same marker two values; disclosure is present but easy to miss

- **Citation:** §2.3 (line 59) vs §3.3 (line 115) and §3.5 (line 129) / Table 2.
- **Evidence:** §3.3 quotes muscle-adjusted values (Cd68 +2.86, Aif1 +2.18 — from
  `etiology_marker_de.csv`), while §3.5 and Table 2 quote the unadjusted per-sample means
  (Cd68 +2.84, Aif1 +2.19 — from `etiology_marker_cross.csv`). Both are correct for their stated
  path and §2.3 explains the split, but a reader can read "Cd68 +2.86" and "Cd68 +2.84" four lines
  apart and suspect an error.
- **Fix:** a one-line footnote under Table 2 ("Bone values differ from §3.3 by ≤0.02 because §3.3
  reports the muscle-adjusted fit and Table 2 the unadjusted per-sample mean, per §2.3").

### m6. Figure 1 caption over-claims provenance

- **Citation:** Figure 1 caption (line 338): "All plotted values are regenerated from
  `results/tables/`."
- **Evidence:** Figure 1 is a schematic diagram (`35_figures.py:34-83`); it plots no data values.
  The statement is true of Fig 2/3 but vacuous/loose for Fig 1.
- **Fix:** scope the sentence to Figures 2-3.

---

## 5. Numbers you independently recomputed

All values below were read from the deposited CSV/JSON with the managed Python and compared to the
manuscript. "match?" = within the manuscript's stated rounding. Source files are under
`results/tables/` unless noted.

| # | Claim (section) | Manuscript | Recomputed | Source | Match? |
|---|---|---|---|---|---|
| 1 | Cd68 Δ vs PTX (Abs/§3.5/Table 2) | +3.42 | +3.4183 | etiology_marker_cross.csv | yes |
| 2 | Aif1 Δ vs PTX (Abs/Table 2) | +3.08 | +3.0757 | etiology_marker_cross.csv | yes |
| 3 | Gal Δ vs PTX (Abs/Table 2) | −2.59 | −2.5921 | etiology_marker_cross.csv | yes |
| 4 | Ccl2 Δ vs OXA (Table 2) | +1.39 | +1.3907 | etiology_marker_cross.csv | yes |
| 5 | Il1b Δ vs PTX (Table 2) | +1.23 | +1.2305 | etiology_marker_cross.csv | yes |
| 6 | Scn10a Δ q vs PTX (Table 2) | 0.026 | 0.0260 | etiology_marker_cross.csv | yes |
| 7 | Calca Δ vs PTX (Table 2) | −0.77 | −0.7698 | etiology_marker_cross.csv | yes |
| 8 | Interferon module t, q (§3.3/Table 1) | +24.9, 6.8e-5 | +24.923, 6.8e-5 | module_level_signal.csv | yes |
| 9 | Glial/immune module t (§3.3/Table 1) | +11.4 | +11.412 | module_level_signal.csv | yes |
| 10 | GABAergic module t (§3.3/Table 1) | −5.71 | −5.708 | module_level_signal.csv | yes |
| 11 | Housekeeping q (§3.3/Table 1) | 0.52 | 0.52371 | module_level_signal.csv | yes |
| 12 | Glial/immune ΔlogFC D7 vs PTX (§3.5/Table 3) | +2.70 | +2.6954 | etiology_module_logfc_cross.csv | yes |
| 13 | Cytokine ΔlogFC D7 vs OXA (§3.5) | +1.75 | +1.7487 | etiology_module_logfc_cross.csv | yes |
| 14 | Neuropeptides Δ vs PTX, q (Table 3) | −1.35, 0.015 | −1.3453, 0.0145 | etiology_module_logfc_cross.csv | yes |
| 15 | Erythrocyte ΔlogFC D7 vs OXA (§3.6/Table 3) | −2.48 | −2.4831 | etiology_module_logfc_cross.csv | yes |
| 16 | Deduplicated KO universe (§3.6) | 132 | 132 | etiology_bootstrap_differences.csv | yes |
| 17 | T1 median \|z\| ctrl vs target (§3.6) | 2.47 vs 1.07 | 2.4702 vs 1.0729 | etiology_bootstrap_differences.csv | yes |
| 18 | T1 Mann–Whitney p (§3.6/Abs/Fig 3) | 0.0019 | 0.001914 | recomputed from CSV | yes |
| 19 | T1 anti-selectivity excl. IFN (§3.6) | p=0.019 | 0.0187 | recomputed | yes |
| 20 | T2 randset-sig controls / total (§3.6) | 14/32 | 14/32 | recomputed | yes |
| 21 | T2 binomial p (§3.6) | 0.0067 | 0.006708 | recomputed | yes |
| 22 | Min q (permutation, 132) (§3.6/Abs) | 0.165 | 0.165447 | etiology_bootstrap_differences.csv | yes |
| 23 | Glial/immune Δ, CI, p (§3.6) | +0.274, 0.203–0.287, 0.0014 | +0.2743, 0.2031–0.2866, 0.0014 | etiology_bootstrap_differences.csv | yes |
| 24 | Largest \|z\| = erythrocyte (§3.6) | +7.39 | +7.388 | etiology_bootstrap_differences.csv | yes |
| 25 | T5 median bootstrap SE (§3.6) | 0.019 | 0.0185 | recomputed | yes |
| 26 | T5 MDE normal / t-based (§3.6) | ≈0.052 / ≈0.10 | 0.0518 / 0.1046 | negative_result_quantification.json | yes |
| 27 | T5 targets > MDE, exact p<0.05 (§3.6) | 41/102, 31/132 | 41/102, 31/132 | recomputed | yes |
| 28 | Median target \|Δ\| (§3.6) | 0.031 | 0.0313 | recomputed | yes |
| 29 | Spearman(network, additive) (§3.4/§3.7/Abs) | 0.991 | 0.99115 | head_to_head_baselines.json | yes |
| 30 | Mean \|reversal\| network/additive/ridge (§3.4) | 0.057/0.033/0.712 | 0.0568/0.0328/0.7119 | head_to_head_baselines.json | yes |
| 31 | Mean % of oracle (§3.4/Table 4) | 20.5 | 20.51 | head_to_head_baselines.json | yes |
| 32 | DE counts adj/unadj (§3.2) | 0/4600/235, 1/2563/958 | identical | run_summary.json | yes |
| 33 | DE median \|logFC\| all / sig D7,D14 (§3.2) | 0.073/0.158/0.119, 0.34, 1.20 | 0.0731/0.1576/0.1187, 0.3398, 1.1994 | DE_D*_vs_Sham_adj.csv | yes |
| 34 | D7 drop: 3,041 genes, Jaccard (§3.2) | 3,041, 0.584 | 3041, 0.584 | de_muscle_drop_sensitivity.csv | yes |
| 35 | Muscle ANOVA F, p (§3.1) | 0.637, 0.605 | 0.637, 0.6054 | run_summary.json | yes |
| 36 | Worst muscle share, sample (§3.1) | 6.26%, GSM4506604 | 6.2573%, GSM4506604 | qc_muscle_scores.csv | yes |
| 37 | Sample correlation range (§3.1) | 0.913–0.994 | 0.9133–0.9935 | qc_sample_correlation.csv | yes |
| 38 | PCA PC1–PC4 (§3.1) | 46.8/17.4/7.2/4.8 | 46.85/17.45/7.23/4.76 | qc_pca_variance.csv | yes |
| 39 | LOO ranges Cd68/Aif1/Gal (Limitations) | [3.28,3.52]/[2.88,3.21]/[−2.69,−2.39] | [3.2795,3.5199]/[2.8821,3.206]/[−2.6902,−2.3865] | etiology_marker_loo.json | yes |
| 40 | Sensitivity min/median ρ (§3.7) | 0.959/0.988, Δ 0.989 | 0.9593/0.9876, 0.9887 | sensitivity_summary.json | yes |
| 41 | DE q<0.05 range across spans (§3.7) | 316–325 | 316–325 | sensitivity_summary.json | yes |
| 42 | Cross-study ρ network vs additive (§3.7) | 0.929 | 0.9286 | sensitivity_summary.json | yes |
| 43 | SC cells; up/down DE (Table 5) | 55,967; 1,517/703 | 55967; 1517/703 | sc_GSE292522_summary.json | yes |
| 44 | Sprr1a logFC, t (Table 5) | +1.55, 54.5 | +1.5539, 54.49 | sc_GSE292522_pseudobulk_DE.csv | yes |
| 45 | Ifit2 / Oas1c q (Table 5) | 0.042 / 0.014 | 0.0417 / 0.0136 | sc_GSE292522_pseudobulk_DE.csv | yes |
| 46 | Gfap logFC, q (Table 5) | +0.004, 0.006 | +0.0036, 0.0061 | sc_GSE292522_pseudobulk_DE.csv | yes |
| 47 | Myeloid fraction Con/B16F10; MWU | 0.12%/0.05%; U=2/14, p=0.114/0.057 | 0.00117/0.000517; 0.1143/0.0571 | sc_GSE292522_summary.json | yes |
| 48 | T3 observed d vs MDE d (Table 5) | 1.27 vs 2.38 | 1.267 vs 2.38 | negative_result_quantification.json | yes |
| 49 | Neuron share 79%→67%; Schwann 19%→29% (§3.8) | both MWU p=0.029 | 78.6%→67.0%; 18.8%→29.3%; p=0.0286 | sc_GSE292522_composition.csv | yes |
| 50 | Compartment Gal / Ifit2 / Isg15 / Ifit1 (§3.8) | +1.47, q<0.001; q=0.046; q=0.007; q=0.004 | 1.472, 1e-4; 0.0459; 0.0071; 0.0037 | sc_GSE292522_celltype_DE.csv | yes |
| 51 | NF-neuron Gal; whole-DRG Gal (§3.8) | +2.54; +1.42 | 2.537; 1.425 | sc_GSE292522_celltype_DE.csv / pseudobulk | yes |
| 52 | 46 myeloid cells, 1 T cell; 7 Gfap UMI (§3.8) | 46, 1, 7 | 46, 1, 7 | composition / gfap_purity_check | yes |
| 53 | Concordance ρ D7/D1/D14 (Table 5) | 0.098/0.112/0.030 | 0.0985/0.1115/0.0299 | sc_GSE292522_concordance_per_day.json | yes |
| 54 | §2.7 floors (§2.7/Limitations) | 3/4,901; 3/71; 2/71 | min paired p=0.0423=3/71; code ⇒ 3/4,901 & 2/71 | marker/module cross + _limma.py:420 | yes |
| 55 | Versions & Python (§2.9) | numpy 2.4.6 … py3.13 | numpy 2.4.6, scipy 1.18.1, pandas 2.3.3, statsmodels 0.15.0, matplotlib 3.11.1, py 3.13.14 | requirements.txt + live env | yes |
| 56 | Test suite "14 groups" (§2.9) | 14 | 14 functions, all pass | scripts/99_unit_tests.py | yes |

**Coverage:** 56 distinct claims checked (≥15 required), spanning Abstract, Methods §2.2–2.9,
Results §3.1–§3.8, Tables 1–5, and figure captions. **0 untraceable, 0 contradicted.** The only
value I could not trace to a *single* artefact line is the §2.7 floor pair `2/71` (shared-control
reversal endpoint); it is defensible from `_limma.exact_two_sided_p` plus the stated
non-antisymmetry argument, but is asserted in the text rather than emitted as a checked number.

---

## 6. Publishability assessment

**Not ready as submitted — Major Revision — but the science and its provenance are sound and the
path to acceptance is short and mechanical.**

Why the provenance dimension is strong:
- The certificate gate runs clean: `287 checks, 0 failed`.
- Table 2 and Table 3 are **genuinely regenerated cell-by-cell** from `etiology_marker_cross.csv`
  and `etiology_module_logfc_cross.csv` (68 cells parsed from the manuscript and compared to the
  CSVs at run time) — this is real regeneration, not pattern-matching. The T1/T2/T5 block likewise
  **recomputes** the Mann–Whitney, binomial and MDE statistics from the deposited CSV and then
  requires the recalculated values to appear in the text.
- 56 independent recomputations all match; no stale or contradicted number survived.
- §2.9 version pins match `requirements.txt` and the live interpreter exactly; the "14-group" test
  suite claim is accurate and the suite passes.
- Figures exist (PNG+PDF), and their captions cite the CSVs their scripts actually read.

Ordered minimum to acceptance:
1. **(Blocking, M1)** Resolve the data-availability statement: publish the repository so the URL
   resolves, or reword to the honest present state and offer on-request access for review. This is
   the only item that must change before I would recommend acceptance.
2. **(Should-fix, m2/m3)** Convert the five under-verified literals (floors `3/4,901`, `3/71`,
   `2/71`; "46 cells"; "eight testable compartments") into run-time-derived assertions, and state
   in the text that the gate anchors ~59% of its checks to a signed-off snapshot (a drift detector)
   rather than recomputing them.
3. **(Nice-to-have)** m4 (compute Fig 3's title p from the CSV), m5 (footnote the adjusted/unadjusted
   double value), m6 (scope the Fig 1 caption).

None of items 2-3 touches any reported number, and item 1 touches none at all. Once the deposit is
honest and reachable, this manuscript's numbers are, to the best of my independent audit, exactly
what it says they are.

---

## Appendix — gate-script audit breakdown (A4, item on `98_check_manuscript_numbers.py`)

| Class of assertion | Runtime count | Share | Nature |
|---|---|---|---|
| Table 2/3 cell regeneration (both sides run-time) | 68 | 24% | fully derived from CSV (strong) |
| T1/T2/T5 recomputation vs JSON | 13 | 5% | expectation recomputed from artefacts (strong) |
| Artefact value vs **hardcoded literal** (`chk`) | ~169 | ~59% | detects artefact drift; expectation frozen |
| Pure string-presence in text (`in_text`) | 37 | 13% | tautological wrt numeric correctness (30 use hardcoded needles; 7 use JSON-derived needles) |
| **Total** | **287** | 100% | matches the gate's own `287 checks, 0 failed` |

Tautologies present: **yes** — the 37 `in_text` checks, 30 of which are hardcoded-needle
string-presence tests. They prove a literal *occurs* in the manuscript, not that the right number
sits in the right place or equals the artefact. Table 2/3 cells, by contrast, are **not**
pattern-matched — they are regenerated from the CSVs.
