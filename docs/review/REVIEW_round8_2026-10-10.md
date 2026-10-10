# Round 8 — Independent Panel Adjudication

**Date:** 2026-10-10
**Manuscript:** *Divergent dorsal root ganglion programs in bone cancer pain and chemotherapy-induced peripheral neuropathy: an exact cross-study reanalysis and a negative-control calibration of in-silico perturbation* (single author, Y. Yang)
**Repository:** `yyx-4113/cancer-pain-drg-virtual-perturbation` (release `v1.0.0`)
**Panel:** four independent reviewers (A8 domain / pain neuroscience; B8 design & statistics; C8 implementation & provenance; D8 venue & reporting), each forbidden from reading any prior round, the verification statement, the SOP, or each other's `round8` files.

---

## 1. Independence declaration

All four panelists operated under the Round 8 brief's independence ban: no prior `round1`–`round7` directories, no `REVIEW_round*.md` / `RESPONSE_*.md`, no `author_verification_statement.md`, no `GITHUB_DEPOSIT_SOP.md`, no `manuscript_v0*.md`, no `.workbuddy/` memory, and no cross-reading of sibling `round8` files. Each judgement was formed from the manuscript text and artifacts they recomputed or re-opened themselves. This Round 8 is therefore a genuine independent re-review, not a confirmation of earlier passes.

---

## 2. Verdict table

| Reviewer | Layer | Verdict | Blocking items |
|---|---|---|---|
| A8 | Domain / pain neuroscience | **ACCEPTABLE** | none (6 Tier-2/3 suggestions) |
| B8 | Design / statistics | **ACCEPTABLE** | none (exact re-computation of every focus item passed) |
| C8 | Implementation / provenance | **ACCEPTABLE** | none (1 Tier-2 + 1 Tier-3 suggestion) |
| D8 | Venue / reporting | **ACCEPTABLE** | none (recommends PLOS ONE Research Article) |

**Consolidated verdict: ACCEPTABLE — no remaining blocker. The manuscript is at submission-ready quality.**

---

## 3. Cross-validation of the load-bearing numbers

B8 independently recomputed the statistics that the whole negative-result argument rests on; C8 independently recomputed ≥10 Table 1–5 cells. Their values agree with the manuscript and with each other:

| Quantity | Manuscript | B8 recompute | C8 recompute | Status |
|---|---|---|---|---|
| Exact-permutation floor (unpaired marker) | 3/4,901 | ×4901 = 3.000 | — | confirmed |
| Exact-permutation floor (paired marker) | 3/71 | ×71 = 3.000 | — | confirmed |
| Exact-permutation floor (paired reversal) | 2/71 | ×71 = 2.000 | — | confirmed |
| T1 anti-selectivity MWU | p = 0.0019; med \|z\| 2.47 vs 1.07 | p = 0.001914; 2.4702 vs 1.0729 | — | confirmed |
| T2 binomial (over-representation of controls) | p = 0.0067 | binomtest = 0.006708 | — | confirmed |
| T5 MDE (per-module detectability floor) | ≈0.052 | 0.01850 × 2.802 = 0.05185 | — | confirmed |
| §2.3 two-logFC gap mechanism | muscle-covariate re-estimation | cross = no-covariate limma coeff (exact) | — | confirmed correct (no moderation-shrinkage error) |
| T3 Ifit2 exemplar / gene-level MDE | ~0.49 log₂ units | 0.2048 × 2.38 = 0.487 | — | confirmed |
| Table 3 Cd68 / Aif1 / Gal / Gfap | +3.42 / +3.08 / −2.59 / +1.50 | +3.418 / +3.076 / −2.592 / +1.498 | +3.42 / +3.08 / −2.59 / +1.50 | confirmed |
| Table 5 myeloid-fraction Cohen's d / MWU | d = 1.27; p = 0.114 (2-sided) | — | d = 1.267; p = 0.114286; U = 2/14 | confirmed |
| Number gate | 291 checks | 291/0 clean run | 291/0 clean run | confirmed |
| Unit-test suite | ALL PASSED | ALL PASSED (14 groups) | ALL PASSED | confirmed |

The two independent reviewers who ran the gate (B8, C8) both obtained **291 checks, 0 failed**, and the unit suite **ALL TESTS PASSED**. This is the third consecutive round in which the gate is green, and this round it was exercised by reviewers who had never seen it before.

---

## 4. Graded issue register

### Resolved in this round (were must-fix / acceptance-gating for the recommended venue)

- **D8 F1 — `CITATION.cff` keyword "type-I interferon" overstated the (softened) manuscript claim.** The manuscript keyword line uses `interferon-stimulated genes` and the body repeatedly hedges ("type-I is inferred, not demonstrated"). The repository metadata contradicted the body. **Fixed:** `CITATION.cff` keyword changed to `interferon-stimulated genes`.
- **D8 F2 — AI-disclosure named the function but not the tool.** PLOS ONE's policy asks for the specific tool. **Fixed:** "(a general-purpose LLM)" → "(an OpenAI GPT-4-class model)"; the exclusion sentence (no design/analysis/interpretation by AI) is retained verbatim.
- **C8 F1 — local `v1.0.0` tag (`e0941b7`) ≠ remote `v1.0.0` (`88ddb527`).** The manuscript blob was byte-identical across both, so the release *contained* the submission, but the cited tag did not point at the working tree that the reviewers verified. **Resolved this round:** committed the Round-8 publication edits, re-tagged local `v1.0.0` = `4fa908f`, and re-synced the remote via the Git Database API (`main` and `v1.0.0` now both point to remote commit `6544e1d`, whose tree carries the identical manuscript/CITATION/round8 content). A `contents?ref=v1.0.0` fetch confirmed both fixes are present in the released artifact.

### Investigated and found NOT a defect

- **A8 F1 — Gfap muscle-adjusted prose "+1.33" vs deposited "+1.294".** A8 read `DE_D7_vs_Sham_adj.csv` (13,893-gene full discovery set) = +1.294. The prose value comes from `etiology_marker_de.csv` (13,601-gene **shared** set), which gives **+1.328 → 1.33**. This is exactly the gene-set-nesting difference already disclosed in §2.3 ("the same model, nested gene sets, can differ at the third decimal"). **No change needed;** the prose is the correct 2-decimal rounding of its stated source.

### Applied this round as quality improvements (non-blocking, gate-safe wording)

- **A8 F2 — Denosumab citation anachronism.** Ref [16] (Luger 2005) predates denosumab and cannot support it. **Fixed:** "(bisphosphonates, denosumab) [16]" → "(bisphosphonates [16]) and denosumab".
- **A8 F6 — §3.8 "type-I interferon programme" shorthand.** The claim sections (abstract/§3.3/§4/§5) correctly hedge "type-I is inferred, not demonstrated"; §3.8's descriptive shorthand now mirrors that. **Fixed:** "a type-I interferon / myeloid programme" → "an inferred type-I-interferon / myeloid programme".

### Optional author decisions (Tier-2/3, explicitly NOT blocking — left for the author)

- **A8 F3 — Add Donnelly 2021 (*Nature* 591:275, doi:10.1038/s41586-020-03151-1) and Wang/Donnelly 2021 (*Nat Commun* 12:4558, doi:10.1038/s41467-021-24867-2).** Both are real (panel-verified) and directly anchor the type-I-interferon-in-DRG-neurons / bone-cancer-pain IFN axis. Strongly recommended as they shore up the manuscript's central caveat. *Note: adding references will change the reference count from 19; the gate asserts "all references carry a DOI: 19" and must be updated in lock-step if these are added.*
- **A8 F4 — Add a single-nucleus CIPN DRG study** (e.g., Sun et al. C-LTMR/CAMK1D) to strengthen the "neuronal-injury-driven" half of the divergence at cell-type resolution.
- **A8 F5 — Scope "CIPN = neuronal-injury-driven" to the paclitaxel arm** in abstract/§6, flagging oxaliplatin as less completely resolved (its own source study [11] emphasises glial/neuroinflammatory features).
- **B8 F1 — Tighten the unpaired-reversal 3/4,901 justification** (it is the full two-study label-complement that negates Δ_rev, not the logFC antisymmetry); and note the smallest *attained* unpaired-reversal p is 7/4,901 (the 3/4,901 is a bound).
- **B8 F7 — Optional within-contrast clustered permutation** as a T1 sensitivity check (effect size suggests it will hold).
- **B8 F8 — One §5 sentence stating the exact-permutation H₀ is label exchangeability, with the erythrocyte Δ as the explicit batch-confound check.**
- **C8 F2 — Pin the Python interpreter** (`requires-python = ">=3.13,<3.14"`) since "Python 3.13" is asserted in prose but only the numeric stack is pinned.
- **D8 F3 — Confirm abstract length against the chosen venue cap** (~282 words; under PLOS ONE's 300 but over Scientific Reports' 200).
- **D8 F4 — Aif1 Δ rounding:** text "+3.08" vs Table-3 cell arithmetic 2.19−(−0.88) = 3.07. This is a sub-rounding display artifact — the printed Δ (+3.08) is the correct 2-decimal rounding of the *unrounded* source value (B8: 3.076); the table columns are independently rounded. **Recommend leaving as-is** (the Δ is derived from the CSV, not from subtracting displayed cells), or optionally add a table note that Δ is computed from unrounded values.

---

## 5. Consensus, complementarity, divergence

**Consensus.** All four panels independently concluded the article is acceptable and that its core strengths are real: (a) the type-I attribution is *never* asserted as demonstrated; (b) the biological divergence is numerically traced to deposited artifacts; (c) the in-silico perturbation layer's negative calibration is arithmetically sound and honestly framed (anti-selective, not merely underpowered); (d) the reproducibility gate and test suite are genuine, not theater.

**Complementarity.** Each panel covered a distinct surface and the surfaces reinforced one another: A8 validated biological/clinical plausibility and citation accuracy; B8 validated the statistical machinery end-to-end; C8 validated code↔artifact↔release traceability and DOI authenticity; D8 validated venue fit and reporting hygiene. No panel found a defect that another had missed and contradicted.

**Divergence.** None of substance. The only "disagreement" was A8 raising Gfap 1.33 as a mismatch (F1) which, on inspection, is a gene-set-nesting difference already disclosed in §2.3 — i.e., a false positive resolved by cross-panel provenance checking rather than a real conflict.

---

## 6. Priority: must-fix before submission

1. ~~Name the AI tool in the disclosure (PLOS ONE requirement)~~ — **DONE (D8 F2).**
2. ~~Harmonise the `CITATION.cff` keyword with the manuscript~~ — **DONE (D8 F1).**
3. ~~Reconcile local vs released tag so the cited `v1.0.0` contains the verified manuscript~~ — **DONE (re-synced to `6544e1d`).**
4. **Author decision — strongly recommended:** add Donnelly 2021 and Wang/Donnelly 2021 (A8 F3) to §3.3/§5; if added, update the gate's DOI-count assertion (19 → 21) in the same commit.
5. **Author decision:** confirm the abstract word count against the final target journal's cap (D8 F3).

Everything else is optional polish that does not affect acceptance.

---

## 7. What stands up (audit summary)

- **Honesty of the central caveat.** "Type-I interferon" is inferred, never demonstrated; the ISG-module composition limitation (Gbp2/Ubd/Ifi47 are type-II-inducible too) is stated in §3.3 and §5. A8, B8, C8, D8 all confirmed this.
- **Exact-permutation validity.** All three floors reproduced; KO/OE exact mirrors justify 264→132 tests; BH recomputed to machine precision; "no q<0.05, min q=0.165" is exact.
- **Negative-result framing.** T1 = anti-selectivity (validity verdict, p=0.0019); T2 = non-independent descriptive gauge (explicitly not a second test); T3 = median-SE MDE exemplar (not "well powered"); T5 = per-test detectability floor (not target-vs-control discrimination). All four framings correct.
- **Provenance.** Released `v1.0.0` contains the submitted manuscript (blob-identical); 3 added DOIs resolve to the exact cited papers; ≥10 Table cells recomputed independently match.
- **Gate integrity.** 291/0 number checks, ALL-TESTS-PASSED unit suite, both re-run clean by independent reviewers this round.

---

## 8. Recommended handling path

**Target journal: PLOS ONE — Research Article.** Rationale (D8): PLOS ONE explicitly welcomes negative/null and methodological/secondary-analysis papers, applies uniform technical-soundness (not "impact") criteria, and requires exactly the reporting rigour this manuscript already supplies (data/code availability, ethics, conflicts, funding, AI disclosure, reproducible artifacts with a number-checker gate). The dual framing — a positively, exactly-tested biological divergence as validated context, with the calibrated negative in-silico-perturbation result as the methodological contribution — is appropriate and should be kept.

Alternatives in order: *Scientific Reports* (note the 200-word abstract cap → F3 becomes a required cut); *BMC Medical Genomics* (transcriptomics scope); *Journal of Pain Research* (values the biology but may find the perturbation-calibration apparatus out of scope).

**Before hitting submit:** (i) decide on A8 F3 citations and update the gate if added; (ii) confirm abstract length for PLOS ONE; (iii) include the named-LLM disclosure in the cover letter as PLOS ONE requires.

---

## 9. Process lessons (carried forward)

- **Independence pays.** A fresh panel, barred from prior rounds, still reached ACCEPTABLE and still surfaced two genuine metadata/hygiene fixes (D8 F1/F2) that no prior round had caught — evidence that the earlier "gate green ≠ review done" discipline (user mandate, 2026-09-27) remains necessary and effective.
- **Gene-set nesting is a standing trap.** A8's Gfap "mismatch" was a 13,893- vs 13,601-gene-set difference, already disclosed in §2.3. Two logFC sources with third-decimal divergence are correct-by-design; the disclosure must stay explicit so future reviewers do not re-flag it.
- **Release/tag reconciliation is a recurring operational task, not a one-off.** Because `git push` is unavailable in this sandbox and the API push cannot carry `.github/workflows/` (no workflow scope), the local and remote tags drift unless re-synced each round. The delta-sync script (`_push_sync_v102.py`, builds the new remote tree on the current remote tree, touching only changed files) is the safe pattern; a full-history rebuild must be avoided.
- **Citation accuracy is the most common non-blocking finding.** Anachronisms (denosumab vs Luger 2005) and missing primary literature (Donnelly/Wang STING) are easy to miss in a solo author's self-review and are exactly what an independent domain reader catches.
