# R3-A5 Review — Journal Editor (venue fit, presentation, submission readiness)

**Manuscript:** `docs/manuscript/manuscript_final.md` (read in full)
**Reviewed as:** first submission to a SCIE venue tolerant of exploratory/computational
evidence (Journal of Pain Research / BMC Medical Genomics / BMC Neurology class).
**Independence:** no prior review files, no other reviewers' outputs, no `.workbuddy/`,
no earlier manuscript versions were opened. All judgements below come from the
manuscript text, the artefact tables in `results/tables/`, the scripts I ran myself,
and web verification of five references.

---

## 1. Verdict

**Minor Revision** (from the venue-fit/presentation perspective; this assessment is
complementary to the scientific review). The paper tells one coherent, internally
consistent story, the numbers I checked all reproduce, and the declarations are
essentially complete. But there are **four desk-screening risks that must be fixed
before submission**: (1) the abstract is 457 words — far over every limit in this
venue class; (2) reference 4 as cited **does not exist**; (3) the GitHub repository
in the data-availability statement returns **HTTP 404** — the code/artefacts are
not actually accessible, so the statement as phrased is currently inaccurate;
(4) a computational paper with **zero figures** is a real desk-rejection risk even
at tolerant venues.

---

## 2. Summary of the manuscript

Single-author computational reanalysis of public rat DRG bulk RNA-seq
(GSE149648, Walker-256 tibial bone cancer pain, n=4/group × 4 groups; GSE160543,
paclitaxel/oxaliplatin/vehicle, n=4/group), with an independent single-cell
validation (GSE292522, mouse B16F10 cancer-nerve-injury DRG, 4+4 animals,
55,967 cells). Two claims: (i) a positive one — bone cancer pain and CIPN drive
qualitatively different DRG programmes (myeloid/type-I-interferon in Walker-256
at day 7; neuronal-injury in CIPN), quantified on canonical markers and 22
predefined modules with an exact exhaustive permutation test over C(8,4)²=4,900
split pairs; (ii) a negative/methodological one — in-silico perturbation
(network-diffusion reversal scores) does not support target ranking on these
data: it correlates 0.991 with an additive no-network baseline, ranks an
erythrocyte negative control above every analgesic target (|z|=7.39), and is
anti-selective (negative controls carry larger |z| than targets, median 2.47 vs
1.07, MWU p=0.0019) while still having adequate per-test power (MDE≈0.105;
largest target Δ=0.274 = 2.6× MDE; 31/132 exact p<0.05). Structure: Intro,
Methods (2.1–2.9), Results (3.1–3.8), Discussion, Limitations, Conclusions,
five tables, 16 references, full declarations.

---

## 3. Major issues

### 3.1 Abstract is 457 words — over every limit in the target venue class

- **Citation:** `manuscript_final.md:17–23` (Abstract block).
- **Evidence (counted precisely, including the `Background./Methods./Results./
  Conclusions.` labels):** **457 words** (453 excluding the four labels). Limits:
  Journal of Pain Research (Dove Press) structured abstract **≤300**; BMC
  Medical Genomics / BMC Neurology structured abstract **≤350**. 457 fails both.
  The Results paragraph alone (~230 words) is a compressed second Discussion.
- **Required fix:** cut to ≤300 words. The strongest cut candidates: the
  muscle-contamination parsing detail (2.2-level), the permutation-floor
  arithmetic (6.1×10⁻⁴ / 3/71 belongs in Methods), and the second half of the
  Results paragraph (the power argument can be one clause: "and the layer is
  adequately powered — its failure is validity, not deafness").

### 3.2 Reference 4 does not exist (fabricated/garbled citation)

- **Citation:** `manuscript_final.md:216` — "Falk S, Dickenson AH. Pain and
  nausea in cancer: mTOR-mediated mechanisms. *Br J Anaesth*. 2015;115(2):169–171."
- **Evidence:** Europe PMC full-text/keyword search for the exact phrase
  "Pain and nausea in cancer" returns no article with that title; an
  author-restricted query (AUTH:"FALK S" AND JOURNAL:"Br J Anaesth" AND
  PUB_YEAR:2015) returns **0 hits**. No Falk & Dickenson paper sits at
  BJA 115(2):169–171. The real, content-equivalent citation is Falk S,
  Dickenson AH. *Pain and nociception: mechanisms of cancer-induced bone pain.*
  J Clin Oncol. 2014;32(16):1647–1654. doi:10.1200/JCO.2013.51.7219 (verified),
  which supports exactly the claim the text hangs on ref [4] (osteoclast
  activation, tumour mediators, innate immune response in bone cancer pain).
- **Required fix:** replace ref 4 with the verified JCO 2014 citation (or
  another verifiable BJA/Expert Opin review). This is a credibility kill-shot
  at any journal if left in; a copyeditor or reviewer who checks it will treat
  the whole reference list as suspect.

### 3.3 Data-availability statement points to a dead URL; "deposited" is currently not true

- **Citation:** `manuscript_final.md:209`.
- **Evidence:** `https://github.com/yyx-4113/cancer-pain-drg-virtual-perturbation`
  returns **404 / "Page not found"** (checked 2026-10-10). The statement reads
  "Analysis code, artefacts and the test suite **are deposited** in the project
  repository…, to be archived at Zenodo upon acceptance." The Zenodo clause is
  honest; the "are deposited" clause is not, since reviewers cannot reach the
  repo. For this venue class (computational reanalysis whose entire value is
  reproducibility — the paper itself advertises a 276-check number gate),
  the code link is checked at submission.
- **Required fix (all of):** (a) create the public GitHub repo and push the
  scripts + artefact tables + test suite **before submission**; (b) deposit a
  Zenodo snapshot and cite the DOI in the statement (many venues now want it at
  submission, not "upon acceptance"); (c) if the repo is currently private,
  the wording is still false for reviewers — make it public or reword to
  "will be made available upon publication", which weakens the paper.

### 3.4 No figures at all — a desk-screening risk at this venue class

- **Citation:** whole manuscript; five tables (`manuscript_final.md:234–319`),
  zero figure blocks.
- **Evidence:** The two headline results are inherently visual: the
  aetiology divergence (Table 2/3 marker/module contrasts) and the
  anti-selectivity result (|z| distribution of targets vs controls, the
  erythrocyte outlier). Both venues named in the brief expect at minimum a
  study design/data-flow figure, a module × time-course heatmap or PCA, and a
  target-vs-control |z| (or reversal-Δ) figure. A computational paper whose
  only visual material is five markdown tables reads as unfinished at
  first screening, and Figure-less submissions are frequently returned
  before review even at tolerant journals.
- **Required fix:** add ≥3 figures: (F1) study overview / two-questions-two-
  answers schematic; (F2) marker and module logFC heat-map across
  Bone-D1/D7/D14 and PTX/OXA (Table 2/3 content); (F3) the anti-selectivity
  evidence — |z| strip-chart or ECDF, targets vs controls, erythrocyte
  highlighted, plus the 0.991 network-vs-additive scatter inset. All are
  one-line renders from the deposited CSVs.

---

## 4. Minor issues

1. **Title is 24 words / ~200 characters** (`manuscript_final.md:1`). Within no
   hard limit but above the ~20-word style norm of Dove and BMC venues; the
   subtitle "an exact cross-study reanalysis and a negative-control calibration
   of in-silico perturbation" is heavy. Consider e.g. "Divergent DRG programs in
   bone cancer pain and CIPN, and a negative-control calibration of in-silico
   perturbation" (17 words).
2. **Overlong paragraphs.** Recomputed paragraph lengths (words):
   §3.8 first paragraph **385**, Limitations **349** (single paragraph for the
   entire section), §2.8 second paragraph **297**, §3.6 third paragraph 253.
   The §3.8 and Limitations blocks each carry 4–5 distinct points; split at
   each "First/Second/Two further caveats" seam. Limitations should be a
   short list or 3 paragraphs (data/batch; perturbation layer; single-cell
   scope + species).
3. **Redundancy between Abstract, §3.6 and Limitations.** The KO/OE mirror
   (r = −1.000), MDE ≈ 0.105, Δ = 0.274 "exceeds it 2.6-fold", and "31 of 132
   tests reach exact permutation p < 0.05" appear in full three times
   (`manuscript_final.md:21, 137–143, 177`), and the anti-selectivity numbers
   (2.47 vs 1.07, p = 0.0019) twice more in Conclusions (`181`). The message
   consistency is good (see §5 table — all copies agree), but state the power
   argument once in §3.6 and cross-reference it from Limitations
   ("not underpowering; see §3.6") instead of repeating all four numbers.
4. **Table captions are terse and two tables are not fully self-explanatory.**
   Table 1 caption does not define "n" (genes per module) or "mean t";
   Table 3 uses "—" cells with no footnote (untested contrast? not
   applicable?) and its q column does not say which contrast each q belongs
   to when several contrasts are shown; Table 5's "U = 2/14" needs "(U statistic
   / n pairs)" gloss. Add one footnote line to each.
5. **CRediT role "investigation" is contestable for a pure secondary analysis**
   (`manuscript_final.md:195`). All 14 CRediT roles are claimed by one author;
   "investigation" normally denotes conducting experiments/recruiting
   subjects. Harmless, but a copyeditor at JPR/BMC may query it; consider
   dropping to the 10 roles that apply. Formatting itself is correct:
   single author, no co-authors, no degree suffixes (MD/PhD absent), ORCID
   and correspondence line present — consistent with the author profile.
6. **Reference-list formatting is uneven.** Only 4 of 16 references carry DOIs
   (refs 9, 10, 15, 16); Dove and BMC both want DOIs on all. Ref 10's author
   list is silently truncated ("Zhai M, Yang S, Lin S, Zhu H." — the real list
   has 7 authors: …, Xu L, Liao H, Song XJ; verified via Europe PMC/PMCID
   PMC8134751) and must end "et al." Ref 15 cites "Baruch EN, et al." while
   ref 10 gives four names — pick one Vancouver convention. 16 references is
   thin for the class (typical 35–60); the in-silico-perturbation background
   in particular is carried by only refs 7–9.
7. **Artefact hygiene (not a manuscript error):**
   `results/tables/negative_result_quantification.json`, field
   `T5_reversal_score_power.interpretation`, says "26/132 exact permutation
   p<0.05" while the same object's `n_exact_permutation_p_lt_0.05` is **31**
   and the manuscript says 31 (correct — I recomputed 31). Fix the
   interpretation string before deposit, or a reviewer who reads the JSON will
   suspect cherry-picking.
8. **Abstract slightly overstates §3.8.** Abstract: "The interferon program is
   **not reproduced** in the B16F10 model." §3.8 itself is careful: the
   compartment test "cannot exclude a myeloid-restricted interferon programme"
   (myeloid ≈ 46 cells, untestable). Add "at pseudobulk resolution" to the
   Abstract sentence.
9. **Em-dash/spacing inconsistency:** Limitations uses "model—a different
   species" (unspaced em dash, line 177) while the rest of the manuscript uses
   spaced dashes. Also "sequencing platform" appears in both §2.1 and the
   Abstract with slightly different framing — fine, but check at copy-edit.

---

## 5. Numbers I independently recomputed

All from `results/tables/etiology_bootstrap_differences.csv` (deduplicated KO
universe: `network_scope=='per_contrast' & perturbation=='knockout'`),
`negative_result_quantification.json`, `etiology_marker_cross.csv`, and
`scripts/98_check_manuscript_numbers.py` (run by me: **"276 checks, 0 failed"**;
its expectations are runtime-derived from the CSV/JSON via `read_csv`/`jload`
and compared to strings in the manuscript file — not hardcoded, see
`scripts/98_check_manuscript_numbers.py:76, 204, 325, 343–376, 405`).

| Claim | Manuscript | Recomputed | Match? |
|---|---|---|---|
| Abstract words | — | 457 (453 excl. labels) | over 300 limit |
| Keywords count | 6 | 6 | OK |
| Running head chars | — | 48 (≤50) | OK |
| Body words to References | (brief said ~5,800) | 5,785 (Intro→Conclusions, excl. headings/tables); +457 abstract = 6,242 | confirmed |
| T1 median \|z\| control vs target | 2.47 vs 1.07 | 2.4702 vs 1.0729 | ✔ |
| T1 MWU two-sided p | 0.0019 | 0.001914 (U = 2102; 102 target, 30 control) | ✔ |
| Largest \|z\| = erythrocyte | 7.39 | 7.3878 (CTRL erythrocyte) | ✔ |
| T2 randset-significant | 32/132, 14 controls (44%) vs base 23%, binomial p = 0.0067 | 32 sig; 14 control (43.75%); base 30/132 = 22.73%; one-sided binomial p = 0.006708 | ✔ |
| T5 median bootstrap SE | 0.019 | 0.0185 | ✔ |
| T5 MDE (unpaired n=2+2, 80%, two-sided) | ≈0.105 | 5.653 × 0.0185 = 0.1046 | ✔ |
| T5 largest target Δ | 0.274, 2.6× MDE | 0.2743 (Glial/immune, Bone_D7 vs PTX); 2.62× | ✔ |
| Exact-permutation p<0.05 count | 31 of 132 | 31 (JSON field; JSON prose "26" is stale) | ✔ |
| Table 2 q for Cd68/Aif1/Gal | 0.012 | 0.012487 (etiology_marker_cross.csv) | ✔ |
| Targets' percentile of effect distribution | 55.7th | 55.7 (JSON T4; spot-checked) | ✔ |
| §3.6 vs Abstract vs Limitations vs Conclusions consistency | — | all four sites quote identical T1/T2/T5 values; no contradictions found | ✔ |

Reference spot-checks (5): **ref 4 — FAIL (nonexistent; see Major 3.2)**;
ref 9 verified (Nat Methods 22(8):1657–1661, 2025, doi correct);
ref 10 verified (Front Mol Neurosci 14:665085, 2021, doi correct; author list
truncated); ref 15 verified (Nature 646:462–473, 2025, doi correct; and the
manuscript's description of Baruch et al.'s neuron-intrinsic type-I IFN finding
is fair, including the acknowledged conflict in §3.8); ref 16 verified
(Nat Neurosci 18:145–153, 2015, doi correct).

---

## 6. What would make this publishable (concrete, ordered)

1. **Replace reference 4** with Falk & Dickenson, J Clin Oncol 2014;32(16):
   1647–1654 (or equivalent verified review). Add DOIs to all 16 references;
   fix ref 10's author list ("… Song XJ, et al." per chosen convention).
2. **Cut the abstract to ≤300 words** (report the new count in the cover
   letter). Keep the three load-bearing numbers: ΔCd68 +3.42/q 0.012,
   anti-selectivity p = 0.0019, and one power clause.
3. **Make the data/code real:** public GitHub repo live at submission
   (scripts + `results/tables/` + test suite + `requirements.txt`), Zenodo DOI
   minted at submission, and reword the statement to cite the DOI. Fix the
   stale "26/132" string in `negative_result_quantification.json` first.
4. **Add 3 figures** (study schema; aetiology divergence heatmap; target-vs-
   control |z| with erythrocyte highlighted + additive-baseline scatter), each
   generated from the deposited CSVs so the number gate can cover them.
5. **Split the two 350+ word paragraphs** (§3.8 opener; Limitations) and
   de-duplicate the power argument (state once in §3.6, cross-refer from
   Limitations).
6. **Footnote Tables 1/3/5** (define n, "—", U = 2/14, and which contrast each
   q belongs to).
7. Optionally shorten the title to ≤20 words and trim CRediT to applicable
   roles. Declarations otherwise complete (Funding, COI, Acknowledgements,
   Generative-AI declaration with appropriate scope and author
   responsibility, Data availability) — nothing fabricated, no invented
   co-authors, no degree suffixes.

With items 1–4 done, this is a credible submission to Journal of Pain Research
or a BMC-series journal of that class; the negative-result layer, exact permutation
discipline, and the 276-check reproducibility gate are genuine differentiators
that fit these venues' stated appetite for neutral/computational evidence.

---
*R3-A5, 2026-10-10.*
