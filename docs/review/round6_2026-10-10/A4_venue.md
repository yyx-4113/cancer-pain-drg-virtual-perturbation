# A4 — Venue / Reporting-Standards Review

**Manuscript:** *Divergent dorsal root ganglion programs in bone cancer pain and chemotherapy-induced peripheral neuropathy: an exact cross-study reanalysis and a negative-control calibration of in-silico perturbation*
**Author:** Yongxin Yang (single author)
**Review round:** 6 (fresh independent review; no prior rounds seen)
**Reviewer role:** Venue / reporting-standards editor
**Files reviewed:** `manuscript_final.md`, `README.md`, `CITATION.cff`, `requirements.txt`, `author_verification_statement.md`, `docs/review/` filename listing only.

---

## Summary verdict (venue layer)

The manuscript is internally well-structured and the reporting is, on the whole, more disciplined than most single-author computational submissions: the exact-permutation design, the negative-control calibration, and the runtime number-checker ("the gate") are genuine strengths. At the venue/reporting layer, however, there are **one hard numeric inconsistency** (the 0.052 vs 0.05 minimum-detectable-Δ restatement), **two moderate cross-reference drifts** (the dual-logFC gap on Gfap/Ccl2 exceeds the paper's own "few hundredths" claim; the 4,900 vs 4,901 permutation denominator is internally contradictory), **one under-disclosed AI statement**, and **one venue-dependent data-availability gap** (GitHub release tag with no DOI-archived snapshot). None is individually fatal, but the AI-disclosure and DOI-archive items will be hard blockers at several of the candidate journals. Specific four-part findings follow.

---

## 1. Version consistency

### V1 — Software version strings are consistent across manuscript, requirements.txt and README (verified, no mismatch)
【Problem】 No mismatch exists between the version strings stated in §2.9 and those in `requirements.txt` and `README.md`; `CITATION.cff` omits the numerical stack, which is normal.
【Evidence】 `manuscript_final.md:101` ("numpy 2.4.6, scipy 1.18.1, pandas 2.3.3, statsmodels 0.15.0 and matplotlib 3.11.1 (`requirements.txt`)"); `requirements.txt:12-16` (identical `==` pins); `README.md:114` (identical pinned block); `README.md:183` ("Python 3.13"); `CITATION.cff` carries no package versions (only `version: 1.0.0`, line 20). Python 3.13 is stated in §2.9 (implied) and README, not in `requirements.txt` — expected, since `requirements.txt` pins libraries not the interpreter.
【Why it matters】 Journals (e.g., Nature Portfolio, PLOS, BMC) routinely diff the Methods environment block against the deposited `requirements.txt`; an exact match here removes a common desk-reject trigger. `CITATION.cff` not repeating the stack is correct practice and is not a defect.
【Specific fix】 No change required. Keep §2.9, `requirements.txt`, and `README.md` in lock-step; optionally add `python==3.13` to `requirements.txt` (or a `RuntimeRequirements` note) so the interpreter version is pinned in one machine-readable place as well.

### V2 — The pinned versions must be confirmed installable from PyPI before submission (soft, venue-relevant)
【Problem】 The pinned stack (scipy 1.18.1, matplotlib 3.11.1, numpy 2.4.6) should be verified to resolve as real published releases; if any are not real PyPI versions the environment is non-reproducible and a reviewer/editor will flag it.
【Evidence】 `manuscript_final.md:101`; `requirements.txt:12-16`. The author-verification statement (item 4) asserts the versions are pinned but does not assert they resolve on PyPI.
【Why it matters】 A fabricated or typo'd version (e.g., a scipy minor that does not exist) is exactly the kind of claim a reproducibility-checking editor verifies, and a non-installable `requirements.txt` invalidates the "every number was produced with exactly these versions" claim that underpins the whole verification package.
【Specific fix】 Paste-ready pre-submission step: `python -m pip install -r requirements.txt` in a clean venv and copy the resolved `pip freeze` into a `requirements.lock` / `environment.yml`, or add a line to `README.md`: "All five pins resolve on PyPI as of 2026-10-10 (verified by `pip install -r requirements.txt`)."

---

## 2. Reporting structure

### R1 — The manuscript names no target journal, which frustrates a venue-fit assessment
【Problem】 Nowhere does the manuscript, README, or `CITATION.cff` state an intended journal; `CITATION.cff:44` and `README.md:6` say only "Manuscript in preparation."
【Evidence】 `CITATION.cff:44` (`journal: "Manuscript in preparation"`); `README.md:6` ("Manuscript in preparation, 2026"); `manuscript_final.md` has no "Prepared for / Submitted to" line.
【Why it matters】 Reporting standards (structured abstract, word limits, data-availability mandate, AI-disclosure format) differ sharply by venue; an editor cannot judge compliance without a named target, and several candidate journals have contradictory requirements (see D1, A1).
【Specific fix】 Add to the cover letter and `CITATION.cff` a `note` field: e.g., "Target journal: PLOS ONE (or Scientific Reports). Article type: Original Research." This also lets the editor confirm the structured abstract and negative-result framing are on-spec.

### R2 — The title is over-long for most candidate venues
【Problem】 The title runs ~200 characters / ~20 words and will breach the hard caps at several pain/genomics journals.
【Evidence】 `manuscript_final.md:1` (full title: "Divergent dorsal root ganglion programs in bone cancer pain and chemotherapy-induced peripheral neuropathy: an exact cross-study reanalysis and a negative-control calibration of in-silico perturbation").
【Why it matters】 *Scientific Reports* caps titles at 20 words / ~200 chars (borderline fail); *Journal of Pain Research* and *BMC Medical Genomics* prefer ≤ 150 chars; *PLOS ONE* "recommends" concise titles. An over-long title is a routine production-stage revision request and signals the piece was not scoped to a venue.
【Specific fix】 Offer a ≤ 150-char alternative, e.g.: "Bone cancer pain and CIPN drive divergent DRG programs, but in-silico perturbation fails negative-control calibration." Keep the current long form only as the running/short title if the venue allows.

### R3 — The running head is within limit but tight; confirm against the chosen venue
【Problem】 The running head is ~48 characters, just under the common ≤ 50-char cap, but several journals use a stricter ≤ 40-char running head.
【Evidence】 `manuscript_final.md:3` ("Divergent DRG programs: bone cancer pain vs CIPN").
【Why it matters】 A 48-char running head passes *PLOS ONE* (≤ 50) but may be truncated at *Scientific Reports* / *BMC* (≤ 40–45). This is a trivial but real production-blocker if the venue is chosen after acceptance.
【Specific fix】 Provide a ≤ 40-char fallback: "Divergent DRG programs in bone vs CIPN pain" (38 chars), and state the chosen venue's head-length rule in the cover letter.

### R4 — Article type / negative-result framing is sound but should be declared explicitly as a "negative-results / methods" original article
【Problem】 The piece is an original research article combining a cross-study reanalysis (positive biological finding) with a calibrated negative result (in-silico perturbation adds nothing). The framing is honest, but the article type is left implicit.
【Evidence】 `manuscript_final.md` IMRaD structure (§1–6), Abstract Background/Methods/Results/Conclusions (lines 15–25), and the explicit negative-control framing in §3.6, §4, §5, §6.
【Why it matters】 The negative-result emphasis is a feature, not a defect, and is best placed at a journal that explicitly welcomes null/validation results: **PLOS ONE** (null results explicitly in scope; method articles accepted), **Scientific Reports** (broad original research), **BMC Medical Genomics** (if framed on the genomic/transcriptomic method), or **Journal of Pain Research** (pain-focused; welcomes validation/negative studies). Naming the type ("Original Research — negative calibration") prevents an editor from misrouting it as a "correspondence" or "letter" and soft-rejecting on scope.
【Specific fix】 In the cover letter state: "Article type: Original Research. We report (i) a positive cross-aetiology DRG finding and (ii) a rigorously calibrated negative result on in-silico perturbation, both of which fall within your journal's original-research scope for computational/translational studies." Emphasise the negative-result calibration in the title/abstract so it is not read as a failed primary hypothesis.

---

## 3. AI disclosure

### A1 — The generative-AI declaration under-discloses relative to standard journal policies
【Problem】 The declaration states an LLM was used for "code debugging and manuscript language editing" but does **not name the tool** and does **not explicitly confirm** that AI was absent from analysis design, statistical interpretation, and result selection — the two points most journals now require.
【Evidence】 `manuscript_final.md:211` ("During the preparation of this work the author used a large language model assistant for code debugging and manuscript language editing. All analysis code was executed and verified locally, every numerical claim was checked against the deposited artefacts, and the author reviewed and edited all content and takes full responsibility for the published article.").
【Why it matters】 *PLOS ONE*, *Nature Portfolio*, *Elsevier*, and *Wiley* AI-disclosure policies require: (a) the **name** of the tool (e.g., ChatGPT, Claude, Gemini, Copilot); (b) a statement that AI was **not** used to generate or interpret the data/results; (c) confirmation that a human takes responsibility. The current text covers (c) and partially (b) but omits (a) entirely, and "manuscript language editing" is vague enough that an editor may request a fuller statement. Given the heavy computational/statistical content, under-disclosure is a real acceptance risk.
【Specific fix】 Paste-ready replacement: "During the preparation of this work the author used a large language model (Claude, Anthropic) solely for code-debugging assistance and for editing the language of the manuscript draft. The LLM was not used to design the analysis, select statistical methods, generate results, or interpret findings. All analysis code was executed and verified locally, every numerical claim was checked against the deposited artefacts by the author, and the author reviewed and edited all content and takes full responsibility for the published article."

---

## 4. Data & code availability

### D1 — GitHub release tag is citable but lacks a DOI-archived snapshot; this is a blocker at several candidate venues
【Problem】 The Data & code availability statement points only to a GitHub repository at "release v1.0.0" with no DOI-archived, version-immutable copy (e.g., Zenodo/figshare/Software Heritage).
【Evidence】 `manuscript_final.md:217` ("Analysis code, artefacts and the test suite are available at https://github.com/yyx-4113/cancer-pain-drg-virtual-perturbation (release v1.0.0), with `requirements.txt` pinning every version used."); `README.md:44` notes "Zenodo 归档已暂缓" (Zenodo archival shelved); `CITATION.cff:17-20` gives the GitHub URL and `version: 1.0.0` but no DOI/archive key.
【Why it matters】 *PLOS ONE*, *Scientific Reports*, *Nature Portfolio*, and many *BMC* journals require that code and data be deposited in a public repository that mints a **persistent DOI** and is archived for the long term; a GitHub release tag alone is frequently judged insufficient because the repo owner can delete or rewrite history. The deliberate omission of Zenodo (per README) is acceptable at *some* venues (e.g., those that accept GitHub + MIT license + CITATION.cff) but is a submission blocker at others. The author should not discover this at technical-review stage.
【Specific fix】 Either (a) mint a Zenodo (or figshare) DOI-archived snapshot of the v1.0.0 tag and add it to the statement: "Code and artefacts: https://github.com/yyx-4113/cancer-pain-drg-virtual-perturbation (release v1.0.0); archived at Zenodo: https://doi.org/10.xxxx/zenodo.XXXXXXX"; or (b) explicitly confirm the chosen journal accepts GitHub release tags as sufficient and state that in the cover letter. Add the DOI to `CITATION.cff` as an `identifiers:`/`doi:` field.

---

## 5. Format hard-fails and internal cross-reference drift

### F1 — HARD: minimum detectable Δ is stated as 0.052 in §3.6 but 0.05 in §5
【Problem】 The same statistic — the minimum detectable cross-aetiology Δ at 80% power with the median bootstrap SE — is reported as "≈0.052" in Results and "≈0.05" in Limitations, a rounding that changes the third significant figure and contradicts the paper's own "every number is gated" claim.
【Evidence】 `manuscript_final.md:143` ("the minimum detectable Δ is ≈0.052 log-units"); `manuscript_final.md:179` ("the median bootstrap standard error the minimum detectable Δ is ≈0.05 log-units").
【Why it matters】 This is the power-floor number that supports the central negative-result conclusion ("the layer resolves large effects, not typical ones"). An editor or the runtime checker should have caught a restatement that disagrees at the third decimal; the inconsistency undermines the "every numerical claim checked against artefacts" assertion (§2.9, §211) and is exactly the class of drift journals' statistical reviewers probe.
【Specific fix】 Use the same value in both places. Recommended: replace both with "≈0.052 log-units" (the more precise figure from §3.6), i.e. line 179 becomes "… the minimum detectable Δ is ≈0.052 log-units — which the largest target-module effect (0.274) clears, but the median target effect (0.031) does not …". Keep the conservative ≈0.10 t-based reference as-is.

### F2 — MODERATE: permutation floor 3/4,901 is given as 6.1 × 10⁻⁴ in §2.7 but 6 × 10⁻⁴ in §5
【Problem】 The exact permutation floor for the four unpaired contrasts is reported as 6.1 × 10⁻⁴ in Methods and as 6 × 10⁻⁴ in Limitations — a rounding that drops a non-trivial digit for a number the paper presents as "exact."
【Evidence】 `manuscript_final.md:87` ("the smallest attainable value is 3/4,901 ≈ 6.1 × 10⁻⁴"); `manuscript_final.md:179` ("a floor of 3/4,901 ≈ 6 × 10⁻⁴").
【Why it matters】 The paper's credibility rests on exact, seed-free permutation p-values; presenting the same exact quantity with two different leading digits (6.1e-4 vs 6e-4) looks like imprecision in a manuscript that elsewhere stresses exactness, and a copyeditor will flag it.
【Specific fix】 Standardise to "3/4,901 ≈ 6.1 × 10⁻⁴" in both §2.7 and §5 (or "≈ 6 × 10⁻⁴" in both if a single significant figure is intended).

### F3 — MODERATE: dual-logFC gap on Gfap and Ccl2 exceeds the paper's own "few hundredths" claim
【Problem】 The manuscript states its two logFC sources "can differ by a few hundredths of a log₂ unit" (§2.3), yet Gfap and Ccl2 differ by ~0.17 and ~0.09 between the §3.3 (muscle-adjusted) and Table 3 (per-sample mean) presentations — far more than "a few hundredths."
【Evidence】 `manuscript_final.md:115` (§3.3 individual-gene values: "Gfap +1.33, Ccl2 +0.81"); `manuscript_final.md:297` (Table 3 per-sample mean: "Gfap | +1.50 | …", "Ccl2 | +0.90 | …"); `manuscript_final.md:59` (§2.3: "they can differ by a few hundredths of a log₂ unit"). For comparison Cd68 2.86→2.84 (0.02) and Aif1 2.18→2.19 (0.01) are within the claimed range, so Gfap/Ccl2 are genuine outliers.
【Why it matters】 A reader (or a statistical reviewer) will see Gfap reported as +1.33 in one place and +1.50 in another and reasonably suspect an error, especially because the paper explicitly claims the two paths agree to within hundredths. The discrepancy is real (different gene sets / adjustment) but is larger than disclosed and is currently unexplained at the point of restatement.
【Specific fix】 Either (a) add a parenthetical at Table 3 / §3.3 noting the source differs and that Gfap/Ccl2 sit at the wide end of the dual-source spread (e.g., "Gfap +1.50 here vs +1.33 in §3.3 because Table 3 uses the unadjusted per-sample mean logFC source, §2.3"); or (b) relax the §2.3 claim from "a few hundredths" to "up to ~0.2 log₂ units for个别 markers such as Gfap," and state that headline markers (Cd68, Aif1, galanin) are quoted from the muscle-adjusted fit.

### F4 — MODERATE: the 4,900 vs 4,901 permutation denominator is internally contradictory
【Problem】 §2.7 first says the enumerated null has 70 × 70 = 4,900 split pairs and that "the observed assignment and its full mirror are both members of the enumerated null," yet it then computes the floor as 3/4,901 — implying a denominator of 4,901 = 4,900 + 1, i.e., the observed split is *outside* the enumerated set. The two statements cannot both be true.
【Evidence】 `manuscript_final.md:87` ("the null over the 70 × 70 = 4,900 (A, B) split pairs is enumerated exhaustively … the observed assignment and its full mirror are both members of the enumerated null … the smallest attainable value is 3/4,901 ≈ 6.1 × 10⁻⁴"). §5 (`manuscript_final.md:179`) repeats "3/4,901."
【Why it matters】 The exactness of the permutation p-value is the paper's methodological backbone; if the observed split is a member of the 4,900 enumerated pairs the correct denominator is 4,900 (floor 3/4,900 = 6.12 × 10⁻⁴), and if it is excluded the denominator is 4,901. The text currently asserts membership *and* uses the excluded denominator, which a statistical reviewer will catch and may read as an arithmetic error rather than a convention.
【Specific fix】 Pick one convention and state it once: e.g., "The 70 × 70 = 4,900 ordered (A,B) split pairs are enumerated; the observed assignment is counted within this set, so the two-sided floor is 3/4,900 ≈ 6.1 × 10⁻⁴." Update both §2.7 and §5 to 4,900, or, if the +1 is intentional (observed counted separately), delete the "both members of the enumerated null" clause.

### F5 — SOFT: reference DOIs and Vancouver formatting should be machine-verified before submission
【Problem】 The 16 references are Vancouver-formatted and carry DOIs, but the set was not independently verified to resolve, and one entry uses "et al." after a single named author (ref 14) which is fine for Vancouver but should be confirmed against the chosen journal's house style.
【Evidence】 `manuscript_final.md:220-235` (refs 1–16, all with `doi:…`); ref 14 "Baruch EN, et al. … Nature. 2025;646:462–473. doi:10.1038/s41586-025-09370-8"; ref 9 "Ahlmann-Eltze C, Huber W, Anders S. … Nat Methods. 2025;22(8):1657–1661."
【Why it matters】 Broken or mismatched DOIs are a routine desk-reject/production flag; several of these (refs 9, 14) are very recent 2025 papers and a wrong page/volume/DOI pair is easy to introduce. A venue editor will run a reference validator.
【Specific fix】 Run each DOI through `https://doi.org/<id>` and confirm a 200/landing page; confirm ref 14's volume/page (Nature 2025;646:462–473) and ref 9's (Nat Methods 2025;22(8):1657–1661) match the publisher records; align author-list truncation (≤ 6 then "et al.") with the target journal's Vancouver variant.

### F6 — SOFT: CITATION.cff and README "Manuscript in preparation" must be updated on acceptance
【Problem】 Both `CITATION.cff` and `README.md` label the work "Manuscript in preparation," which is correct pre-submission but will be stale metadata at the deposited repo once a journal is chosen.
【Evidence】 `CITATION.cff:44` (`journal: "Manuscript in preparation"`); `README.md:6` ("Manuscript in preparation, 2026").
【Why it matters】 `CITATION.cff` is machine-readable and indexed by Zenodo/GitHub; leaving "Manuscript in preparation" after publication breaks citation resolution and the "preferred-citation" block.
【Specific fix】 On acceptance, set `CITATION.cff` `preferred-citation.journal` to the real journal name, year, volume, article number, and add `doi:`; update `README.md:6` to the citation. This is a post-acceptance task, not a blocker.

### F7 — SOFT: single-author CRediT list is acceptable but should be flagged as self-attested
【Problem】 A single author claims all 11 CRediT roles (conceptualization through visualization), which is permissible for a solo computational study but is unusual enough that some editors request a confirming sentence.
【Evidence】 `manuscript_final.md:199` ("Yongxin Yang: conceptualization, methodology, software, validation, formal analysis, investigation, data curation, writing — original draft, writing — review and editing, visualization.").
【Why it matters】 Credibility, not compliance: a one-person study genuinely spanning all roles is fine, but the absence of any "these roles were performed solely by the author" note can prompt a query. The `author_verification_statement.md:5` (B.M. only, no research degree) makes the solo-execution claim more salient.
【Specific fix】 Add one sentence to the Author contributions section: "As the sole author, Y.Y. performed all roles listed above; no contributions were made by others." (Optional, not required.)

---

## § Stands up (things the manuscript does right, at the venue layer)

1. **Exact, seed-free permutation design is reported with full internal consistency on the headline markers.** Cd68 Δ = +3.42 (q = 0.012), Aif1 Δ = +3.08 (q = 0.012), galanin Δ = −2.59 (q = 0.012) appear identically in §3.5, Table 3, and the Abstract; the 0.991 network-vs-additive correlation is stated identically in §3.4, §3.7, §4, Table 2, and the Abstract; 2.47-vs-1.07 / p = 0.0019 anti-selectivity is identical in §3.6, §6, and Figure 3 caption. These core numbers are clean.
2. **Software versions are locked across manuscript, `requirements.txt`, and `README.md`** (V1) — a real reproducibility asset that most submissions lack, and the `statsmodels`-is-mandatory guard (§2.9, README) is a credible safeguard against silent prior degradation.
3. **Negative-result framing is honest and well-structured.** The paper does not dress the perturbation failure up as a positive finding; it uses negative-control modules, degree-matched nulls, and an additive baseline exactly as a venue that welcomes null results (PLOS ONE, Scientific Reports) expects. The IMRaD + structured abstract + Keywords + Running head are all present and on-spec for a generic computational-biology/pain journal.
4. **Ethics, funding, conflicts, and data-source statements are complete and correctly scoped** to a secondary reanalysis (no new animals/humans), with original approvals cited ([10,11,14]).

---

## § Questions for the authors

1. Which journal are you targeting? (This determines whether the missing Zenodo DOI — D1 — and the AI-tool-name gap — A1 — are blockers or non-issues, and whether the title length — R2 — must be cut.)
2. For the minimum-detectable-Δ: is the true value 0.052 or 0.05 (F1)? Please confirm against the deposited artefact and make both §3.6 and §5 agree; also confirm the permutation denominator is 4,900 or 4,901 (F4).
3. Why do Gfap (+1.33 vs +1.50) and Ccl2 (+0.81 vs +0.90) differ by ~0.1–0.17 between §3.3 and Table 3 when §2.3 claims the two logFC sources agree to "a few hundredths" (F3)? Is the gap driven by the 13,893- vs 13,601-gene nested fit, the muscle covariate, or the per-sample-mean source, and should the §2.3 claim be widened?
4. Will you name the specific LLM and confirm it was absent from analysis design/interpretation (A1)? Which tool, and can you attest it did not select methods or results?
5. Can you produce a `pip freeze` / lockfile proving the five pins resolve on PyPI today (V2)?
6. Are all 16 reference DOIs live, and do refs 9 and 14 match the publisher's current volume/page (F5)?

---

## § What I actually checked

**Files read in full:**
- `docs/manuscript/manuscript_final.md` (348 lines, full).
- `README.md` (204 lines, full).
- `CITATION.cff` (45 lines, full).
- `requirements.txt` (17 lines, full).
- `author_verification_statement.md` (58 lines, full) — permitted for this reviewer.
- `docs/review/` directory listing (filename conventions only; did NOT open REVIEW_round*, round1–5_*, R3_editorial_review.md, GITHUB_DEPOSIT_SOP.md, or other reviewers' round6 files).

**Greps run (on `manuscript_final.md`):**
- Version strings: `numpy 2.4.6`, `scipy 1.18.1`, `pandas 2.3.3`, `statsmodels 0.15.0`, `matplotlib 3.11.1`, `Python 3.13` → matched §2.9, `requirements.txt:12-16`, `README.md:114,183`; absent from `CITATION.cff` (expected).
- Key statistics: `0.0019`, `0.0067`, `44.7`, `0.052`, `0.05`, `2.47`, `1.07`, `1.27`, `2.38`, `0.991`, `0.929`, `4,900`, `4,901`, `3/71`, `2/71`, `1,517`, `703`, `55,967`.
  - `0.0019` → §3.6 (line 141) and §6 (line 185): **identical, consistent.**
  - `2.47` vs `1.07` → §3.6 (141), §6 (185), Fig 3 caption (347): **identical, consistent.**
  - `1.27` vs `2.38` → §3.8 (155) and Table 5 (328): **identical (1.27), consistent** — no `1.267` form exists anywhere; the suspected 1.27/1.267 drift is **not present**.
  - `0.991` → §3.4 (125), §3.7 (151), §4 (169), Table 2 (281), Abstract (21): **identical, consistent**; `0.929` (§3.7) is a distinct cross-study quantity, correctly separated.
  - `4,900` → §2.7 (87), Table 3 caption (283): consistent as "enumerated pairs"; `4,901` → §2.7 (87), §5 (179): the +1 denominator is **internally contradictory with the "observed is a member" clause (F4).**
  - `3/71` (§2.7 0.042; §5) and `2/71` (§2.7 0.028; §5): **consistent.**
  - `1,517`/`703` → §3.8 (155) and Table 5 (324): **identical, consistent.**
  - `55,967` → §2.8 (95), §3.8 (155), Table 5 (323): **identical, consistent.**
- Numeric restatements manually cross-checked: Gfap +1.33 (§3.3, 115) vs +1.50 (Table 3, 297) → **drift of 0.17, exceeds disclosed tolerance (F3)**; Ccl2 +0.81 (115) vs +0.90 (297) → **drift of 0.09 (F3)**; Cd68/Aif1/Il1b within 0.02.
- `0.052` (§3.6, 143) vs `0.05` (§5, 179) → **HARD mismatch (F1).**
- `6.1 × 10⁻⁴` (§2.7, 87) vs `6 × 10⁻⁴` (§5, 179) → **rounding drift (F2).**

**Mismatches / inconsistencies found (summary):**
- F1 (hard): 0.052 vs 0.05 minimum-detectable Δ.
- F2: 6.1e-4 vs 6e-4 permutation floor.
- F3: Gfap/Ccl2 dual-source gap > disclosed "few hundredths."
- F4: 4,900 enumerated vs 4,901 denominator internally contradictory.
- A1: AI disclosure omits tool name + analysis/interpretation attestation.
- D1: no DOI-archived code snapshot (GitHub tag only).
- V2/F5/F6/F7: soft/venue-dependent verification items.
- Verified consistent: all version strings (manuscript/requirements/README); 0.0019, 2.47/1.07, 1.27/2.38, 0.991, 4,900-pairs wording, 3/71, 2/71, 1,517/703, 55,967; Fig 1–3 and Tables 1–5 all present and referenced; ethics/funding/conflicts complete.

**Forbidden files (not opened, per instruction):** `docs/review/REVIEW_round*`, `docs/review/round1..5_*`, `docs/review/R3_editorial_review.md`, `GITHUB_DEPOSIT_SOP.md`, `manuscript_v0.*.md`, and all other reviewers' `round6_2026-10-10/` files (A1_domain.md, A2_design.md only noted by filename, not read for content).

— End of A4 venue/reporting-standards review —
