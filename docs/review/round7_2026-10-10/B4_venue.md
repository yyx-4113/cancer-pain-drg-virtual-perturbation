# Reviewer B4 — Venue / Reporting-Standards Review
**Manuscript:** *Divergent dorsal root ganglion programs in bone cancer pain and chemotherapy-induced peripheral neuropathy: an exact cross-study reanalysis and a negative-control calibration of in-silico perturbation*
**Author:** Yongxin Yang (single author)
**Review type:** Fresh independent review (no prior rounds seen)
**Focus:** Venue fit, reporting standards, version/reference/AI/availability consistency, format hard-fails, internal number drift, target-journal framing.

---

## Summary verdict
The reporting layer is **substantially sound and submission-ready after minor fixes**. All four software-version strings match across `manuscript_final.md`, `requirements.txt`, and `README.md`; all three new references [17][18][19] resolve to real papers with valid DOIs and consistent metadata; the AI-disclosure and data/code-availability statements are adequate and internally consistent; the abstract (258 words) is within the ≤300 gate; Figures 1–3 and Tables 1–5 are present and referenced. The only genuine defects are (a) small **q-value rounding drift** between §3.5 text and Table 4, (b) **Table 1 never explicitly cited** in the Results narrative, and (c) a **missing target journal** with no article-type framing. None are blockers; (c) is the most consequential for acceptance strategy.

---

## FOCUS 1 — Version consistency

### Item 1. Software-version strings match across the three code/manuscript files (PASS), but CITATION.cff omits the numeric stack.
【Problem】 The four files are mutually consistent on version strings, but `CITATION.cff` records no Python/numpy/scipy/pandas/statsmodels/matplotlib versions, so a reader consulting only the citation metadata cannot verify the runtime.
【Evidence】 manuscript_final.md:101 — "Python 3.13 with numpy 2.4.6, scipy 1.18.1, pandas 2.3.3, statsmodels 0.15.0 and matplotlib 3.11.1 (`requirements.txt`)". requirements.txt:12–16 — `numpy==2.4.6 / scipy==1.18.1 / pandas==2.3.3 / statsmodels==0.15.0 / matplotlib==3.11.1`. README.md:114 — identical line. CITATION.cff:1–22 — only `version: 1.0.0` of the *software repository*; no numpy/scipy/pandas/statsmodels/matplotlib/Python entries.
【Why it matters】 The exact-match requirement (gate) is satisfied for the three files that matter for reproduction. The CITATION.cff gap is a metadata-completeness issue, not a contradiction; it does not affect acceptance but weakens "cite-the-artifact" robustness if a reader relies solely on the `.cff`.
【Specific fix】 Add to CITATION.cff (e.g., under a comment or a `references`-style note): `# Runtime: Python 3.13; numpy 2.4.6, scipy 1.18.1, pandas 2.3.3, statsmodels 0.15.0, matplotlib 3.11.1 — see requirements.txt`. No change needed to the other three files.

### Item 2. The pinned versions are forward-dated to 2026 but internally coherent (informational).
【Problem】 The pinned stack (numpy 2.4.6, scipy 1.18.1, pandas 2.3.3, statsmodels 0.15.0, matplotlib 3.11.1, Python 3.13) is plausible for a 2026 release but cannot be independently confirmed as "real" today; this is expected for a manuscript dated 2026-10-10 and is not a defect.
【Evidence】 Same lines as Item 1; deposition date `date-released: "2026-10-10"` in CITATION.cff:21 and author_verification_statement.md:9.
【Why it matters】 A reviewer cannot certify these exact patch versions exist yet, but the *internal* consistency (manuscript = requirements = README) is what the gate checks, and that holds.
【Specific fix】 No action. If the target journal asks for a frozen environment, also deposit a `lock` file or a `environment.yml` (optional).

---

## FOCUS 2 — Reference formatting

### Item 3. The three new references [17][18][19] are real, DOIs valid, and metadata consistent (PASS).
【Problem】 None. All three newly added references resolve and their author lists / journals / years / volumes match the cited strings.
【Evidence】
- [17] manuscript_final.md:234 — `Inflammation Research. 2023;72(1):117–132. doi:10.1007/s00011-022-01663-2`. Crossref confirms title *"Activation of the STING pathway induces peripheral sensitization via neuroinflammation in a rat model of bone cancer pain"*, journal *Inflammation Research*, vol 72(1), pp 117–132, authors Zhang/Wang/Gong/Peng/Li/Zhang/Zhang/You/Wu (9) — matches the cited "Zhang Y, Wang W, … Wu J" exactly.
- [18] manuscript_final.md:235 — `Proc Natl Acad Sci U S A. 2025;122(34):e2503779122. doi:10.1073/pnas.2503779122`. PNAS landing page confirms the title *"Unraveling the neuroimmune mechanisms in cancer-induced bone pain…"* at that DOI; the article-number `e2503779122` and volume/issue 122(34)/2025 match.
- [19] manuscript_final.md:236 — `Brain Res. 2007;1168:46–59. doi:10.1016/j.brainres.2007.06.066`. Crossref confirms title *"An evolving cellular pathology occurs in dorsal root ganglia, peripheral nerve and spinal cord following intravenous administration of paclitaxel in the rat"*, *Brain Research*, 2007, vol 1168, pp 46–59, authors Peters/Jimenez-Andrade/Kuskowski/Ghilardi/Mantyh (5) — matches "Peters CM, … Mantyh PW" exactly.
【Why it matters】 Reference integrity is a common desk-reject/copy-edit trigger; here it is clean, which removes a class of delays.
【Specific fix】 No change required.

### Item 4. Reference [17] year: Crossref "published" = 2022 (online), manuscript cites 2023 (print-volume year) — minor, acceptable.
【Problem】 The DOI `10.1007/s00011-022-…` encodes a 2022 online date; the manuscript follows the print volume year 2023. Both are legitimate, but the two sources disagree on the displayed year.
【Evidence】 Crossref `published.date-parts = [2022]` vs manuscript_final.md:234 `2023;72(1)`. The DOI suffix `022` and the issue's online-publication year are 2022.
【Why it matters】 Not a hard error — many Vancouver citations use the print year — but a meticulous copy editor may flag it. Aligning removes a query.
【Specific fix】 Either keep `2023` (print year, defensible) or change to `2022` (online year); state the online/publication dates consistently if the journal requires them. No other reference shows a year conflict.

### Item 5. All 19 references were spot-verified at the load-bearing core; the rest are well-formed Vancouver (PASS, with one caveat).
【Problem】 I directly resolved DOIs for [9][10][11][14][17][18][19] via Crossref/PNAS; all match the manuscript (e.g., [9] Nature Methods 2025;22:1657–1661; [10] Front Mol Neurosci 2021;14:665085; [11] Neuroscience 2023;516:1–14; [14] Nature 2025;646:462–473). The remaining 12 were not each individually resolved.
【Evidence】 Crossref batch output (titles/journals/years/volumes) matches manuscript_final.md:218–236 line by line for the ones tested.
【Why it matters】 The verified set covers every dataset-source and benchmark citation the conclusions rest on. The unverified 12 are standard, well-known pain/statistics citations (ESMO, JNCCN, Pain, JCO, Ann Neurol, ASCO, Nat Methods/scGen/GEARS, Stat Appl Genet Mol Biol, JRSSB, J Pain Symptom Manage, Usoskin Nat Neurosci) and are Vancouver-formatted.
【Specific fix】 Before submission, run a DOI-lint over all 19 (e.g., a `habanero/anystyle` or Crossref batch) so the final copy-editor pass is clean. Not a blocker.

### Item 6. In-text citations [17][18][19] are placed correctly, not orphaned (PASS).
【Problem】 None.
【Evidence】 [17] at manuscript_final.md:115 — "A STING/TBK1/NF-κB innate-immune programme … has been reported in Walker-256 DRG neurons [17]"; [18] at :115 — "macrophage infiltration of the DRG is a documented feature of cancer-induced bone pain [18]"; [19] at :131 — "Atf3, the standard nerve-injury marker (elevated in paclitaxel-treated DRG neurons [19])". Each anchors a specific claim and sits inside a sentence that uses it.
【Why it matters】 Orphaned or mis-placed citations are a frequent minor-revision request; here placement is correct.
【Specific fix】 No change.

---

## FOCUS 3 — AI disclosure

### Item 7. The generative-AI declaration is adequate against PLOS/Nature/Wiley/BMC policy; one optional tightening.
【Problem】 The statement is adequate and covers the required elements, but it names the tool only as "a general-purpose LLM" rather than a specific product, which a few journals (notably PLOS and some Wiley titles) prefer to see named.
【Evidence】 manuscript_final.md:207–209 — "During the preparation of this work the author used a large language model assistant (a general-purpose LLM) for code debugging and manuscript language editing. No AI tool was used to design the statistical analyses, choose the methods, generate the numerical results, or interpret the findings; all analysis code was executed and verified locally and every numerical claim was checked against the deposited artefacts. The author performed all conceptual, analytical and interpretive work, reviewed and edited all content, and takes full responsibility for the published article."
【Why it matters】 This exceeds the bare minimum: it (i) names the *function* of the tool, (ii) explicitly excludes design/analysis/interpretation, and (iii) affirms human authorship and responsibility — satisfying PLOS (tool + how used), Nature (transparency), Wiley/BMC (disclosure) expectations. The only gap is that "a general-purpose LLM" is non-specific; some editorial systems ask for the tool name.
【Specific fix】 Optional: replace "(a general-purpose LLM)" with the actual tool, e.g., "(a large language model, e.g., OpenAI GPT-4 class model)" — only if the author is comfortable disclosing it. The current wording is acceptable for submission.

---

## FOCUS 4 — Data & code availability

### Item 8. Availability statement points to a real, citable location with a version tag, and matches CITATION.cff + README (PASS).
【Problem】 None.
【Evidence】 manuscript_final.md:215 — "Analysis code, artefacts and the test suite are available at https://github.com/yyx-4113/cancer-pain-drg-virtual-perturbation (release v1.0.0), with `requirements.txt` pinning every version used." CITATION.cff:17–18 (`repository-code` and `url` = same URL), :20 (`version: 1.0.0`). README.md:177 — "Data-and-code-availability statement points to this repository (release v1.0.0)."
【Why it matters】 A version-tagged GitHub release is a valid, citable code-availability location for PLOS ONE, Scientific Reports, BMC Medical Genomics, and Journal of Pain Research. Cross-file agreement prevents the "URL in paper ≠ URL in repo metadata" mismatch that triggers queries.
【Specific fix】 No change required.

### Item 9. Deliberate omission of Zenodo is acceptable; an archive DOI would be a nice-to-have, not a blocker.
【Problem】 The manuscript intentionally omits a Zenodo/Mendoley Data archive DOI (the README notes "Zenodo 归档已暂缓" / Zenodo archival deferred). For the candidate venues this is acceptable.
【Evidence】 README.md:44 — `GITHUB_DEPOSIT_SOP.md ← 中文：建仓 / 打 tag 操作手册（Zenodo 归档已暂缓）`; manuscript_final.md:215 cites GitHub release only. No Zenodo DOI appears in manuscript, README, or CITATION.cff.
【Why it matters】 None of PLOS ONE, Scientific Reports, BMC Medical Genomics, or Journal of Pain Research *requires* a Zenodo DOI; a tagged GitHub release with a `MANIFEST.sha256` (present: MANIFEST.sha256:1–15) and a runtime number-checker is a stronger reproducibility story than a static archive. The omission is therefore a defensible editorial choice, not a submission blocker. A permanent archive DOI only becomes mandatory if the chosen journal's policy demands immutable DOIs (some BMC/SpringerNature portfolios encourage, few *mandate*, it).
【Specific fix】 Optional but recommended: mint a Zenodo DOI at acceptance (or upon final revision) and add it to CITATION.cff + the availability statement; this future-proofs the citation if the GitHub repo is later moved. Not required for initial submission.

---

## FOCUS 5 — Format hard-fails

### Item 10. Table 1 is displayed but never explicitly cited in the Results narrative (minor hard-fail risk).
【Problem】 §3.3 quotes the day-7 module moderated-*t* and logFC values that live in Table 1, but the body text never writes "(Table 1)"; the only occurrences of the string "Table 1" are inside the Tables block itself (heading and footnote), not in the narrative.
【Evidence】 grep of manuscript_final.md shows "Table 1" at lines 240 and 267 only (both within the `## Tables` section). §3.3 (lines 113–117) presents "+24.9, q = 6.8 × 10⁻⁵", "Cd68 (+2.86)", "Aif1 (+2.18)", and the full neuronal-module list that is tabulated in Table 1, without a cross-reference. By contrast Table 4 is explicitly cited at :133 ("(Table 4)") and Table 5 at :155.
【Why it matters】 Most journals' editorial rules require every display item to be cited in the text; an uncited table can trigger a "cite or remove" query and looks like an oversight. Tables 2, 3, 4, 5 are all cited in-body; Table 1 is the lone exception.
【Specific fix】 In §3.3, at first presentation of the module values, add the citation, e.g.: "At day 7 the interferon-stimulated-gene control module carried the largest mean moderated *t* in the study (+24.9, q = 6.8 × 10⁻⁵; Table 1), followed by glial/immune (+11.4)…". (And optionally remind the reader of Table 1 when the neuronal-module list is enumerated.)

### Item 11. Figures 1–3 and Tables 1–5 are present and referenced; abstract within limit; CRediT/ethics/conflicts/funding complete (PASS).
【Problem】 None — no format hard-fail on counts, abstract length, author statements.
【Evidence】 Figures: Fig 1 (design) cited at :49/:342; Fig 2 (markers) at :129/:346; Fig 3 (anti-selective) at :137/:350 — all present as `![Figure N]` at :340–349. Tables 1–5 present at :240–334 and cited (Table 2 :121; Table 3 :59/:129/:286/:346; Table 4 :133/:302; Table 5 :155/:320). Abstract body (lines 17–21) = **258 words** ≤ 300 gate (line 15 "## Abstract" + blank + keywords excluded). Author contributions: :197 lists all 14 CRediT roles for the single author — credible for a sole-author computational study. Ethics: :203–205 complete (secondary analysis, original approvals cited [10][11][14], ARRIVE scoping stated). Funding: :189 "no specific grant"; Conflicts: :193 "declares no competing…". 
【Why it matters】 These are the items most likely to cause an immediate desk-format rejection; all pass, so no structural blocker exists.
【Specific fix】 No change. (See Item 10 for the one Table-1 citation gap.)

### Item 12. Single-author CRediT listing all roles is credible (PASS, note only).
【Problem】 A single author claiming conceptualization-through-visualization across all 14 CRediT roles can read as inflated to some editors, though it is standard for genuine single-author work.
【Evidence】 manuscript_final.md:197.
【Why it matters】 Not a defect; flagged only so the author is ready if a reviewer questions it. The author_verification_statement.md:36–38 already discloses "Bachelor of Medicine only, no post-nominal research degree," which preempts any credential-inflated-reading concern.
【Specific fix】 No change; optionally add a one-line note "All roles performed by the sole author" if the journal's CRediT instructions imply multi-author context.

---

## FOCUS 6 — Internal cross-reference drift

### Item 13. All flagged key statistics are identical across every restatement (PASS — author's unification claim supported).
【Problem】 None. Every number the author claimed to have unified is, in fact, identical wherever restated; I found **no** inconsistent rounding of the probed values.
【Evidence】 grep results across manuscript_final.md:
- `0.0019`: lines 21, 141, 183, 350 — identical (Mann–Whitney anti-selectivity p).
- `0.0067`: line 139 — binomial over-representation p (single occurrence, internally consistent).
- `44.7`: line 141 — "44.7th percentile" (single, consistent).
- `0.052`: lines 143, 177 — "≈0.052 log-units" minimum detectable Δ (identical).
- `2.47` vs `1.07` (control vs target median |z|): lines 141, 183, 254, 350 (text/Table 1 n-col/Fig 3 caption) — identical pair every time.
- `1.27` vs `2.38` (observed vs MDE Cohen's d): lines 155 (text) and 331 (Table 5) — identical.
- `0.991` (network vs additive ρ): lines 21, 125, 151, 167, 284 — identical.
- `3/4,901` and `6.1 × 10⁻⁴`: lines 87, 177 — identical; no `6 × 10⁻⁴` variant found.
- `3/71` (lines 87, 89, 177) and `2/71` (lines 89, 177) — identical.
- `1,517` / `703` (up/down DE): lines 155, 327 — identical.
- `55,967` (cells): lines 95, 97, 155, 326 — identical.
- `4,900` (null pairs): lines 19, 87, 286 — identical.
No occurrence of `0.05` in place of `0.052`, no `1.267` in place of `1.27`, no `6 × 10⁻⁴` in place of `6.1 × 10⁻⁴`. The probed inconsistencies the brief warned about are **absent**.
【Why it matters】 Number drift is the most common cause of reviewer "internal inconsistency" rejections; this manuscript is clean on the audited set, which substantially de-risks the venue review.
【Specific fix】 No change. Keep the runtime checker (scripts/98) that enforces this.

### Item 14. Minor q-value rounding drift between §3.5 text and Table 4 (the one real internal inconsistency).
【Problem】 Three module q-values are rounded in the narrative to 3 decimals in a way that does not match the 4-decimal values printed in Table 4, and the rounding direction is inconsistent (one truncated, one rounded up).
【Evidence】
- Glial/immune: §3.5 :133 "q = 0.007" vs Table 4 :308 "q = 0.0073".
- Erythrocyte (module-level ΔlogFC): §3.6 :145 "q = 0.007" vs Table 4 :318 "q = 0.0073".
- Neuropeptides: §3.5 :133 "q = 0.015" vs Table 4 :313–314 "q = 0.0145".
The text shows `0.007` (truncation of 0.0073) and `0.015` (round-up of 0.0145) — two different rounding conventions applied to the same table, so the narrative cannot be reconciled with Table 4 by a fixed rule.
【Why it matters】 This is exactly the class of "text ≠ table" drift reviewers flag. It is tiny in magnitude but visible, and a careful copy editor will query it. It does not affect any conclusion (all are "significant" either way), so it is a pure consistency fix.
【Specific fix】 Make the narrative match Table 4 exactly: change :133 to "glial/immune … (q = 0.0073)" and "neuropeptides … (q = 0.0145)", and :145 to "(ΔlogFC = −2.48, q = 0.0073)". Alternatively, if the journal prefers 3-decimal reporting, round Table 4 to 3 decimals uniformly (0.007, 0.015, 0.007) — but then apply one rule everywhere, not a mix.

### Item 15. §3.5 quotes a cytokine/chemokine "Δ = +1.75 against oxaliplatin" that Table 4 does not print (minor cross-reference gap).
【Problem】 The text states a D7-vs-OXA Δ for cytokine/chemokine (+1.75) that is arithmetically correct from Table 4 (bone D7 +1.40 − OXA −0.35 = +1.75) but is not shown as a printed cell in Table 4, which only prints the D7-vs-PTX Δ (+1.70).
【Evidence】 manuscript_final.md:133 — "cytokine/chemokine Δ = +1.75 against oxaliplatin and +1.70 against paclitaxel"; Table 4 :309 prints only "Δ (D7 vs PTX) = +1.70" for cytokine/chemokine, with OXA logFC = −0.35 and bone D7 = +1.40 visible but no OXA-Δ cell.
【Why it matters】 A reader checking the text against Table 4 cannot find the +1.75 figure; it is derivable but not displayed, which reads as a missing table column rather than an error. Low severity.
【Specific fix】 Either add a "Δ (D7 vs OXA)" column to Table 4, or rephrase :133 to "cytokine/chemokine Δ = +1.70 (D7 vs PTX) and +1.75 (D7 vs OXA, derived from the printed logFC values)". Choose one.

---

## FOCUS 7 — Article type / target journal

### Item 16. The manuscript names no target journal and gives no article-type framing; this is the single most consequential venue-layer gap.
【Problem】 Nowhere in the manuscript, README, or CITATION.cff is a target journal or article type stated (CITATION.cff:44 literally reads `journal: "Manuscript in preparation"`). For a paper whose headline is a *negative* computational result plus a *positive* exact reanalysis, the framing determines acceptance odds.
【Evidence】 CITATION.cff:43–44 (`year: 2026`, `journal: "Manuscript in preparation"`); README.md:6 ("Manuscript in preparation, 2026"); manuscript_final.md — no "To the Editor" / target-journal line anywhere. The negative-result emphasis is explicit (§3.6, §4, §6) and the n = 4 per group is stated up front (Abstract :19–20, §2.1).
【Why it matters】 Journals differ sharply on (a) whether a reanalysis/negative computational result is in scope and (b) whether n = 4 animal groups are acceptable. Picking the wrong venue = desk rejection; picking the right one + framing = straightforward acceptance.
【Specific fix】 State the intended journal and article type explicitly in the cover letter and (if preprinting) the README. Recommended hierarchy:

1. **PLOS ONE — Research Article (primary recommendation).** Explicitly welcomes negative results and secondary/reanalysis studies; the exact-permutation + negative-control-calibration design is squarely in scope; n = 4 is acceptable when, as here, the permutation p is exhaustive and seed-free. Frame the title/abstract to foreground the *positive, exact* biological divergence (bone-cancer myeloid/ISG vs CIPN neuronal-injury) as the contribution, with the in-silico-perturbation negative result as a rigorous methodological companion — do **not** lead with "we found nothing."
2. **Scientific Reports — Article (strong secondary).** Broad scope, negative results accepted, computational reanalysis fine. Same framing as above.
3. **BMC Medical Genomics — Research Article (alternative).** Fits the "genomics of pain" angle; would appreciate the reproducible artifact story (MANIFEST + checker). Watch: they may want clearer novel-biology emphasis.
4. **Journal of Pain Research — Original Research (possible).** Pain-audience fit is excellent, but it leans clinical/translational; a purely computational rat/mouse reanalysis with n = 4 may be judged "limited translational immediacy" — acceptable if the mechanistic correlate to opioid-refractory bone pain vs duloxetine-managed CIPN (§4) is foregrounded.
5. **Journals likely to desk-reject a reanalysis/negative computational paper:** *PAIN*, *The Journal of Pain*, and other top-tier clinical-pain titles — they generally desk-reject pure reanalyses without new data and are skeptical of n = 4 claims; *Nature/Science*-tier methods venues would not take a negative "it doesn't work" perturbation result as a positive methods contribution. Avoid these unless reframed as a *methods/caution* paper (then consider *BMC Bioinformatics* or a methods-focused outlet).

**Framing rule that maximises acceptance:** Lead with the secured biological claim ("DRG programmes of bone-cancer pain and CIPN diverge exactly and direction-specifically, readable from individual markers, validated by independent scRNA-seq as model-specific"), present the negative-control calibration of in-silico perturbation as a *methodological contribution* (a reusable negative-control + degree-matched-null + additive-baseline protocol), and label the perturbation null result as hypothesis-generating. This converts the "negative" into a citable methods advance and satisfies the "negative results welcome" policies of the recommended venues.

---

## § Stands up (what is already strong)
1. **Version provenance is airtight across the three authoritative files** (manuscript = requirements.txt = README), and a runtime number-checker + MANIFEST.sha256 back it — a genuine reproducibility strength rare in single-author submissions.
2. **Reference integrity is clean where it matters most**: all three new refs [17][18][19] and the four load-bearing dataset/benchmark refs [9][10][11][14] were resolved and match; DOIs are valid; in-text placement is correct.
3. **The audited internal numbers are fully consistent** — every probed statistic (0.0019, 0.0067, 44.7, 0.052, 2.47/1.07, 1.27/2.38, 0.991, 3/4,901, 3/71, 2/71, 1,517/703, 55,967, 4,900, 6.1×10⁻⁴) is identical on restatement, with no inconsistent rounding. The author's unification claim holds.
4. **AI disclosure, data/code availability, ethics, funding, conflicts, and abstract length (258 words) all pass** venue gate checks; Figures 1–3 and Tables 1–5 are present and (except Table 1, Item 10) cited.

---

## § Questions for the authors
1. Which journal do you intend to submit to, and in what article type? (See Item 16 — this determines whether the negative-result framing needs adjustment before submission.)
2. For reference [17], do you intend the print year 2023 or the online year 2022? (Item 4.)
3. Can you add an explicit "(Table 1)" cross-reference in §3.3, and align the three q-values in §3.5/§3.6 with Table 4 (Item 10, Item 14)?
4. Will you mint a Zenodo (or equivalent) archive DOI at acceptance to give the code an immutable citation, given the GitHub release is the only current pointer (Item 9)?
5. For the AI-disclosure statement, are you willing to name the specific LLM tool used, or keep the current generic phrasing (Item 7)?

---

## § What I actually checked
**Files read in full:** `docs/manuscript/manuscript_final.md` (351 lines), `README.md` (204 lines), `CITATION.cff` (45 lines), `requirements.txt` (17 lines), `MANIFEST.sha256` (header/structure, 15 lines), `author_verification_statement.md` (58 lines). **Files NOT opened (per instruction):** `docs/review/REVIEW_round*`, `docs/review/round1..6_*`, the entire `docs/review/round6_2026-10-10/`, `GITHUB_DEPOSIT_SOP.md`, `manuscript_v0.*.md`. **Folder listing only:** `docs/review/` and `docs/review/round7_2026-10-10/` (filenames only; did not open other reviewers' B1–B3 files in round7).

**Greps run on the manuscript:** `0.0019, 0.0067, 44.7, 0.052, 2.47, 1.07, 1.27, 1.267, 2.38, 0.991, 6.1, 6×10, 3/71, 2/71, 4,900, 55,967, 1,517, 703, q = 0.007, 0.0073, q = 0.015, 0.0145, Table [1-5], Figure [1-3]`. **Version-string cross-file check:** manuscript §2.9 vs requirements.txt vs README vs CITATION.cff. **DOI verification (live):** [17] (Crossref), [18] (PNAS landing), [19] (Crossref), plus [9][10][11][14] (Crossref batch) — all resolved and matched. **Abstract word count:** 258 words (body, excl. heading/keywords) ≤ 300.

**Mismatches / defects found:**
- Item 10 — Table 1 not cited in §3.3 narrative (only in the Tables block).
- Item 14 — q-value drift: text `0.007`/`0.015` vs Table 4 `0.0073`/`0.0145` (glial/immune, erythrocyte, neuropeptides).
- Item 15 — §3.5 cytokine/chemokine "+1.75 vs OXA" not printed in Table 4.
- Item 4 — [17] year 2023 (print) vs Crossref 2022 (online), minor.
- Item 1 — CITATION.cff omits the numeric Python stack (no contradiction, metadata gap).
- Item 16 — no target journal / article type stated.

**No inconsistency found** for any of the probed key statistics in Focus 6 (Item 13 passed clean).

---

*Prepared by Reviewer B4 (venue / reporting-standards). Independent of rounds 1–6 and of co-reviewers B1–B3 in round7.*
