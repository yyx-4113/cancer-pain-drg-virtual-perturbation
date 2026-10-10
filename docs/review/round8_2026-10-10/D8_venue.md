# Reviewer D8 — Independent venue / reporting-standards review

**Manuscript:** *Divergent dorsal root ganglion programs in bone cancer pain and chemotherapy-induced peripheral neuropathy: an exact cross-study reanalysis and a negative-control calibration of in-silico perturbation* (single author, Y. Yang)
**Mode:** Fresh, independent first-submission review. No prior review rounds, author statements, or other reviewers' files were read.
**Scope of this layer:** version-string consistency, reference formatting/DOI resolution, AI-disclosure adequacy, data/code availability, format hard-fails, internal number drift, and article-type / target-journal framing.

---

## Findings

### F1 — CITATION.cff keyword "type-I interferon" overstates the (deliberately softened) manuscript claim
【Problem】 The repository metadata asserts a stronger claim than the manuscript body supports.
【Evidence】 `CITATION.cff:26` lists the keyword `type-I interferon`. The manuscript keyword line (`manuscript_final.md:25`) instead uses the softer `interferon-stimulated genes`, and the body repeatedly hedges: `manuscript_final.md:115` ("we did not distinguish a type-I from a broader … innate-immune programme"; "'type-I' is inferred, not demonstrated"), `:165` ("an interferon-stimulated-gene programme (type-I not distinguished from type-II)"), and `:179` (states bulk rat DRG "could not determine whether it is type-I, type-II, or mixed").
【Why it matters】 This is exactly the kind of metadata/title-body drift a screening editor or indexer will flag. The manuscript's own careful softening is undermined by a keyword that commits to "type-I."
【Specific fix】 Change the `CITATION.cff` keyword to `interferon-stimulated genes` (matching the manuscript keyword list), or drop it. Keep "type-I" only in the hypothesis-discussion prose where it is explicitly qualified.

### F2 — AI-disclosure names the function but not the specific tool (required by the recommended venue)
【Problem】 The declaration is well-structured and properly scoped, but does not name the model/service used.
【Evidence】 `manuscript_final.md:209-211`: "the author used a large language model assistant (**a general-purpose LLM**) for code debugging and manuscript language editing. No AI tool was used to design the statistical analyses, choose the methods, generate the numerical results, or interpret the findings."
【Why it matters】 The framing excludes design/analysis/interpretation correctly (good), but PLOS ONE's AI-disclosure policy asks authors to name the specific tool used. As I recommend PLOS ONE as the primary target, this becomes a submission-required item.
【Specific fix】 Replace "a general-purpose LLM" with the concrete product/model name (e.g., "ChatGPT (OpenAI)" or the specific assistant). Keep the exclusion sentence verbatim — it already satisfies the "no design/analysis/interpretation" requirement.

### F3 — Abstract length is compliant but tight against the PLOS ONE cap
【Problem】 The structured abstract sits close to the word limit of the recommended venue.
【Evidence】 Word count of the abstract body (excluding the keyword line) is **≈282 words** (regex token count 299 incl. hyphenated forms and the "×"/"–" symbols; whitespace count 282 incl. the four bolded section labels). PLOS ONE caps abstracts at **300 words**; Scientific Reports at **200**.
【Why it matters】 282 is under 300, so this is *not* a hard fail, but it leaves little margin and a different word-counting convention (e.g., counting "Background./Methods./Results./Conclusions." labels or the structured sub-headings) could tip it over at editorial copy-editing.
【Specific fix】 No change strictly required for PLOS ONE. If targeting Scientific Reports (200-word cap), the abstract must be cut by ~80 words. For PLOS ONE, a light trim of one or two subordinate clauses is advisable as a safety margin; verify the final count with the journal's own counter at submission.

### F4 — Minor rounding drift in the Aif1 cross-aetiology Δ between text and Table 3
【Problem】 The text Δ for Aif1 does not exactly equal the arithmetic of the cited table cells.
【Evidence】 `manuscript_final.md:129` states Aif1 "Δ = +3.08 (q = 0.012)" for day-7 vs paclitaxel, but Table 3 cells are +2.19 (D7) − (−0.88) (PTX) = **3.07**; and the day-14 contrast is given as "+2.84" while the table cells +1.95 − (−0.88) = **2.83**. Cd68 is internally exact (2.84−(−0.58)=3.42 ✓; 2.49−(−0.58)=3.07 ✓).
【Why it matters】 A ~0.01 offset is within rounding tolerance and does not change any conclusion, but it is the one place where a text figure and a table cell disagree by more than rounding symmetry, and a number-checker (scripts/98) keyed to deposited CSVs would be the place to confirm which literal is authoritative.
【Specific fix】 Reconcile the Aif1 Δ literals to the Table 3 cells (use 3.07 and 2.83, or confirm the table cells themselves are the rounded values and adjust the text to match). Verify against `etiology_marker_cross.csv`.

### F5 — Reference [17] carries a print-year (2023) that differs from its online year (2022)
【Problem】 The DOI suffix signals an online-first date one year earlier than the cited print year.
【Evidence】 `manuscript_final.md:236` cites Zhang et al., *Inflammation Research* **2023**;72(1):117–132, doi:`10.1007/s00011-022-01663-2` (the "022" prefix is the 2022 online-publication identifier).
【Why it matters】 Not a defect — Vancouver permits citing the print year, and 2023 is what the volume/issue carries. Flagged only for completeness so the author is aware the same article may appear as "2022 (Epub)" in some databases.
【Specific fix】 No change required. Optionally add "Epub 2022" if the target journal's style favours it; otherwise leave as is.

---

## What stands up (≥3)

1. **Version-string consistency is exact.** `numpy 2.4.6 / scipy 1.18.1 / pandas 2.3.3 / statsmodels 0.15.0 / matplotlib 3.11.1` appear identically in `manuscript_final.md:101`, `README.md:114`, and `requirements.txt:12-16`; Python 3.13 is stated in both the manuscript and README. `CITATION.cff` correctly carries no package versions (it is software metadata with its own `version: 1.0.0`, matching the `release v1.0.0` referenced in text and README). No drift.
2. **All Tables 1–5 and Figures 1–3 are cited in the text.** Table 1 (§3.3, §3.4), Table 2 (§3.4), Table 3 & Figure 2 (§3.5), Table 4 (§3.5), Table 5 (§3.8), Figure 1 (§2.1), Figure 3 (§3.6). No orphaned or uncited display items.
3. **References are well-formed Vancouver with resolvable DOIs.** All 19 entries ([1]–[19]) follow a consistent Vancouver template and are cited in the body. I spot-verified the six requested DOIs plus [19] by resolving them at `doi.org`; all returned HTTP 302 to a valid publisher landing page: [9] Nature Methods, [10] Front Mol Neurosci, [11] Neuroscience, [14] Nature, [17] Inflammation Research, [18] PNAS, [19] Brain Res.
4. **Data/code availability points to a plausible, internally consistent repository.** `manuscript_final.md:217` gives `https://github.com/yyx-4113/cancer-pain-drg-virtual-perturbation` (release v1.0.0); the same URL/version appears in `CITATION.cff:17-21` and README. Format is a syntactically valid GitHub release URL (no network verification of the live tag was performed, per instructions).
5. **AI disclosure is substantively adequate** — it names the *function* (code debugging, language editing) and explicitly *excludes* design, method choice, result generation, and interpretation. This satisfies the core editorial requirement; only the specific-tool name is missing (F2).
6. **Mandatory sections are complete.** Funding (`:189`), Conflicts of interest (`:193`), Author contributions (`:197`), Ethics statement (`:205`) are all present. The single-author CRediT list (conceptualization through visualization) is credible for a solo computational reanalysis; "data curation" is light (public data) but defensible.
7. **Internal numbers are largely self-consistent.** §3.3 module *t* values, §3.4 reversal/z/oracle/ρ, §3.5 marker and module Δ/q, and §3.6 anti-selectivity statistics (median |z| 2.47 vs 1.07, p = 0.0019) all match their tables/figure captions. The two logFC sources (muscle-adjusted vs raw per-sample mean) are explicitly disclosed in §2.3 and the Cd68 2.86↔2.84 / Aif1 2.18↔2.19 difference is explained rather than contradicted.

---

## Questions (for the author, not blockers)

- Q1. Will the `release v1.0.0` GitHub tag be live and immutable at the time of submission? Some venues (and the number-checker gate) assume the deposited artefacts are retrievable; confirm the tag exists and is not a branch.
- Q2. For PLOS ONE, the AI-disclosure must also appear in the cover letter — is the author prepared to supply that, and to name the tool (F2)?
- Q3. The README describes the study as a "negative result, by design" for the perturbation layer, while the title/abstract lead with the biological divergence. Is the intended *primary* contribution the methodological caution (favours PLOS ONE / Scientific Reports) or the biological reanalysis (favours a pain/neuro journal)? The current dual framing is acceptable, but the authors should be deliberate about which audience they lead with.

---

## What I actually checked

- Read the full `manuscript_final.md`, `README.md`, `CITATION.cff`, `requirements.txt` (the four files permitted for this layer).
- Ran a script to extract/count the abstract (≈282 words), enumerate citations used in the body (all of [1]–[19] present, including the grouped `[1,2]`/`[3,4]`/`[5,6]`/`[7,8]` forms), and confirm 19 reference entries.
- Ran a script comparing package version strings and the Python-3.13 mention across the four files.
- Grepped for the `type-I interferon` phrase across manuscript, CITATION.cff, and README to localise the drift (F1).
- Resolved 7 spotlight DOIs (HEAD at doi.org) — all valid.
- Manually cross-checked text-vs-table numbers in §3.3–§3.6 and Tables 1–5 / Figure 3 captions; found one sub-rounding Aif1 discrepancy (F4).
- Confirmed every Table/Figure is referenced in-text.
- Did **not** read any `docs/review/REVIEW_round*.md`, `RESPONSE_*.md`, `round1_*`–`round7_*` directories, `author_verification_statement.md`, `GITHUB_DEPOSIT_SOP.md`, `manuscript_v0*.md`, or other `round8_*` reviewers' files, per the independence rules.

---

## Article-type / target-journal framing assessment

**Recommended primary target: PLOS ONE — Research Article.**
Rationale: PLOS ONE explicitly welcomes negative/null results and methodological/secondary-analyses papers, applies uniform technical-soundness (not "impact") criteria, and requires exactly the reporting rigour this manuscript already supplies (data/code availability, ethics, conflicts, funding, AI disclosure, reproducible artefacts with a number-checker gate). The manuscript's honest "two conclusions of unequal weight" structure and its explicit labelling of the perturbation layer as hypothesis-generating are a good fit for that venue's negative-result policy.

**Acceptable alternatives, in order:** *Scientific Reports* (also negative-result-friendly; note the 200-word abstract cap — F3 would then become a required cut); *BMC Medical Genomics* (genomics scope fits the transcriptomics reanalysis, but the pain/neuro framing is secondary); *Journal of Pain Research* (would value the biological divergence but may view the heavy perturbation-calibration apparatus as out of scope).

**Is the current framing appropriate, or does it need reframing?** The current framing — a positive, exactly-tested biological divergence as the validated context, with the negative in-silico-perturbation calibration presented as the methodological contribution — is **appropriate for acceptance at PLOS ONE / Scientific Reports and needs no structural reframing.** The only adjustments are cosmetic/housekeeping: harmonise the `type-I interferon` keyword (F1), name the LLM for PLOS compliance (F2), and confirm the abstract length against the specific venue cap (F3). I would *not* recommend leading the paper as a pure pain-biology positive result, because the manuscript's genuine novelty is the calibrated negative finding; the dual framing already preserves this and should be kept.

---

**VERDICT: ACCEPTABLE — no remaining venue/reporting blocker.** (Minor, non-blocking items to address at submission: F1 keyword harmonisation, F2 name the LLM, F3 confirm abstract length for the chosen venue, F4 Aif1 Δ rounding. None is a blocker.) Recommended target journal + article type: **PLOS ONE, Research Article**.
