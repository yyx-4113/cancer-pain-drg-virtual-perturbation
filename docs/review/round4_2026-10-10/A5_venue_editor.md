# R4-A5 Review — Venue Editor (Submission Readiness)

*Independent Round-4 panel member. Read only `_PANEL_BRIEF.md` and
`docs/manuscript/manuscript_final.md` (plus public artefacts I cite, and public
web/API resources for DOI and repository verification). No other file under
`docs/review/`, no `manuscript_v0.*`, no `.workbuddy/` file, no
`SUBMISSION_MANIFEST.md` was opened. Manuscript treated as a first submission.*

---

## 1. Verdict

**Minor Revision.**

The manuscript is **not publishable as-is**, but every remaining item on my
(submission-readiness) list is a bounded, non-scientific fix that a single
editing pass resolves. My seat is editorial; the scientific soundness verdict
belongs to the statistics/methods and domain reviewers. On *readiness grounds
alone* I would not let it through triage unchanged, because one of the fixes —
the Data & code availability link — bears directly on the reproducibility claim
that is the spine of a negative-result reanalysis.

Nothing I found requires new data, new analysis, or re-interpretation.

---

## 2. Summary of the manuscript (shows I read it)

A single-author computational reanalysis (title, line 1; author + single
affiliation, lines 5–7; ORCID line 9) of two public bulk rat DRG RNA-seq
studies on the same platform — GSE149648 (Walker-256 tibia bone cancer pain;
sham/d1/d7/d14) and GSE160543 (paclitaxel/oxaliplatin neuropathy) — plus an
independent mouse single-cell validation (GSE292522, B16F10 nerve-injury).
It makes two claims of unequal weight: (i) a *positive* biological claim that
bone cancer pain and CIPN drive qualitatively different DRG programmes, tested
by exhaustive label permutation with the exact p floor worked out per endpoint
(§2.7, lines 87–89); and (ii) a *negative/methodological* claim that in-silico
perturbation of these bulk transcriptomes provides no valid analgesic-target
ranking, failing its own negative controls and additive baseline (§3.4, §3.6,
lines 119–145). The single-cell model is used to show the Walker-256
interferon programme is *model-specific*, not a cancer-pain universal (§3.8,
lines 153–161). Results are stated to be reproducible via 14 test groups and a
number-checking gate (§2.9, line 101). Five tables and three figures support it.

---

## 3. Major issues (numbered, each with citation, evidence, required fix)

### M1. The Data & code availability URL does not resolve (HTTP 404) — must fix

**Citation.** "Data and code availability", line 213:
`... deposited in the project repository (https://github.com/yyx-4113/cancer-pain-drg-virtual-perturbation), to be archived at Zenodo upon acceptance`.
**Evidence.** I checked both the page and the API:

```
$ curl -s -o /dev/null -w "%{http_code}" -L \
   https://github.com/yyx-4113/cancer-pain-drg-virtual-perturbation
404
$ curl -s ~/api.github.com/repos/yyx-4113/cancer-pain-drg-virtual-perturbation
{"message": "Not Found", ... "status": "404"}
$ curl -s -o /dev/null -w "%{http_code}" ~/api.github.com/users/yyx-4113
200
```

The **user account exists (200) but the repository does not (404)** — so this is
a placeholder, not a private-repo artefact. The statement is worded as if the
code is live ("are deposited in the project repository"), which is not
transparent about its placeholder status.
**Why it matters for this paper specifically.** §2.9 (line 101) and §2.2–2.6
put the reproducibility machinery front-and-centre; a reanalysis whose sole
code/data pointer is dead at the moment of submission undercuts the paper's own
central methodological claim and is exactly what a soundness-oriented venue's
data-availability check is built to catch at triage.
**Required fix (one of).** (a) Make the repository public before submission and
confirm the link returns 200; or (b) if it must stay private until acceptance,
rewrite the statement so it is honest about status, e.g. "analysis code and
artefacts are available from the corresponding author on reasonable request and
will be released publicly (GitHub, then archived at Zenodo) upon acceptance" —
and provide the code to the editor at submission. Either way, do not assert
"are deposited" against a 404. Recommend also minting the Zenodo DOI at
acceptance (current "to be archived upon acceptance" wording is otherwise
standard and acceptable).

### M2. Table 4 and Table 5 are never cited in the text — must fix

**Citation.** Table 4 heading (line 299) and Table 5 heading (line 316) exist,
but a full-text scan of §1–§6 finds **no in-text call-out of "Table 4" or
"Table 5"** anywhere. The only in-text table citations are Table 1 (§3.3,
line 115), Table 2 and Figure 2 (§3.5, line 129) and Table 3 (§3.5, line 133).
**Evidence.** Regex scan over the body (everything before `## References`):
in-text mentions are `Table 1`, `Table 2`, `Table 3` only; `Table 4` and
`Table 5` appear solely as their own headings.
**Why it matters.** Vancouver/ICMJE venues require every table to be cited in
the text, in numerical order. Two orphaned tables are a structural defect that
at minimum triggers a formatting return and at worst signals to a reviewer that
the calibration table (§3.4) and the validation table (§3.8) were bolted on.
**Required fix.** Add call-outs at the natural seams: cite **Table 4** in §3.4
(the day-7 reversal/calibration paragraph, line 121–125) and **Table 5** in
§3.8 or the Discussion summary of the single-cell result (line 155–161).

### M3. Abstract length: recomputed 325 words (body) / 329 (with headings), not 295 — must fix

**Citation.** Abstract, lines 17–23; panel brief states "295-word abstract".
**Evidence.** Exact count of the four abstract paragraphs (markdown bold
stripped):

```
total tokens including section labels: 329
excluding the 4 label words (Background./Methods./Results./Conclusions.): 325
```

Word boundaries: `Background.`=token 0, `Methods.`=52, `Results.`=156,
`Conclusions.`=285. So the body is **325 words**; the stated **295 is not
reproducible** under a plain word-split (and is ~30 words low).
**Why it matters.** A 325-word structured abstract exceeds the common
250–300-word cap at most SCIE venues; over-length abstracts are frequently
desk-returned *before* review, and the mismatch with the declared count is a
small reporting-integrity problem at submission. (Some venues — e.g. a 350-word
cap — would tolerate it, so the fix is venue-dependent.)
**Required fix.** Confirm the target venue's abstract cap; if ≤300, trim ~25–35
words (the Methods and Results sentences in the abstract are the natural
targets) and correct any stated word count to the recomputed value.

---

## 4. Minor issues (numbered, same format)

### m1. Reference [14] is cited out of numerical order

**Citation.** References list, line 230 (ref 14, Luger 2005) vs the in-text
first appearances of refs 15 (line 95, §2.8) and 16 (line 97, §2.8).
**Evidence.** Scan of first-citation order over the body:
`1,2,3,4,5,6,7,8,9,10,11,12,13,15,16,…,14` — i.e. [14] first appears in the
Discussion (line 167), **after** [15] and [16] first appear in Methods §2.8.
**Fix.** Renumber by first appearance (Baruch→[14], Usoskin→[15], Luger→[16]),
or move the [14] citation earlier. All 16 references are otherwise cited and
present.

### m2. No standalone Ethics / animal-reporting (ARRIVE-style) declaration item

**Citation.** Declarations block, lines 189–213: Funding, COI, Author
contributions, Acknowledgements, Generative-AI, Data & code availability are
all present; there is **no** separate "Ethics approval" / animal-ethics item.
The relevant content exists but is buried in §2.1 (line 49): "no new animal
experiments were performed… the original collections were approved by the
respective institutional animal-care committees."
**Fix.** Add an explicit "Ethics statement" declaration (secondary analysis of
publicly deposited data; original studies' approvals). If the target venue
expects ARRIVE 2.0 items for animal work, state that ARRIVE reporting is
inherited from the depositing studies and/or supply the item list — the brief
names ARRIVE-style reporting as the standard in play.

### m3. Figure 1 caption contains an inaccurate boilerplate clause

**Citation.** Figure 1 legend, line 338, final sentence: "All plotted values
are regenerated from `results/tables/`."
**Evidence.** Figure 1 (I inspected `figures/fig1.png`) is a study-design
schematic — boxes and arrows, no plotted data values, and it does not read from
`results/tables/`. The clause is copied from the data-figure captions.
**Fix.** Delete the clause from the Figure 1 legend (Figures 2 and 3, which *do*
read from CSVs, may keep their provenance sentences).

### m4. Title is long (202 characters)

**Citation.** Line 1. Many SCIE venues cap titles (commonly ~120–150 chars, or
~15–20 words). Consider trimming the post-colon qualifier.
**Fix.** Optional but advisable: shorten to the discovery claim plus the
methodology tag.

### m5. Figures are supplied inline as `.png`; provide discrete files + a legends list

**Citation.** Figures section, lines 336–346 embed `figures/fig1-3.png` with
bold captions adjacent.
**Evidence.** `figures/` contains both `.png` (4278×2118, 2704×1995, 2709×1984
px, all 600 dpi) and matching `.pdf` (vector) — i.e. the assets *are*
publication-quality in principle and meet the usual ≥300 dpi rule. The issue is
packaging, not quality.
**Fix.** For submission, upload figures as separate high-resolution files
(PDF/EPS/TIFF preferred) and place legends in a dedicated "Figure legends"
section or at the end of the manuscript, in numbering order. No panel labels
are needed (all three figures are single-panel).

---

## 5. Numbers I independently recomputed

| Claim | Manuscript / brief | Recomputed (my run) | Match? |
|---|---|---|---|
| Abstract length | 295 words (brief) | **325** body / **329** with headings | **No** |
| Keywords | — | **6** (cancer pain; DRG; CIPN; type-I interferon; in-silico perturbation; negative control) | ok (3–10 typical) |
| Running head | "Divergent DRG programs: bone cancer pain vs CIPN" | **48 characters** | ok (≤50) |
| Title length | — | **202 characters** | long |
| Tables cited in text, in order | all tables cited | Tables 1,2,3 cited; **Table 4 & 5 NOT cited** | **No** |
| Figures cited in order | — | Fig 1 (L49) → Fig 2 (L129) → Fig 3 (L137) in order | yes |
| References cited / present | 16 | 16 present, all 16 cited | yes |
| Reference first-appearance order | sequential | **[14] appears after [15],[16]** | **No** |
| Data & code availability URL | repo link | **HTTP 404** (GitHub API also 404; user 200) | **No** |
| Author degree suffix (MD/PhD/etc.) | none | grep for MD/PhD/MD/MBBS/MSc/… → **no match** | **pass** |
| DOI 10.1038/s41592-025-02772-6 | Ahlmann-Eltze, *Nat Methods* 2025;22:1657–1661 | Crossref: title/journal/year/vol/page/authors all match | yes |
| DOI 10.1038/s41586-025-09370-8 | Baruch et al., *Nature* 2025;646:462–473 | Crossref: matches (Baruch, Gleber-Netto, Nagarajan) | yes |
| DOI 10.1016/j.neuroscience.2023.02.007 | Sun et al., *Neuroscience* 2023;516:1–14 | Crossref: matches | yes |
| DOI 10.3389/fnmol.2021.665085 | Zhai et al., *Front Mol Neurosci* 2021;14:665085 | Crossref: matches (vol 14) | yes |
| Table 2 Cd68 Δ (D7 vs PTX) | +3.42, q = 0.012 | `etiology_marker_cross.csv` Δ = +3.4183, q = 0.012487 | yes |
| Table 3 Glial/immune Δ (D7 vs PTX) | +2.70, q = 0.0073 | `etiology_module_logfc_cross.csv` Δ = +2.6954, q = 0.007345 | yes |
| Table 1 CTRL interferon, day 7 | t = +24.92, q = 6.8 × 10⁻⁵ | `module_level_signal.csv` mean_t = 24.923, q = 0.000068 | yes |
| Figure resolution | "publication-quality" | 600 dpi PNG + vector PDF for all three | yes |
| Main text length | ~6,150 words (brief) | ~6,965 words start→Data availability | differs (context only) |

DOI verification: `curl -L https://doi.org/<doi>` returned HTTP 200 (or 403 at
a bot-blocked LWW page for ref 3) for every DOI spot-checked, and the Crossref
records match the manuscript's titles, journals, years, volumes, pages and
author surnames. The 4-reference spot-check is clean.

---

## 6. Publishability assessment

**Is it ready?** No — not as-is. It is a **Minor Revision** on readiness
grounds: the science is untouched by any of my findings, but three items must
be fixed before an editor would pass it to review (M1 dead data link, M2 two
uncited tables, M3 abstract over-length/under-counted).

**Ordered minimum to reach submission-ready:**

1. **Make the code/data public** (or rewrite the Data-availability statement to
   be honest and supply code to the editor) so the URL resolves — M1. This is
   the highest-priority item because it defends the paper's core claim.
2. **Add in-text call-outs for Table 4 and Table 5** in numerical order — M2.
3. **Bring the abstract within the target venue's cap** and correct the stated
   count — M3.
4. Renumber reference [14] to first-appearance order — m1.
5. Add a standalone Ethics/ARRIVE declaration item — m2.
6. Delete the inaccurate Figure-1 caption clause; supply figures as discrete
   files with a legends list — m3, m5.
7. Optional: shorten the title — m4.

**Desk-reject risk at a soundness-oriented SCIE venue.** Moderate. The paper is
unusually well-defended for a negative result (exact permutation floors,
degree-matched nulls, additive/ridge baselines, a number-checking gate, an
explicit limitations section). The risk is *not* soundness; it is **remit and
framing**. A negative/methodological result about a specific in-silico pipeline
can be desk-rejected as "not of sufficient general interest" unless the venue is
one that explicitly publishes negative, confirmatory or computational-soundness
work, and unless the *transferable* contribution is the headline. The
biological claim is deliberately partial (the interferon programme is shown to
be model-specific), so the durable contribution is the **calibration recipe**:
negative-control modules inside the perturbation itself, degree-matched
permutation nulls, an additive baseline, and force-inclusion of target genes
(§4, line 173).

**Single change that would most reduce desk-reject risk.** Reposition the
submission so the generalisable calibration recipe — not the pipeline's failure
per se — is the lead contribution, and target it at a venue whose scope
explicitly welcomes negative/confirmatory/computational results. Concretely:
state in the cover letter that the paper is a negative-result methodological
calibration and name the venue's relevant article type, and mirror that framing
in the first two sentences of the Abstract and in the title so the editor reads
a reusable standard, not a single failed pipeline. (This is a framing change,
not new analysis.)

---

### Verification honesty notes
- Repository URL: checked live on 2026-10-10; **404 at both the HTML page and
  the GitHub REST API** while the owning user returns 200 — reported as
  unreachable/uncreated, not assumed.
- All four spot-checked DOIs resolved and matched their Crossref metadata; I did
  not use any journal metric (impact factor, CiteScore, acceptance rate) in this
  report, so no metric needed online verification.
