# A4 — Venue / Reporting-Standard Audit

**Reviewer role:** Independent manuscript-review panel — venue & reporting-standard auditor
**Manuscript:** "Divergent dorsal root ganglion programs in bone cancer pain and chemotherapy-induced peripheral neuropathy: a calibrated in-silico perturbation analysis" (single-author original research, first submission)
**Target standard:** Rigorous computational-biology / pain-neuroscience journal
**Audit basis:** `docs/manuscript/manuscript_final.md` (read in full) plus its embedded tables. No other review-round file was consulted.

---

## Overall verdict on ARTICLE TYPE

The "original research" framing is **defensible but under-justified**. Strand A (the bone-cancer-pain vs CIPN DRG divergence, with exact-permutation marker evidence, plus the model-specificity of the interferon programme) is a genuine, novel *comparative* finding. Strand B (the honest-null perturbation result and its three-part calibration prescription) is the more generally portable contribution. The combination justifies original-research status — but only if the authors stop implying they generated primary data and instead state plainly that this is a **re-analysis / comparative re-analysis of two published GEO series plus an independent public single-cell set**. As written, the paper reads at points as if the biology were a fresh experiment. See item A1.

The manuscript does **not** headline its weakest claim (the failed perturbation ranking) as a positive — it is correctly presented as a null. Its strongest claims (model-specificity + the calibration prescription) are present and prominent in the abstract, Discussion and Conclusions. The one real emphasis defect is that the *running head* and the subtitle centre the **method** rather than the **biological divergence**, which risks an editor filing the paper as a methods note. See item A2.

---

## Findings (Problem / Evidence / Why it matters / Specific fix)

### A1. "Original research" is acceptable only if the re-analysis status is declared up front
- **【Problem】** The manuscript never states that both discovery datasets are *previously published* and that no primary data were generated, so the "original" claim is implicitly about novel experiments rather than novel analysis.
- **【Evidence】** `manuscript_final.md:43` "**GSE149648** — rat DRG, Walker-256 mammary carcinoma cells implanted into the tibia… [10]" and `:45` "**GSE160543** — rat DRG, paclitaxel, oxaliplatin or vehicle… [11]" — both cite primary-publication references [10],[11]; the single-cell set GSE292522 (`:92`) is likewise public. Nowhere does the text say "we re-analysed" or "secondary analysis of public data."
- **【Why it matters】** For a rigorous venue, "original research" on re-analysed public data is publishable but must be self-described as such; otherwise the contribution looks overstated and the novelty is ambiguous (is it the data, the comparison, or the calibration?).
- **【Specific fix】** Add one sentence to the end of §2.1 (or the start of §1): *"Both discovery series and the validation set are publicly deposited; this work is a comparative re-analysis of those published datasets, and its originality lies in the cross-aetiology comparison and the calibrated perturbation layer, not in new data generation."*

### A2. The running head and subtitle over-weight the method and under-weight the biological anchor
- **【Problem】** The running head leads with the perturbation method, while the durable contribution is the aetiology divergence; this mismatches the (correct) biological emphasis of the title main clause and Conclusions.
- **【Evidence】** `:3` Running head: "Calibrated in-silico perturbation of cancer-pain DRG transcriptomes"; `:1` title main clause: "Divergent dorsal root ganglion programs in bone cancer pain and chemotherapy-induced peripheral neuropathy" (biology-led) + subtitle "a calibrated in-silico perturbation analysis" (method-led). The Conclusions (`:169`) and abstract (`:23`) correctly lead with biology.
- **【Why it matters】** A method-led running head invites classification as a technical/methods note and can bury the strand-A finding in editorial triage; the two halves of the paper should be visually consistent in emphasis.
- **【Specific fix】** Change the running head to foreground the divergence and the calibration equally, e.g.: *"Divergent DRG programs in bone-cancer pain and CIPN: a calibrated perturbation audit."* Keep the title main clause as-is.

### A3. Code / artefact availability is a HARD FAIL — no resolvable location is given
- **【Problem】** The availability statement points to a "project repository" with no URL, archival DOI, or accession, which fails the mandatory code-availability requirement for a computational paper.
- **【Evidence】** `:175` "Analysis code, artefacts and the test suite are deposited at the project repository, with `requirements.txt` pinning every version used." No link, DOI, or repository name anywhere in the manuscript. (`:97` references a "deposited CSV and JSON artefacts" and a "14-group" test suite, but gives no way to obtain them.)
- **【Why it matters】** A rigorous computational-biology venue will desk-reject or mandate a Machinable/Archived code artefact (GitHub + Zenodo/Figshare DOI, or an accepted repository). "Project repository" is not resolvable and is not archival.
- **【Specific fix】** Replace with a specific, persistent identifier: *"Analysis code and deposited artefacts (CSV/JSON outputs and the 14-group test suite) are available at [GitHub URL], archived with DOI [Zenodo/Figshare DOI]; `requirements.txt` pins numpy 2.4.6, scipy 1.18.1, pandas 2.3.3, statsmodels 0.15.0, matplotlib 3.11.1, Python 3.13."*

### A4. The Conclusions over-generalise "bone cancer pain" and over-claim "model-specific" beyond §3.8
- **【Problem】** The opening conclusion line attributes the interferon programme to "bone cancer pain" generally, while §3.8 explicitly limits it to the Walker-256 tibia model; and "model-specific" is used where "not reproduced in the one alternative (mouse, different-injury) model tested" would be proportionate.
- **【Evidence】** `:169` "Bone cancer pain and chemotherapy-induced peripheral neuropathy drive qualitatively different programmes in the DRG — myeloid and interferon-driven in the first…" vs `:147` "the interferon / myeloid programme is a feature of the Walker-256 tibia model at its day-7 peak, not a universal signature of cancer-pain DRG." Abstract `:23` similarly says "myeloid/interferon-driven in one" before the later "model-specific" qualifier. The evidence for model-specificity is weak concordance: `:145` "Spearman ρ = 0.098 between pseudobulk log₂ fold changes" plus a different dominant programme in B16F10 (a *different species and injury paradigm*, `:92`,`:147`).
- **【Why it matters】** Calling the interferon signal "bone-cancer-pain-driven" in the conclusion while §3.8 calls it model-specific is internally inconsistent and invites the very over-generalisation the Limitations section warns against ("a single-model transcriptional signature must not be generalised to 'cancer pain' as a monolith," `:147`). "Model-specific" also cannot be cleanly separated from a species effect (rat vs mouse) given only one alternative model.
- **【Specific fix】** Scope the claim at first mention: *"In the Walker-256 tibia model, bone cancer pain is myeloid- and interferon-driven, whereas paclitaxel neuropathy is neuronal-injury-driven…"* and replace "model-specific" with *"not a universal or cross-species cancer-pain signature; it was not reproduced in an independent mouse B16F10 model."*

### A5. Zero figures — atypical and weakens presentation for this venue
- **【Problem】** The manuscript contains five tables and **no figures at all**; a computational-biology / pain-neuroscience paper normally carries at least a PCA/sample-correlation, a volcano/DE overview, a module-`t` heatmap, a network-propagation schematic, and a concordance scatter.
- **【Evidence】** Tables 1–5 are present (`:198`–`:280`); a full-text scan finds no `Figure` reference and no figure caption. The body refers to "PCA gave 46.8%, 17.4%…" (`:103`) and "Spearman ρ = 0.991" network-vs-additive (`:121`) without any visual.
- **【Why it matters】** For the stated venue, all-evidence-in-tables reads as a technical report, not a research article, and obscures the very patterns (divergence direction, redundancy of propagation) the authors most want readers to see.
- **【Specific fix】** Add at minimum: (i) a sample-PCA / correlation heatmap (§3.1); (ii) a volcano or bar of day-7 module mean `t` (Table 1 → Fig); (iii) a scatter of network vs additive reversal scores showing ρ=0.991 (§3.4); (iv) a concordance scatter for the B16F10 vs Walker-256 ρ=0.098 (§3.8). Each needs a caption stating what it shows and the panel/statistic.

### A6. Reference [11] is malformed (Vancouver) and DOI presence is inconsistent
- **【Problem】** Reference 11 lists four authors then "et al.", which violates Vancouver (list all if ≤6, or first 3 + et al. if >6); and DOIs are present on [9],[10],[11] but absent on [12],[13],[14], with no stated policy.
- **【Evidence】** `:189` "11. Sun W, Hao Y, Li R, Ho IHT, et al. Comparative transcriptome… *Neuroscience*. 2023;516:1–14. doi:…" — four named + et al. Compare `:184` ref 6 (Loprinzi CL, Lacchetti C, Bleeker J, et al. — correct, >6) and `:190`–`:192` refs [12],[13],[14] with no DOI.
- **【Why it matters】** Inconsistent reference formatting and missing DOIs are a common desk-reject/format-revision trigger at Vancouver-style journals and signal lack of copy-editing care.
- **【Specific fix】** Either list all authors of [11] (if ≤6) or use "Sun W, Hao Y, Li R, et al." (first 3 + et al.); and add DOIs to [12],[13],[14] or state a journal-specific "DOI not assigned" policy. Audit all 14 references for the same first-N + et al. rule.

### A7. No power statement and the pseudobulk @ n=4 limitation is under-discussed
- **【Problem】** There is no statement of achieved statistical power / minimum detectable effect at n=4, and the single-cell *pseudobulk* validation's low-power implication for the "interferon flat" reading is not addressed.
- **【Evidence】** §5 (`:165`) discusses the permutation floor (1/4,901) and that "reversal scores are sensitive to individual animals" but never states power or MDE. For the B16F10 validation, `:145` reports myeloid fraction "~0.1%" and "every interferon-stimulated gene tested was flat (|logFC| < 0.11, q > 0.15)," and `:279` myeloid 0.12% vs 0.05%. A pseudobulk of 4-vs-4 samples with a near-absent myeloid compartment has very low power to detect a real (small) interferon signal — the null there may be a power artefact, not a true absence.
- **【Why it matters】** The §3.8 conclusion that interferon is "model-specific / not universal" rests partly on a *null* in a low-power pseudobulk; if that null is under-powered, the model-specificity claim is weaker than stated. Transparency on power is required for proportionality.
- **【Specific fix】** Add to §5: a one-line power/MDE statement for n=4 permutation tests (e.g., the smallest detectable marker Δ given the exhaustive null), and explicitly state that the B16F10 interferon "flat" result is constrained by n=4 pseudobulk power and a ~0.1% myeloid compartment, so absence of signal is not proof of absence.

### A8. Abstract Results is over-long/dense and the target journal's limits are unstated
- **【Problem】** The abstract is structured correctly (Background/Methods/Results/Conclusions) but the Results subsection is extremely dense (~250 words of nested statistics), risking word/section limits at the (unstated) target venue.
- **【Evidence】** `:21` Results packs: stage-dependent DEG counts, module `t`=+24.9, Cd68/Aif1/Gfap log₂, cross-aetiology Δ values, the 264→132 dedup, minimum q=0.165, 32/132 vs 14 negative-controls, erythrocyte z=−5.78, B16F10 neuronal-injury, ρ=0.991 — all in one paragraph. No target journal or abstract word cap is named anywhere.
- **【Why it matters】** Without a named target venue, the authors cannot self-check against that journal's abstract word/structure rules; an over-long Results paragraph is a frequent format-revision item.
- **【Specific fix】** Name the target journal in the cover letter and trim the Results abstract to its three load-bearing facts (divergence is exact; perturbation fails calibration; interferon is model-specific), moving secondary numbers to the body. Verify the abstract against that journal's word cap (commonly 250–350 words).

### A9. Within-contrast significance (Table 4) vs the "perturbation fails" headline needs explicit reconciliation
- **【Problem】** Table 4 reports several *significant* reversal scores (q=0.011, 0.016, 0.020) against the degree-matched null, which can be misread as "perturbation works," contradicting the paper's headline null — the two are different tests but the manuscript does not flag the distinction prominently.
- **【Evidence】** `:259`–`:266` Table 4: interferon KO q=0.011, glial/immune q=0.011, pan-neuronal q=0.011, K⁺/HCN q=0.011, Na⁺ q=0.020, myelin q=0.020, erythrocyte q=0.016, neuropeptides q=0.016 — all "significant" vs degree-matched null at day 7. Contrast `:133` "no module's reversal Δ reached q < 0.05 on the exact permutation p; the smallest was… q = 0.165" (cross-aetiology) and abstract `:21` "No module's *perturbation* reversal difference survived multiplicity correction."
- **【Why it matters】** A skimming reader (or a methods-averse editor) could cite Table 4 as evidence the method "works," undermining the paper's central methodological point; the distinction (within-contrast vs cross-aetiology; significance vs a null that still can't separate targets from erythrocyte controls) must be unmissable.
- **【Specific fix】** Add a one-sentence bridge at the top of §3.6 (or end of §3.4): *"Table 4 significance is against the degree-matched null within a single contrast; it does not survive the cross-aetiology comparison (§3.6, minimum q=0.165) and, critically, cannot separate targets from negative-control/erythrocyte modules — both facts are required to read Table 4 correctly."*

### A10. The methodological prescription is slightly under-sold as the memorable takeaway
- **【Problem】** The single most generalisable contribution — the "negative-control + degree-matched null + additive baseline (+ force-include targets)" prescription — is buried in Discussion §4 and not elevated to a headline takeaway.
- **【Evidence】** `:159` "Three practices separated a publishable-looking hierarchy from an interpretable result… A fourth is not optional either: force the target genes onto the graph." This is the field-level contribution but the title/running head/abstract-Results do not foreground it; abstract Conclusions (`:23`) mentions it only in the final clause.
- **【Why it matters】** For a computational-biology venue, the calibrated-perturbation prescription is what will be cited; under-selling it lets the paper read as "we tried perturbation and it failed" rather than "here is how to do perturbation honestly."
- **【Specific fix】** Elevate the prescription in the abstract Conclusions to a standalone sentence, e.g.: *"We propose a minimal calibration standard for in-silico perturbation of small-n bulk transcriptomes: negative-control modules carried through the perturbation, degree-matched permutation nulls, a reported additive baseline, and forced inclusion of target genes."*

---

## § Stands up (what is already at standard)

1. **Data accession codes are complete and correct.** GSE149648, GSE160543 (discovery) and GSE292522 (validation) are all named with platform IDs (`:43`,`:45`,`:92`) — a clean reporting win versus many computational submissions.
2. **Software versions are fully pinned.** `:97` gives Python 3.13, numpy 2.4.6, scipy 1.18.1, pandas 2.3.3, statsmodels 0.15.0, matplotlib 3.11.1, and ties a hard import-failure guard to the limma-trend prior — exactly the dependency transparency a methods reviewer wants.
3. **Cross-references and in-text citations are internally consistent.** Every `§3.x` pointer (3.3, 3.5, 3.8, etc.) resolves to an existing section, and every in-text citation [1]–[14] maps to a reference that is also used in text (no orphans, no missing). No broken section citations were found.
4. **The honest-null is reported without softening.** The abstract (`:21`,`:23`) and Conclusions (`:169`) state plainly that the data "do not support module-level analgesic target ranking" and that the reversal score "tracks cellular composition rather than target pharmacology" — proportionality is, on the whole, well-held (see A4 for the one scoping exception).

---

## § Questions for the authors

1. Which **specific journal** is the target? This determines the abstract word cap, reference style细则, and whether a figures-only presentation is acceptable — please state it so the format items above can be checked against the right rules.
2. Will you provide a **resolvable code/artefact DOI** (GitHub + Zenodo/Figshare)? Without it the paper cannot clear a computational-venue mandatory-availability bar (A3).
3. In §3.8, the B16F10 "interferon flat" null comes from an n=4 *pseudobulk* with a ~0.1% myeloid compartment — do you consider that a true absence or a power-limited null, and should the model-specificity claim be softened accordingly (A7)?
4. Table 4 shows several q≈0.01 reversal scores versus the degree-matched null; can you add the explicit within-contrast-vs-cross-aetiology reconciliation sentence (A9) so the null is not misread as a positive?
5. Given both discovery sets are published, would you accept a title/subtitle tweak that labels this a **comparative re-analysis** (A1/A2), or do you intend to claim primary-data novelty?

---

## § What I actually checked

- Read `docs/manuscript/manuscript_final.md` in full (281 lines), including all five tables and the reference list.
- Verified every `§3.x` cross-reference resolves to an existing section (§3.1–§3.8, §4, §5, §6 all present).
- Verified every in-text citation [1]–[14] is present in the reference list and used at least once in text; flagged ref [11] Vancouver malformation and DOI inconsistency (A6).
- Cross-checked numeric consistency between Abstract and body for the headline statistics: DEG counts (0/4,600/235, `:21` vs `:107`), Cd68/Aif1 Δ (+3.42/+3.08, `:21` vs `:125`), 264→132 dedup and minimum q=0.165 (`:21` vs `:133`), ρ=0.991 network-vs-additive (`:21` vs `:121`), B16F10 ρ=0.098 (`:21` vs `:145`) — all consistent.
- Confirmed presence of GSE accessions and pinned software versions (A3 stand-up vs A3 fail on repository URL).
- Confirmed absence of any figure or figure caption (A5) and the method-led running head (A2).
- Assessed article-type fit (A1), conclusion proportionality (A4), abstract coherence vs body (A8), and the headline-vs-evidence emphasis (A10).
- Did **not** read any other file in `docs/review/round1_2026-10-10/` (independence discipline); all judgements are from the manuscript text and tables alone.
