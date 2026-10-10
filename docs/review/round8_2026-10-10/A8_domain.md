# Reviewer A8 — Domain review (pain neuroscience / anaesthesiology)

**Manuscript:** "Divergent dorsal root ganglion programs in bone cancer pain and chemotherapy-induced peripheral neuropathy: an exact cross-study reanalysis and a negative-control calibration of in-silico perturbation"
**Role:** Clinician-scientist in pain neuroscience / anaesthesiology; FRESH, INDEPENDENT read.
**Independence statement:** I did not read any prior review/response files, `author_verification_statement.md`, `GITHUB_DEPOSIT_SOP.md`, `manuscript_v0*.md`, any `.workbuddy/` memory, or any other reviewer's `round8` file. Every judgement below comes from the manuscript text and the source/artefact files I opened myself, plus the four web checks I ran by DOI/title.

---

## Findings

### Finding 1 — Gfap muscle-adjusted value does not match the deposited artefact (Tier 3, non-blocking)
【Problem】 §4 states Gfap was "mildly up in the rat bulk DRG (log₂ +1.33, muscle-adjusted limma-trend; unadjusted +1.50, Table 3)". The unadjusted limb reconciles perfectly, but the muscle-adjusted limb does not match the deposited D7 adjusted DE table.
【Evidence】 `results/tables/DE_D7_vs_Sham_adj.csv` (muscle-adjusted) gives Gfap logFC = +1.294 (t = 7.05, q = 0.0015). `results/tables/etiology_marker_cross.csv` gives Gfap Bone-day-7 = +1.498 (unadjusted), which matches the "+1.50" in Table 3 and the "+1.50" in §4 prose. So the unadjusted reconciliation is exact; only the "muscle-adjusted +1.33" is ~0.04 above the deposited +1.294.
【Why it matters】 It is a trivial numeric inconsistency, but for a reviewer layer whose entire credibility rests on "every number is traced to the artefact," the one prose number that does not match the artefact is the Gfap muscle-adjusted value. It also slightly undercuts the manuscript's own strongest selling point (numerical traceability).
【Specific fix】 "We note Gfap was mildly up in the rat bulk DRG (log₂ +1.29, muscle-adjusted limma-trend; unadjusted +1.50, Table 3); if the +1.33 derives from a different gene set (e.g., the 13,601-gene cross-study filter versus the 13,893-gene discovery set), state that explicitly so the two adjusted values are not read as contradictory."

### Finding 2 — Clinical citation is anachronistic for denosumab (Tier 3, non-blocking)
【Problem】 §4 attributes the anti-resorptive responsiveness of bone-metastases pain to "(bisphosphonates, denosumab) [16]", but ref [16] (Luger NM et al., *J Pain Symptom Manage* 2005) predates denosumab.
【Evidence】 Manuscript §4, sentence: "…responsive to anti-resorptive agents (bisphosphonates, denosumab) [16]…". Reference [16] in the reference list is Luger 2005. Denosumab for bone metastases was approved in the 2010s, so a 2005 review cannot support it.
【Why it matters】 The underlying clinical claim is correct (bone-metastases pain is partially opioid-refractory and does respond to anti-resorptives, including denosumab), so this is a citation-accuracy problem rather than a factual one; but an anaesthesiology/pain readership will notice the anachronism immediately.
【Specific fix】 "…responsive to anti-resorptive agents (bisphosphonates [16]; denosumab)…" (i.e., move denosumab out of the [16] bracket and cite a source that actually discusses denosumab in bone-metastases pain, e.g., the relevant ASCO/NCCN bone-health guidance), or simply drop denosumab from the [16] bracket.

### Finding 3 — The flagship STING → type-I-interferon-in-DRG-neurons literature is missing (Tier 2, should cite, not blocking)
【Problem】 The manuscript's single most important biological caveat is "type-I interferon was inferred, not demonstrated." The natural literature to anchor that inference — and to show it is biologically plausible in bone-cancer-pain DRG — is the Duke/Ji-lab STING work, none of which is cited. Only the secondary Chinese STING paper [17] is cited for the STING/type-I point.
【Evidence】 Donnelly CR, Jiang C, et al., "STING controls nociception via type I interferon signalling in sensory neurons," *Nature* 2021;591:275–280 (doi:10.1038/s41586-020-03151-1) — establishes STING→type-I IFN (IFN-α/β) autocrine/paracrine signalling in DRG nociceptors, with behavioural/electrophysiology validation in a bone-cancer-pain model. Wang K, Donnelly CR, et al., "STING suppresses bone cancer pain via immune and neuronal modulation," *Nat Commun* 2021;12:4558 (doi:10.1038/s41467-021-24867-2) — shows STING agonists suppress bone-cancer-pain and bone destruction in multiple mouse models in an IFN-I/IFNAR-dependent manner. Both are directly on point for the day-7 programme and for the bone-cancer-pain IFN axis specifically.
【Why it matters】 Citing [17] alone leaves the "could this be type-I?" discussion resting on a single secondary paper that itself mainly reports STING/TBK1/NF-κB + IL-1β/IL-6/TNF-α (not type-I IFN directly). The Donnelly/Wang papers are the primary evidence that STING-driven type-I IFN signalling operates in DRG neurons and in bone-cancer pain, so they belong in §3.3 and §5 whether the authors ultimately conclude type-I, type-II, or mixed.
【Specific fix】 Add to §3.3 (and §5): "A STING→type-I-interferon axis in DRG nociceptors has been demonstrated directly in nerve-injury and bone-cancer-pain models (Donnelly et al., 2021; Wang et al., 2021), which makes a type-I component of the day-7 programme biologically plausible even though the present bulk rat data cannot resolve interferon subclasses."

### Finding 4 — Single-nucleus CIPN DRG literature that corroborates the neuronal-injury axis is missing (Tier 2, should cite, not blocking)
【Problem】 The CIPN side of the divergence ("paclitaxel neuropathy is neuronal-injury-driven") is presented only at bulk-marker level. Recent single-nucleus DRG studies of paclitaxel CIPN provide cell-type-resolved evidence for exactly this claim and should be cited to strengthen it.
【Evidence】 Single-nucleus RNA-seq of paclitaxel-treated mouse DRG shows the most pronounced transcriptomic changes in type-C low-threshold mechanoreceptor (C-LTMR) subtypes with suppressed potassium-current/microtubule/mitochondrial programmes, plus universal CAMK1D up-regulation (e.g., Sun W et al., bioRxiv 2023.12.31.573703; the peer-reviewed C-LTMR/CAMK1D report, 2024–2025). These are independent, higher-resolution confirmations that paclitaxel "directly impacts somatosensory neurons" (a conclusion the manuscript's own source study [11] also reaches at bulk level).
【Why it matters】 The manuscript is arguing for a qualitative divergence using bulk data; citing single-nucleus work that localises the paclitaxel signature to neuronal subtypes (and, separately, shows oxaliplatin skews more glial — see Finding 5) makes the "neuronal-injury-driven" half of the divergence considerably more robust and shows the authors know the cell-type-resolution literature.
【Specific fix】 Add to §3.5/§4: "The neuronal-injury axis of paclitaxel neuropathy is supported at single-nucleus resolution, where C-LTMR subtypes and CAMK1D show the largest paclitaxel-driven transcriptomic shifts (Sun et al., 2023/2024), consistent with the bulk-marker picture reported here."

### Finding 5 — "CIPN = neuronal-injury-driven" slightly over-generalises across both chemotherapy agents (Tier 3, precision, non-blocking)
【Problem】 The abstract and §6 frame the second aetiology as a single "neuronal-injury-driven" programme, subsuming both paclitaxel and oxaliplatin under that label. The clean neuronal-injury signature in the data is the paclitaxel arm; oxaliplatin is less completely resolved and its own source study emphasises glial/neuroinflammatory features.
【Evidence】 Manuscript abstract: "neuronal-injury-driven in the second"; §6: "neuronal-injury-driven in the second". Table 3/§3.5 anchor the divergence on Bone-day-7 vs paclitaxel (Cd68 Δ +3.42, galanin Δ −2.59). Oxaliplatin myeloid markers are roughly flat (Cd68 +0.40, Aif1 −0.46). The source study [11] (Sun 2023) concludes "paclitaxel might directly affect somatosensory neurons while oxaliplatin primarily targets dividing cells and immune cells." 
【Why it matters】 Not a contradiction — oxaliplatin is genuinely part of "CIPN" — but the headline implies a single mechanism for both agents when the within-CIPN evidence (both the present data and [11]) is that paclitaxel is neuronal and oxaliplatin is more glial/inflammatory. Precision here prevents a reviewer or reader from over-reading the divergence.
【Specific fix】 "…neuronal-injury-driven in the paclitaxel arm of the second aetiology (oxaliplatin is less completely resolved here and, per the source study [11], skews more glial/neuroinflammatory), and the difference is large, direction-specific and statistically exact."

### Finding 6 — §3.8 shorthand reads "type-I interferon programme" for an inferred attribution (Tier 3, clarity, non-blocking)
【Problem】 §3.8 repeatedly labels the §3.3 finding as "a type-I interferon / myeloid programme at the Walker-256 day-7 peak" and contrasts it with B16F10 being "not type-I interferon," even though §3.3 and §5 explicitly state type-I is inferred, not demonstrated.
【Evidence】 §3.8 (manuscript line ~155): "The Section 3.3 finding, a type-I interferon / myeloid programme at the Walker-256 day-7 peak…"; "…not type-I interferon." §3.3/§5: "type-I is inferred, not demonstrated" / "we could not determine whether it is type-I, type-II, or mixed." The shorthand is internally consistent with the caveat but a reader skimming §3.8 could take "type-I" as established in the claim.
【Why it matters】 This is the exact attribution the brief asked me to police. The claim sections (abstract/§3.3/§4/§6) are clean, but §3.8's descriptive shorthand should mirror the hedged wording so the inferred status is never lost in the results narrative.
【Specific fix】 "The Section 3.3 finding — the inferred type-I / interferon-stimulated-gene programme at the Walker-256 day-7 peak…" and "…the B16F10 model shows no such inferred type-I signal."

---

## What stands up (≥3, with evidence)

1. **The type-I attribution is honestly NOT asserted as demonstrated anywhere in the claims.** Abstract describes the day-7 programme as "an interferon-stimulated-gene (innate-immune) programme" with no type-I claim; §3.3 states "type-I is inferred, not demonstrated (see §5)" and explicitly notes the ISG module contains Gbp2/Ubd/Ifi47, which are also robustly type-II-inducible; §4 writes "an interferon-stimulated-gene programme (type-I not distinguished from type-II)"; §5 repeats "we could not determine whether it is type-I, type-II, or mixed"; §6 says "myeloid and interferon-associated" with no type-I wording. The ISG-module-composition limitation is therefore honestly noted, as required.

2. **The paclitaxel "no myeloid up-regulation" wording is correctly time-point-specific and cites ref [19].** §4: "CD68⁺ macrophage hyperplasia has been documented in paclitaxel DRG at later time points (day 6, ref [19]), so our 'no myeloid up-regulation' statement is time-point-specific and does not claim paclitaxel never recruits macrophages." I verified ref [19] (Peters 2007, Brain Res 1168:46–59): its abstract states ATF3 is up in DRG neurons from day 1 and CD68⁺ macrophage hyperplasia in the DRG is "not observed until day 6 post-infusion." The hedge is exactly correct.

3. **The headline divergence is numerically traced and internally consistent.** Deposited `DE_D7_vs_Sham_adj.csv` gives Cd68 +2.862 and Aif1 +2.177 (muscle-adjusted), matching §3.3's "+2.86"/"+2.18" to the stated precision; `etiology_marker_cross.csv` gives Cd68 Bone-D7 +2.839 / PTX −0.579 / Δ +3.418 / q 0.012 and Aif1 Δ +2.838 / q 0.012, matching Table 3 and the abstract (Cd68 Δ +3.42, Aif1 Δ +3.08, q 0.012). The §2.3 reconciliation note ("Cd68 2.86 versus 2.84; Aif1 2.18 versus 2.19") is borne out. This is a genuinely reproducible cross-study result, not a tuning artefact.

4. **All four added references [17][18][19][14] are real, correctly attributed, and used accurately.** [17] (Zhang 2023, Inflamm Res, doi:10.1007/s00011-022-01663-2) genuinely used Walker-256 intratibial rats and reported STING/TBK1/NF-κB activation in DRG neurons; [18] (Wang 2025, PNAS 122(34):e2503779122) genuinely reports macrophages infiltrating DRG/sciatic nerve in a cancer-induced bone pain model; [19] (Peters 2007) genuinely reports ATF3-in-neurons (day 1+) and day-6 CD68⁺ hyperplasia; [14] (Baruch 2025, Nature 646:462–473) genuinely reports neuron-intrinsic type-I-interferon inflammation in injured neurons. No misattribution.

5. **The §3.8 Baruch [14] reframing is honest.** The manuscript correctly notes agreement with [14] on a neuronal-injury programme (Sprr1a/Atf3) and a genuine, narrower discrepancy on neuronal ISG up-regulation, and it lists two explicit caveats (ATF3-downstream injury response; pseudobulk sensitivity / CD45-depleted preparation). It does not over-claim resolution of the disagreement.

6. **The §4 clinical sentence is hypothesis-generating, not over-claimed.** It states the split is "consonant with, but does not by itself establish" the clinical bone-metastases-vs-CIPN distinction, explicitly flags the translation as "hypothesis-generating," and notes "no human validation exists." That is the correct posture for a rat n = 4 bulk reanalysis.

---

## Questions for the authors

1. §4 cites [16] (Luger 2005) for denosumab in bone-metastases pain. Did you intend a different citation for denosumab, or should it be removed from the [16] bracket (see Finding 2)?
2. The muscle-adjusted Gfap in §4 prose is +1.33, but the deposited D7 adjusted DE gives +1.294. Which value does your number-checker assert, and should the prose be aligned to the artefact (Finding 1)?
3. Given oxaliplatin's own source study [11] emphasises glial/neuroinflammatory features for oxaliplatin (vs neuronal for paclitaxel), should the "CIPN = neuronal-injury-driven" headline in the abstract/§6 be scoped to the paclitaxel arm, with oxaliplatin flagged as less completely resolved (Finding 5)?
4. Would you consider adding Donnelly 2021 (*Nature*) and Wang/Donnelly 2021 (*Nat Commun*) to §3.3/§5, since they directly bear on whether the day-7 programme could be type-I and on the bone-cancer-pain IFN axis specifically (Finding 3)?
5. The erythrocyte-control Δ (−2.48) is interpreted as "a direct measurement of how much the two DRG preparations differ in blood content." The two studies also differ in dissection/perfusion protocol — is erythrocyte transcript the only plausible driver, or could part of that Δ reflect non-erythrocyte dissection differences? (Minor; the compositional argument does not depend on it.)

---

## What I actually checked

**Files read (independently, this round):**
- `docs/manuscript/manuscript_final.md` — full text, including the lines truncated in the primary display (recovered via PowerShell: the §3.8 tail "…a single-model transcriptional signature must not be generalised to 'cancer pain' as a monolith…" and the §4 tail "…partially refractory to opioids and responsive to anti-resorptive agents (bisphosphonates, denosumab) [16]… hypothesis-generating").

**Data artefacts verified with the managed Python (`C:/Users/Administrator/.workbuddy/binaries/python/envs/default/Scripts/python.exe`), `PYTHONIOENCODING=utf-8`:**
- `results/tables/DE_D7_vs_Sham_adj.csv`: Cd68 +2.862 (t 18.05, q 1.6e-5), Aif1 +2.177 (t 14.05, q 5.2e-5), Gfap +1.294 (t 7.05, q 0.0015), Scn10a −0.387, Gal −0.692, Calca −0.201.
- `results/tables/etiology_marker_cross.csv`: Cd68 Bone-D7 +2.839 / PTX −0.579 / Δ +3.418 / q 0.012; Aif1 Δ +2.838 / q 0.012; Gal Bone-D7 −0.672 / PTX +1.920 / Δ −2.592 / q 0.012; Gfap Bone-D7 +1.498; Atf3 PTX +1.122.

**Web checks run (by DOI / title):**
- [17] doi:10.1007/s00011-022-01663-2 — REAL; Walker-256 intratibial rat; STING/TBK1/NF-κB in DRG neurons; cytokines IL-1β/IL-6/TNF-α; cites Donnelly 2021 for type-I IFN. Accurate.
- [18] doi:10.1073/pnas.2503779122 — REAL (title confirmed via PNAS + PubMed); mouse CIBP; macrophages infiltrate DRG/sciatic nerve, M1 augmented; two-phase SP/CCL2–CCL3 paradigm. Accurate.
- [19] doi:10.1016/j.brainres.2007.06.066 — REAL; ATF3 up in DRG neurons from day 1 through day 10; CD68⁺ macrophage hyperplasia in DRG not until day 6. Both manuscript uses accurate.
- [14] doi:10.1038/s41586-025-09370-8 — REAL; injured neurons autonomously initiate IL-6- and type-I-interferon-mediated inflammation; ATF3/Ifnar1. Supports §3.8 reframing. Accurate.
- Donnelly 2021 doi:10.1038/s41586-020-03151-1 — REAL (verified for the must-cite in Finding 3).
- Wang/Donnelly 2021 doi:10.1038/s41467-021-24867-2 — REAL (verified for the must-cite in Finding 3).

**Discrepancies found:** Only the Gfap muscle-adjusted prose value (+1.33) versus the deposited artefact (+1.294); every other checked number (Cd68, Aif1, Gal, Atf3, Gfap unadjusted, the cross-aetiology Δs and qs) reconciles exactly. No reference misattribution. No clinically or mechanistically false statement identified.

---

## VERDICT

ACCEPTABLE — no remaining domain blocker.

(The six findings are Tier 2/3: two should-cite additions (Donnelly 2021; Wang/Donnelly 2021; and a single-nucleus CIPN DRG study) and four minor precision/citation/clarity fixes. None is a domain blocker — there is no false biological claim, the type-I attribution is properly hedged as inferred-not-demonstrated throughout the claims, the paclitaxel myeloid statement is correctly time-point-specific and cited, the divergence is numerically traced to the deposited artefacts, and the references are accurate. Per the review brief, an honestly hedged, numerically traced reanalysis should not be failed for the absence of wet-lab validation.)
