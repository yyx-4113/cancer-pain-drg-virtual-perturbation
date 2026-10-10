# Reviewer B1 — Domain review (clinician-scientist, pain neuroscience / anaesthesiology)

**Manuscript:** *Divergent dorsal root ganglion programs in bone cancer pain and chemotherapy-induced peripheral neuropathy: an exact cross-study reanalysis and a negative-control calibration of in-silico perturbation* (single author, Y. Yang)

**Charge:** Fresh, independent first-pass review. Verify the four claimed biological fixes and look for residual or newly-introduced biological over-statements. I read the manuscript in full and web-verified the four disputed references ([14], [17], [18], [19]).

---

## Overall verdict

The methodological contribution (negative-control calibration of in-silico perturbation; exact-permutation cross-aetiology comparison) is strong and well-defended. The *biological* claims are much better than the prior-round description you relayed — the §3.3 reframing away from "interferon event, not ion-channel event," the §3.8 Baruch [14] reframing, the §4 Gfap caution, and the §4 clinical sentence are each, in their own section, honestly written. **However, one over-statement was not just left in place but sharpened into a more specific and still-unsupported claim: "type-I interferon."** The abstract, §3.3 (by forward-reference), §4 and §6 now label the day-7 programme "type-I interferon," but the manuscript never distinguishes type-I from type-II interferon, the ISG module it scores actually *contains type-II-associated genes*, and the §3.3→§3.8 "probe" of this question points to a section that does not perform it. This is the single most important item below (F1) and should block acceptance until corrected.

---

## Findings (each item: Problem / Evidence / Why it matters / Specific fix)

### F1 — The "type-I interferon" attribution is never validated, and the ISG module itself contains type-II-associated genes

【Problem】 The manuscript repeatedly calls the day-7 programme a *type-I* interferon-stimulated-gene programme, but it never distinguishes type-I (IFN-α/β) from type-II (IFN-γ) interferon, and the very gene set used to score the programme includes canonical type-II-IFN-inducible genes.

【Evidence】 Abstract (line 21): "Day 7 is dominated by a **type-I** interferon-stimulated-gene programme (module t = +24.9, q = 6.8 × 10⁻⁵)." §3.3 heading (line 113) and text (line 115) frame it as "interferon-stimulated-gene programme," but §3.3 then states "whether the ISG signature reflects a **type-I** interferon programme specifically (versus a broader cytokine-driven innate-immune response) is probed in §3.8." §4 (line 165) asserts "**type-I interferon** signalling," and §6 (line 183) calls BCP "myeloid and **interferon-driven**." The actual ISG control module is defined in Table 1 footnote as *Irf1, Ifitm3, Ifi47, Bst2, Rsad2, Gbp2, Ubd, Ly6e*. Of these, **Gbp2** (guanylate-binding protein 2), **Ubd** (ubiquitin D / FAT10) and **Ifi47** (IRG/immunity-related GTPase) are classically and strongly induced by **type-II IFN (IFN-γ)** at least as much as by type-I; several ISG panels use Gbp2/Ubd as *selective type-II markers*. No interferon gene itself (Ifnb1, Ifna family, Ifng) and no IRF7/STAT1 phosphorylation/IRF9 evidence is presented anywhere. §3.8 does **not** contain a type-I-vs-type-II probe — it compares Walker-256 (rat, bulk) with B16F10 (mouse, scRNA) for *model specificity* and discusses the Baruch neuron-intrinsic discrepancy. So the question §3.3 says is "probed in §3.8" is not answered there.

【Why it matters】 This is precisely the class of claim the prior round flagged (over-stating an "interferon event"). The author has converted a vague "interferon" into a specific "type-I interferon" without adding the evidence a specific claim requires. A reader taking the abstract at face value will cite this as showing type-I IFN in Walker-256 DRG neurons, which the data cannot support. It also creates an internal contradiction: §3.3 (line 115) says "whether … type-I … specifically … is probed in §3.8," yet §3.8 (lines 153–159) never performs that probe. Credibility of the headline biological claim is at stake.

【Specific fix】
- Abstract (line 21): replace "Day 7 is dominated by a type-I interferon-stimulated-gene programme (module t = +24.9, q = 6.8 × 10⁻⁵)" with:
  > "Day 7 is dominated by an interferon-stimulated-gene (ISG) programme (module t = +24.9, q = 6.8 × 10⁻⁵); we infer a likely type-I signature from the canonical-ISG composition but did not distinguish type-I from type-II interferon in these bulk data."
- §3.3 (line 115): replace "whether the ISG signature reflects a type-I interferon programme specifically (versus a broader cytokine-driven innate-immune response) is probed in §3.8" with:
  > "we did not distinguish a type-I from a broader (e.g., type-II/interferon-γ or cytokine-driven) innate-immune programme; this is a limitation because the interferon-stimulated-gene control module (Irf1, Ifitm3, Ifi47, Bst2, Rsad2, Gbp2, Ubd, Ly6e) contains genes (Gbp2, Ubd, Ifi47) that are also robustly induced by type-II interferon, so 'type-I' is inferred, not demonstrated."
- §4 (line 165): replace "type-I interferon signalling and myeloid marker up-regulation (Cd68, Aif1)" with:
  > "an interferon-stimulated-gene programme (most likely type-I, but type-I was not distinguished from type-II interferon) and myeloid marker up-regulation (Cd68, Aif1)."
- §6 (line 183): replace "myeloid and interferon-driven" with "myeloid and interferon-associated."

---

### F2 — Abstract says the interferon programme is "absent" from B16F10, but §3.8 reports ISGs are down-regulated, not absent

【Problem】 The abstract over-states the B16F10 validation as "interferon absent," contradicting its own §3.8, which reports several ISGs significantly *down*-regulated.

【Evidence】 Abstract (line 21): "The interferon program is absent from the B16F10 model, whose response is neuronal-injury-driven." §3.8 (line 155, Table 5 line 329): "Interferon-stimulated genes were not uniformly flat: **Ifit2 and Oas1c were significantly down-regulated (q = 0.042 and 0.014) and no ISG was up-regulated**"; Table 5 also lists "Ifit2↓ neurons (q = 0.046), Isg15↓ satellite glia (q = 0.007)." So the correct statement is "no ISG up-regulation; some ISGs down-regulated," not "absent."

【Why it matters】 "Absent" implies zero ISG signal, which is false and could be seized on by a sceptical reviewer as a factual error in the abstract that is contradicted by the body. It also slightly inflates the contrast the author is trying to draw (Walker-256 ISG-up versus B16F10 ISG-flat/down).

【Specific fix】 Abstract (line 21): replace "The interferon program is absent from the B16F10 model" with:
> "The interferon up-regulation seen at Walker-256 day 7 is absent from the B16F10 model (its ISGs are flat or down-regulated)."

---

### F3 — Gfap log₂ value is internally inconsistent between §4 prose (+1.33) and Table 3 (+1.50)

【Problem】 §4 states Gfap is up +1.33 in the rat bulk DRG, but Table 3 (the reader-checkable source) lists Gfap +1.50 at bone day 7; the discrepancy is unexplained.

【Evidence】 §4 (line 165): "We note Gfap was mildly up in the rat bulk DRG (log₂ **+1.33**)". Table 3 (line 300): "Gfap | **+1.50** | +0.59 | +0.09 | +0.05 | — | ns". No footnote reconciles the two. For comparison, §2.3 (line 59) explains that muscle-adjusted and unadjusted logFC differ only by hundredths for well-expressed markers (Cd68 2.86 vs 2.84; Aif1 2.18 vs 2.19), so a +1.33-vs-+1.50 gap of 0.17 is far outside that documented range and cannot be the muscle-adjusted value. The source of +1.33 is unidentified.

【Why it matters】 A numeric mismatch between a prose claim and its own table undermines trust in every other number. It also matters biologically: the entire §4 Gfap caution ("we do not claim glial activation as established") hinges on the reader accepting the quoted magnitude; a wrong magnitude in the caution sentence is self-defeating.

【Specific fix】 §4 (line 165): replace "(log₂ +1.33)" with "(log₂ +1.50, Table 3)" — or, if +1.33 is a deliberately muscle-adjusted value, add an explicit qualifier: "(log₂ +1.33, muscle-adjusted; unadjusted +1.50, Table 3)."

---

### F4 — §4 labels the biological conclusion "secure" while embedding a component the same section admits is model-specific

【Problem】 The §4 biological conclusion opens as "secure," yet it asserts "type-I interferon signalling" as part of that secure BCP programme, while the very next sentences concede the interferon signal is model-specific to Walker-256 and the type-I attribution is unvalidated (see F1).

【Evidence】 §4 (lines 163–165): "The biological one is **secure**. … bone cancer pain at its transcriptional peak is an innate-immune event: **type-I interferon signalling** and myeloid marker up-regulation (Cd68, Aif1) …. The independent B16F10 single-cell model … reproduces the neuronal-injury programme … but not the interferon programme, so **the interferon signal is model-specific**." The robust, cross-validated claim is the *divergence* (myeloid-up in BCP vs neuronal-injury in CIPN). The interferon programme is explicitly model-specific and, per F1, its "type-I" qualifier is unsupported. Calling the package "secure" while it contains a model-specific, type-unvalidated component is internally inconsistent.

【Why it matters】 Reviewers will read "secure" as endorsing the type-I interferon programme as established biology, which the manuscript itself qualifies two sentences later. This over-reach is what earlier rounds objected to, and it survives in the framing word.

【Specific fix】 §4 (line 163): replace "The biological one is secure. Bone cancer pain and paclitaxel neuropathy drive qualitatively different programmes in the DRG. In the Walker-256 tibia model, bone cancer pain at its transcriptional peak is an innate-immune event: type-I interferon signalling and myeloid marker up-regulation (Cd68, Aif1), with the neuronal compartment diluted rather than specifically remodelled." with:
> "The biological one is secure in its broad shape: bone cancer pain and paclitaxel neuropathy drive qualitatively different programmes in the DRG. In the Walker-256 tibia model, bone cancer pain at its transcriptional peak is an innate-immune event — an interferon-stimulated-gene programme (type-I not distinguished from type-II) and myeloid marker up-regulation (Cd68, Aif1), with the neuronal compartment diluted rather than specifically remodelled."

---

### F5 — Paclitaxel "no myeloid infiltration at all" is over-stated relative to the author's own cited literature [19]

【Problem】 §4 claims paclitaxel neuropathy shows "no myeloid infiltration at all," but the very reference cited two sentences later for ATF3 ([19] Peters 2007) documents CD68⁺ macrophage hyperplasia/hypertrophy in the DRG by day 6, and the claim is time-point-dependent on the author's own bulk sampling.

【Evidence】 §4 (line 165): "Paclitaxel neuropathy is a neuronal injury event, with galanin, CGRP and Atf3 up and **no myeloid infiltration at all** — Aif1 and Cd68 are if anything slightly down." The supporting data are Table 3 (line 293–294): Aif1 −0.88, Cd68 −0.58 in paclitaxel. But [19] (PubMed 17698044; verified) states: "hypertrophy and hyperplasia of CD68(+) macrophages in the DRG and sciatic nerve … were not observed until **day 6** post-infusion" after ATF3 up-regulation in neurons from day 1. So macrophage recruitment in paclitaxel DRG is a *later* event, exactly the kind of time-course the manuscript's single bulk snapshot cannot rule out. Citing [19] for ATF3 while ignoring its macrophage finding is selective.

【Why it matters】 "No myeloid infiltration at all" will collide with the established paclitaxel-CIPN literature (macrophage/monocyte involvement is documented) and with [19] itself. It also undercuts the clean BCP-vs-CIPN binary the paper's headline rests on, because the difference may be one of *timing* (early neuronal injury in both; myeloid influx later and perhaps larger in BCP) rather than *kind*.

【Specific fix】 §4 (line 165): replace "Paclitaxel neuropathy is a neuronal injury event, with galanin, CGRP and Atf3 up and no myeloid infiltration at all — Aif1 and Cd68 are if anything slightly down" with:
> "At the time point profiled here, paclitaxel neuropathy is dominated by a neuronal-injury signal — galanin, CGRP and Atf3 up — with myeloid markers (Aif1, Cd68) if anything slightly down; we note that CD68⁺ macrophage hyperplasia has been documented in paclitaxel DRG at later time points (day 6, ref [19]), so our 'no myeloid up-regulation' statement is time-point-specific and does not claim paclitaxel never recruits macrophages."

---

### F6 — The §3.3 "probe" of type-I specificity points to §3.8, which does not perform it (focus item #1, explicit)

【Problem】 This is the explicit focus item #1: the §3.3 reframing is *otherwise* defensible (see Stands-up S1), but its promise that type-I-vs-broader-innate specificity "is probed in §3.8" is false — §3.8 addresses model specificity and the Baruch discrepancy, not type-I specificity.

【Evidence】 §3.3 (line 115): "whether the ISG signature reflects a type-I interferon programme specifically (versus a broader cytokine-driven innate-immune response) is probed in §3.8." §3.8 (lines 153–159): content is (i) B16F10 shows neuronal-injury not type-I IFN; (ii) per-compartment ISGs flat/down; (iii) the Baruch [14] neuron-intrinsic IFN-I discrepancy. None of these tests type-I vs type-II. The "probe" therefore does not exist; the sentence is a dangling forward-reference.

【Why it matters】 A reader sent to §3.8 to resolve the type-I question finds only a model-comparison and a citation discussion. This is the mechanism by which the unsupported "type-I" claim persists (it is asserted in the abstract and §4, "deferred" to §3.8, and never actually adjudicated). Fixing F1's wording without removing this forward-reference leaves the illusion of a test that was not run.

【Specific fix】 Same replacement as F1's §3.3 fix (delete the "is probed in §3.8" sentence and convert it into an explicit limitation noting the ISG module contains type-II-associated genes). Optionally add to §5 Limitations a sentence: "Because bulk rat DRG does not resolve interferon subclasses, we could not determine whether the day-7 ISG programme is type-I, type-II, or mixed; a targeted assay of Ifnb1/Ifng and IRF7/STAT1 signalling would be required."

---

## § Stands up (things I suspected were still broken but found correctly fixed)

**S1 — §3.3 reframing away from "interferon event, not ion-channel event" is now scientifically defensible.** The new heading (line 113) "prominent interferon-stimulated-gene programme and myeloid marker up-regulation, on a background of neuronal-transcript dilution" is supported by the numbers: ISG module t = +24.9 (Table 1, line 244); Cd68 +2.86 and Aif1 +2.18 rising (§3.3, line 115); and the neuronal modules all falling by near-identical amounts with the compositional controls inside that range (GABAergic −5.71 … pan-neuronal control −5.20 … erythrocyte −4.80, §3.3 line 117; Table 1). The dilution argument (myeloid influx lowers every neuronal gene's transcript share) is a legitimate inference from the compositional controls, and the author no longer claims selective ion-channel remodelling. This is a genuine, correct fix.

**S2 — §3.8 Baruch [14] reframing is honest and accurate.** I web-verified the paper (Baruch EN et al., *Nature* 2025;646:462–473, doi:10.1038/s41586-025-09370-8). Its abstract states: "The injured neurons respond by autonomously initiating IL-6- and type I interferon-mediated inflammation to promote nerve healing and regeneration," and the mechanism is "conditional knockout of … Atf3 … knockout of interferon-α receptor signalling (Ifnar1)." So the manuscript's representation (line 159: "ref [14] reports a neuron-intrinsic type-I-interferon (IFN-I) response … downstream of ATF3 and part of a nerve-healing/regeneration programme") is faithful. The manuscript's own B16F10 data show Atf3 +0.84 and Sprr1a +1.55 (Table 5, line 328) — agreeing on the neuronal-injury programme — while Ifit2 is down in neurons (q = 0.046), which genuinely *contradicts* a neuronal type-I IFN transcriptional up-regulation (if IFN-I were up in neurons, ISGs such as Ifit2 would be up, not down). The new wording "the studies agree on a neuronal-injury programme and differ on whether ISGs are transcriptionally up-regulated in neurons at our resolution" is a fair, correctly hedged statement, and the author honestly notes the caveat that pseudobulk averaging across all neurons limits sensitivity to a sub-neuronal population and that the myeloid compartment (where [14]'s programme could also reside) is untestable due to CD45 depletion. This is a clean correction.

**S3 — §4 Gfap caution is now consistent with §2.8/§3.8.** §4 (line 165) drops the prior "(Gfap up)" from the secure conclusions and adds: "in the companion mouse dataset the satellite-glia compartment carried almost no Gfap transcript, consistent with GFAP being a weak marker of DRG glia, so we do not claim glial activation as established." This matches §2.8 (line 97: "the satellite-glia bin carries almost no Gfap transcript") and §3.8 (line 159: "7 unique-molecule counts across all eight animals, zero in five of them"). The caution is adequate and internally consistent. (Only the magnitude +1.33 is wrong — see F3.)

**S4 — §4 clinical sentence is now hypothesis-generating and the citations are appropriate.** §4 (line 165): "The split is consonant with, but does not by itself establish, the clinical distinction that bone-metastases pain is partially refractory to opioids and responsive to anti-resorptive agents (bisphosphonates, denosumab) [16] whereas CIPN is managed with agents such as duloxetine [1,2]; our data suggest a mechanistic correlate … that future human or ex-vivo DRG work should test directly." The "does not by itself establish" and "future work should test directly" framing is appropriately non-over-claimed. [16] Luger 2005 (*J Pain Symptom Manage*) is a bone-cancer-pain mechanism/therapy review supporting anti-resorptive/bisphosphonate therapy; [1] (ESMO 2018) and [2] (NCCN 2019) are the standard cancer-pain guidelines that list duloxetine for CIPN. Both citation pairs are reasonable for what is actually stated.

**S5 — All three new references [17][18][19] are REAL, correctly formatted, and support their claims.** (a) [17] Zhang Y et al., *Inflammation Research* 2023;72(1):117–132, doi:10.1007/s00011-022-01663-2 — verified via PubMed (PMID 36346430) and the Springer PDF: it uses *intratibial Walker 256* in rats and shows "STING/TBK1/NF-κB pathway activation … in the DRGs," with type-I interferon in its vocabulary (it cites Donnelly 2021, "STING controls nociception via type I interferon signalling"). So "a STING/TBK1/NF-κB innate-immune programme upstream of type-I interferon … in Walker-256 DRG neurons" (§3.3, line 115) is accurately supported. (b) [18] Wang C et al., *PNAS* 2025;122(34):e2503779122, doi:10.1073/pnas.2503779122 — verified (PMID 40828016): "Macrophages were found to infiltrate the DRGs and the sciatic nerves … in mice with CIBP." So "macrophage infiltration of the DRG is a documented feature of cancer-induced bone pain" (§3.3, line 115) is accurate. (c) [19] Peters CM et al., *Brain Res* 2007;1168:46–59, doi:10.1016/j.brainres.2007.06.066 — verified (PMID 17698044): "up-regulation of ATF3 in a subpopulation of large and small DRG neurons … present through day 10." So "Atf3 … elevated in paclitaxel-treated DRG neurons [19]" (§3.5, line 131) is accurate. One minor note (non-blocking): [18] is a *mouse* CIBP model, not the rat Walker-256 discovery model, so it corroborates the general phenomenon rather than replicating the exact dataset — acceptable as supporting literature, but the manuscript could say "in a mouse CIBP model" to be precise.

---

## § Questions for the authors (please answer; do not leave to inference)

1. **Type-I vs type-II:** Did you examine the interferon genes themselves (Ifnb1, Ifna family, Ifng) or any type-I-selective readout (e.g., IRF7, Mx1, Isg15 kinetics, STAT1/IRF9) in the Walker-256 day-7 bulk data? If not, what justifies the "type-I" label beyond module composition? (See F1.)
2. **Gfap magnitude:** Is the +1.33 in §4 (line 165) a muscle-adjusted or otherwise derived value, or a typo for the Table 3 +1.50? Please state the exact source so the text and table can be reconciled (F3).
3. **Paclitaxel time point:** What post-injection day was the paclitaxel group in GSE160543 sampled at? [19] shows CD68⁺ macrophage hyperplasia only from day 6, so the "no myeloid infiltration" claim's validity depends entirely on your sampling time (F5).
4. **Myeloid fraction in Walker-256 bulk:** Your bulk Cd68 +2.84 / Aif1 +2.19 (Table 3) is interpreted as myeloid influx, and you cite [18] (mouse CIBP macrophage infiltration). Do you have any direct cellular-resolution or IHC evidence for macrophage infiltration in the *rat Walker-256* model itself, or is the myeloid claim purely inferential from bulk marker up-regulation? (This bears on how "secure" the innate-immune conclusion can be.)
5. **B16F10 ISG down-regulation mechanism:** In §3.8 you report Ifit2/Oas1c down-regulated and "no ISG up-regulated," yet Baruch [14] reports neuron-intrinsic IFN-I *up* in the same B16F10 model. Before concluding "model-specific," have you considered that the IFN-I signal in [14] may reside in a rare neuronal subpopulation (e.g., non-peptidergic nociceptors) that your pseudobulk averaging (and CD45 depletion) would mask? Your current text acknowledges the resolution limit but frames the difference as definitive; please clarify the strength of that inference.

---

## § What I actually checked

- **Read in full:** `manuscript_final.md` (all 351 lines, including Abstract, §1–§6, References [1]–[19], Tables 1–5, Figures 1–3 captions).
- **Marker-direction internal consistency (focus #1 / mandatory recompute #1):** Confirmed Cd68 (+2.86 muscle-adj §3.3 / +2.84 Table 3 / +2.84 §3.5), Aif1 (+2.18 / +2.19 / +2.19), Gal/Galanin (−0.67 / −0.67, +1.92 / +1.92, Δ −2.59), Scn10a (−0.46 / −0.46, +0.24 / +0.24, Δ −0.70) are mutually consistent across §3.3, §3.5 and Table 3. **Only Gfap is inconsistent: +1.33 (§4) vs +1.50 (Table 3) — see F3.**
- **Ref [14] web-verification (mandatory recompute #2):** Confirmed Baruch EN et al., *Nature* 2025;646:462–473, doi:10.1038/s41586-025-09370-8, reports neuron-intrinsic type-I IFN downstream of ATF3 for nerve healing; confirms the manuscript's §3.8 representation is accurate and that the manuscript's own Ifit2-down-in-neurons (q = 0.046) genuinely contradicts a neuronal IFN up-regulation.
- **Refs [17][18][19] DOI spot-checks (mandatory recompute #3):** All three DOIs resolve to the papers described — [17] Inflammation Research 2023;72:117–132 (Walker-256 rat STING/TBK1/NF-κB in DRG); [18] PNAS 2025;122(34):e2503779122 (macrophage infiltration of DRG in mouse CIBP); [19] Brain Res 2007;1168:46–59 (ATF3 up-regulation in paclitaxel DRG neurons, plus CD68⁺ macrophage hyperplasia by day 6). All support the claims they are cited for, with the caveats noted in S5 and F5.
- **Biological-mechanism re-scan for new over-statements:** Identified the sharpened-but-unsupported "type-I interferon" labelling (F1), the "absent" vs "down-regulated" mismatch (F2), the "secure" framing enclosing a model-specific component (F4), and the paclitaxel myeloid over-statement (F5). The abstract's "type-I interferon-stimulated-gene programme" does **not** honestly match §3.3/§3.8, which at most show a general ISG programme whose type is unvalidated — this is the headline mismatch the prior round's spirit required fixing and which remains.

---

## Bottom line

Two of the four claimed fixes are clean and verified (§3.8 Baruch reframing; §4 Gfap caution). The §3.3 reframing is correct in substance but tainted by a dangling "type-I specificity probed in §3.8" promise that §3.8 does not fulfil. The §4 clinical sentence is appropriately hypothesis-generating. **The blocker is the unsupported "type-I interferon" claim (F1/F2/F4/F6):** it is asserted in the abstract and conclusions, "deferred" to a section that does not test it, and is contradicted by the composition of the author's own ISG module (which contains type-II-associated genes). Correct F1–F6 and this is a solid domain contribution; as written, the headline biological label over-reaches.
