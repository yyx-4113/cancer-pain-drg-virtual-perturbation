# Independent review — Domain layer (pain neuroscience / DRG biology)

**Reviewer role:** independent domain expert (dorsal-root-ganglion biology, cancer pain, CIPN)
**Manuscript:** "Divergent dorsal root ganglion programs in bone cancer pain and chemotherapy-induced peripheral neuropathy: a calibrated in-silico perturbation analysis"
**Manuscript file reviewed:** `docs/manuscript/manuscript_final.md`
**Status of independence:** treated as a first submission; no other reviews, response letters, manifests, or project-status files were read.
**Verdict (domain layer):** *Major revision.* The biological backbone — that bone-cancer pain and CIPN drive qualitatively different DRG programmes, and that the in-silico perturbation layer is a composition/topo artefact — is, on the whole, sound and in several places genuinely rigorous. But the manuscript over-states a "satellite gliosis" claim from a single marker, conflates myeloid-support with interferon-support, labels a mouse-vs-rat confound as "model-specific," and contains a concrete transcription error in Table 4 that misrepresents two control modules. These are fixable but must be fixed before acceptance.

---

## Issues

### Issue 1 — Table 4 misreports two control modules (knockout reversal copied from overexpression)
【Problem】 Two rows of Table 4 ("Reversal (knockout)" column) report the *overexpression* values, not the knockout values, flipping the sign of the erythrocyte and neuropeptide controls.

【Evidence】
- `virtual_perturbation_scores.csv`, contrast D7, perturbation `knockout`:
  - CTRL erythrocyte → `reversal_network = -0.0744`, `z_vs_perm = -2.568`, `pct_of_oracle = 32.98`, `p_empirical = 0.00599`, `q_empirical = 0.0195`.
  - Neuropeptides → `reversal_network = -0.0675`, `z_vs_perm = -2.449`, `pct_of_oracle = 26.85`, `p_empirical = 0.00799`, `q_empirical = 0.0195`.
- The manuscript Table 4 reports for these rows `+0.074 / +2.52 / 33.0 / 0.0040 / 0.016` (erythrocyte) and `+0.067 / +2.56 / 26.9 / 0.0040 / 0.016` (neuropeptides). These numbers match the **overexpression** rows in the CSV exactly (`+0.0744 / +2.519 / 32.98 / 0.0040 / 0.016` and `+0.0675 / +2.557 / 26.85 / 0.0040 / 0.016`), not the knockout rows. The other six rows of Table 4 match the knockout column correctly, so this is isolated to two control rows.
- Biological consistency check: the erythrocyte module is a *down* module at day 7 (`module_level_signal.csv`: mean logFC `-1.147`), so a knockout (further down) must give a *negative* reversal — the CSV knockout value is correct, the Table-4 printed value is wrong.

【Why it matters】 This is not cosmetic. The erythrocyte and neuropeptide controls are the linchpin of the compositional-artifact argument in §3.4/§3.6: the authors need them to behave like the neuronal "down" modules. Table 4 as printed makes both look like *positive* (interferon-like) reversals, which is internally inconsistent with Table 1 (where both are negative-logFC modules) and undercuts the very narrative the table is meant to support. A reader checking Table 4 against Table 1 will find an unresolvable contradiction.

【Specific fix】 Regenerate Table 4 from the knockout rows of `virtual_perturbation_scores.csv`. The correct entries are:
- CTRL erythrocyte (knockout, D7): `-0.074`, `-2.57`, `33.0`, `p = 0.0060`, `q = 0.020`.
- Neuropeptides (knockout, D7): `-0.068`, `-2.45`, `26.9`, `p = 0.0080`, `q = 0.020`.
Add a one-line note that the table was regenerated from the deposited artefacts.

---

### Issue 2 — "Satellite gliosis" is over-read from a single marker and is conflated with the myeloid signal
【Problem】 The manuscript headlines "satellite gliosis (Gfap +1.33)" as a co-equal component of the day-7 programme alongside "myeloid infiltration" and "type-I interferon," but the glial/immune module is dominated by myeloid markers and Gfap is the only glia marker in it.

【Evidence】
- `module_level_signal.csv` D7, module "Glial / immune", genes = `Gfap;Aif1;Itgam;Cd68;S100a9`. Four of five are myeloid/microglial (Aif1, Itgam, Cd68, S100a9); only Gfap is glial. Module mean logFC = `+2.01`, driven by Cd68 (`+2.86`, `etiology_marker_de.csv`) and Aif1 (`+2.18`) — Gfap is the smallest of the three at `+1.33`.
- The abstract/§3.3/§4 present "satellite gliosis" as a definitive programme. The Limitations section (line 165) already concedes it "cannot be separated from Schwann-cell change," yet the headline claims are not correspondingly softened.
- In the independent sc validation, Gfap is *significantly* up in B16F10 (`sc_GSE292522_pseudobulk_DE.csv`: logFC `+0.0036`, t `+5.83`, q `0.0061`) but with a trivial magnitude (~1.003-fold). So even where a glial signal exists, it is biologically marginal.

【Why it matters】 In bulk DRG, a single marker at +1.33 log2 is consistent with — but not proof of — satellite-glia activation; it could reflect perivascular/endoneurial glia, Schwann-cell change, or low-abundance contamination. Asserting "satellite gliosis" as a named programme oversells the evidence and invites a reviewer (and a glial-biased reader) to reject the whole day-7 conclusion. It also blurs the more defensible claim: an innate-immune/myeloid event with a modest glial component.

【Specific fix】 Replace the phrasing. Suggested:
> "Day 7 was dominated by a type-I interferon / myeloid programme (Cd68 +2.86, Aif1 +2.18 log2; interferon-stimulated module mean t = +24.9), accompanied by a modest glial component (Gfap +1.33 log2) that is consistent with, but not by itself conclusive for, satellite-glia activation and cannot be separated from Schwann-cell change."

---

### Issue 3 — The interferon claim is supported by the module, but the cited "individual genes" are not type-I ISGs
【Problem】 §3.3 states "Individual genes agree: Cd68 +2.86, Aif1 +2.18, Gfap +1.33, Ccl2 +0.81, Il1b +0.80," implying these support the type-I interferon programme. None of them are classical type-I interferon-stimulated genes.

【Evidence】 The five genes listed are myeloid/glia/chemokine markers (Cd68, Aif1, Gfap = myeloid/glia; Ccl2, Il1b = chemokines/cytokines). The only genuinely type-I ISG evidence is the 8-gene "CTRL interferon" module (`Irf1;Ifitm3;Ifi47;Bst2;Rsad2;Gbp2;Ubd;Ly6e`), whose mean moderated t = `+24.9`, q = `6.8×10⁻⁵` (`module_level_signal.csv`, verified). The text never reports any individual ISG logFC.

【Why it matters】 The sentence reads as if five independent genes confirm a type-I IFN event; a careful reader will notice they confirm *myeloid/innate immune* infiltration, not interferon signalling per se. The interferon conclusion is still valid (the module is bona fide), but the write-up misattributes the support and weakens trust.

【Specific fix】 Split the claim:
> "The type-I interferon programme is supported by the interferon-stimulated module itself (mean moderated t = +24.9, q = 6.8×10⁻⁵; genes Irf1, Ifitm3, Ifi47, Bst2, Rsad2, Gbp2, Ubd, Ly6e). Separately, the myeloid/innate-immune compartment is strongly up (Cd68 +2.86, Aif1 +2.18 log2) and is accompanied by glial (Gfap +1.33) and chemokine (Ccl2 +0.81, Il1b +0.80) upregulation."

---

### Issue 4 — "Model-specific" cannot be distinguished from "species-specific"
【Problem】 The conclusion that the interferon programme is "model-specific … not a universal cancer-pain signature" conflates model identity with species identity. Walker-256 is *rat* tibia implantation; GSE292522 (B16F10) is *mouse* sciatic-nerve-injury. The two differ in species, injury paradigm, and quantification simultaneously.

【Evidence】
- `sc_GSE292522_summary.json`: dataset "GSE292522 (10x scRNA-seq mouse DRG, B16F10 melanoma model)" — i.e., mouse, not rat.
- The manuscript itself notes (§3.8, Limitations) that B16F10 is "a different species, injury paradigm and quantification," yet the title, abstract (line 21), and §3.8 repeatedly call the finding "model-specific." The data can support "not generalizable beyond the rat Walker-256 model"; they cannot support attributing the difference to *model* rather than *species*.

【Why it matters】 A reviewer will flag that a rat-vs-mouse comparison cannot isolate model from species. Over-claiming "model-specific" invites the objection that the absence of interferon in B16F10 is simply a mouse-vs-rat effect (or a nerve-injury-vs-tibia-tumour effect), which would shrink the substantive conclusion. The honest, defensible claim is narrower and should be stated as such.

【Specific fix】 Hedged wording:
> "The interferon/myeloid programme seen at Walker-256 day 7 is not reproduced in an independent mouse B16F10 cancer-nerve-injury model, which instead shows a neuronal-injury programme (Sprr1a, Atf3). Because the two models also differ in species and injury paradigm, we can conclude the interferon signature is not a universal cancer-pain programme, but we cannot, from this comparison alone, separate model-specific from species-specific effects."

---

### Issue 5 — The single-cell ISG statement is factually wrong ("q > 0.15")
【Problem】 §3.8 states "Every interferon-stimulated gene tested was flat (|logFC| < 0.11, q > 0.15)." The "|logFC| < 0.11" half holds, but the "q > 0.15" half is contradicted by at least one gene.

【Evidence】 `sc_GSE292522_focus_markers.csv`: `Ifit2` → logFC `-0.048`, t `-3.857`, **q = 0.042** (significant, and it is *down*, not up). `sc_GSE292522_summary.json` correctly lists `ISG_up_in_B16F10 = []` (no *up*-regulated ISG), so the spirit — "no interferon UP programme" — survives; but the literal "every ISG flat / q > 0.15" is false, and the gene that violates it is significantly (if trivially) regulated.

【Why it matters】 Precise quantitative claims are exactly what gives this paper its credibility; a one-line overstatement of a q-value threshold is the kind of error a methods reviewer will catch and use to question the rest. It does not change the conclusion, but it must be corrected.

【Specific fix】
> "No interferon-stimulated gene was significantly up-regulated in B16F10 (none met q < 0.05 with positive logFC); the few ISGs that reached significance were marginally down (e.g., Ifit2, logFC −0.05, q = 0.042), and all had |logFC| < 0.11."

---

### Issue 6 — The B16F10 validation is not purely "neuronal-injury only"; a glial signal is present
【Problem】 §3.8 frames B16F10 as showing "a neuronal-injury programme (Sprr1a, Atf3) rather than the interferon programme" and omits that Gfap is significantly up there too.

【Evidence】 `sc_GSE292522_pseudobulk_DE.csv`: Gfap logFC `+0.0036`, t `+5.83`, q `0.0061` (significant, though trivial magnitude — see Issue 2). Sprr1a `+1.55` (t `54.5`) and Atf3 `+0.84` (t `13.1`, q `6.3×10⁻⁵`) verified. (Galanin is also strongly up in B16F10, `+1.42`, consistent with an injury programme, but the authors do not highlight it — acceptable, just note it reinforces their injury narrative.)

【Why it matters】 Minor, but it tightens the contrast. The authors' "neuronal-injury vs interferon" dichotomy is slightly too clean; acknowledging the small glial component makes the comparison more credible and pre-empts the critique that they cherry-picked the B16F10 readout.

【Specific fix】 Add one sentence in §3.8 noting Gfap is significantly but marginally up in B16F10 (q = 0.006, |logFC| ≈ 0), confirming the glial component is trivial there compared with the myeloid/IFN component in Walker-256.

---

### Issue 7 — Missing must-cite literature (a reviewer will demand these)
【Problem】 Several citation gaps weaken the biological framing and the cell-type inference.

【Evidence / required additions】
1. **GSE292522 source publication.** The single-cell validation dataset is used as a primary result but its originating paper is never cited in the 14-reference list. This is a硬性 omission — the dataset cannot be presented without citing the study that produced it.
2. **DRG single-cell / single-nucleus atlases that *define* the marker sets.** The authors infer cell-type composition (neurons, satellite glia, Schwann, myeloid, T cells, endothelial, fibroblast, erythrocyte) from marker genes without citing the atlases that validated those markers in DRG. At minimum cite a mouse DRG classification reference (e.g., Usoskin et al., *Nat Neurosci* 2015) and a recent DRG snRNA-seq/nuclei atlas. Without this, the composition call (and the "myeloid fraction ~0.1%") rests on asserted markers.
3. **Bone-cancer-pain immune/glial DRG transcriptome literature beyond Zhai 2021 [10].** There is a substantial literature on innate-immune and glial activation (IL-1β, TNF, CCL2, p38 MAPK) in bone-cancer-pain DRG/sciatic that should be cited to position the day-7 finding.
4. **CIPN DRG transcriptomics / snRNA-seq.** Recent single-nucleus studies of paclitaxel/oxaliplatin DRG showing neuronal-injury and regeneration programmes (Atf3, Sprr1a, Gap43) should be cited to contextualise the §3.5/§3.8 CIPN findings.
5. **Galanin in pain.** The "injury-associated neuropeptide" label for galanin should cite the galanin-in-neuropathic-pain literature, not just assert it.
6. **CIPN clinical-pharmacology RCTs.** The statement that CIPN "responds partly to agents acting on neuronal excitability" should cite the actual trials — duloxetine (Smith et al., CALGB 170601) and the gabapentin/pregabalin evidence base — rather than only the ESMO/NCCN guidelines [1,2].
7. **Bone-metastasis-pain RCTs for anti-resorptive/anti-inflammatory strategy.** Cite the denosumab and bisphosphonate bone-pain trials (e.g., Wong/Conte-level evidence) alongside Mantyh [3,14].

【Why it matters】 Citation completeness is part of biological credibility. Items 1–2 are close to mandatory for a transcriptomics paper that infers cell types and validates on an external dataset.

【Specific fix】 Add the cited works to the References and a sentence in §1/§2.8 anchoring the marker-based composition call to the relevant atlas.

---

### Issue 8 — Clinical-pharmacology statements are accurate but should state evidence level
【Problem】 The clinical read (bone pain → anti-resorptive/anti-inflammatory; CIPN → neuronal-excitability agents) is broadly correct, but the equivalence implied between agents is not supported by equal evidence.

【Evidence】 Bone: anti-resorptives (bisphosphonates, denosumab) and NSAIDs/steroids are indeed first-line for bone-metastases pain (Mantyh 2013 [3], Luger 2005 [14]) — accurate. CIPN: duloxetine has the strongest RCT support; gabapentin/pregabalin evidence in CIPN is modest and mixed. The manuscript's phrasing ("responds partly to agents acting on neuronal excitability") is defensible but flattened.

【Why it matters】 A clinical-pharmacology claim in a mechanistic paper should not overstate analgesic efficacy of gabapentinoids in CIPN, where the evidence is weaker than for duloxetine.

【Specific fix】 Specify: "…CIPN responds, at best modestly, to gabapentin/pregabalin and more robustly to the SNRI duloxetine, agents that act on neuronal excitability."

---

### Issue 9 — The very low absolute myeloid fraction in B16F10 warrants a methods caveat
【Problem】 The estimated myeloid cell fraction is ~0.12% (control) and ~0.05% (B16F10). A ~0.1% myeloid fraction in DRG is implausibly low for a tissue that normally contains several-percent resident macrophages/microglia-like cells, suggesting the marker-based cell-type caller may be under-counting.

【Evidence】 `sc_GSE292522_summary.json`: myeloid fraction Con `0.00117`, B16F10 `0.000517`, MWU p `0.971`. The relative conclusion (no increase in B16F10) is unaffected, but the absolute number is suspicious.

【Why it matters】 If the myeloid caller is insensitive, the "no myeloid infiltration in B16F10" claim is weaker than presented (it could be a detection floor, not a true absence). It does not overturn the IFN conclusion but should be flagged.

【Specific fix】 Add a caveat in §3.8 / Limitations that the marker-based myeloid fraction is a lower bound and that the key comparison is the *relative* Con-vs-B16F10 change (non-significant), not the absolute percentage.

---

## What stands up (verified, and I had suspected otherwise)

1. **The day-7 type-I interferon / myeloid event is real and strong.** Module-level interferon mean moderated t = `+24.9`, q = `6.8×10⁻⁵` (recomputed from `module_level_signal.csv`). The 8-gene ISG set (Irf1, Ifitm3, Ifi47, Bst2, Rsad2, Gbp2, Ubd, Ly6e) is genuinely a type-I signature. Cd68 (+2.86) and Aif1 (+2.18) on bulk DRG at day 7 are large, biologically coherent innate-immune signals. I had wondered whether this was an artefact of the one contaminated muscle sample (GSM4506604), but the day-7 signal survives its removal (run_summary / de_muscle_drop_sensitivity), so it is robust.

2. **The cross-aetiology marker directions are sensible and exact-permutation supported.** Cd68 and Aif1 rise in bone-cancer pain and fall in paclitaxel neuropathy (Δ = +3.42 and +3.08, q = 0.012, recomputed from `etiology_marker_cross.csv`); galanin moves oppositely (Δ = −2.59, q = 0.012); Atf3 is higher in paclitaxel (+1.13) than at bone day 7 (+0.22). This is exactly the direction a pain neurobiologist would predict and it is verified.

3. **The honest-null perturbation result is the strongest part of the paper.** Network vs additive Spearman ρ = `0.991` (recomputed from `head_to_head_baselines.json`); of 132 degree-matched-null tests, 32 were "significant" of which 14 were negative controls including the erythrocyte set; minimum exact-permutation q = `0.165` (recomputed from `etiology_summary.json`); KO and OE are exact mirrors. The claim that in-silico perturbation here tracks composition/topology, not pharmacology, is well supported and rigorously demonstrated.

4. **The exact-permutation cross-aetiology design (4,900 enumerable null pairs, floor 1/4,901) is a genuine methodological strength**, not a claim I could fault on the numbers I checked (Table 2 and Table 3 values recomputed and matched: e.g., Glial/immune Δ = +2.70, q = 0.0073; Cytokine/chemokine Δ = +1.70, q = 0.0073; Neuropeptides Δ = −1.35, q = 0.0145).

---

## Questions for the authors (do not guess; answer explicitly)

1. **GSE292522 provenance.** Which publication produced GSE292522, and will you cite it? (Required — see Issue 7.1.) Is the B16F10 "cancer-nerve-injury" model a tumour inoculated onto/around the sciatic nerve, and what is the exact injury paradigm?
2. **Disentangling model from species.** Given Walker-256 is rat and B16F10 is mouse, do you have any same-species two-model or same-model two-species comparison that would let you claim "model-specific" rather than merely "not universal"? If not, will you rephrase to "not generalizable beyond the rat Walker-256 model"?
3. **Cell-type calling.** How exactly were the eight cell-type fractions estimated from the 10x data (which marker sets, which assignment rule), and do you have a positive control showing the myeloid caller detects known DRG macrophages? The ~0.1% myeloid fraction is far below expected resident macrophage levels (Issue 9).
4. **Individual ISGs at day 7.** Will you report the individual interferon-stimulated gene logFCs at Walker-256 day 7 (not just the module mean), so the interferon claim is transparent at gene level (Issue 3)?
5. **Table 4 regeneration.** Will you regenerate Table 4 from the knockout columns and confirm the erythrocyte/neuropeptide rows (Issue 1)? I verified the other six rows match the knockout data exactly; only these two are wrong.
6. **Galanin heterogeneity.** Galanin is labelled "injury-associated" yet is *down* in bone-cancer pain and *up* in paclitaxel — do you intend to make a point about injury-type-specific neuropeptide programmes, or is galanin being used purely as a directional marker?

---

## What I actually checked

**Manuscript:** `docs/manuscript/manuscript_final.md` (full read, all 281 lines including Tables 1–5).

**Source artefacts recomputed (numbers above are mine, not transcribed from the manuscript):**
- `module_level_signal.csv` — D7 interferon module (mean t `+24.923`, mean logFC `+2.614`, q `6.76×10⁻⁵`); erythrocyte D7 module (mean logFC `-1.147`).
- `etiology_marker_de.csv` — Bone_D7 Cd68 `+2.858`, Aif1 `+2.180`, Gfap `+1.328`; PTX Atf3 `+1.126`.
- `etiology_marker_cross.csv` — Cd68 D7-vs-PTX Δ `+3.418` (q `0.0125`); Aif1 `+3.076` (q `0.0125`); Gal `−2.592` (q `0.0125`); Scn10a `−0.702` (q `0.0260`); Calca `−0.770` (q `0.0272`).
- `etiology_module_logfc_cross.csv` — Glial/immune D7-vs-PTX Δ `+2.695` (q `0.0073`); Cytokine/chemokine `+1.702` (q `0.0073`); Neuropeptides `−1.345` (q `0.0145`).
- `virtual_perturbation_scores.csv` — D7 knockout rows for CTRL erythrocyte (`−0.0744`, z `−2.568`, q `0.0195`) and Neuropeptides (`−0.0675`, z `−2.449`, q `0.0195`); confirmed these match the *overexpression* rows the manuscript printed in Table 4.
- `head_to_head_baselines.json` — Spearman network-vs-additive `0.9911`; mean |network| `0.0568`; mean |additive| `0.0328`; mean % oracle `20.5`.
- `etiology_summary.json` — minimum exact-permutation q `0.165`; 32/132 degree-matched-null "significant" (18 target / 14 control); KO–OE exact-mirror design.
- `sc_GSE292522_summary.json` — 55,967 cells; myeloid fraction Con `0.00117` / B16F10 `0.000517`; MWU p `0.971`; concordance ρ(D7) `0.098`; ISG_up_in_B16F10 = `[]`.
- `sc_GSE292522_focus_markers.csv` and `sc_GSE292522_pseudobulk_DE.csv` — Sprr1a `+1.554` (t `54.49`), Atf3 `+0.835` (q `6.3×10⁻⁵`), Gfap `+0.0036` (q `0.0061`), Cd68/Aif1 unchanged; **Ifit2 q `0.042`** (contradicts manuscript "q > 0.15").
- `run_summary.json` — n_sig_q05_adj D7 = `4600` (matches §3.2); muscle ANOVA F `0.637`, p `0.605`.

**Not read (per independence discipline):** any other reviewer file, response/revision letters, `SUBMISSION_MANIFEST.md`, deposit SOP, author-verification statement, and everything else under `docs/review/round1_2026-10-10/`.

---

## Summary for the editor (domain layer)
Accept in principle after major revision. The core biology — divergence of bone-cancer-pain (innate-immune/myeloid) and CIPN (neuronal-injury) DRG programmes, plus the rigorous demonstration that in-silico perturbation here is a compositional artefact — is sound and, on the numbers I recomputed, reproducible. The blocking issues are: (i) the Table 4 transcription error (Issue 1), which must be regenerated from the knockout columns; (ii) over-claiming "satellite gliosis" from a single marker and conflating myeloid with interferon support (Issues 2–3); (iii) the unjustified "model-specific" label that confounds species (Issue 4); (iv) the incorrect "q > 0.15" ISG statement (Issue 5); and (v) mandatory missing citations, chiefly the GSE292522 source paper and the DRG cell-type atlases (Issue 7). None require new experiments; all require careful rewriting and one table regeneration.
