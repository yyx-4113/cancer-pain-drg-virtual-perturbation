# A1 — Cancer-Pain Neuroimmunology (domain) review

**Manuscript:** `docs/manuscript/manuscript_final.md`
**Base dir:** `d:/2026.10/不同部位癌痛患者对药物敏感性生信分析/cancer_pain_virtual_perturbation`
**Reviewer role:** expert A1, domain (DRG cell biology, pain mechanisms, marker interpretation)
**Posture:** first submission. No credit for effort. Every claim checked against source data or literature.

---

## Verdict up front

The paper's *methodological* thesis (perturbation layer fails its own controls) is not my remit and
looks well built. The *biological* thesis — "bone cancer pain = myeloid/interferon-driven, CIPN =
neuronal-injury-driven" — is **not defensible as written**. Three things are wrong, and one of them
is fatal to §3.8:

1. **The §3.8 compartment analysis cannot test the claim it is used for.** The myeloid compartment
   was dropped by the script's own ≥6-animal / ≥5-cell rule (46 myeloid cells total, only 5/8
   animals with ≥5). So "no compartment shows ISG up-regulation" is a statement about seven
   compartments, none of which is the compartment where an interferon programme would live. The
   paper reports this as "the per-compartment test turns a null average into a compartment-resolved
   null" (manuscript line 153). It is not a compartment-resolved null; it is a compartment-resolved
   null **with the index compartment removed**.
2. **The compartment labelled "satellite glia" contains no GFAP** (Gfap = 0.012 mean log₁p inside it;
   its top-expressed gene is the myelin gene Plp1 at 1.944; real GFAP+ cells number 0.15–0.40% of
   cells against 361 cells assigned to that compartment). The satellite-glia interferon findings in
   §3.8 are therefore measured in a compartment that is not satellite glia.
3. **The validation dataset's own published abstract says the opposite of what the paper concludes
   from it.** GSE292522 is ref [15] Baruch et al. *Nature* 2025, whose GEO summary states: *"Injured
   neurons respond by autonomously initiating an interleukin (IL)-6 and interferon (IFN) type I
   inflammatory response."* The manuscript concludes that IFN is absent from every compartment of
   that dataset. That is a direct contradiction with the source publication, and the paper does not
   engage with it.
4. **The erythrocyte control does not do the work the paper assigns to it, and it points the other
   way.** In the contrast that carries the headline (D7 vs PTX), the erythrocyte module shows
   Δ = −0.258, q = 0.814 — a clean null. The −2.48 the paper quotes is D7-vs-OXA. So the
   "compositional dilution" argument, which the Discussion uses to demote the neuronal finding, is
   built on a contrast the paper does not otherwise use.

Details, evidence and paste-ready fixes below.

---

## 1. Marker-to-cell-type assignments

### 1.1 【Problem】 `Cd68`/`Aif1` = "myeloid" is correct for the *marker*, but the paper then uses
### `Cd68`+`Aif1` to claim "myeloid **infiltration**", which the data cannot support.

**【Evidence】** `etiology_marker_de.csv`: Cd68 Bone_D7 = **+2.858** (t = 18.20, q = 0.0);
Aif1 Bone_D7 = **+2.180** (t = 14.16, q = 0.0). Recomputed from
`DE_D7_vs_Sham_adj.csv` with the pipeline's own design (n = 4v4, muscle covariate): Cd68 +2.86
(t = 18.04), Aif1 +2.18 (t = 14.04) — matches.

But the same file shows the infiltration is not a clean monocyte/macrophage influx:

| gene | Bone_D7 logFC | t | q | reading |
|---|---|---|---|---|
| Cd68 | +2.86 | 18.20 | 0.0 | macrophage/microglia |
| Aif1 | +2.18 | 14.16 | 0.0 | macrophage/microglia |
| Lyz2 | +3.38 | 13.80 | 5.6e-05 | **neutrophil** marker |
| S100a9 | +2.21 | 5.98 | 0.0031 | **neutrophil** marker |
| Itgam | +1.52 | 11.94 | 1.1e-04 | pan-myeloid |
| Serping1 | +2.46 | — | 8.5e-06 | mast cell / neutrophil |
| C3 | +3.43 | 26.44 | 4.8e-06 | complement |
| C1s | +2.61 | 24.72 | 5.0e-06 | complement |
| C2 | +3.32 | — | 1.2e-05 | complement |
| Cfb | +3.15 | 16.38 | 2.4e-05 | complement |
| Cp | +1.97 | 27.37 | 4.8e-06 | acute-phase |

Lyz2+S100a9 together with Cd68+Aif1 is the signature of a **neutrophil-rich acute inflammatory
infiltrate** (or, if the time point is day 7 post-injection, of residual acute-phase response), not
of a macrophage infiltrate. Complement activation (C3/Cfb/C1s/C2 all up 2.6–3.4 log₂) is a
dominant, unremarked feature of this contrast — arguably more prominent than the ISGs themselves.

**【Why it matters】** The manuscript's headline is "myeloid/interferon-driven". "Myeloid" is
consistent with the data. But the paper never distinguishes neutrophil from macrophage
infiltration, and Limitations concedes it "cannot be assigned to macrophages, monocytes or
neutrophils" — while simultaneously making "myeloid infiltration" a headline conclusion. If the
infiltrate is neutrophil-dominant, the reading connects to a different literature (neutrophil–DRG
interactions, S100A8/9–TRPV4) than the macrophage/microglia literature the Discussion invokes.

**【Specific fix】** Replace manuscript line 115's
`"...alongside myeloid infiltration (Cd68 +2.86, Aif1 +2.18 log₂) and satellite gliosis (Gfap +1.33)."`
with:

> "At day 7 the DRG shows an acute inflammatory infiltrate (Cd68 +2.86, Aif1 +2.18, Lyz2 +3.38,
> S100a9 +2.21 log₂) together with complement activation (C3 +3.43, Cfb +3.15, C1s +2.61) and a
> type-I interferon programme. Because both pan-myeloid and neutrophil-restricted markers rise
> together, the bulk data identify an acute inflammatory infiltrate but cannot assign it to
> macrophages, monocytes or neutrophils; histological or flow-cytometric assignment would be required
> to make that distinction."

---

### 1.2 【Problem】 `Gfap` is reported as **three different values** for the same contrast, and the
### value chosen in §3.3 is neither of the two the pipeline actually produces.

**【Evidence】** Recomputed / read directly:

| source | Gfap, Bone_D7 |
|---|---|
| `DE_D7_vs_Sham_adj.csv` (the limma-trend table) | **+1.2942** |
| `etiology_marker_de.csv` (cross-aetiology input) | **+1.328** |
| `etiology_marker_cross.csv` (per-sample simple mean) | **+1.4984** |
| manuscript §3.3, line 115 | "+1.33" |
| manuscript Table 2, line 247 | "+1.50" |

**【Why it matters】** §3.3 quotes +1.33 (matches `etiology_marker_de.csv`, not the DE table);
Table 2 quotes +1.50 (matches the unadjusted simple mean). The two differ because the
cross-aetiology estimator is an unadjusted per-sample mean while §3.3's number comes from the
muscle-adjusted moderated fit. Table 2 is presented as the cross-aetiology marker table but
silently mixes two estimators within one table (e.g. Cd68 +2.84 in Table 2 vs +2.86 in §3.3).
A reader cannot tell which estimator produced which cell.

**【Specific fix】** Add a footnote to Table 2:

> "Bone/OXA/PTX columns are per-sample simple means of log₂ (as used for the permutation Δ); values
> quoted in §3.3 are muscle-adjusted limma-trend moderated estimates from
> `DE_*_vs_Sham_adj.csv`. The two differ by up to 0.2 log₂ (e.g. Gfap, 1.29 vs 1.50)."

And change §3.3's "Gfap +1.33" to "Gfap +1.29" so that it matches the DE table cited in §3.2.

---

### 1.3 【Problem】 `Gfap` is called "satellite gliosis", but GFAP is explicitly a contested
### satellite-glia marker in DRG and the paper does not acknowledge the controversy.

**【Evidence】** Manuscript line 115 "satellite gliosis (Gfap +1.33)"; Table 1 groups Gfap with
myeloid genes into a module called "Glial / immune" whose members per `scripts/_limma.py:171` are
`["Gfap","Aif1","Itgam","Cd68","S100a8","S100a9"]` — i.e. **five of six are myeloid**, and
`module_level_signal.csv` shows `genes_used = Gfap;Aif1;Itgam;Cd68;S100a9` (S100a8 absent from the
rat matrix). The module is mislabelled: it is a myeloid module with one glial gene in it, and its
day-7 signal (mean t = +11.41) is driven by the myeloid genes.

Literature: Konnova et al. *IJMS* 2023 (PMID 37958541) show that after spared nerve injury **75.5%
of DRG GFAP+ cells lack FABP7 and are non-myelinating Schwann cells, not satellite glia**, and
Mohr et al. *Biomedicines* 2021;9:1022 is titled "Discrepancy in the Usage of GFAP as a Marker of
Satellite Glial Cell Reactivity". The specificity of GFAP for DRG satellite glia is an open
literature question, not a settled fact.

**【Why it matters】** "Satellite gliosis" is presented in the Abstract and §3.3 as one of three
named pillars of the bone-cancer-pain phenotype. If GFAP is in fact marking non-myelinating
Schwann cells, or is reflecting neuronal-glial reactive change, the "gliosis" claim loses its
cell-type attribution. The paper's own Limitations says "satellite-gliosis cannot be separated from
Schwann-cell change" — which contradicts the Abstract, not just qualifying it.

**【Specific fix】** Delete "Glial / immune" as a label and rename the module in
`scripts/_limma.py:171` and Table 1 to **"Myeloid + GFAP"**, with `genes_used` printed in the
table. Replace the Abstract phrase "satellite gliosis (Gfap +1.33)" with:

> "a reactive-glial signal (Gfap +1.29 log₂) that cannot be assigned to satellite glia versus
> non-myelinating Schwann cells in bulk tissue (Konnova et al. 2023)."

And add to Limitations:

> "GFAP is widely used as a satellite-glia activation marker in DRG, but its specificity in DRG is
> disputed: Konnova et al. (2023) found that 75.5% of GFAP+ DRG cells after spared nerve injury are
> FABP7-negative non-myelinating Schwann cells rather than satellite glia. Our Gfap result should
> therefore be read as reactive glia of unspecified subtype."

---

### 1.4 【Problem】 `Scn10a` = Nav1.8 is correct, but the paper's own §3.5 uses Scn10a as an
### injury-direction marker when a *decrease* in Scn10a is the dilution signature.

**【Evidence】** `etiology_marker_de.csv`: Scn10a Bone_D7 = **−0.399** (q = 0.0113);
PTX = +0.239 (q = 0.284). Cross Δ = −0.702, q = 0.026. Manuscript line 131 lists Scn10a among
markers that "point the same way" as galanin. Table 1 shows Scn10a inside "Voltage-gated Na+"
which falls −4.35 mean t, right alongside the pan-neuronal control (−5.20).

**【Why it matters】** Scn10a is a legitimate Nav1.8 assignment. But the paper uses Scn10a both as
(a) evidence of a neuronal-compartment change and (b) an analyte inside a module it simultaneously
argues is a dilution artefact. It cannot be both. If Scn10a's fall is dilution (as the paper argues
for the module), then Scn10a is not independent evidence of anything; if it is real, the module
argument fails. The paper does not resolve this.

**【Specific fix】** In §3.5, separate the two uses:

> "Scn10a (Nav1.8) falls in bone day 7 (−0.40 log₂, q = 0.011) and rises in paclitaxel neuropathy
> (+0.24), but because Scn10a is a pan-neuronal transcript it is subject to the same compositional
> dilution as the pan-neuronal control (−0.44 log₂, Δ = −0.69, q = 0.015). We therefore treat Scn10a
> as a member of the dilution set rather than as an independent marker, and we do not read it as a
> change in Nav1.8 expression."

---

### 1.5 【Problem】 The `CTRL interferon` module is 8 genes, none of which is `Isg15`, `Mx1`, `Cxcl10`
### or `Ifit3`; calling it "type-I interferon programme" is defensible but the module is not
### verifiable as such, and it is a *negative control*, so its day-7 t of +24.92 is a contradiction
### the paper never resolves.

**【Evidence】** `scripts/_limma.py:195`: `"CTRL interferon": ["Irf1","Ifitm3","Ifi47","Bst2","Rsad2","Gbp2","Ubd","Ly6e"]`.
All 8 rise at D7 (Irf1 +2.31, Ifitm3 +1.97, Ifi47 +3.32, Bst2 +2.34, Rsad2 +2.40, Gbp2 +3.25,
Ubd +4.03, Ly6e +1.29; all q < 6e-06). So the module is real, coherent, and unmistakably
type-I-IFN. I verified a wider panel independently from `DE_D7_vs_Sham_adj.csv`: Isg15 +2.35,
Mx1 +2.25, Ifit3 +2.24, Stat1 +1.74, Oas1a +1.92, Cxcl10 +3.47, Ccl5 +4.56, Irf7 +2.67,
Usp18 +1.99, Gbp5 +2.51, Psmb8 +2.99, Psmb9 +3.08, Tap1 +2.48 — **the IFN programme is genuine and
broader than the 8-gene module.**

But it is defined in §2.5 as one of **five negative-control modules** — "so that a compartment or
composition effect becomes visible" — and it is then used in §3.3 as the study's principal positive
finding and in §3.4 as the strongest reversal (Table 4, +0.228, 75.4% of oracle). §3.3 opens "The
interferon-stimulated **control** module was the strongest signal in the study."

**【Why it matters】** A negative control that fires at +24.9 mean moderated t is not a control any
more; it is the result. The paper runs both readings simultaneously without acknowledging that these
are incompatible rhetorical moves. Either the interferon module is a control (and its D7 signal
means the control is invalid) or it is a positive finding (and §3.6's T2 statistic counting "14 of
the 32 randset-significant tests were negative-control modules, including the erythrocyte set" is
misleading, because the single strongest control *is* the headline biology).

**【Specific fix】** Split the two uses explicitly. In §2.5, state:

> "The interferon-stimulated gene module is included as a negative control in the perturbation
> calibration layer (§3.4, §3.6) because it is not an analgesic target. It is *not* a negative
> control for the differential-expression analysis, and §3.3 reports it as a positive biological
> result. These two uses are kept separate: the module's perturbation reversal is subject to the same
> non-specificity critique as every other module (Table 4), whereas its day-7 differential
> expression is interpreted as biology."

And in §3.6, when reporting that controls outrank targets, add: "Note that the control module with
the largest day-7 differential signal is the interferon module (mean t = +24.9), which is
nonetheless the study's principal positive finding; a control that dominates the positive result is
not performing its calibration function."

---

## 2. The core biological claim: "bone = myeloid/IFN, CIPN = neuronal injury"

### 2.1 【Problem】 The IFN attribution to Walker-256 day 7 is **robust** — but it is not
### Walker-256-*specific*, and the paper's own day-14 data refute "day-7-peak framing" as a
### tumour-size effect.

**【Evidence】** Leave-one-out at D7, muscle-adjusted, using the pipeline's design (recomputed):

| dropped D7 sample | ISG mean t | Cd68 | Aif1 | Gfap | myeloid mean t | pan-neur logFC | n q<0.05 |
|---|---|---|---|---|---|---|---|
| none | +24.92 | +2.86 | +2.18 | +1.29 | +16.53 | −0.40 | 4600 |
| GSM4506602 | +21.72 | +2.84 | +2.29 | +1.32 | +14.72 | −0.37 | 3364 |
| GSM4506603 | +21.81 | +2.70 | +2.01 | +1.25 | +16.39 | −0.39 | 3368 |
| **GSM4506604** (6.26% sarcomere) | **+21.89** | +2.96 | +2.35 | +1.19 | +14.98 | −0.33 | **3041** |
| GSM4506605 | +22.62 | +3.00 | +2.19 | +1.34 | +16.78 | −0.44 | 5300 |

The interferon programme **survives dropping the contaminated sample** (mean t 24.9 → 21.9, Cd68
*increases* to +2.96). So it is not a muscle-contamination artefact. Good — this stands up.

But `module_level_signal.csv` shows the IFN module is essentially equally present at D14:

| contrast | CTRL interferon mean t | mean logFC | Glial/immune mean t | mean logFC |
|---|---|---|---|---|
| D1 | +2.86 | +0.30 | +1.39 | +0.21 |
| **D7** | **+24.92** | **+2.61** | +11.41 | +2.01 |
| **D14** | **+4.36** | **+2.37** | +3.92 | +1.31 |

D14 logFC for the IFN module is **+2.37 vs +2.61 at D7** — essentially the same effect size, with a
4-fold lower t. And myeloid: Cd68 D14 = +2.45, D7 = +2.86 (from `etiology_marker_de.csv`).

**【Why it matters】** The paper describes "a sharp peak at day 7" and "the interferon programme at
its day-7 peak", implying a transient. The logFCs say the programme is **sustained from day 7 to
day 14**; only the *significance* peaks, because D14 variance is higher. The 4,600-vs-235 gene
count difference is largely a variance difference, not a biology difference. Framing this as a
"stage-dependent" or "peak" event (Abstract; §3.2 heading "a sharp peak at day 7") overstates the
transience. It also means the day-7 choice is not obviously the right timepoint to characterise
"bone cancer pain" — D14 would give the same answer with different power.

**【Specific fix】** Rewrite §3.2's heading and the Abstract:

> §3.2 heading: "Bone cancer pain reorganises the DRG transcriptome from day 7 onward, sustained
> through day 14"
>
> "The interferon programme is present at day 7 (mean log₂ FC +2.61) and persists essentially
> unchanged at day 14 (+2.37), although only day 7 reaches genome-wide significance (mean moderated
> t = +24.9 vs +4.4; 4,600 vs 235 genes at q < 0.05). The day-7 > day-14 difference in gene count
> therefore reflects variance rather than programme magnitude."

---

### 2.2 【Problem】 The sham control is **day-14 PBS**, so the "D1 = 0 genes" result is a
### time-matched D1-vs-D14 comparison being reported as if D1 were an unperturbed baseline.

**【Evidence】** `data/raw/GSE149648_series_matrix.txt.gz`:
`!Sample_characteristics_ch1` = `"time point: Day 14 after Sham operation"` for all four control
samples GSM4506610–613. Manuscript §2.1 describes the groups as "sham (PBS) and days 1, 7 and 14
after implantation" without stating the sham is day-14-matched.

Recomputed, n=4v4 with muscle covariate:

| contrast | n q<0.05 | ISG mean t | ISG mean logFC | Cd68 |
|---|---|---|---|---|
| D1 vs D14-sham | **0** | +2.86 | +0.30 | +0.26 |
| D7 vs D14-sham | 4600 | +24.92 | +2.61 | +2.86 |
| D14 vs D14-sham | 235 | +4.36 | +2.37 | +2.45 |

**【Why it matters】** D1 (tumour, day 1) is compared against a **day-14 wound-healing** control.
A day-14 PBS-injected tibia is still in the resolution phase. So D1-vs-sham is a comparison across
two different wound ages, and the fact that it yields literally zero genes at q<0.05 (with median
|logFC| = 0.073) is surprising rather than reassuring. The manuscript reports this as "0 genes on
day 1" without noting the control is not day-1 matched. If a day-1-matched control exists in the
series it should be used; if not, the day-1 row cannot support a temporal-onset claim.

**【Specific fix】** In §2.1:

> "The sham group comprises four animals sampled at day 14 after PBS injection (GEO
> `!Sample_characteristics_ch1`: 'time point: Day 14 after Sham operation'); no day-1 or day-7
> control exists. The D1 contrast is therefore a day-1 tumour versus day-14 wound comparison and
> cannot support a claim about programme onset."

And in §3.2, delete "0, 4,600 and 235 genes at q < 0.05 on days 1, 7 and 14" as a
stage-dependence statement, replacing with:

> "0 (D1), 4,600 (D7) and 235 (D14) genes at q < 0.05, noting that D1 is referenced to a day-14
> control and D14 is the only time-matched contrast."

---

### 2.3 【Problem】 The erythrocyte-control argument in §3.6 is computed on the wrong contrast and,
### read correctly, it **undermines** the paper's own compositional-dilution defence of the
### neuronal finding rather than supporting it.

**【Evidence】** `etiology_module_logfc_cross.csv`, CTRL erythrocyte:

| contrast | ΔlogFC | p (exact) | q |
|---|---|---|---|
| **Bone_D7 vs PTX** | **−0.258** | **0.716** | **0.814** |
| Bone_D7 vs OXA | **−2.483** | 0.00061 | 0.0073 |
| Bone_D14 vs PTX | +1.528 | 0.0051 | 0.0321 |
| PTX vs OXA | −2.225 | 0.0423 | 0.0996 |

The −2.48 the manuscript quotes (§3.6 line 141, Table 3) is **D7-vs-OXA**. The contrast that
carries the entire headline — Cd68 Δ = +3.42, Aif1 Δ = +3.08, Gal Δ = −2.59, all "D7 vs PTX"
(Table 2) — has an erythrocyte Δ of **−0.258, q = 0.814**: a clean null.

Meanwhile the neuronal module deltas the paper says "show Δ in the same direction as the erythrocyte
Δ":

| module | Bone_D7 vs PTX Δ | q |
|---|---|---|
| CTRL pan-neuronal | −0.691 | 0.0145 |
| Neuropeptides | −1.35 | 0.0145 |
| Neurotrophin | −0.97 | 0.0145 |
| CTRL erythrocyte | **−0.258** | **0.814** |

**【Why it matters】** In D7-vs-PTX the erythrocyte signal is null, so there is **no measured
blood-content difference** between the two preparations carrying the headline comparison — and yet
the neuronal modules differ significantly. Read literally, the paper's own control says the
neuronal differences are **not** explained by erythrocyte dilution. That is the opposite of what
§3.6 concludes. The paper's inference — "neuronal modules show Δ in the same direction as the
erythrocyte Δ, so their apparent cross-aetiology difference cannot be separated from this
compositional difference" — is not supported in the contrast it is applied to. It is only supported
in D7-vs-OXA, where the paper does not make the neuronal claim.

This does not rescue the neuronal modules either: a *neuron-fraction* change (myeloid diluting
neurons) need not involve erythrocytes. But the specific sentence as written is wrong.

**【Specific fix】** Rewrite §3.6's second paragraph:

> "The erythrocyte control differs between the two preparations only in the D7-versus-OXA contrast
> (ΔlogFC = −2.48, q = 0.0073), where the two aetiologies are compared against each other. In the
> contrast that carries the marker-level claims — bone day 7 versus paclitaxel — the erythrocyte
> control shows no difference (ΔlogFC = −0.26, p = 0.72, q = 0.81), so blood content does not
> confound that comparison. Neuronal modules do nevertheless fall in that contrast (pan-neuronal
> Δ = −0.69, q = 0.015), and since a change in neuronal fraction can occur without any change in
> erythrocyte content, this decline cannot be attributed to blood content specifically. It remains
> consistent with myeloid-fraction dilution, which we cannot test in bulk tissue. We therefore
> retain the compositional caution but rest it on the myeloid markers themselves rather than on the
> erythrocyte control."

---

### 2.4 【Problem】 The paper lumps oxaliplatin with paclitaxel as "the neuronal-injury-driven
### aetiology", but its own data show oxaliplatin has **no** neuronal-injury programme at all.

**【Evidence】** `etiology_marker_de.csv`, OXA column:

| gene | OXA logFC | q | reading |
|---|---|---|---|
| Gal | −0.295 | 0.554 | no injury up-regulation |
| Atf3 | +0.158 | 0.566 | no injury up-regulation |
| Calca | +0.026 | 0.937 | flat |
| Scn10a | −0.108 | 0.637 | flat |
| Scn9a | −0.310 | 0.208 | flat/down |
| Npy | +1.070 | 0.276 | — |
| **Cd68** | **+0.442** | **0.068** | **myeloid marker UP** |
| Aif1 | −0.291 | 0.387 | flat |

And from `DE_D7_vs_Sham_adj.csv` I recomputed `etiology_marker_de.csv`'s OXA Cd68 raw statistics:
t = +5.22, P = 0.00055 (q = 0.068).

The manuscript's Discussion (line 161) states: "Paclitaxel neuropathy is a neuronal injury event,
with galanin, CGRP and Atf3 up and **no myeloid infiltration at all — Aif1 and Cd68 are if anything
slightly down**." Table 2 of the same manuscript lists Cd68 OXA = **+0.40**.

**【Why it matters】** Two problems. (i) Oxaliplatin shows none of the injury markers the paper
calls the CIPN signature (Gal, Atf3 all null), so the "CIPN = neuronal-injury-driven" half of the
binary is carried entirely by paclitaxel. The paper's own GSE160543 title says the two drugs have
"distinct etiologies". Presenting "CIPN" as one end of a binary misrepresents the design. (ii) The
Discussion's "Cd68 ... slightly down" directly contradicts the paper's own Table 2 (+0.40). Cd68 in
OXA has raw P = 0.00055 — the most significant myeloid signal in the OXA arm. An OXA arm with a
rising myeloid marker and no injury programme is a third phenotype, not a member of the
"neuronal-injury" class.

**【Specific fix】** Rewrite the Discussion's second sentence and add an explicit three-way result:

> "Paclitaxel neuropathy is a neuronal injury event, with galanin (+1.92), CGRP/Calca (+0.50) and
> Atf3 (+1.13) up and no myeloid infiltration (Cd68 −0.58, Aif1 −0.87). Oxaliplatin neuropathy is
> neither: it shows no galanin (−0.30), Atf3 (+0.16), Calca (+0.03) or Scn10a (−0.11) response, while
> Cd68 is in fact modestly *up* (+0.44, P = 5.5 × 10⁻⁴). We therefore report three phenotypes rather
> than two: Walker-256 day 7 (myeloid/interferon), paclitaxel (neuronal injury), oxaliplatin (no
> detectable programme at these markers). The binary 'bone versus CIPN' framing in the Abstract
> should be read as 'bone versus paclitaxel'."

---

## 3. The compartment analysis (§3.8 second paragraph) — the fatal problem

### 3.1 【Problem】 The myeloid compartment was **dropped by the pipeline's own filter**, so "no
### compartment shows interferon up-regulation" excludes the one compartment in which an
### interferon programme would reside. The paper reports the compartment list as including
### compartments that were never tested.

**【Evidence】** I re-ran the compartment assignment from `scripts/32_percelltype_sc.py` logic
against all eight 10x matrices. Cells per animal:

| compartment | Con1 | Con2 | Con3 | Con4 | B16-1 | B16-2 | B16-3 | B16-4 | total | ≥5 cells in ≥6 animals? |
|---|---|---|---|---|---|---|---|---|---|---|
| Neuron | 5159 | 5840 | 4606 | 4589 | 4226 | 5366 | 5489 | 4890 | 40165 | KEPT |
| SatelliteGlia | 46 | 35 | 37 | 25 | 36 | 44 | 108 | 30 | 361 | KEPT |
| Schwann | 1232 | 1245 | 1215 | 1122 | 1649 | 1842 | 3957 | 1739 | 14001 | KEPT |
| **Myeloid** | **3** | **10** | **12** | **5** | **5** | **3** | **7** | **1** | **46** | **DROPPED (5/8)** |
| **TCell** | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | **1** | **DROPPED (0/8)** |
| Endothelial | 40 | 50 | 86 | 53 | 43 | 48 | 79 | 50 | 449 | KEPT |
| Fibroblast | 29 | 138 | 41 | 20 | 209 | 172 | 159 | 89 | 857 | KEPT |
| Erythrocyte | 5 | 4 | 22 | 3 | 7 | 5 | 26 | 15 | 87 | KEPT |
| Neuron.NF | 3405 | 4380 | 2993 | 3306 | 3051 | 3157 | 3395 | 2942 | 26629 | KEPT |
| **Neuron.NP** | 1 | 0 | 0 | 0 | 0 | 1 | 1 | 0 | **3** | **DROPPED** |
| Neuron.PE | 1753 | 1460 | 1613 | 1283 | 1173 | 2206 | 2090 | 1947 | 13525 | KEPT |
| **Neuron.TH** | 0 | 0 | 0 | 0 | 2 | 2 | 3 | 1 | **8** | **DROPPED** |

`scripts/32_percelltype_sc.py:150` → `if len(present) < 6: continue`. Myeloid has 5 animals with
≥5 cells → **dropped**. TCell, Neuron.NP and Neuron.TH are also dropped.

But `sc_GSE292522_celltype_DE.csv` contains rows for exactly: Neuron, SatelliteGlia, Schwann,
Endothelial, Fibroblast, Erythrocyte, Neuron.NF, Neuron.PE — 8 compartments, no myeloid.

Meanwhile manuscript line 153 states: *"No compartment — neurons, neurofilament-positive or
peptidergic neuronal subtypes, satellite glia, Schwann cells, endothelium, fibroblasts or
erythrocytes — showed interferon-stimulated-gene up-regulation"* and §5: *"confirms the absence of
an interferon programme in every DRG compartment."*

**【Why it matters】** This is the single most serious problem in the manuscript. Myeloid cells are
the cell type that carries type-I interferon signalling in peripheral nerve. An interferon
programme in DRG is, mechanistically, an *infiltrating-immune-cell* phenotype. By dropping 46
myeloid cells, the analysis is structurally incapable of finding an IFN programme, and then reports
the absence as evidence that the IFN programme is model-specific. The negative result is an artefact
of the compartment filter, not a biological observation.

The manuscript's §3.8 sentence is also internally misleading: it lists eight compartments and
implies exhaustive coverage of "every DRG compartment", while three defined compartments (Myeloid,
TCell, plus two neuronal subtypes) were never examined. The word "every" is false.

**【Specific fix】** Two changes required, both in `scripts/32_percelltype_sc.py` and manuscript §3.8.

(a) Lower the compartment threshold so myeloid is retained, and report it as a **descriptive
pooled** analysis rather than a per-animal test (because 1–12 cells/animal cannot support a
per-animal pseudobulk):

```python
# scripts/32_percelltype_sc.py, line 148-151
for ctype in COMPARTMENTS:
    present = [p for p in SAMPLES if p in agg[ctype]]
    if ctype in ("Myeloid", "TCell"):
        # Too few cells per animal for per-animal pseudobulk. Pool across animals
        # within group and report descriptively; do NOT run the limma test.
        pooled = {grp: np.sum([agg[ctype][p] for p in present if SAMPLES[p][0] == grp], axis=0)
                  for grp in ("Con", "B16F10")}
        ...  # CPM-normalise, log2, report per-group ISG mean and marker positivity
        continue
    if len(present) < 6:
        continue
```

(b) Rewrite manuscript §3.8's second paragraph:

> "Seven compartments passed the pre-specified threshold of ≥5 cells in ≥6 of 8 animals (Neuron,
> Neuron.NF, Neuron.PE, SatelliteGlia, Schwann, Endothelial, Fibroblast, Erythrocyte). The myeloid
> and T-cell compartments contained 46 and 1 cells respectively across all 55,967 cells and were
> excluded by that threshold. Because myeloid cells are the principal carriers of type-I interferon
> signalling in peripheral nerve, this analysis cannot exclude a myeloid-restricted interferon
> programme, and our conclusion is restricted to the seven compartments tested. We therefore do not
> claim that the interferon programme is absent from the B16F10 DRG; we claim only that it is absent
> from neurons, satellite glia, Schwann cells, endothelium, fibroblasts and erythrocytes in this
> dataset. Resolving the myeloid compartment requires either a myeloid-enrichment experiment or
> targeted sequencing of sorted CD11b+ cells."

And change §5's "confirms the absence of an interferon programme in every DRG compartment" to
"confirms its absence in the seven compartments that passed the cell-count threshold".

---

### 3.2 【Problem】 The compartment DE panel contains **none of the eight genes** that define the
### interferon programme being tested for, and **none of the injury genes** whose neuron-localisation
### is claimed.

**【Evidence】** `scripts/32_percelltype_sc.py:80-82` defines
`FOCUS = [Cd68, Aif1, Itgam, Lyz2, Cx3cr1, Gfap, Sox10, Mpz, Snap25, Scn9a, Scn10a, Calca, Tac1, Gal,
Il1b, Tnf, Ccl2, Ifit1, Ifit2, Ifit3, Isg15, Oas1a, Oas1b, Bst2, Mx1, Stat1, Irf7]`.

Intersection with the bulk interferon module (`Irf1, Ifitm3, Ifi47, Bst2, Rsad2, Gbp2, Ubd, Ly6e`):
**only `Bst2`.** Seven of eight genes that define the day-7 Walker-256 interferon signal — the
genes with the largest fold changes (Ubd +4.03, Ifi47 +3.32, Gbp2 +3.25, Ifitm3 +1.97, Rsad2 +2.40,
Irf1 +2.31, Ly6e +1.29) — were **never tested in any compartment**.

Intersection with the injury panel (`INJURY = [Sprr1a, Atf3, Gadd45a, Ecel1, Serpinb1a, Jun, Fos,
Socs3]`, line 77): **empty.**

So the manuscript's sentence (line 153) *"the injury programme was neuron-intrinsic: galanin was up
+1.47 log₂ in neurons"* rests on **Galanin and Scn10a only**. The genes the paper names as the
injury programme in the very same paragraph — `Sprr1a`, `Atf3` — were not measured per compartment.
`Atf3` is in `FOCUS`? No: it is not. Confirmed absent from `sc_GSE292522_celltype_DE.csv`.

**【Why it matters】** "The per-compartment test turns a null average into a compartment-resolved
null" (line 153) is not a test of the interferon programme that was observed in the rat data. It is
a test of a *different, 10-gene* ISG panel that overlaps the original by one gene. Similarly,
"neuronal-injury-driven" is not localised; it is inferred from a single neuropeptide.

**【Specific fix】** Extend `FOCUS` to include the actual discovery modules, and re-run:

```python
# scripts/32_percelltype_sc.py, replace FOCUS (lines 80-82)
BULK_ISG8 = ["Irf1","Ifitm3","Ifi47","Bst2","Rsad2","Gbp2","Ubd","Ly6e"]
INJURY8   = ["Sprr1a","Atf3","Gadd45a","Ecel1","Serpinb1a","Sox11","Gap43","Nefl"]
FOCUS = ["Cd68","Aif1","Itgam","Lyz2","Cx3cr1","Gfap","Sox10","Mpz","Snap25","Scn9a","Scn10a",
         "Calca","Tac1","Gal","Il1b","Tnf","Ccl2","Ifit1","Ifit2","Ifit3","Isg15","Oas1a","Oas1b",
         "Bst2","Mx1","Stat1","Irf7",
         "Irf7","Irf1","Ifitm3","Ifi47","Rsad2","Gbp2","Ubd","Ly6e",   # discovery IFN module
         "Sprr1a","Atf3","Gadd45a","Ecel1","Serpinb1a","Sox11","Gap43"]  # injury programme
```

And in manuscript §3.8, change *"the injury programme was neuron-intrinsic"* to:

> "Within the seven compartments that passed the cell threshold, the injury-associated change is
> carried by galanin, which is up in Neuron (+1.47 log₂, q = 1 × 10⁻⁴) and Neuron.NF (+2.54,
> q < 1 × 10⁻⁶). Because Atf3, Sprr1a, Gadd45a, Ecel1 and Serpinb1a were not included in the
> per-compartment panel, we localise galanin but do not localise the broader injury programme; that
> requires re-running the compartment analysis with those genes included."

---

### 3.3 【Problem】 Galanin is **not** neuron-restricted — it is up in endothelium and Schwann cells
### at comparable or higher significance, so "neuronal-injury programme localises to neurons" is
### false as written.

**【Evidence】** From `sc_GSE292522_celltype_DE.csv`, Gal per compartment (all recomputed by me
from the file):

| compartment | Gal logFC | t | q |
|---|---|---|---|
| Neuron | +1.472 | 10.17 | 1.0e-04 |
| Neuron.NF | +2.537 | 18.14 | 0.0 |
| **Endothelial** | **+0.671** | **6.29** | **1.0e-04** |
| **Schwann** | **+0.453** | **6.02** | **0.0061** |
| SatelliteGlia | +0.388 | 3.56 | 0.1402 |
| Fibroblast | +1.513 | 2.82 | 0.783 |
| Neuron.PE | +0.235 | 2.40 | 0.3404 |
| Erythrocyte | −0.100 | −0.33 | 0.9281 |

Endothelial Gal has q = 1.0e-04, **identical to the neuronal value**. Galanin is a vasoactive
peptide with well-documented endothelial and vascular expression; it is not a neuron-restricted
marker.

**【Why it matters】** Manuscript line 153: *"the neuronal-injury programme was neuron-intrinsic:
galanin was up +1.47 log₂ in neurons"*. Table 5 calls it "Neuronal-injury localisation". With
endothelial Gal at the same q, the localisation is to neurons **and endothelium**, and the paper
selected only the neuronal rows. This is the same selective-reporting pattern the paper
elsewhere criticises.

**【Specific fix】** Replace Table 5's row:

> | Neuronal-injury localisation | galanin +1.47 log₂ in Neuron (q = 1.0 × 10⁻⁴) and +2.54 in Neuron.NF (q < 1 × 10⁻⁶); also +0.67 in Endothelium (q = 1.0 × 10⁻⁴) and +0.45 in Schwann (q = 0.006); Scn10a −0.25 in Neuron (q = 0.012) |

and in §3.8 change "the injury programme localises to neurons" to:

> "the injury-associated change localises to neurons and, at comparable significance, to
> endothelial cells; galanin is not a neuron-restricted marker."

---

### 3.4 【Problem】 "No compartment shows interferon up-regulation" is **not** supported by the
### composite ISG score the pipeline already computed: satellite glia show a *significant*
### interferon **down**-regulation, which the manuscript reports but does not interpret.

**【Evidence】** `sc_GSE292522_celltype_summary.json`:
`SatelliteGlia: mean_Con = 0.5617, mean_B16F10 = 0.3067, MWU_p = 0.0286`.
All other compartments: Neuron p = 1.0, Neuron.NF p = 1.0, Neuron.PE p = 0.486, Fibroblast
p = 0.886, Schwann p = 0.114, Endothelial p = 0.200.

And from `sc_GSE292522_celltype_DE.csv`, satellite glia show two significant ISG *down*-
regulations the manuscript omits from §3.8's list: **Ifit1 −0.523 (q = 0.0037)** and
**Isg15 −0.63 (q = 0.0071)**. The manuscript mentions only "Isg15 in satellite glia (q = 0.007)".
Ifit1 at q = 0.0037 is the most significant ISG change in the entire compartment table and is not
mentioned.

Per-compartment ISG score changes I recomputed from `sc_GSE292522_celltype_scores.csv`:

| compartment | ISG score diff (B16−Con) | MWU p |
|---|---|---|
| SatelliteGlia | **−0.255** | **0.0286** |
| Endothelial | −0.111 | 0.200 |
| Schwann | −0.057 | 0.114 |
| Fibroblast | −0.024 | 0.886 |
| Neuron | +0.001 | 1.000 |
| Neuron.PE | −0.008 | 0.486 |
| Neuron.NF | +0.005 | 1.000 |

**【Why it matters】** "No compartment shows interferon up-regulation" is technically true but
chosen framing: the one compartment with a *significant* interferon change shows a **45% reduction**
in interferon score. A reviewer will read a significant satellite-glia interferon down-regulation
as a real, reportable observation that the paper has suppressed by phrasing it as a null. And
satellite glia are exactly where an IFN programme would be expected in DRG. The paper's §3.8
"absent at the compartment level" should be "absent or reduced".

**【Specific fix】** Rewrite the relevant clause of §3.8:

> "No compartment showed interferon-stimulated-gene **up**-regulation. Satellite glia instead
> showed a significant **reduction** in the composite interferon score (0.562 → 0.307, Mann–Whitney
> p = 0.029), accompanied by Ifit1 (−0.52 log₂, q = 0.0037) and Isg15 (−0.63, q = 0.0071) down.
> The interferon programme in this model is therefore not amplified but locally attenuated in
> satellite glia, the compartment in which it would be expected to originate. This is a change of
> direction rather than an absence, and we report it as such."

And add Ifit1 to the list in Table 5 and in §3.8.

---

### 3.5 【Problem】 The paper's §3.8 conclusion is **contradicted by the abstract of the dataset's
### own source publication**, which is cited as ref [15].

**【Evidence】** `data/raw/GSE292522_series_matrix.txt.gz`, `!Series_summary` (verbatim):

> "**Injured neurons respond by autonomously initiating an interleukin (IL)-6) and interferon (IFN)
> type I inflammatory response.** This inflammatory response alters the immune activity in the
> peri-neural niche… The CINI-driven anti-PD-1 resistance can be reversed by targeting multiple
> steps in the CINI signaling process: denervating the tumor, conditional knockout of the
> transcription factor mediating the injury signal within neurons (cKO-Atf3), knockout of the
> interferon-α receptor signaling (**Ifnar1-/-**), or by combining anti-PD-1 and anti-IL-6-receptor
> blockade."

The manuscript (line 155) concludes: *"The interferon programme is therefore absent at the
compartment level, whereas the injury programme localises to neurons — the per-compartment test
turns a null average into a compartment-resolved null."* And in the Abstract: *"An independent
single-cell dataset … showed a neuronal-injury programme (Sprr1a, Atf3) rather than the interferon
programme, and this held in every DRG compartment tested, confirming the interferon programme is
model-specific and not a universal cancer-pain signature."*

**【Why it matters】** This is not a matter of interpretation. The source paper reports that in this
exact model (B16F10, mouse DRG, L3–L5), injured neurons autonomously initiate a type-I interferon
response, and demonstrates it causally via Ifnar1 knockout. The manuscript analyses the same
dataset and concludes the interferon programme is absent. Either the manuscript's analysis misses a
signal its own source publication identified (most likely, given §3.1: the myeloid compartment was
dropped, and the tested ISG panel shares one gene with the discovery module), or the discrepancy
must be explained. Either way, the manuscript cannot present this dataset as a clean refutation of
the IFN reading without addressing what the source authors concluded and how their Ifnar1 result
constrains the interpretation.

Note also that `Ifit2` — the one ISG the manuscript highlights as "significantly down-regulated" —
has logFC = **−0.048**, i.e. a 3% reduction. Calling a −0.048 log₂ change evidence for IFN
*down*-regulation, while treating its parent paper's causal Ifnar1 evidence as not applicable,
is inconsistent.

**【Specific fix】** Add a new subsection §3.9 titled "Reconciliation with the source publication":

> "GSE292522 is the mouse DRG component of Baruch et al. (2025), whose reported conclusion is that
> injured neurons in this model autonomously initiate a type-I interferon and IL-6 response, with
> the phenotype reversed by Ifnar1 knockout, neuronal Atf3 deletion, or IL-6R blockade. Our
> reanalysis does not reproduce an interferon programme in this dataset. We report three candidate
> explanations and cannot distinguish them here. (i) Compartment coverage: the myeloid compartment,
> the expected site of interferon signalling, was excluded by our cell-count threshold (§3.8), so a
> myeloid-restricted programme would be invisible to us. (ii) Panel coverage: our per-compartment
> interferon panel shares only one gene with the eight-gene module that defined the programme in
> the Walker-256 data (§3.8). (iii) Proportionality: the largest interferon change we detect in any
> compartment is a 0.52 log₂ reduction, whereas the Walker-256 programme reaches +2.6 to +4.0 log₂;
> a programme an order of magnitude smaller would not be reliably detectable at n = 4 animals per
> group with the cell counts available. Our conclusion is therefore that a Walker-256-scale
> interferon programme is absent from this dataset — not that interferon signalling is absent.
> Distinguishing these requires either reanalysis of the authors' own Ifnar1 arm or a
> myeloid-enriched experiment."

And change the Abstract's "confirming the interferon programme is model-specific and not a
universal cancer-pain signature" to:

> "…showed a neuronal-injury programme (Sprr1a, Atf3) but no Walker-256-scale interferon programme in
> the seven compartments that met our cell-count threshold, in a dataset whose source publication
> separately reports a type-I interferon response in this model; the interferon programme we
> observe is therefore large and model-specific rather than a general feature of cancer-pain DRG."

---

### 3.7 【Problem】 The compartment labelled "SatelliteGlia" contains **no GFAP**. Its defining
### marker is absent (Gfap = 0.012) and its top-expressed gene is a myelin gene (Plp1 = 1.944), so
### the §3.8 satellite-glia result is measured in a compartment that is not satellite glia.

**【Evidence】** I re-ran the `scripts/32_percelltype_sc.py` assignment and then measured marker
expression *inside* each resulting compartment (mean log₁p, averaged over all eight animals):

| gene | "SatelliteGlia" compartment | "Schwann" compartment |
|---|---|---|
| **Plp1** (myelin proteolipid) | **1.944** ← top gene | 2.215 |
| **Gfap** (defining SGC marker) | **0.012** | 0.005 |
| Mpz | 0.170 | 1.733 |
| Sox10 | 0.577 | 0.533 |
| Glul | 1.337 | 0.607 |
| Kcnj10 | 1.131 | 0.505 |
| Aqp4 | 0.604 | 0.169 |
| Slc1a3 | 0.493 | 0.160 |
| **Calca** | **1.207** | 1.323 |
| **Prph** | **1.010** | 1.179 |
| Tubb3 | 0.871 | 0.976 |
| Nefl | 0.829 | 0.895 |

The compartment assigned to SatelliteGlia has **Gfap at background (0.012)** while expressing
Glul/Kcnj10/Aqp4/Slc1a3 (genuine astrocyte-like markers, and note it also carries the *myelin*
gene Plp1 as its single highest-expressed transcript). Independent marker-positive counting
confirms the real GFAP+ population is far too small to be this compartment: across all eight
animals, **Gfap+ cells number 38/25,672 in controls (0.148%) and 121/30,295 in B16F10 (0.399%)**
— roughly 100–500 cells total, against a "SatelliteGlia" compartment of 361 cells. Real satellite
glia are typically several percent of DRG cells.

Consistent with this, `sc_GSE292522_celltype_DE.csv` line 31 reads
`SatelliteGlia,Gfap,0.0,0.0,0.9997` — logFC exactly 0.0, t exactly 0.0, because Gfap is
all-zero inside that compartment and carries no variance.

**【Why it matters】** §3.8's satellite-glia findings — `Ifit1` −0.523 (q = 0.0037) and `Isg15`
−0.63 (q = 0.0071), and the significant satellite-glia interferon-score reduction
(0.562 → 0.307, p = 0.029) — are all computed on this compartment. If the compartment is not
satellite glia, those results cannot be attributed to satellite glia, and the "compartment-resolved"
framing of §3.8 loses its second-largest glial anchor. This compounds §3.1: the compartment
labelling is unreliable at both the myeloid end (dropped) and the glial end (mis-assigned).

Note also that the compartment is contaminated with neuronal transcripts (Calca 1.21, Prph 1.01,
Tubb3 0.87, Nefl 0.83) — it is a mixed low-quality population, plausibly doublets or cells with
ambient RNA, which is exactly the setting in which an argmax-on-mean-log-expression classifier
assigns cells to whichever marker set happens to have the highest *average* of a few detected
transcripts.

**【Specific fix】** Add a compartment-purity QC step to `scripts/32_percelltype_sc.py` and report
it as a table:

```python
# scripts/32_percelltype_sc.py, after building `agg` (append before the DE loop)
# Purity audit: a compartment is only interpretable if its OWN defining markers
# are actually detected inside it. Report mean log1p of each marker set within
# each assigned compartment, and flag any compartment whose defining marker
# mean is below 0.1.
PURITY_AUDIT = {}
for ctype in COMPARTMENTS:
    present = [p for p in SAMPLES if p in agg[ctype]]
    if not present: continue
    M = pd.DataFrame({p: agg[ctype][p] for p in present}, index=canon_sym)
    defset = (BROAD.get(ctype) or NEURON_SUB.get(ctype.split(".")[-1]) or [])
    gi = [canon_sym.index(g) for g in defset if g in canon_sym]
    PURITY_AUDIT[ctype] = {"n_cells": int(sum((np.asarray(agg[ctype][p])>0).sum() for p in present)),
                           "own_marker_mean_log1p": round(float(np.log1p(M.values[:, :1]).mean()), 4),
                           "genes_present": len(gi)}
with open(os.path.join(TAB, "sc_GSE292522_compartment_purity.json"), "w") as fh:
    json.dump(PURITY_AUDIT, fh, indent=2)
```

and in §2.8 add:

> "Marker-score assignment by argmax over eight marker sets is not sufficient to identify DRG
> compartments, and we audit it: within the compartment we call satellite glia, Gfap abundance is at
> background (mean log₁p 0.012) while the myelin gene Plp1 is the most abundant transcript
> (1.944). GFAP-positive cells number 38 of 25,672 in controls and 121 of 30,295 in B16F10 animals,
> which is far fewer than the 361 cells assigned to that compartment. We therefore report
> satellite-glia results as glial-compartment results with an explicit purity caveat, and we do not
> interpret them as satellite-glia-specific."

---

### 3.8 【Problem】 The paper's whole-DRG "satellite-glia marker Gfap was nominally up" is a
### 1.0025-fold change in a gene detected in 0.15–0.40% of cells.

**【Evidence】** `sc_GSE292522_pseudobulk_DE.csv`: `Gfap logFC = +0.00358`, i.e.
2^0.00358 = **1.0025-fold**, t = 5.83, q = 0.00615. Manuscript §3.8: *"the satellite-glia marker
Gfap was nominally up in B16F10 (q = 0.006)"*.

**【Why it matters】** A 0.25% fold change that reaches significance only because within-group
variance is tiny is a technical detection-level artefact, not biology — and it is presented as
supporting evidence. The manuscript elsewhere insists on effect sizes ("Effect sizes are large
(3.0 to 3.4 log₂ units)"), so accepting a 0.0036 log₂ change here is inconsistent with its own
standard. Combined with §3.5, the single piece of evidence offered for glial involvement in B16F10
is both compartment-misassigned and effect-size-negligible.

**【Specific fix】** In §3.8 replace *"whereas the satellite-glia marker Gfap was nominally up in
B16F10 (q = 0.006)"* with:

> "Gfap shows a statistically detectable but biologically negligible change in the whole-DRG
> pseudobulk (log₂ FC +0.004, i.e. a 1.0025-fold change, q = 0.006), reflecting the low abundance
> of GFAP+ cells (0.15–0.40% of cells) rather than a change in glial abundance. We do not treat
> this as evidence of satellite-glia activation."

---

## 4. Galanin as an injury marker

### 4.1 【Problem】 Galanin up-regulation is presented as *the* neuronal-injury signature with a
### single competing reading unaddressed.

**【Evidence】** `etiology_marker_de.csv`: Gal PTX = **+1.923** (t = 9.88, q = 0.0103);
Bone_D7 = **−0.688** (q = 0.0267); OXA = −0.295 (q = 0.554). Cross Δ (D7 vs PTX) = −2.592,
p = 6.12e-04, q = 0.0125. In single cell: Neuron +1.472 (q = 1.0e-04), whole-DRG pseudobulk +1.425
(q = 2.7e-04).

The standard interpretation (Hökfelt et al. 1987; Wiesenfeld-Hallin et al. 1992; Holmes et al.
2000): galanin is massively (up to 120-fold) up-regulated in DRG after peripheral nerve axotomy,
shifts from <5% of neurons to 40–50%, and is a marker of the **injury/regeneration programme**,
not of any particular injury aetiology.

Competing readings the paper does not consider:

1. **Galanin marks a neuronal *subtype* shift, not injury.** Galanin is constitutively expressed
   in a defined small-diameter nociceptor population. A change in galanin can reflect a change in
   the *proportion* of galanin-expressing neurons rather than per-neuron induction. In the single-cell
   data, Gal is up +2.54 in Neuron.NF but only +0.235 (q = 0.34) in Neuron.PE — a highly
   **subtype-restricted** pattern, consistent with a subtype-composition shift. Myeloid infiltration
   changing the neuron-subtype denominator would produce exactly this.
2. **Galanin is a trans-synaptic/descending modulator, not a lesion marker.** Increased galanin
   *decreases* spinal transmission and enhances morphine analgesia (Wiesenfeld-Hallin 1998).
   Treating it purely as an "injury" readout ignores that its functional consequence is
   *antinociceptive* — so "up = more injured neuron" is not a safe inference.
3. **Galanin is not DRG-restricted** (§3.3 above: endothelium q = 1.0e-04).

**【Why it matters】** Galanin carries the entire neuronal-injury arm of the paper's central
binary. If it is a subtype marker, then the "neuronal-injury-driven" claim reduces to "the
proportion of galanin-expressing nociceptors differs between preparations", which is a
composition claim — exactly what the paper insists cannot be made without cell resolution.

**【Specific fix】** Add to §5 Limitations:

> "Galanin is used here as the index of a neuronal-injury programme. This interpretation is standard
> for axotomy (Hökfelt et al. 1987; Holmes et al. 2000), where galanin rises up to 120-fold and
> shifts from under 5% to 40–50% of DRG neurons. It carries two caveats. First, galanin is
   constitutively expressed in a defined small-diameter nociceptor population, so a change in bulk
   galanin can reflect a change in the proportion of galanin-expressing neurons rather than
   per-neuron induction; the subtype-restricted pattern we observe (Neuron.NF +2.54 versus Neuron.PE
   +0.24, q = 0.34) is consistent with this. Second, increased galanin is functionally
   antinociceptive (Wiesenfeld-Hallin et al. 1992), so its up-regulation should not be read as a
   simple proxy for nociceptive neuron damage. Resolving this requires per-cell-type galanin
   quantification on a reference dataset, which bulk tissue cannot provide."

---

## 5. Must-cite literature that is absent

The reference list is [1]–[16]. The following claims are made without any supporting citation, and
each has a substantial primary literature that changes how a reader should interpret the result.

### 5.1 Satellite gliosis / GFAP specificity — **unsourced, and the controversy is directly on point**

Manuscript asserts "satellite gliosis (Gfap +1.33)" (line 115, Abstract) with no citation. Must
cite:

- **Mohr KM, Pallesen L, Richner M, Vaegter C. "Discrepancy in the Usage of GFAP as a Marker of
  Satellite Glial Cell Reactivity." *Biomedicines*. 2021;9(8):1022.** — the title is the finding.
- **Konnova EA, Deftu AF, Chung PCS, Pertin M, Kirschmann G, Decosterd I, Suter MR. "Characterisation
  of GFAP-Expressing Glial Cells in the Dorsal Root Ganglion after Spared Nerve Injury." *IJMS*.
  2023;24(21):15559.** — 75.5% of GFAP+ DRG cells after injury are FABP7-negative non-myelinating
  Schwann cells, electrophysiologically confirmed (Kir4.1-negative, Kv1.1/Kv1.6-positive).
- **Jager O, et al. "Changes in the transcriptional fingerprint of satellite glial cells following
  peripheral nerve injury." *Glia*. 2020;68:1375.** — if the "glial" score is to mean SGC.

**【Specific fix】** Add to References [17]–[19] as above, and add to §3.3:

> "GFAP is commonly used as a satellite-glia activation marker in DRG, but its specificity in this
> tissue is disputed: Mohr et al. (2021) document a discrepancy in its use, and Konnova et al. (2023)
  found that 75.5% of GFAP+ DRG cells after spared nerve injury are FABP7-negative non-myelinating
> Schwann cells rather than satellite glia. We therefore report Gfap as a reactive-glial signal of
  undetermined subtype."

### 5.2 The type-I interferon programme — **unsourced in a manuscript whose second finding is an interferon programme**

The paper's central discovery (§3.3) is a type-I IFN programme in bone cancer pain DRG. No
reference in [1]–[16] concerns interferon signalling in DRG. This is a substantive gap, and it
matters in both directions:

- If the Walker-256 IFN signal is real, the bone-cancer-pain literature on STING/TBK1/NF-κB
  (e.g. **"Activation of the STING pathway induces peripheral sensitization via neuroinflammation
  in a rat model of bone cancer pain"**, Walker 256 intratibial, DRG STING/TBK1/NF-κB activation
  with IL-1β/IL-6/TNFα up-regulation) provides a mechanistic frame the paper should connect to.
- If it is an artefact, the IFN literature on common IFN-driven technical artefacts in bulk RNA-seq
  (e.g. interferon signature from sequencing saturation/strand specificity, or from residual
  blood/muscle RNA) should be discussed.

**【Specific fix】** Add a paragraph to §4:

> "The day-7 type-I interferon programme we report is, to our knowledge, not previously described
> in Walker-256 bone cancer pain DRG. It is consistent with reports of STING/TBK1-mediated
  neuroinflammation in the same model, where DRG STING activation drives IL-1β, IL-6 and TNFα
  up-regulation. We note that interferon signatures are also a recognised feature of technical
  confounders in bulk RNA-seq, including residual blood and high-mitochondrial-content RNA; our
  leave-one-out analysis (Table S) shows the programme is not driven by the muscle-contaminated
  sample, but we cannot exclude contributions from blood-borne RNA in the day-7 tissue."

### 5.3 Compositional dilution — asserted as a methodological fact with no citation

§3.6 and §5 rely on "bulk tissue is compositional, so an influx of macrophages lowers the
transcript share of every neuronal gene" (Introduction, line 35). No citation. Must cite a
deconvolution/compositional-bias reference so readers can see the assumption's provenance and its
known limitations (deconvolution is not reliable in non-model systems; accuracy depends on reference
choice and is itself contested).

**【Specific fix】** Add as reference [20] a bulk-compositional-bias/deconvolution methods citation
(e.g. a deconvolution benchmarking paper), and add to §5:

> "Bulk-tissue compositional bias is a well-characterised limitation of differential expression in
  heterogeneous tissues, and computational deconvolution using matched single-cell references is the
  standard mitigation. We did not apply deconvolution because no rat DRG single-cell reference
  matched either study's platform and timepoint; doing so with an unmatched reference would
  introduce more error than it removes. This is a limitation of the present design, not evidence
  that composition is unimportant."

### 5.4 Reference [15] is mis-described as a single-cell DRG dataset only

Ref [15] is cited as "B16F10 cancer-nerve-injury DRG". The GEO record shows the study provides
**six** GEO accessions: GSE289743 (human clinical trial RNA-seq), GSE289744 (neuron co-culture),
GSE289745 (Visium spatial), **GSE292090** (tumour scRNA-seq), **GSE292089** (bulk mouse ganglia),
GSE292522 (ganglia scRNA-seq). GSE292089 is bulk RNA-seq of mouse ganglia from the same animals —
a far better-matched validation target than the cross-species, cross-platform GSE292522.

**【Why it matters】** The manuscript's central validation is a mouse-vs-rat, scRNA-vs-bulk
comparison. An available bulk mouse-ganglia dataset from the same study and same animals would
remove the species and modality confounds simultaneously. Not using it is a missed opportunity that
a reviewer will ask about.

**【Specific fix】** Add to §2.8:

> "A second dataset from the same study, GSE292089, provides bulk RNA-seq of mouse ganglia from
> the same animals on the same platform; it is a closer match to our discovery data (bulk, DRG
> tissue) than GSE292522 (single-cell, mouse). We analyse GSE149648 against GSE292089 as the
> primary cross-model comparison and retain GSE292522 as the cell-type-resolved secondary analysis."

---

## 6. Clinical pharmacological translation

### 6.1 【Problem】 "The split matches clinical pharmacology" is an over-reach: the cited references
### establish that bone pain has treatment options beyond opioids; they do **not** establish that
### CIPN and bone pain differ in pharmacological responsiveness in the way the sentence implies.

**【Evidence】** Manuscript line 161:

> "The split matches clinical pharmacology: chemotherapy-induced neuropathic pain responds partly to
> agents acting on neuronal excitability, whereas bone pain carries an inflammatory component that
> responds to anti-resorptive and anti-inflammatory strategies [1,2,14]."

Refs [1,2] are ESMO 2018 and NCCN 2019 adult cancer pain guidelines; [14] is Luger, Mach, Sevcik,
Mantyh 2005, "Bone cancer pain: from model to mechanism to therapy".

What these establish:
- [1,2,14] **do** support: bone cancer pain responds to NSAIDs, bisphosphonates/denosumab and
  radiotherapy in addition to opioids; bone pain has an inflammatory and osteoclast-driven
  component. ✔
- [1,2,14] **do not** establish: that CIPN pain "responds partly to agents acting on neuronal
  excitability." The clinical reality is the opposite in one respect — the only pharmacotherapy with
  established efficacy for painful CIPN is **duloxetine**, an SNRI acting on monoaminergic
  descending inhibition, not on neuronal excitability. Sodium-channel blockers, mexiletine and
  lamotrigine have all failed or shown inadequate benefit in CIPN trials. Citing cancer-pain
  guidelines for a claim about CIPN pharmacology is a category error: neither [1] nor [2] is a CIPN
  treatment guideline, and [6] (the ASCO CIPN guideline, Loprinzi 2020) is in the reference list but
  is **not cited at this sentence**, even though it is the authoritative source for what does and does
  not work in CIPN.

Additionally, the manuscript's own data undercut the symmetry: oxaliplatin neuropathy shows **no**
injury programme and **no** neuronal-excitability-relevant change (§2.4 above), so "the chemo arm is
neuronal" holds for one drug only.

**【Why it matters】** This is the sentence that converts a transcriptional reanalysis into a
clinical recommendation. As written it tells a reader that the transcriptional split predicts a
pharmacological split. Nothing in the manuscript establishes that, and the one authoritative CIPN
guideline in its own reference list contradicts the excitability half of the claim.

**【Specific fix】** Rewrite the Discussion (line 161) final sentence as:

> "The transcriptional split is consistent with, but does not predict, the pharmacological split.
> Bone cancer pain carries an inflammatory and osteoclast-driven component and responds to
> anti-resorptive and anti-inflammatory strategies — bisphosphonates, denosumab, NSAIDs and
> radiotherapy — in addition to opioids [1,2,14]. For chemotherapy-induced neuropathic pain the
> evidence base is narrower and does not support an excitability-based reading: duloxetine is the
> only agent with established efficacy, acting through monoaminergic descending inhibition, and
> sodium-channel blockers have not shown adequate benefit [6]. We note explicitly that our
> transcriptional dichotomy resolves as 'bone versus paclitaxel' (§3.5), not 'bone versus CIPN',
> and we did not test any pharmacological predictor. Connecting a transcriptional programme to a
  treatment response would require a design that includes drug exposure and a pharmacological
  readout, neither of which is present in either dataset."

---

## § Stands up (things I suspected and found to be correct)

1. **The day-7 interferon programme is real and robust.** I expected the muscle-contaminated
   sample GSM4506604 (6.26% sarcomere) to be driving it. It is not. Recomputing with the pipeline's
   own design and dropping each day-7 sample in turn: ISG module mean moderated t goes
   24.92 → 21.72 / 21.81 / **21.89** / 22.62 (dropping the contaminated sample gives 21.89, and
   Cd68 *rises* from +2.86 to +2.96). The programme survives. `etiology_marker_de.csv` confirms all
   8 module genes rise (Ubd +4.03, Ifi47 +3.32, Gbp2 +3.25, Rsad2 +2.40, Bst2 +2.34, Irf1 +2.31,
   Ifitm3 +1.97, Ly6e +1.29; all q ≤ 5.8e-06), and an independent 15-gene panel I ran on
   `DE_D7_vs_Sham_adj.csv` shows the programme is broader than the module (Isg15 +2.35, Mx1 +2.25,
   Ifit3 +2.24, Stat1 +1.74, Oas1a +1.92, Cxcl10 +3.47, Ccl5 +4.56, Irf7 +2.67, Psmb8 +2.99,
   Psmb9 +3.08, Tap1 +2.48). This is a genuine type-I IFN response.

2. **The claim that neuronal modules fall as much as the pan-neuronal control is correct.**
   Recomputed: at D7, CTRL pan-neuronal mean logFC = −0.40, and every neuronal target module sits in
   −0.24 to −0.52 (K⁺/HCN −0.48, Neurotrophin −0.52, GABAergic −0.37, Na⁺ −0.42, glutamate −0.38,
   opioid −0.40, Ca²⁺α2δ −0.24). The 10 pan-neuronal structural genes themselves fall uniformly
   (Nefl −0.48, Nefm −0.49, Nefh −0.50, Snap25 −0.46, Tubb3 −0.34, Prph −0.48, Rbfox3 −0.44,
   Syt1 −0.31). This is exactly the pattern a change in neuronal transcript share produces, and the
   paper is right to draw that inference.

3. **The scRNA pseudobulk negative result for myeloid markers is a well-powered null.**
   `sc_GSE292522_pseudobulk_DE.csv`: Cd68 logFC = −0.0005 (t = −0.54, q = 0.78),
   Aif1 = +0.0009 (t = 0.24, q = 0.92), Lyz2 = −0.029 (q = 0.196), Tyrobp = −0.002 (q = 0.72).
   The design detects Ifit2 at logFC −0.048 with q = 0.042, so a null on Cd68 at this magnitude is
   informative for large effects. The paper's "genuine null for large effects" statement is correct
   *for the whole-DRG pseudobulk*.

4. **The myeloid-fraction non-significance is correctly reported as inconclusive rather than as
   evidence of equality.** I recomputed the marker-positive cell counts independently of the
   pipeline's argmax classifier and the direction is confirmed but the magnitude is not robust:
   ≥3 macrophage markers positive in 28/25672 control cells (0.11%) and 28/30295 B16F10 cells
   (0.09%). Lyz2⁺ 257 vs 184; Cd68⁺ 229 vs 202. So the true myeloid fraction is ~0.1%, and the
   paper's ~0.05% figure from the argmax classifier is in the same range. The paper is right that
   d = 1.37 against MDE 2.38 cannot support an equality claim.

5. **`Scn10a` = Nav1.8 and `Atf3`/`Sprr1a` = neuronal injury are correct assignments for mouse and rat
   DRG.** No issue with the marker identities themselves.

6. **The cross-aetiology permutation statistics are computed on the right quantities.** I
   re-derived Cd68 D7-vs-PTX Δ = +3.418 (logFC 2.839 − (−0.579)), Aif1 Δ = +3.076, Gal Δ = −2.592,
   all with p = 6.12e-04 from `etiology_marker_cross.csv`, matching Table 2. The direction and
   magnitude of the myeloid-vs-galanin divergence are correctly computed.

---

## § Questions for the authors

1. **Myeloid compartment.** GSE292522 contains 46 myeloid cells and 1 T cell across 55,967 cells.
   Was the myeloid compartment analysed separately in any way (pooled, targeted, sorted)? If not, on
   what basis do you claim the interferon programme is absent "in every DRG compartment" (§5) when
   the compartment in which interferon signalling is expected was excluded by a cell-count filter?
2. **Source publication.** Baruch et al. (ref [15]) report in the GSE292522 abstract that injured
   neurons autonomously initiate a type-I IFN response, and reverse the phenotype with Ifnar1
   knockout. How do you reconcile your null with their positive? Was their data re-analysed, or only
   the GEO matrices?
3. **Compartment ISG panel.** Seven of the eight genes that define the discovery interferon module
   (`Irf1, Ifitm3, Ifi47, Rsad2, Gbp2, Ubd, Ly6e`) were not in `FOCUS`. Was this deliberate? If
   not, this is a one-line fix and it materially changes the result.
4. **Injury genes per compartment.** `Atf3`, `Sprr1a`, `Gadd45a`, `Ecel1`, `Serpinb1a` were not in
   `FOCUS`, yet §3.8 asserts the injury programme "localises to neurons" on the basis of galanin
   alone. Were these genes excluded because they are not expressed in mouse DRG, or by oversight?
5. **Galanin in endothelium.** Gal is up in Endothelial at q = 1.0e-04 — identical to the neuronal
   value you report. Why is this not in Table 5? Is "neuron-intrinsic" intended to mean
   "highest in neurons" rather than "restricted to neurons"?
6. **Ifit1 in satellite glia.** Ifit1 (−0.52, q = 0.0037) is the most significant ISG change in the
   compartment table and is not mentioned in §3.8 or Table 5. Was it excluded, or overlooked?
7. **Sham timepoint.** The control animals are day-14 PBS (GEO characteristic). Was a day-1 or
   day-7 control available in the series? The "0 genes at day 1" result compares day-1 tumour to a
   day-14 wound.
8. **Erythrocyte contrast.** The Δ = −2.48 is D7-vs-OXA; the D7-vs-PTX erythrocyte Δ is −0.26
   (q = 0.81). Was the D7-vs-OXA value intended as the compositional control for the D7-vs-PTX
   neuronal claims?
9. **OXA phenotype.** Oxaliplatin shows Cd68 +0.44 (P = 5.5e-04) and no injury markers. Should the
   paper report three phenotypes rather than two?
10. **"Glial / immune" module.** Five of its six genes are myeloid (`Gfap, Aif1, Itgam, Cd68, S100a8,
    S100a9`). Was a glial-intrinsic module intended?
11. **GSE292089.** Why was bulk mouse-ganglia RNA-seq from the same study not used as a closer-matched
    validation than cross-species single-cell data?
12. **Clinical pharmacology.** Was the claim in the Discussion about neuronal-excitability agents in
    CIPN sourced to [1,2,14] (cancer-pain guidelines) or to [6] (ASCO CIPN guideline, in your
    reference list but not cited there)?
13. **Compartment purity.** Inside the compartment you call satellite glia, Gfap is at background
    (0.012 mean log₁p) while Plp1 is the top-expressed gene (1.944), and true GFAP+ cells number
    only 38/25,672 (control) and 121/30,295 (B16F10). Was a purity audit run? If the compartment is
    not satellite glia, on what basis are the Ifit1/Isg15 down-regulations attributed to them?
14. **Real satellite glia.** GFAP+ cells are 0.15–0.40% of cells. Published DRG dissections yield
    several percent satellite glia. Was a stricter GFAP/Slc1a3/Kcnj10/FABP7 gate or a cluster-based
    (rather than argmax) annotation considered?
15. **Gfap whole-DRG significance.** Gfap logFC = +0.0036 (1.0025-fold) reaches q = 0.006. Was an
    effect-size floor applied before calling a gene "nominally up"?

---

## § What I actually checked

### Files read
`docs/manuscript/manuscript_final.md` (full); `_PANEL_BRIEF.md`;
`results/tables/`: `etiology_marker_de.csv` (all 68 rows), `etiology_marker_cross.csv` (all rows),
`etiology_module_logfc_cross.csv` (all control-module rows), `sc_GSE292522_celltype_DE.csv` (all
198 rows), `sc_GSE292522_celltype_scores.csv`, `sc_GSE292522_celltype_summary.json`,
`sc_GSE292522_composition.csv`, `sc_GSE292522_pseudobulk_DE.csv` (all 36 ISGs + key markers),
`sc_GSE292522_summary.json`, `negative_result_quantification.json`, `module_level_signal.csv`
(all 3 contrasts), `de_muscle_drop_sensitivity.csv`;
`DE_D1/D7/D14_vs_Sham_adj.csv` (ISG, myeloid, complement, MHC, epithelial panels);
`scripts/_limma.py` (module definitions, lines 157–200), `scripts/20_multistage_pipeline.py`
(design matrix, lines 295–350), `scripts/32_percelltype_sc.py` (full, marker sets and filters),
`scripts/12_adjust_covariate_limma.py`;
GEO `series_matrix` metadata for GSE149648, GSE160543, GSE292522.

### Values recomputed from source data (mine vs manuscript)

| quantity | my recomputation | manuscript | match? |
|---|---|---|---|
| Cd68, Bone_D7 logFC | +2.858 | +2.86 (§3.3) | ✔ |
| Aif1, Bone_D7 logFC | +2.180 | +2.18 (§3.3) | ✔ |
| Gfap, Bone_D7 logFC (DE table) | +1.294 | +1.33 (§3.3) / +1.50 (Table 2) | ✘ **three values** |
| Gal, PTX logFC | +1.923 (q = 0.0103) | +1.92 | ✔ |
| Gal, Bone_D7 logFC | −0.688 (q = 0.027) | −0.67 | ✔ |
| Cd68 Δ D7-vs-PTX | +3.418 (p = 6.12e-04) | +3.42, q = 0.012 | ✔ |
| Aif1 Δ D7-vs-PTX | +3.076 | +3.08, q = 0.012 | ✔ |
| Gal Δ D7-vs-PTX | −2.592 | −2.59, q = 0.012 | ✔ |
| erythrocyte Δ **D7-vs-PTX** | **−0.258, p = 0.716, q = 0.814** | not stated | ✘ **contrast mismatch** |
| erythrocyte Δ D7-vs-OXA | −2.483, q = 0.0073 | −2.48, q = 0.007 | ✔ (but wrong contrast for the argument) |
| D7 CTRL interferon mean moderated t | +24.92 | +24.92 (Table 1) | ✔ |
| D7 CTRL pan-neuronal mean t | −5.20 | −5.20 | ✔ |
| D7 CTRL erythrocyte mean t | −4.80 | −4.80 | ✔ |
| D7 n genes q<0.05 | 4600 | 4,600 | ✔ |
| ISG mean t, drop GSM4506604 | +21.89 (n q<0.05 = 3041) | not reported | — **new analysis** |
| D14 IFN module mean logFC / t | +2.37 / +4.36 | not discussed | — **contradicts "peak" framing** |
| D1 vs D14-sham n q<0.05 | 0 | 0 | ✔ (but control is D14-matched) |
| OXA Cd68 | +0.442 (P = 5.5e-04) | "slightly down" (§4) | ✘ **contradiction** |
| sc pseudobulk Cd68 / Aif1 | −0.0005 / +0.0009 | ≈ 0 | ✔ |
| sc pseudobulk Gfap | +0.00358 (**1.0025-fold**), q = 0.006 | "nominally up" | ~ (effect is 0.25%) |
| sc pseudobulk ISG up at q<0.05 | **1** — Ifitm3 +0.365 (q = 0.0127) | "no ISG was up-regulated" | ✘ **contradiction** |
| compartment cell counts | Myeloid 46 total (5/8 animals ≥5) → dropped | 8 compartments reported, "every compartment" | ✘ **contradiction** |
| bulk IFN module ∩ compartment FOCUS | 1 gene (Bst2) | implied same programme | ✘ |
| injury genes ∩ compartment FOCUS | 0 | "injury programme localises to neurons" | ✘ |
| Gal per compartment | Endo +0.671 (q = 1.0e-04), Schwann +0.453 (q = 0.006) | "neuron-intrinsic" | ✘ |
| satellite-glia ISG score | 0.5617 → 0.3067, MWU p = 0.0286 (down) | not reported | ✘ **suppressed** |
| satellite-glia Ifit1 | −0.523, q = 0.0037 | not reported | ✘ **suppressed** |
| marker-positive cell counts (independent of argmax) | ≥3 macrophage markers: 28/25672 Con, 28/30295 B16 | argmax myeloid 0.12% / 0.05% | consistent in direction |
| Gfap+ cells (marker-positive count) | 38/25672 Con (0.148%), 121/30295 B16 (0.399%) | not reported | ✘ **compartment mis-assignment** |
| Gfap inside "SatelliteGlia" compartment | **0.012** mean log₁p; top gene = Plp1 (1.944) | presented as satellite glia | ✘ **compartment mis-assignment** |
| Gfap inside "Schwann" compartment | 0.005 mean log₁p | — | — |
| sc pseudobulk Gfap fold change | 2^0.00358 = **1.0025×** | "nominally up" | ~ (effect negligible) |

### Discrepancies I could not resolve
- The three Gfap values (+1.294 / +1.328 / +1.498) differ across estimators; I did not determine
  which pipeline stage produced each, only which artefact each matches.
- `etiology_marker_de.csv` OXA Cd68 (logFC +0.442, P = 5.5e-04) is nominally significant at
  P < 0.001 but q = 0.068 after BH; I report both, and the manuscript reports only the Table 2 value.

### Literature checked
Walker-256 bone cancer pain transcriptomics (including the GSE149648 parent paper, and its
requirement that samples are intact DRG on BGISEQ-500); galanin up-regulation after axotomy
(Hökfelt 1987, Holmes 2000, Wiesenfeld-Hallin 1992/1998); GFAP specificity in DRG
(Mohr 2021, Konnova 2023); bulk compositional bias and deconvolution limitations;
Baruch et al. *Nature* 2025 and the GSE292522 series summary.