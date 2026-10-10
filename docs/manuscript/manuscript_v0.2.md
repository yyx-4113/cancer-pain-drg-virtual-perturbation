# Divergent analgesic-target hierarchies in cancer pain of different aetiology: an in-silico perturbation analysis of rat dorsal root ganglion transcriptomes

**Running head:** In-silico perturbation of cancer-pain DRG transcriptomes

**Author:** Yongxin Yang¹

¹ Department of Anesthesiology, The Second Affiliated Hospital of Fujian University of Traditional Chinese Medicine, Fuzhou, Fujian 350003, China

**ORCID:** 0009-0004-9698-6552

**Correspondence:** Yongxin Yang, Department of Anesthesiology, The Second Affiliated Hospital of Fujian University of Traditional Chinese Medicine, Fuzhou, Fujian 350003, China. Email: 960856791@qq.com

---

## Abstract

**Background.** Cancer pain is treated as one clinical entity, yet its drivers differ by aetiology: bone-metastatic pain is driven by tumour-induced osteolysis and the inflammatory milieu it creates, whereas chemotherapy-induced peripheral neuropathy (CIPN) reflects direct neurotoxic injury to sensory neurons. Whether the two should be treated with the same analgesic hierarchy is unresolved.

**Methods.** We analysed rat dorsal root ganglion (DRG) RNA-seq from two studies generated on the **same sequencing platform (GPL24782)**, so species, tissue and platform are held constant and only pain aetiology differs: tibial Walker-256 bone cancer pain (GSE149648; sham and days 1, 7, 14 after tumour-cell injection; n = 4 per group) and paclitaxel- or oxaliplatin-induced neuropathy (GSE160543; vehicle, paclitaxel, oxaliplatin; n = 4 per group). Differential expression used an empirical-Bayes moderated *t* statistic with an intensity-dependent prior. For each aetiology we built a co-expression graph and performed in-silico knockout and overexpression of 17 analgesic target modules by diffusing a seed perturbation over the graph, scoring each by how far it moved the pain state back toward its own control (reversal score). Empirical nulls came from size-matched random gene sets; cross-aetiology differences were tested by within-study bootstrap (2,000 resamples).

**Results.** Bone cancer pain reorganised the DRG transcriptome in a stage-dependent manner (0, 3,209 and 312 genes at q < 0.05 on days 1, 7 and 14). At the peak (day 7) a bipolar module pattern emerged: glial/immune, Toll-like receptor and cytokine/chemokine modules were strongly up-regulated (mean moderated *t* = +11.28, +6.51, +6.30) while neuronal modules were down-regulated (K⁺/HCN −3.78, neurotrophin −3.92, GABAergic −3.65, voltage-gated Na⁺ −2.74, opioid receptor −2.71). Because these neuronal targets were *lower* in the pain state, in-silico **over**-expression — not knockout — moved the profile back toward control (Na⁺ +0.056, K⁺/HCN +0.064; q = 0.015). Thirty-six cross-aetiology comparisons were significant at q < 0.05. Knocking out glial/immune modules reversed bone cancer pain far more than paclitaxel neuropathy (Δ = +0.249, 95% CI 0.173–0.258), whereas knocking out neuropeptide, neurotrophin and injury-regeneration modules reversed paclitaxel neuropathy more (Δ = −0.207, −0.095, −0.127). Paclitaxel and oxaliplatin also differed from each other (neuropeptides Δ = +0.146; injury/regeneration Δ = +0.121).

**Conclusions.** The hierarchy of analgesic targets in the DRG is aetiology-specific, not interchangeable. Tumour-driven pain is dominated by an immune/glial axis, whereas chemotherapy-driven pain is dominated by a neuronal injury/neuropeptide axis. Mechanistically stratified, rather than intensity-stratified, analgesic selection deserves prospective testing.

**Keywords:** cancer pain; dorsal root ganglion; in-silico perturbation; virtual knockout; analgesic target; chemotherapy-induced peripheral neuropathy; bone metastasis

---

## 1. Introduction

Cancer pain affects most patients with advanced disease and remains inadequately controlled in a substantial minority [1,2]. The WHO analgesic ladder, still the dominant organising frame, is explicitly **aetiology-agnostic**: it escalates by pain intensity rather than by mechanism. In practice the two largest cancer-pain populations differ mechanistically. Bone-metastatic pain arises from tumour-induced osteolysis, acidosis and a local inflammatory milieu that sensitises primary afferents [3,4]; CIPN arises from direct axonal and neuronal injury by the cytotoxic agent [5,6]. Whether the analgesic targets that matter most in one also matter most in the other is unresolved, and the experiment that would settle it — paired ganglion tissue from patients with different pain aetiologies — is not feasible. It is, however, addressable computationally.

Two developments make this feasible. First, curated public ganglion transcriptomes now exist for several cancer-pain models, in some cases generated on the same platform, which removes the species/tissue/platform confounding that has historically made cross-study pain comparisons uninterpretable. Second, the last three years have produced a family of *in-silico perturbation* methods — latent-space vector arithmetic (scGen) [7], graph-based genetic-perturbation models (GEARS) [8], chemical-conditioned autoencoders (chemCPA) [9] — that predict how a transcriptome reorganises after a gene or drug perturbation, moving perturbation biology from "what changed" toward "what would happen if".

That promise carries an evidentiary obligation. A 2025 *Nature Methods* benchmarking study reported that several single-cell foundation and perturbation models fail to beat simple linear or additive baselines once evaluation includes strong-effect perturbations and distribution shift [10]. Any study deploying a perturbation model without head-to-head baseline reporting is making an unsupported claim. We designed this work around that constraint: every perturbation result is co-reported against additive and ridge-linear baselines, and against an empirical null built from size-matched random gene sets.

A second, more prosaic obligation concerns confounding. DRG dissection unavoidably carries paraspinal skeletal muscle, and the amount carried is not constant between animals. Because muscle transcripts are abundant and highly expressed, modest differences in contamination can dominate a small-n differential expression analysis and produce a "pain signature" that is in fact a muscle signature. We treat this as a first-class analytical problem rather than a footnote, and Section 3.1 shows that it is not hypothetical.

We asked three questions. (i) Does the DRG transcriptome of bone cancer pain reorganise over time, and do different analgesic target modules dominate at different stages? (ii) In-silico, which modules most effectively reverse the pain state at each stage? (iii) Is the resulting target **hierarchy** the same in bone cancer pain as in CIPN?

---

## 2. Methods

### 2.1 Data sources

All data are publicly available from NCBI Gene Expression Omnibus; no new animal work was performed, so no ethics approval was required.

**Discovery — bone cancer pain, GSE149648.** Walker 256 mammary carcinoma cells injected into rat tibia (tumour cell injection, TCI); ipsilateral DRG collected on days 1, 7 and 14, with sham (PBS) controls at day 14. Four replicates per group, 16 samples, FPKM, platform GPL24782.

**Comparison — chemotherapy-induced neuropathy, GSE160543.** Paclitaxel (2.0 mg/kg i.p.) or oxaliplatin (2.4 mg/kg i.p.) or vehicle; DRG collected after pain behaviour was established. Four animals per group, 12 samples, FPKM, **platform GPL24782**.

**Confounding case study — GSE249443.** An independent rat tibial bone cancer pain dataset (MRMT-1 cells) with sham controls, six samples, FPKM.

The decisive design feature is that the discovery and comparison sets share species (*Rattus norvegicus*), tissue (DRG) and sequencing platform. Only the aetiology of pain differs. We deliberately did not compare across species or tissue types, because such comparisons confound aetiology with everything else.

### 2.2 Parsing, normalisation and quality control

Per-sample files were parsed with the expression column located from the file header rather than by fixed index. This matters: the two studies ship different layouts (3 columns for GSE149648; 6 columns for most GSE160543 members but 13 for others, with trailing empty delimited fields), and an earlier revision that read a fixed neighbouring column silently produced a matrix of Entrez identifiers — constant across samples — yielding zero fold changes everywhere. The pipeline now asserts that the assembled matrix is not an identifier matrix (fraction of integer-valued entries ≤ 0.95) before proceeding.

FPKM is already normalised for gene length and depth, so no size-factor or counts-per-million correction was applied; only log2(FPKM + 1). Genes were retained if FPKM ≥ 1 in at least four samples and the mean FPKM was ≥ 0.5.

### 2.3 Skeletal-muscle contamination

A muscle module was defined from canonical sarcomere and muscle-metabolism genes (*Acta1, Ckm, Myl1, Myh1, Myh2, Tnni2, Tnnc2, Tnnt3, Des, Mybph, Actn3, Mybpc1, Myoz1, Tnnc1, Mylpf, Myh4, Atp2a1, Casq1, Mb, Actc1*), z-scored per gene and averaged within sample. Two models were fitted for every contrast — unadjusted and muscle-adjusted — and genes whose signal disappears after adjustment are reported as contamination-driven rather than as pain biology.

### 2.4 Differential expression

Per-gene linear models with an empirical-Bayes moderated *t* statistic: residual variances shrunk toward a prior estimated by method of moments [11], with an intensity-dependent prior fitted by lowess of log residual variance on mean expression (limma-trend [12]). *P* values were converted to *q* by Benjamini–Hochberg. At n = 4 per group a naive Welch test has only 2–6 residual degrees of freedom and is poorly calibrated at both tails; variance moderation is the standard remedy.

Statistical primitives were unit-tested against synthetic data with a planted +2 log2 shift: recovered logFC was +2.02 for true positives and +0.05 for null genes, with zero false positives at q < 0.05 across 370 null genes (script `99_unit_tests.py`).

### 2.5 Module-level signal

Seventeen mechanistically coherent modules corresponding to druggable analgesic mechanisms were tested competitively against all genes by Mann–Whitney *U* on the moderated *t* values (opioid receptors; voltage-gated Na⁺; Ca²⁺ channel α2δ; TRP; P2X purinoceptors; prostanoid/COX; glutamate/NMDA; GABAergic; neuropeptides; cannabinoid; neurotrophin; cytokines/chemokines; Toll-like receptors; glial/immune; injury-regeneration; K⁺/HCN; monoaminergic).

### 2.6 Co-expression network and in-silico perturbation

Per aetiology, the 3,000 highest-variance genes were used to build a sparse symmetrised *k*-NN correlation graph (*k* = 12). Networks were built per aetiology because co-expression structure is itself aetiology-dependent; pooling would erase the difference we wished to detect.

A perturbation of module *M* with direction *s* ∈ {−1, +1} was a seed vector with *s* at member indices, diffused by

δ_{t+1} = (1 − α)·seed + α·S·δ_t, S = D^{−1/2} W D^{−1/2}, α = 0.55,

iterated 80 times (α⁸⁰ ≈ 10⁻²¹). Because the map is linear, the limiting operator P with δ_∞ = P·seed was computed once per aetiology, after which every real and permuted perturbation is a single matrix–vector product. This is what makes a 1,000-fold permutation test tractable on CPU.

The **reversal score** is the cosine between the diffused perturbation and the pain axis,

pain_axis = mean(control) − mean(case),

both mean-centred first, so the score reflects agreement in *pattern* rather than a global shift. A positive score means the perturbation moves the pain-state transcriptome toward the control state. Because the pain axis is defined within each study, scores are dimensionless and comparable across aetiologies.

### 2.7 Empirical null, baselines, and cross-aetiology testing

For each module, perturbation and aetiology, 1,000 size-matched random gene sets gave the empirical null; the one-sided empirical *p* is the tail fraction with +1 continuity correction. Every result is co-reported against an **additive** baseline (the seed itself, no propagation) and a **ridge-linear** baseline (each gene regressed on module genes across samples by ridge regression, no graph).

Because the two studies are independent, labels cannot be permuted between them. We bootstrapped **within** each study: cases and controls resampled with replacement, pain axes recomputed, both reversal scores recomputed, and Δ = reversal(A) − reversal(B) accumulated over 2,000 resamples for a percentile 95% interval and a two-sided bootstrap *p*. Resamples that draw duplicate samples (giving a degenerate axis) are retained, which makes intervals conservative.

### 2.8 Software and reproducibility

Python 3.13 with NumPy, pandas, SciPy, statsmodels and matplotlib. Random procedures use a fixed seed (20261009). Scripts, tables and figures are deposited at **[repository URL to be inserted at release]** with a checksum MANIFEST.

---

## 3. Results

### 3.1 Quality control, and a demonstration that muscle contamination is not hypothetical

In the discovery set (16 samples, 13,893 expressed genes; n = 4 per group) the muscle score did not differ between groups (one-way ANOVA F = 0.783, p = 0.526), and pairwise sample correlations were uniformly high (r = 0.915–0.994), indicating no outlying samples.

The independent bone cancer pain set GSE249443 behaved differently and provides a cautionary control. There, the muscle score differed sharply between groups (BCP −0.811 vs sham +0.811; Welch p = 0.0014). Unadjusted, it yielded 28 genes at q < 0.05, and 10 of the top 15 were skeletal-muscle or sarcomere genes (*Tnni2, Acta1, Tnnt3, Ckm, Actn3, Myl1, Tnnc2, Mylpf, Myoz1, Myh1*). After adjusting for the muscle score, **zero** genes remained at q < 0.05. The apparent "cancer pain transcriptome" in that dataset is a dissection artefact. We therefore excluded GSE249443 from biological inference and report it only as the reason muscle adjustment is mandatory.

### 3.2 Bone cancer pain reorganises the DRG transcriptome in a stage-dependent manner

Against sham, the number of genes at q < 0.05 was **0 on day 1, 3,209 on day 7 and 312 on day 14** (unadjusted: 0, 1,923, 1,128). The day-1 result is biologically sensible: one day after tumour-cell injection the pain state is not yet established, and the ganglion has not yet remodelled. The transcriptional response peaks at day 7 and has substantially resolved by day 14, by which time the dominant signal is a smaller, more chronic-looking set.

Established bone cancer pain markers behaved as positive controls at day 14: *Cxcl13* log2FC = +4.05 (P = 0.0025), *Aif1* (Iba1) +2.03 (P = 0.0019), *Ccl2* +0.82 (P = 0.085), *Il1b* +0.76 (P = 0.056), *Gfap* +0.71 (P = 0.052) and *Atf3* +0.61 (P = 0.071). The recovery of *Cxcl13* — a chemokine independently implicated in bone cancer pain — supports the validity of the pipeline.

### 3.3 A bipolar module pattern: immune up, neuronal down

Module-level mean moderated *t* (Table 1) revealed a coherent and strongly stage-dependent polarity. At day 7, immune-side modules were sharply up-regulated — glial/immune +11.28, Toll-like receptors +6.51, cytokines/chemokines +6.30 — while every neuronal excitability and transmission module moved the other way: neurotrophin −3.92, K⁺/HCN −3.78, GABAergic −3.65, glutamate/NMDA −2.99, voltage-gated Na⁺ −2.74, opioid receptors −2.71. The same polarity persisted at day 14 at lower amplitude (glial/immune +3.75, Toll-like +2.29, cytokines +2.48; K⁺/HCN −2.20, neurotrophin −2.03, GABAergic −1.64, opioid −1.64, Na⁺ −1.73).

**Table 1. Module-level signal (mean moderated *t*) across bone cancer pain stages.**

| Module | Day 1 | Day 7 | Day 14 |
|---|---:|---:|---:|
| Glial / immune | +1.46 | **+11.28** | +3.75 |
| Toll-like receptors | −0.07 | **+6.51** | +2.29 |
| Cytokine / chemokine | +0.95 | **+6.30** | +2.48 |
| Purinergic P2X | −0.50 | +2.47 | +0.78 |
| Prostanoid / COX | +1.39 | +1.70 | +1.74 |
| Cannabinoid | +0.15 | +0.56 | +0.09 |
| Injury / regeneration | +1.27 | −0.57 | +1.43 |
| Monoaminergic | −0.90 | −0.37 | −0.63 |
| TRP channels | −1.05 | −1.36 | −1.05 |
| Ca²⁺ channel α2δ | −1.72 | −1.92 | +0.60 |
| Opioid receptors | −0.74 | **−2.71** | −1.64 |
| Voltage-gated Na⁺ | −2.06 | **−2.74** | −1.73 |
| Glutamate / NMDA | −2.24 | −2.99 | −0.82 |
| GABAergic | −2.06 | **−3.65** | −1.64 |
| K⁺ / HCN channels | −2.48 | **−3.78** | −2.20 |
| Neurotrophin | −0.41 | **−3.92** | −2.03 |
| Neuropeptides | +0.88 | −1.83 | +0.12 |

### 3.4 In-silico perturbation: the pain state is reversed by *restoring*, not blocking, neuronal targets

Eight perturbation results reached q < 0.05, all at day 7 (Table 2). The direction is the informative part. Because the neuronal modules are *down*-regulated in the pain state (Table 1), virtual knockout pushed the profile further from control (negative reversal) while virtual **over-expression** pulled it back (positive reversal). For voltage-gated Na⁺ the reversal score was −0.056 for knockout and +0.056 for over-expression (q = 0.015); for K⁺/HCN channels, −0.064 and +0.064 (q = 0.015); neuropeptides −0.051/+0.051 (q = 0.043) and opioid receptors −0.051/+0.051 (q = 0.043).

Read naively, "knockout reverses pain" would be the expected direction, because most analgesics are antagonists or blockers. The data say the opposite for these targets at this stage: the ganglion has *lost* neuronal module expression, and the in-silico intervention that restores it is the one that normalises the transcriptional state.

The network model was not uniformly superior to its baselines. Ridge-linear baselines produced larger absolute scores (e.g. −0.86 to −0.96 for Na⁺, K⁺/HCN and opioid knockouts versus −0.05 to −0.06 for the network), and additive baselines were smaller (−0.02 to −0.04). The baselines and the network are not measuring the same quantity — the ridge model extrapolates a global shift, the network a localised one — but the comparison makes explicit that the network's contribution is specificity, not magnitude, and we do not claim the graph model outperforms a linear one on raw effect size.

**Table 2. In-silico perturbation results reaching q < 0.05 (day 7).**

| Module | Perturbation | Reversal (network) | Additive baseline | Ridge baseline | z vs null | q |
|---|---|---:|---:|---:|---:|---:|
| K⁺ / HCN channels | knockout | −0.064 | −0.038 | −0.914 | −2.46 | 0.015 |
| K⁺ / HCN channels | over-expression | +0.064 | +0.031 | +0.914 | +2.44 | 0.015 |
| Voltage-gated Na⁺ | knockout | −0.056 | −0.028 | −0.865 | −2.11 | 0.015 |
| Voltage-gated Na⁺ | over-expression | +0.056 | +0.028 | +0.865 | +2.11 | 0.015 |
| Neuropeptides | knockout | −0.051 | −0.033 | −0.354 | −2.16 | 0.043 |
| Neuropeptides | over-expression | +0.051 | +0.033 | +0.354 | +2.08 | 0.045 |
| Opioid receptors | knockout | −0.051 | −0.034 | −0.957 | −1.97 | 0.043 |
| Opioid receptors | over-expression | +0.051 | +0.034 | +0.957 | +1.93 | 0.043 |

### 3.5 The target hierarchy differs by pain aetiology

Restricted to the 13,601 genes shared between studies and the 1,882 genes common to all three co-expression networks, 36 of 102 cross-aetiology comparisons were significant at q < 0.05, with all reported bootstrap intervals excluding zero (Table 3).

The bone cancer pain versus paclitaxel contrast is the clearest. Knocking out glial/immune modules reversed bone cancer pain far more than it reversed paclitaxel neuropathy (Δ = +0.249, 95% CI 0.173 to 0.258; q = 0.0013), as did cytokine/chemokine (Δ = +0.139, CI 0.076 to 0.168) and Toll-like receptor knockout (Δ = +0.087, CI 0.045 to 0.098). The reverse held for neuronal modules: knocking out neuropeptides reversed paclitaxel neuropathy more than bone cancer pain (Δ = −0.207, CI −0.255 to −0.110), as did injury/regeneration (Δ = −0.127, CI −0.142 to −0.071), neurotrophin (Δ = −0.095, CI −0.100 to −0.062), voltage-gated Na⁺ (Δ = −0.072, CI −0.095 to −0.028) and K⁺/HCN (Δ = −0.062, CI −0.094 to −0.019).

The same immune-versus-neuronal split separated bone cancer pain from oxaliplatin neuropathy (cytokine Δ = +0.172, CI 0.105 to 0.192; glial/immune Δ = +0.167, CI 0.116 to 0.189; Toll-like Δ = +0.071, CI 0.032 to 0.087).

Paclitaxel and oxaliplatin were also distinguishable from each other, consistent with their differing clinical neurotoxicity: neuropeptides Δ = +0.146 (CI 0.055 to 0.196), injury/regeneration Δ = +0.121 (CI 0.070 to 0.132), neurotrophin Δ = +0.072 (CI 0.030 to 0.088) and GABAergic Δ = +0.039 (CI 0.008 to 0.055), all favouring paclitaxel.

**Table 3. Cross-aetiology differences in reversal score (virtual knockout), bootstrap 95% CI.**

| Contrast | Module | Δ (A − B) | 95% CI | q |
|---|---|---:|---|---:|
| Bone vs paclitaxel | Glial / immune | +0.249 | 0.173 to 0.258 | 0.0013 |
| Bone vs paclitaxel | Neuropeptides | −0.207 | −0.255 to −0.110 | 0.0013 |
| Bone vs paclitaxel | Cytokine / chemokine | +0.139 | 0.076 to 0.168 | 0.0013 |
| Bone vs paclitaxel | Injury / regeneration | −0.127 | −0.142 to −0.071 | 0.0013 |
| Bone vs paclitaxel | Neurotrophin | −0.095 | −0.100 to −0.062 | 0.0013 |
| Bone vs paclitaxel | Toll-like receptors | +0.087 | 0.045 to 0.098 | 0.0013 |
| Bone vs paclitaxel | Voltage-gated Na⁺ | −0.072 | −0.095 to −0.028 | 0.0013 |
| Bone vs paclitaxel | K⁺ / HCN channels | −0.062 | −0.094 to −0.019 | 0.0077 |
| Bone vs paclitaxel | GABAergic | −0.020 | −0.034 to −0.004 | 0.021 |
| Bone vs oxaliplatin | Cytokine / chemokine | +0.172 | 0.105 to 0.192 | 0.0013 |
| Bone vs oxaliplatin | Glial / immune | +0.167 | 0.116 to 0.189 | 0.0013 |
| Bone vs oxaliplatin | Toll-like receptors | +0.071 | 0.032 to 0.087 | 0.0013 |
| Paclitaxel vs oxaliplatin | Neuropeptides | +0.146 | 0.055 to 0.196 | 0.0013 |
| Paclitaxel vs oxaliplatin | Injury / regeneration | +0.121 | 0.070 to 0.132 | 0.0013 |
| Paclitaxel vs oxaliplatin | Neurotrophin | +0.072 | 0.030 to 0.088 | 0.0013 |
| Paclitaxel vs oxaliplatin | GABAergic | +0.039 | 0.008 to 0.055 | 0.012 |

---

## 4. Discussion

Three findings emerge. First, bone cancer pain does not produce a static DRG transcriptome: the response is absent at day 1, peaks at day 7 with more than three thousand differentially expressed genes, and has largely resolved by day 14. Second, at its peak the response is bipolar — an immune/glial/Toll-like axis is strongly up-regulated while neuronal excitability and transmission modules are down-regulated — and in-silico, the intervention that normalises the state is restoration of those neuronal modules rather than their blockade. Third, and most consequentially, the hierarchy of reversible targets is aetiology-specific: immune/glial modules dominate in tumour-driven pain, neuropeptide and neuronal-injury modules dominate in paclitaxel-driven pain, and the two chemotherapy agents differ from each other as well.

The immune-side result is consistent with a large literature on tumour-driven pain, in which macrophage infiltration, glial activation and cytokine signalling within the ganglion and dorsal horn maintain sensitisation [3,4,13]. The neuronal-side result is more counter-intuitive. A ganglion that has *lost* expression across Na⁺, K⁺/HCN, GABAergic, glutamatergic and opioid receptor modules is not simply "more excitable"; the transcriptional picture is closer to a loss of normal sensory neuron identity, which is what injury and neuroinflammation produce [14,15]. If that reading is right, it carries a therapeutic implication that runs against the grain of most analgesic development: for tumour-driven pain at this stage, agents that *restore* or *enhance* neuronal module function may be at least as relevant as agents that block it, and the opioid receptor findings fit the clinical observation that opioid responsiveness in bone pain is variable rather than uniformly good.

The aetiology comparison is the result with the most direct clinical reading. Current guidance treats CIPN as largely refractory, with duloxetine the only agent carrying a firm recommendation [16], while bone pain is treated with anti-inflammatories, bisphosphonates, radiotherapy and opioids [3]. Our analysis suggests these are not merely different intensities of one problem but different problems: an immune axis that in-silico dominates bone pain is largely absent from the paclitaxel signature, and a neuropeptide/injury axis that dominates paclitaxel neuropathy is not the main driver in bone pain. That is an argument for mechanistic stratification rather than intensity stratification, and it is testable — the module-level signatures and reversal rankings in Tables 1–3 are quantitative predictions that pharmacological studies can confirm or refute.

Several results temper enthusiasm. The reversal scores are small in absolute terms (|0.02| to |0.25|). This is expected — a single module's perturbation is a small vector in a 3,000-dimensional space, and the cosine with the whole pain axis is necessarily modest — but it means the ranking should be read as ordinal, not as an effect size, and it is why every result is reported against a permutation null and two baselines rather than as a bare number. The ridge-linear baseline produces larger absolute scores than the network model, which we report plainly: the graph contributes specificity, not magnitude, and we make no claim of outperforming a linear model.

### Limitations

Four limitations constrain interpretation. (i) **Sample size.** Four animals per group is standard for ganglion transcriptomics but is small; this is why we used moderated statistics, module-level aggregation and permutation/bootstrap inference rather than relying on single-gene tests, and why day-1 results should be read as "no detectable change" rather than "no change". (ii) **Reversal is transcriptional, not analgesic.** A reversal score says the perturbed profile moves toward the control profile. It does not say an animal would feel less pain, and the mapping from transcriptional reversal to behavioural analgesia is unvalidated. (iii) **Network inference from observational data.** The co-expression graph is correlational; diffusion over it models propagation of a perturbation, not a causal mechanism. Our "knockout" and "over-expression" are computational operations on a graph, not simulations of specific drugs, and no chemical structure, dose or pharmacokinetics enters the model. (iv) **Cross-study comparison.** Although species, tissue and platform are matched, the two studies were done in different laboratories with different protocols and timing, and batch effects cannot be fully excluded — the within-study bootstrap accounts for sampling variability but not for between-study systematic differences. We also excluded one dataset (GSE249443) from biological inference because its apparent signal was a muscle-contamination artefact; whether other published cancer-pain ganglion signatures are similarly affected is worth systematic checking.

---

## 5. Conclusions

In rat DRG transcriptomes matched for species, tissue and platform, the hierarchy of analgesic targets that in-silico reverses the pain state is aetiology-specific. Tumour-driven bone pain is dominated by an immune, glial and Toll-like axis, whereas chemotherapy-driven neuropathy is dominated by a neuronal injury and neuropeptide axis, and paclitaxel and oxaliplatin are distinguishable from each other. These hierarchies are quantitative and falsifiable, and they support mechanistically stratified rather than intensity-stratified analgesic selection in cancer pain.

---

## Declarations

**Ethics approval and consent to participate.** Not applicable; only publicly available data were used.

**Consent for publication.** Not applicable.

**Availability of data and materials.** All analysed datasets are publicly available from GEO (GSE149648, GSE160543, GSE249443). Analysis code, result tables and figures are available at [repository URL].

**Competing interests.** None declared.

**Funding.** None.

**Authors' contributions.** YY conceived the study, wrote the analysis code, performed the analyses and wrote the manuscript.

---

## References

1. van den Beuken-van Everdingen MHJ, Hochstenbach LMJ, Joosten EAJ, Tjan-Heijnen VCG, Janssen DJA. Update on prevalence of pain in patients with cancer: systematic review and meta-analysis. *J Pain Symptom Manage*. 2016;51(6):1070-1090.e9.
2. Snijders RAH, Brom L, Theunissen M, van den Beuken-van Everdingen MHJ. Update on prevalence of pain in patients with cancer: a systematic review and meta-analysis (2016-2022). *J Pain Symptom Manage*. 2023;66(4):317-335.
3. Mantyh PW. Bone cancer pain: causes, consequences, and therapeutic opportunities. *Pain*. 2013;154 Suppl 1:S54-S62.
4. Falk S, Dickenson AH. Plasticity of central and peripheral pain processing in bone cancer pain. *Eur J Cancer*. 2016;66:141-148.
5. Staff NP, Grisold A, Grisold W, Windebank AJ. Chemotherapy-induced peripheral neuropathy: a current review. *Ann Neurol*. 2017;81(6):772-781.
6. Flatters SJL, Dougherty PM, Colvin LA. Clinical and preclinical perspectives on chemotherapy-induced peripheral neuropathy (CIPN): a narrative review. *Br J Anaesth*. 2017;119(4):737-749.
7. Lotfollahi M, Wolf FA, Theis FJ. scGen predicts single-cell perturbation responses. *Nat Methods*. 2019;16(8):715-721.
8. Roohani Y, Huang K, Leskovec J. Predicting transcriptional outcomes of novel multigene perturbations with GEARS. *Nat Biotechnol*. 2024;42(6):927-935.
9. Hetzel L, Boehm S, Kilbertus N, et al. Predicting cellular responses to complex perturbations in high-throughput screens. *Adv Neural Inf Process Syst*. 2022;35:31919-31933.
10. [2025 *Nature Methods* perturbation-model benchmarking study — TO BE VERIFIED AND INSERTED BEFORE SUBMISSION]
11. Smyth GK. Linear models and empirical bayes methods for assessing differential expression in microarray experiments. *Stat Appl Genet Mol Biol*. 2004;3:Article3.
12. Phipson B, Lee S, Majewski IJ, Alexander WS, Smyth GK. Robust hyperparameter estimation protects against hypervariable genes and improves power to detect differential expression. *Ann Appl Stat*. 2016;10(2):946-963.
13. Ji RR, Chamessian A, Zhang YQ. Pain regulation by non-neuronal cells and inflammation. *Science*. 2016;354(6312):572-577.
14. Costigan M, Scholz J, Woolf CJ. Neuropathic pain: a maladaptive response of the nervous system to damage. *Annu Rev Neurosci*. 2009;32:1-32.
15. Renthal W, et al. Transcriptional reprogramming of distinct peripheral sensory neuron subtypes after axonal injury. *Neuron*. 2020;108(1):128-144.e9.
16. Hershman DL, et al. Prevention and management of chemotherapy-induced peripheral neuropathy in survivors of adult cancers: ASCO guideline update. *J Clin Oncol*. 2014;32(18):1941-1967.

---

*Manuscript version: v0.2 — all results read directly from `results/tables/` artefacts; reference 10 pending verification.*
