# Divergent analgesic-target hierarchies in cancer pain of different aetiology: an in-silico perturbation analysis of rat dorsal root ganglion transcriptomes

**Running head:** In-silico perturbation of cancer-pain DRG transcriptomes

**Author:** Yongxin Yang¹

¹ Department of Anesthesiology, The Second Affiliated Hospital of Fujian University of Traditional Chinese Medicine, Fuzhou, Fujian 350003, China

**ORCID:** 0009-0004-9698-6552

**Correspondence:** Yongxin Yang, Department of Anesthesiology, The Second Affiliated Hospital of Fujian University of Traditional Chinese Medicine, Fuzhou, Fujian 350003, China. Email: 960856791@qq.com

---

## Abstract

**Background.** Cancer pain is treated as a single clinical entity, yet its dominant drivers differ by aetiology: bone-metastatic pain is driven by tumour–bone destruction and the inflammatory milieu it creates, whereas chemotherapy-induced peripheral neuropathy (CIPN) reflects direct neurotoxic injury to sensory neurons. Whether the two should be treated with the same analgesic hierarchy is unresolved, and no study has asked this question at the level of the sensory ganglion transcriptome.

**Methods.** We analysed rat dorsal root ganglion (DRG) RNA-seq from two studies generated on the **same sequencing platform (GPL24782)**, so that species, tissue and platform are held constant and only the aetiology of pain differs: tibial Walker-256 bone cancer pain (GSE149648; sham and days 1, 7 and 14 after tumour cell injection; n = 4 per group) and paclitaxel- or oxaliplatin-induced neuropathy (GSE160543; vehicle, paclitaxel, oxaliplatin; n = 4 per group). After confounding-controlled normalisation, differential expression used an empirical-Bayes moderated *t* statistic with an intensity-dependent prior. We then built an aetiology-specific co-expression graph and performed in-silico knockout and overexpression of 17 analgesic target modules by diffusing a seed perturbation over the graph, scoring each perturbation by how far it moved the pain state back toward its own control (reversal score). Empirical null distributions came from size-matched random gene sets; cross-aetiology differences were tested by within-study bootstrap.

**Results.** [TO BE FILLED FROM results/tables]

**Conclusions.** [TO BE FILLED]

**Keywords:** cancer pain; dorsal root ganglion; in-silico perturbation; virtual knockout; analgesic target; chemotherapy-induced peripheral neuropathy; bone metastasis

---

## 1. Introduction

Cancer pain affects the majority of patients with advanced disease and remains inadequately controlled in a substantial minority [1,2]. The WHO analgesic ladder, still the dominant organising frame, is explicitly **aetiology-agnostic**: it escalates by pain intensity rather than by mechanism. In practice the two largest cancer-pain populations differ mechanistically. Bone-metastatic pain arises from tumour-induced osteolysis, acidosis and a local inflammatory milieu that sensitises primary afferents [3,4]; chemotherapy-induced peripheral neuropathy arises from direct axonal and neuronal injury by the cytotoxic agent itself [5,6]. It is therefore reasonable to ask whether the analgesic targets that matter most in one also matter most in the other. That question is difficult to address experimentally — it would require paired ganglion tissue from patients with different pain aetiologies — but it is addressable computationally.

Two recent developments make this feasible. First, curated public ganglion transcriptomes now exist for several cancer-pain models, in some cases generated on the same platform, which removes the species/tissue/platform confounding that has historically made cross-study pain comparisons uninterpretable. Second, the last three years have produced a family of *in-silico perturbation* methods — latent-space vector arithmetic (scGen) [7], graph-based genetic-perturbation models (GEARS) [8], and chemical-conditioned autoencoders (chemCPA) [9] — that predict how a cell's transcriptome reorganises after a gene or drug perturbation. These methods have moved perturbation biology from "what changed" toward "what would happen if".

That promise comes with an evidentiary obligation. A 2025 benchmarking study in *Nature Methods* reported that several single-cell foundation and perturbation models fail to beat simple linear or additive baselines once the evaluation includes strong-effect perturbations and distribution shift [10]. Any manuscript that deploys a perturbation model without head-to-head baseline reporting is therefore making an unsupported claim. We designed this study around that constraint: every perturbation result is co-reported against additive and ridge-linear baselines, and against an empirical null built from size-matched random gene sets.

A second, more prosaic obligation concerns confounding. DRG dissection unavoidably carries paraspinal skeletal muscle, and the amount carried is not constant between animals. Because skeletal muscle transcripts are abundant and highly expressed, even modest differences in contamination can dominate a small-n differential expression analysis and produce a "pain signature" that is in fact a muscle signature. We treat this as a first-class analytical problem rather than a footnote.

We asked three questions. (i) Does the DRG transcriptome of bone cancer pain reorganise over time, and if so do different analgesic target modules dominate at different stages? (ii) In-silico, which target modules most effectively reverse the pain state at each stage? (iii) Is the resulting target **hierarchy** the same in bone cancer pain as in chemotherapy-induced neuropathy?

---

## 2. Methods

### 2.1 Data sources

All data are publicly available from the NCBI Gene Expression Omnibus. No new animal work was performed, so no ethics approval was required.

**Discovery set — bone cancer pain, GSE149648.** Walker 256 mammary carcinoma cells were injected into the rat tibia (tumour cell injection, TCI); ipsilateral DRG were collected on day 1, day 7 and day 14, with sham (PBS injection) controls at day 14. Four biological replicates per group, 16 samples total, raw gene-level counts, platform GPL24782.

**Comparison set — chemotherapy-induced neuropathy, GSE160543.** Rats received intraperitoneal paclitaxel (2.0 mg/kg) or oxaliplatin (2.4 mg/kg), or vehicle; DRG were collected after pain behaviour was established. Four animals per group, 12 samples total, raw gene-level counts, **platform GPL24782**.

**Replication set — bone cancer pain, GSE249443.** An independent rat tibial bone cancer pain dataset (MRMT-1 cells) with sham controls, six samples total, FPKM values.

The decisive design feature is that the discovery and comparison sets share species (*Rattus norvegicus*), tissue (DRG) and sequencing platform (GPL24782). Only the aetiology of pain differs. We deliberately did **not** compare across species or across tissue types, because such comparisons confound aetiology with everything else.

### 2.2 Normalisation and quality control

Per-sample count files were assembled into a gene × sample matrix. Size factors were estimated by the median-of-ratios method, and counts converted to log2 counts-per-million, log2(CPM + 1). Genes were retained if total counts across samples were ≥ 10 and the gene was detected in at least four samples.

Sample-level quality was assessed by pairwise Pearson correlation and by a contamination score (Section 2.3).

### 2.3 Skeletal-muscle contamination: measurement and control

DRG dissection carries adjacent paraspinal skeletal muscle. To quantify this, we defined a muscle module from canonical sarcomere and muscle-metabolism genes (*Acta1, Ckm, Myl1, Myh1, Myh2, Tnni2, Tnnc2, Tnnt3, Des, Mybph, Actn3, Mybpc1, Myoz1, Tnnc1, Mylpf, Myh4, Atp2a1, Casq1, Mb, Actc1*), z-scored each gene across samples and averaged within sample to give a per-sample muscle score.

Two models were then fitted for every contrast: an unadjusted model and a model including the (mean-centred) muscle score as a covariate. Genes whose apparent differential signal disappears after adjustment are reported as **contamination-driven**, not as pain biology. This is a conservative rule and it is applied uniformly; in the replication set it is decisive (Section 3.1).

### 2.4 Differential expression

Differential expression used a linear model per gene with an empirical-Bayes moderated *t* statistic. Briefly, residual variances were shrunk toward a prior estimated by method of moments across genes [11], and — because RNA-seq variance depends strongly on mean intensity — the prior was made intensity-dependent by fitting a lowess curve of log residual variance on mean log-expression and moderating the *scaled* variances (the limma-trend procedure [12]). *P* values were converted to *q* values by Benjamini–Hochberg.

This choice matters at n = 4 per group. A naive Welch *t*-test at this sample size has only 2–6 residual degrees of freedom and is poorly calibrated at both tails; moderated variance borrows strength across genes and is the standard remedy.

### 2.5 Module-level signal

Single-gene tests are under-powered at this sample size, so we additionally tested **analgesic target modules** — mechanistically coherent gene sets corresponding to druggable analgesic mechanisms (opioid receptors; voltage-gated sodium channels; Ca²⁺ channel α2δ subunits; TRP channels; P2X purinoceptors; cyclo-oxygenase/prostanoid receptors; glutamate/NMDA; GABAergic; neuropeptides; cannabinoid; neurotrophin; cytokines/chemokines; Toll-like receptors; glial/immune; injury-regeneration; K⁺/HCN channels; monoaminergic). For each module and contrast, module-level signal was the mean moderated *t* across member genes, tested competitively against all genes by Mann–Whitney *U*.

### 2.6 Co-expression network

For each aetiology separately, the highest-variance genes (3,000 by default) were used to build a sparse *k*-nearest-neighbour correlation graph (*k* = 12): for each gene, the 12 strongest absolute correlations were retained, the matrix symmetrised as (W + Wᵀ)/2. Networks were built **per aetiology** because co-expression structure is itself aetiology-dependent; pooling would erase exactly the difference we wish to detect.

### 2.7 In-silico perturbation

A perturbation of module *M* with direction *s* ∈ {−1, +1} (knockout, overexpression) was encoded as a seed vector with *s* at the indices of genes in *M* and zero elsewhere. The seed was diffused over the graph by

δ_{t+1} = (1 − α)·seed + α·S·δ_t,  S = D^{−1/2} W D^{−1/2},  α = 0.55,

iterated to convergence (80 iterations; α⁸⁰ ≈ 10⁻²¹). Because the map is linear, the limiting operator P with δ_∞ = P·seed was computed **once** per aetiology, after which every real and permuted perturbation reduces to a single matrix–vector product. This is what makes a 1,000-fold permutation test computationally tractable on CPU.

The perturbation was scored by the **reversal score**:

reversal = cosine( P·seed , pain_axis ),  pain_axis = mean(control) − mean(case),

with both vectors mean-centred before the cosine so that the score reflects agreement in *pattern* rather than a global shift. A positive score means the in-silico perturbation moves the pain-state transcriptome back toward the control state. Because pain_axis is defined within each study, reversal scores are dimensionless and comparable across aetiologies.

### 2.8 Empirical null and baseline comparison

For each module, perturbation and aetiology, 1,000 size-matched random gene sets were drawn and their reversal scores computed, giving an empirical null; the one-sided empirical *p* value is the tail fraction with +1 continuity correction. *q* values were computed within aetiology.

Every perturbation result is co-reported against two baselines:

- **Additive (no propagation):** the seed itself, scored directly. This isolates whatever the network contributes.
- **Ridge-linear:** each gene regressed on the module's genes across samples by ridge regression, then the response to shifting the module by the same amount is predicted — no graph, no propagation.

This head-to-head is not optional decoration. Benchmarking has shown perturbation models frequently fail to outperform such baselines [10]; if the network adds nothing, the additive result will say so.

### 2.9 Cross-aetiology comparison

Because the two studies are independent, labels cannot be permuted between them. We therefore bootstrapped **within** each study: cases and controls were resampled with replacement, both pain axes recomputed, and both reversal scores recomputed; the difference Δ = reversal(A) − reversal(B) was accumulated over 2,000 resamples to give a percentile 95% interval and a two-sided bootstrap *p* value (the smaller tail fraction, doubled, floored at 1/n_boot). Note that bootstrap resampling can draw duplicate samples; the resulting degenerate axes are retained rather than filtered, which makes the interval conservative.

### 2.10 Software and reproducibility

Python 3.13 with NumPy, pandas, SciPy, statsmodels (lowess) and matplotlib. All random procedures use a fixed seed (20261009). Analysis scripts, result tables and figures are deposited at **[repository URL to be inserted on release]** with a checksum MANIFEST.

---

## 3. Results

[TO BE FILLED FROM results/tables — every number below must be read from an artefact file]

### 3.1 Data quality and the muscle-contamination problem

### 3.2 Bone cancer pain reorganises the DRG transcriptome over time

### 3.3 Analgesic target modules are engaged at different stages

### 3.4 In-silico perturbation: a stage-dependent target hierarchy

### 3.5 The target hierarchy differs by pain aetiology

### 3.6 Independent replication

### 3.7 Head-to-head against baselines

---

## 4. Discussion

[TO BE FILLED]

---

## 5. Limitations

[TO BE FILLED]

---

## 6. Conclusions

[TO BE FILLED]

---

## Declarations

**Ethics approval and consent to participate.** Not applicable; only publicly available data were used.

**Consent for publication.** Not applicable.

**Availability of data and materials.** All analysed datasets are publicly available from GEO (GSE149648, GSE160543, GSE249443). Analysis code and result tables are available at [repository URL].

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
10. [2025 Nature Methods perturbation-benchmarking study — TO BE VERIFIED AND INSERTED]
11. Smyth GK. Linear models and empirical bayes methods for assessing differential expression in microarray experiments. *Stat Appl Genet Mol Biol*. 2004;3:Article3.
12. Phipson B, Lee S, Majewski IJ, Alexander WS, Smyth GK. Robust hyperparameter estimation protects against hypervariable genes and improves power to detect differential expression. *Ann Appl Stat*. 2016;10(2):946-963.

---

*Manuscript version: v0.1 (skeleton; Results/Discussion pending artefact-derived numbers)*
