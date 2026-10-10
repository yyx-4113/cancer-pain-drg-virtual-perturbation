# Divergent dorsal root ganglion programs in bone cancer pain and chemotherapy-induced peripheral neuropathy: a calibrated in-silico perturbation analysis

**Running head:** Calibrated in-silico perturbation of cancer-pain DRG transcriptomes

**Author:** Yongxin Yang¹

¹ Department of Anesthesiology, The Second Affiliated Hospital of Fujian University of Traditional Chinese Medicine, Fuzhou, Fujian 350003, China

**ORCID:** 0009-0004-9698-6552

**Correspondence:** Yongxin Yang, Department of Anesthesiology, The Second Affiliated Hospital of Fujian University of Traditional Chinese Medicine, Fuzhou, Fujian 350003, China. Email: 960856791@qq.com

---

## Abstract

**Background.** Bone cancer pain and chemotherapy-induced peripheral neuropathy (CIPN) are both treated as "cancer pain", but their drivers are assumed to differ—tumour-associated inflammation versus direct neurotoxicity. Whether the dorsal root ganglion (DRG) transcriptional programmes actually diverge, and whether in-silico perturbation can rank analgesic targets from bulk transcriptomes, are separate questions that are often conflated.

**Methods.** We analysed rat DRG RNA-seq from two studies generated on the same platform (GPL24782), so species, tissue and platform are held constant: tibial Walker-256 bone cancer pain (GSE149648; sham and days 1, 7, 14; n = 4 per group) and paclitaxel- or oxaliplatin-induced neuropathy (GSE160543; vehicle, paclitaxel, oxaliplatin; n = 4 per group). Differential expression used an empirical-Bayes moderated *t* statistic with an intensity-dependent (limma-trend) prior and a skeletal-muscle contamination covariate. We performed in-silico knockout and overexpression of 17 analgesic target modules by diffusing a seed perturbation over a co-expression graph, scoring each by its cosine with the pain-to-control axis (reversal score). The perturbation layer was calibrated against five negative-control modules, a degree-matched permutation null, an oracle anchor, and additive and ridge-linear baselines. Cross-aetiology differences were tested by exhaustive label permutation (C(8,4)² = 4,900 null pairs; exact p, no resampling error).

**Results.** Bone cancer pain reorganised the DRG transcriptome in a stage-dependent manner (0, 4,600 and 235 genes at q < 0.05 on days 1, 7 and 14). In the Walker-256 tibia model, day 7 was dominated by a type-I interferon programme (module mean moderated *t* = +24.9, q = 6.8 × 10⁻⁵) with myeloid infiltration (Cd68 +2.86, Aif1 +2.18 log₂) and satellite gliosis (Gfap +1.33). Neuronal target modules were down-regulated by essentially the same amount as the pan-neuronal (−5.20) and myelin (−4.10) negative controls (Na⁺ −4.35, K⁺/HCN −5.64, GABAergic −5.71, opioid −3.76), indicating a compartment effect rather than five independent mechanisms. Cross-aetiology, the divergence was large and direction-specific: myeloid markers were up in bone cancer pain and down in paclitaxel neuropathy (Cd68 Δ = +3.42, Aif1 Δ = +3.08, both q = 0.012), whereas the injury-associated neuropeptide galanin moved in the opposite direction (Δ = −2.59, q = 0.012). The in-silico perturbation layer did not survive calibration. Knockout and overexpression are exact mirror images (identical exact-permutation p to eight decimals), so the 264 contrast-level tests reduce to 132 independent hypotheses per network scope; even with the corrected count, no module's reversal Δ reached q < 0.05 (minimum q = 0.165). Against the degree-matched random-gene-set null, 32/132 tests were "significant" but 14 were negative-control modules—including the erythrocyte set with the largest standardised effect—so the signal is non-specific. An independent single-cell dataset (B16F10 cancer-nerve-injury DRG, n = 4 per group) showed a neuronal-injury programme (Sprr1a, Atf3) rather than the interferon programme, confirming the latter is model-specific and not a universal cancer-pain signature. Network propagation was almost perfectly redundant with an additive model (Spearman ρ = 0.991, n = 132).

**Conclusions.** The two aetiologies produce qualitatively different DRG programmes—myeloid/interferon-driven versus neuronal-injury-driven—and this is demonstrable at the level of individual canonical markers with exact permutation tests. An independent single-cell model shows the interferon programme is model-specific rather than a universal cancer-pain signature, so cross-study generalisation requires the negative controls and exact tests used here. The same data do not support module-specific analgesic target ranking by in-silico perturbation: at n = 4 per group the reversal score is dominated by cellular-composition differences between datasets. In-silico perturbation should be reported with negative controls, degree-matched nulls and an additive baseline, or not at all.

**Keywords:** cancer pain; dorsal root ganglion; in-silico perturbation; virtual knockout; chemotherapy-induced peripheral neuropathy; type-I interferon; satellite glia; negative control

---

## 1. Introduction

Cancer pain is treated as a single clinical entity and managed with an escalation ladder that is largely intensity-driven rather than mechanism-driven [1,2]. Yet the biology at the first-order sensory neuron is unlikely to be the same when the driver is a tumour growing inside bone and when it is a neurotoxic chemotherapeutic agent. Bone cancer pain models consistently show osteoclast activation, tumour-derived mediators, and a pronounced innate immune response in the DRG [3,4]; chemotherapy-induced peripheral neuropathy (CIPN) models show axonal degeneration, altered ion-channel trafficking and a neuronal injury/regeneration transcriptional programme [5,6].

Two developments make a direct comparison feasible. First, multiple DRG RNA-seq datasets from different pain aetiologies are now public, and some share species, tissue and sequencing platform, which removes the three largest sources of technical heterogeneity in a cross-study comparison. Second, "virtual experiment" methods—in-silico knockout and overexpression, in which a perturbation is seeded on a gene set and propagated over a data-derived network—have become a standard way to move from differential expression towards a ranked, mechanistic hypothesis about which targets to modulate [7,8].

The enthusiasm for these methods has run ahead of their validation. A 2025 benchmark in *Nature Methods* reported that deep-learning perturbation models do not yet outperform simple linear or additive baselines [9]. That result was established in cell-line perturbation assays; it has not, to our knowledge, been extended to the case that concerns translational pain research, namely perturbation of predefined target modules in small-n bulk tissue transcriptomes. Two features of that setting make it fragile. Bulk tissue measurements are compositional: an influx of macrophages mechanically lowers the transcript share of every neuronal gene, so a "neuronal module is down-regulated" finding can be a statement about cell-type proportions rather than about gene regulation. And network propagation is degree-dependent: a seed on a hub gene diffuses further than a seed on a peripheral gene, so a module's score reflects its topological position as well as its biology.

We therefore set out to do two things and to keep them separate. The first is a biological question: do the DRG transcriptional programmes of bone cancer pain and CIPN diverge, and if so, along which axes? The second is a methodological question: after calibration against negative controls, degree-matched nulls, an oracle anchor and a simple additive baseline, does in-silico perturbation add anything to the answer?

## 2. Methods

### 2.1 Datasets

**GSE149648** — rat DRG, Walker-256 mammary carcinoma cells implanted into the tibia (TCI); sham (PBS) and days 1, 7 and 14 after implantation; n = 4 per group; FPKM quantified on platform GPL24782 [10]. Files were downloaded from GEO as per-sample supplementary tables.

**GSE160543** — rat DRG, paclitaxel, oxaliplatin or vehicle; n = 4 per group; FPKM on the same platform GPL24782 [11].

Both series provide bulk DRG, not single-cell data, and both are rat. Species, tissue and platform are therefore matched; laboratory, animal batch and dissection protocol are not, and this is treated as a limitation rather than adjusted away.

### 2.2 Parsing and normalisation

Per-sample files were parsed by locating the FPKM column from the header with a right-to-left numeric fallback, because the two series ship different column layouts and GSE160543 is not internally uniform (six columns for most members, thirteen for some because of trailing empty tab-delimited fields). Reading a fixed column index selects the Entrez identifier in both series and yields a matrix of identifiers—a silent failure that produces zero fold change everywhere while reporting success. Each parsed matrix is therefore asserted to be non-integer-valued (fraction of integer-valued entries ≤ 0.95) before use.

FPKM is nominally length- and depth-normalised, but per-sample FPKM totals are not equal across groups (a 6–9% global offset in GSE149648), so a median-of-expressed-genes scaling factor was applied before the log₂ transform. Genes were retained if FPKM ≥ 1 in at least four samples and mean FPKM ≥ 0.5: 13,893 genes in GSE149648 and 13,601 genes shared between the two studies.

### 2.3 Confounding control: skeletal-muscle contamination

DRG dissection systematically carries paraspinal skeletal muscle. A muscle score was computed as the mean z-score of 20 sarcomere and muscle-specific genes (Acta1, Ckm, Myl1, Myh1, Myh2, Tnni2, Tnnc2, Tnnt3, Des, Mybph, Actn3, Mybpc1, Myoz1, Tnnc1, Mylpf, Myh4, Atp2a1, Casq1, Mb, Actc1), and every differential-expression model was fitted both unadjusted and with the muscle score as a covariate. Contamination was also quantified absolutely as the share of total FPKM contributed by those genes, and samples above a substantive threshold of 1% were flagged and carried through a leave-one-out sensitivity analysis. A MAD-based robust z was computed on log₁₀ share rather than on the raw share, because the raw distribution is extremely right-skewed and its MAD collapses towards zero, which inflates the z statistic by orders of magnitude and flags the cleanest samples as outliers.

### 2.4 Differential expression

Moderated *t* statistics were computed with an empirical-Bayes variance shrinkage (Smyth 2004) and a limma-trend intensity-dependent prior variance obtained by lowess regression of log residual variance on mean intensity (span 0.4) [12]. At n = 4 per group a naive *t*-test has 2–6 residual degrees of freedom and is badly calibrated at both tails; moderation borrows strength across the ~14,000 genes. *q* values are Benjamini–Hochberg [13].

### 2.5 Target and control modules

Seventeen analgesic target modules were defined a priori from pharmacology, not from the data: opioid receptors, voltage-gated Na⁺, Ca²⁺ channel α2δ, TRP channels, purinergic P2X, prostanoid/COX, glutamate/NMDA, GABAergic, neuropeptides, cannabinoid, neurotrophin, cytokine/chemokine, Toll-like receptors, glial/immune, injury/regeneration, K⁺/HCN channels, monoaminergic.

Five negative-control modules were defined to test specificity: pan-neuronal structural genes, myelin/Schwann genes, housekeeping genes, interferon-stimulated genes and erythrocyte genes. The first four are chosen so that a compartment or composition effect would be visible; the erythrocyte set is a pure technical control, since red-cell content varies with perfusion and dissection and carries no pain biology.

### 2.6 Co-expression graph and in-silico perturbation

For each contrast a k-NN co-expression graph was built on the highest-variance genes (k = 12, symmetrised). Two properties of this construction are load-bearing and were each the site of a bug that changed the results:

*No self-loops.* The diagonal of the correlation matrix is set to zero before neighbour selection. Setting it to −1 instead, on the reasoning that this "excludes self", has the opposite effect: selection is on |C|, |−1| = 1 is the row maximum, so the self-pair is always selected and the value written is the signed −1, producing a negative self-loop at every node.

*Target genes are force-included.* Selecting nodes by variance alone drops most analgesic target genes, because ion channels and opioid receptors are moderately expressed with unremarkable across-sample variance. Without force-inclusion the GABAergic module is represented by one node and the P2X module by one node, and modules are perturbed on a graph that does not contain them. Every gene of every module is therefore forced onto the graph.

A knockout or overexpression is a seed vector in gene space, diffused as δ_{t+1} = (1−α)·seed + α·S·δ_t with S = D^{−1/2} W D^{−1/2}, α = 0.55. Because the map is linear, the limiting operator P with δ_∞ = P·seed was precomputed once per contrast; |α·S| ≤ α < 1 guarantees convergence, and the closed form P = (1−α)(I − αS)⁻¹ agrees with 80 iterations to 2.2 × 10⁻¹⁶.

The reversal score is cosine(P·seed, pain_axis), with pain_axis = mean(control) − mean(case), both mean-centred. Positive means the perturbation moves the pain state toward its own control. Three calibrations are reported alongside it:

- **Oracle anchor.** The best achievable reversal for a seed of the same size, obtained by perturbing the *k* genes already furthest along the pain axis with the correct sign. A cosine of 0.05 is uninterpretable without this reference; expressed as a percentage of the oracle it is.
- **Degree-matched permutation null.** Random gene sets are drawn so that each gene is replaced by a random gene from the same decile of network degree. Uniform random sets are not a valid null for a diffusion score, because hub seeds propagate further and a hub-rich module would look significant for a purely topological reason.
- **Baselines.** An additive model (no propagation) and a ridge-linear model that regresses every gene on the target genes across samples, fitted leave-one-sample-out so it cannot fit the pain axis in-sample.

### 2.7 Cross-aetiology comparison

Two endpoints were used, and they are not interchangeable.

*Module mean log₂ fold change.* For each module, the mean logFC (case versus control) within each study; the cross-aetiology difference Δ = logFC(A) − logFC(B). Because logFC does not involve the empirical-Bayes prior, it can be recomputed for every one of the C(8,4) = 70 case/control assignments of eight samples, so the null over the 70 × 70 = 4,900 (A, B) split pairs is enumerated exhaustively. This gives an exact permutation p with no Monte-Carlo error and no seed dependence, and a hard floor of 1/4,901 rather than an arbitrary 1/2,000. Where the two contrasts share their control samples (both bone-cancer days versus sham; paclitaxel and oxaliplatin versus the same vehicle), the split index is paired so that both sides use the same control subset. The same test was applied to individual canonical marker genes.

*Perturbation reversal difference.* Δ = reversal(A) − reversal(B), tested with the same exhaustive permutation scheme. Because a permutation test gives no interval, a within-study bootstrap (2,000 resamples) provides the 95% interval for Δ; cases and controls are resampled separately and, where the control group is shared, jointly. A second null asks a different question—is *this module's* Δ larger than that of arbitrary gene sets of the same size and degree (500 degree-matched draws)?

Networks were built per contrast (each on its own eight samples), because co-expression structure is itself treatment-dependent and a shared network removes part of the difference being tested. A pooled-network variant is reported as a sensitivity analysis.

### 2.8 Independent single-cell validation

To test whether the Walker-256 interferon programme generalises beyond that single model, an independent single-cell RNA-seq dataset of mouse DRG was analysed: GSE292522 (10x Genomics 3′, GPL34328; four DRG from B16F10 cancer-nerve-injury animals and four littermate controls, 55,967 cells after quality filtering). The eight per-sample 10x matrices were read directly and collapsed to a gene × sample pseudobulk (sum of counts per gene), converted to log₂(CPM + 1), and compared between the four tumour and four control samples with the same limma-trend model used in the discovery pipeline. Cell-type composition was estimated per sample from canonical marker sets (neurons, satellite glia, Schwann, myeloid, T cells, endothelial, fibroblast, erythrocyte) by assigning each cell to the marker set with the highest mean log-expression; the myeloid fraction was compared across groups with a Mann–Whitney U test. Convergent validation correlated the B16F10-versus-control pseudobulk log₂ fold change against the Walker-256 day-1/day-7/day-14 bulk logFC over shared genes (Spearman ρ).

### 2.9 Software and reproducibility

Python 3.13 with numpy 2.4.6, scipy 1.18.1, pandas 2.3.3, statsmodels 0.15.0 and matplotlib 3.11.1 (`requirements.txt`). statsmodels is imported at the top level and its absence raises a fatal error rather than silently degrading limma-trend to a constant prior. All statistical and perturbation primitives live in one module so that the discovery and cross-aetiology scripts cannot drift. A test suite of 14 groups asserts the properties that past bugs violated—no self-loops, knockout and overexpression as exact mirror images, identifier-matrix rejection, D14 not collapsing into D1, degree-matched draws staying in their bin, exhaustive splits containing the true assignment, and per-column normalisation of the batched cosine. A checker script asserts every number reported here against the deposited CSV and JSON artefacts.

## 3. Results

### 3.1 Data quality and confounding

All 16 samples of GSE149648 and all 12 of GSE160543 were parsed and verified. Sample-to-sample correlation ranged 0.913–0.994 in the discovery set; PCA gave 46.8%, 17.4%, 7.2% and 4.8% for PC1–PC4. The muscle score did not differ by group (ANOVA F = 0.637, p = 0.605). One sample, GSM4506604 (day 7), carried 6.26% of its total FPKM in sarcomere genes, an order of magnitude above every other sample (next highest 0.16%); it is the only sample above the 1% threshold and was carried through a leave-one-out sensitivity analysis.

### 3.2 Bone cancer pain reorganises the DRG transcriptome with a sharp day-7 peak

After adjustment for muscle contamination, 0, 4,600 and 235 genes reached q < 0.05 on days 1, 7 and 14 respectively (5 and 5,867 at q < 0.10), with median |logFC| of 0.073, 0.158 and 0.119. Unadjusted counts were 1, 2,563 and 958; the adjusted count is *higher* at day 7 because the muscle covariate absorbs residual variance and raises the prior degrees of freedom, which increases rather than decreases sensitivity. Excluding GSM4506604 reduced the day-7 count to 3,041 genes (Jaccard 0.584) while leaving the fold-change ranking essentially unchanged (Spearman ρ = 0.912), so the day-7 peak is robust but its exact size is driven partly by one contaminated sample.

### 3.3 Day 7 is dominated by an interferon programme, not by ion channels

Module-level mean moderated *t* at day 7 (Table 1) separates the modules cleanly. The interferon-stimulated control module was the single strongest signal (+24.9, q = 6.8 × 10⁻⁵), followed by glial/immune (+11.4), Toll-like receptors (+6.4) and cytokine/chemokine (+5.9); the housekeeping control was flat (+0.65, q = 0.52). Individual genes agree: Cd68 +2.86, Aif1 +2.18, Gfap +1.33, Ccl2 +0.81, Il1b +0.80 log₂.

Every neuronal module was down-regulated, and by amounts statistically indistinguishable from the negative controls that index neuronal transcript share: GABAergic −5.71, K⁺/HCN −5.64, neurotrophin −5.29, pan-neuronal control −5.20, erythrocyte control −4.80, glutamate/NMDA −4.35, voltage-gated Na⁺ −4.35, myelin/Schwann control −4.10, opioid receptor −3.76, Ca²⁺ α2δ −3.51. The target modules span −5.71 to −3.51 and the two compositional controls sit in the middle of that range. Five neuronal modules being down by the same amount as a pan-neuronal structural control is one finding about cellular composition, not five findings about five analgesic targets.

### 3.4 In-silico perturbation reproduces the compartment effect and adds nothing over an additive baseline

At day 7, the largest reversal was produced by knocking out the interferon module (+0.228, z = +4.61 versus the degree-matched null, 75% of the oracle, p = 0.002, q = 0.011), followed by glial/immune knockout (+0.155, z = +3.27, 61% of the oracle, q = 0.011). Because interferon and immune genes are *up* in the pain state, knocking them down moves the profile back—as expected.

The neuronal modules ran the other way, and again in lockstep with the compositional controls: pan-neuronal −0.109, K⁺/HCN −0.092, Na⁺ −0.090, myelin −0.090 (all z ≈ −2.6 to −3.6). Because these genes are down in the pain state, knockout pushes the profile further from control; overexpression pushes it back. The effect is real in the sense that it is far out in the null tail, and uninterpretable in the sense that a pan-neuronal structural gene set does the same thing.

Network propagation was almost perfectly redundant with the additive baseline: Spearman ρ = 0.991 across all 132 module × contrast × perturbation conditions. Mean absolute reversal was 0.057 for the network model and 0.033 for the additive model—the graph inflates magnitude by about 70% but leaves the ranking essentially untouched. On average the modules reached 20.5% of their oracle value.

### 3.5 The two aetiologies diverge, and the divergence is direction-specific

At the level of individual canonical markers, the divergence is large and survives multiplicity correction (Table 2, exact permutation). Myeloid markers were up in bone cancer pain and *down* in paclitaxel neuropathy: Cd68 +2.84 versus −0.58 (Δ = +3.42, q = 0.012), Aif1 +2.19 versus −0.88 (Δ = +3.08, q = 0.012); the day-14 contrast gives the same picture (Cd68 Δ = +3.07, Aif1 Δ = +2.84). Chemokine and cytokine markers follow: Ccl2 Δ = +1.39 versus oxaliplatin (q = 0.017), Il1b Δ = +1.23 versus paclitaxel (q = 0.018).

The injury-associated neuropeptide galanin moves in the opposite direction: −0.67 in bone cancer pain versus +1.92 in paclitaxel neuropathy (Δ = −2.59, q = 0.012). Calca (−0.27 versus +0.50, Δ = −0.77, q = 0.027) and Scn10a (−0.46 versus +0.24, Δ = −0.70, q = 0.026) point the same way. Atf3, the canonical injury marker, is +1.13 in paclitaxel neuropathy and +0.22 at day 7 of bone cancer pain.

At module level the same conclusion holds with greater power (Table 3): glial/immune ΔlogFC = +2.70 (bone day 7 versus paclitaxel, q = 0.007), cytokine/chemokine Δ = +1.81 (versus oxaliplatin) and +1.70 (versus paclitaxel), Toll-like receptors Δ = +1.28, with neurotrophin (Δ = −0.97), Ca²⁺ α2δ (Δ = −0.60) and neuropeptides (Δ = −1.35, q = 0.015) in the opposite direction.

### 3.6 The perturbation layer does not survive calibration

The same comparison done with reversal scores does not hold up. We evaluated 22 modules × 6 aetiology pairs × 2 perturbation directions = 264 contrast-level tests. Knockout and overexpression are exact mirror images: for every (module, aetiology pair), reversal(KO) = −reversal(OE) and the exact-permutation p is identical to eight decimal places (verified across all 132 module-pair combinations, max|Δ| = 0). They are therefore a single independent hypothesis, and the correct multiple-testing count is 132 independent tests per network scope (264 across the per-contrast and pooled-by-study scopes). Even with the corrected denominator, no module's reversal Δ reached q < 0.05 on the exact permutation p; the smallest was glial/immune (Δ = +0.274, 95% interval 0.203–0.287, p = 0.0014, q = 0.165 under the per-contrast count). Sixty-five of the 264 (or 32/132 deduplicated) tests were "significant" against the degree-matched random-gene-set null, but 14 of those 32 were negative-control modules—including the erythrocyte set with the single largest standardised effect (z = −5.78 for paclitaxel versus oxaliplatin) and the pan-neuronal control (z = −3.74); the remaining 18 were analgesic targets. A test that ranks an erythrocyte gene set above every analgesic target, and that cannot separate targets from controls, is measuring cellular composition and graph topology rather than analgesic target biology.

The erythrocyte control also differs between the two datasets at module level (ΔlogFC = −2.48, q = 0.007), which is a direct measurement of the fact that the two DRG preparations differ in blood content. The neuronal modules' Δ has the same sign as the erythrocyte Δ, so their apparent cross-aetiology difference cannot be separated from this compositional difference. The immune modules' Δ has the *opposite* sign, which is why the myeloid finding survives and the ion-channel finding does not.

### 3.7 The findings are robust to the analysis choices; the perturbation ranking is not informative

The pipeline has four free parameters with no canonical value—diffusion restart probability α, graph degree k, number of high-variance genes, and lowess span. We swept α ∈ {0.30, 0.45, 0.55, 0.70, 0.85}, k ∈ {5, 8, 12, 20} and 1,000–3,000 high-variance genes (60 settings), and separately α × k for the cross-aetiology Δ (15 settings). The module ranking was stable throughout: Spearman ρ against the reference configuration had a minimum of 0.959 and a median of 0.988 for the discovery ranking, and a minimum of 0.989 for the cross-aetiology Δ ranking; the top five modules were identical in every one of the 60 settings. Differential expression was insensitive to the lowess span (316–325 genes at q < 0.05 across spans 0.2–0.6 and with the trend disabled entirely; Spearman ρ of moderated *t* versus the reference ≥ 0.9999). So the results in Sections 3.2, 3.3 and 3.5 are not artefacts of tuning.

The robustness is exactly why the null result in Section 3.6 matters: the perturbation ranking is reproducible, and it is reproducibly uninformative. On the cross-study gene space the network model correlated 0.929 with the additive model, confirming the 0.991 obtained in the discovery set.

### 3.8 An independent single-cell dataset shows the interferon programme is model-specific

The Section 3.3 finding—a type-I interferon / myeloid programme at the Walker-256 day-7 peak—could be a property of that model or a general feature of cancer-pain DRG. An independent single-cell dataset separates the two. GSE292522 profiled mouse DRG (10x Genomics; four B16F10 cancer-nerve-injury animals versus four controls, 55,967 cells). Pseudobulk differential expression between the tumour and control samples (limma-trend, n = 4 per group) identified 1,517 up- and 703 down-regulated genes at q < 0.05; the dominant programme was neuronal injury and regeneration—Sprr1a (+1.55 log₂, t = 54.5), Atf3 (+0.84), Gadd45a, Serpinb1a, Ecel1—not type-I interferon. Every interferon-stimulated gene tested was flat (|logFC| < 0.11, q > 0.15), the myeloid markers Cd68 and Aif1 were unchanged, and the estimated myeloid cell fraction was ~0.1% and if anything lower in B16F10 (Mann–Whitney p = 0.97). Concordance with the Walker-256 day-7 programme was weak (Spearman ρ = 0.098 between pseudobulk log₂ fold changes; highest at day 1, ρ = 0.112, and vanishing by day 14, ρ = 0.030).

This qualifies the Section 3.3 conclusion: the interferon / myeloid programme is a feature of the Walker-256 tibia model at its day-7 peak, not a universal signature of cancer-pain DRG. An independent cancer-associated nerve-injury model instead recapitulates the neuronal-injury programme that, in the cross-aetiology comparison of Section 3.5, characterised CIPN. The divergence is therefore between models as much as between aetiologies. It sharpens the methodological point of the paper: a single-model transcriptional signature must not be generalised to "cancer pain" as a monolith, and cross-study claims require the negative controls and exact permutation tests applied here.

## 4. Discussion

Two conclusions, of different strength, come out of the same two datasets.

The first is biological and, we think, solid. In the Walker-256 tibia model, bone cancer pain at its transcriptional peak is an innate-immune event: type-I interferon signalling, myeloid infiltration and satellite-gliosis, with the neuronal compartment diluted rather than specifically remodelled. Paclitaxel neuropathy is a neuronal injury event: galanin, calcitonin-gene-related peptide and Atf3 up, with no myeloid infiltration at all—Aif1 and Cd68 are if anything slightly *down*. The effect sizes are large (3–3.4 log₂ units for Cd68 and Aif1 between aetiologies), direction-specific, and confirmed by exact permutation tests that have no resampling error. The independent B16F10 single-cell model (Section 3.8) reproduces the neuronal-injury programme (Sprr1a, Atf3) but not the interferon programme, so the interferon signal is model-specific; the neuronal-injury programme is the feature shared across nerve-injury aetiologies. This is consistent with the distinct clinical pharmacology of the two conditions: neuropathic pain from chemotherapy responds partly to agents acting on neuronal excitability, whereas bone pain has a substantial inflammatory component that responds to anti-resorptive and anti-inflammatory strategies [1,2,14].

The second is methodological and is a caution. In-silico perturbation, applied as it usually is—a target module, a co-expression graph, a score, a p-value against random gene sets—produced here a hierarchy that does not survive its own controls. Neuronal target modules behaved exactly like a pan-neuronal structural control. An erythrocyte control out-ranked every real target. The network was 0.991-correlated with a model that does not use a network at all. Each of these is a statement that the score is tracking cellular composition and graph topology rather than target pharmacology.

This is the same verdict the 2025 *Nature Methods* benchmark reached for far more sophisticated models [9], arrived at in a different and, we suspect, more普遍 setting. Bulk tissue perturbation studies inherit the compositional problem that cell-line perturbation assays largely avoid, and the problem is worst exactly where the clinical question is most pressing: small-n animal pain models where an immune influx is part of the phenotype.

Three specific practices made the difference between a publishable-looking hierarchy and an interpretable result, and none of them is expensive. Negative-control modules must be carried through the perturbation, not just the differential expression. Permutation nulls must be matched on network degree, because diffusion is degree-driven. And an additive baseline must be reported; if the network does not change the ranking, the graph is decoration. A fourth—force-including the target genes on the graph—is not a nicety: without it most ion-channel and opioid-receptor genes were not nodes at all, and modules were being perturbed on a network that did not contain them.

We also note two places where we were wrong earlier in this project and where the correction changed the conclusion. A one-sided p-value chosen from the perturbation label (knockout → lower tail) silently deleted the strongest reversal in the study, because it was positive. And a batched cosine that normalised a matrix of permutation axes by the grand mean and the Frobenius norm, instead of column by column, made every cross-aetiology Δ wrong by up to a factor of five while leaving the point estimates correct—a discrepancy visible only if the two code paths are compared against each other. Both are recorded in the deposited test suite.

## 5. Limitations

Both discovery datasets are rat, bulk DRG, n = 4 per group. Cellular resolution is absent, so "myeloid infiltration" is inferred from marker genes and cannot be assigned to macrophages, monocytes or neutrophils, and "satellite gliosis" cannot be separated from Schwann-cell change. The two series come from different laboratories with different dissection and perfusion protocols; platform, species and tissue are matched, batch is not, and the erythrocyte-control difference (ΔlogFC = −2.48) is a direct measurement of how much that matters. One day-7 sample carries 6.26% sarcomere transcript; removing it lowers the day-7 differential-expression count from 4,600 to 3,041 (Jaccard 0.584), so the peak is robust but its magnitude is not. At n = 4 per group the exact permutation p has a floor of 1/4,901, and the reversal score is sensitive to single animals; leave-one-out ranges are reported in the deposited artefacts. The single-cell validation (GSE292522) is a mouse B16F10 cancer-nerve-injury model—a different species, injury paradigm and quantification (10x scRNA-seq) from the rat Walker-256 discovery—so its role is to show that the interferon programme does not generalise, not to confirm it, and the weak cross-model concordance (ρ = 0.098) itself argues against over-generalising a single model. Finally, nothing here is a human validation, and the translation from a rat tibia tumour model to patients with bone metastases involves a species extrapolation that we have not bridged.

## 6. Conclusions

Bone cancer pain and chemotherapy-induced peripheral neuropathy drive qualitatively different programmes in the DRG—myeloid and interferon-driven in the first, neuronal-injury-driven in the second—and that difference is large, direction-specific and statistically exact. The same data do not support ranking analgesic targets by in-silico perturbation: once negative-control modules, degree-matched nulls and an additive baseline are included, the apparent target hierarchy dissolves into a cellular-composition effect. For translational pain research the implication is that mechanism stratification should be built on marker-level immunecompartment evidence, which these datasets do support, and that in-silico perturbation results should be reported with the calibrations described here or treated as hypothesis generation only.

---

## Data and code availability

All data are public: GSE149648 and GSE160543 (NCBI Gene Expression Omnibus) for discovery, and GSE292522 (NCBI GEO) for single-cell validation. Analysis code, artefacts and the test suite are deposited at the project repository, with `requirements.txt` pinning every version used.

## References

1. Fallon M, Giusti R, Aielli F, et al. Management of cancer pain in adult patients: ESMO Clinical Practice Guidelines. *Ann Oncol*. 2018;29(Suppl 4):iv166–iv191.
2. Swarm RA, Paice JA, Bedard A, et al. Adult cancer pain, Version 3.2019. *J Natl Compr Canc Netw*. 2019;17(8):977–1006.
3. Mantyh PW. Bone cancer pain: causes, consequences, and therapeutic opportunities. *Pain*. 2013;154(Suppl 1):S54–S62.
4. Falk S, Dickenson AH. Pain and nausea in cancer: mTOR-mediated mechanisms. *Br J Anaesth*. 2015;115(2):169–171.
5. Staff NP, Grisold A, Grisold W, Windebank AJ. Chemotherapy-induced peripheral neuropathy: a current review. *Ann Neurol*. 2017;81(6):772–781.
6. Loprinzi CL, Lacchetti C, Bleeker J, et al. Prevention and management of chemotherapy-induced peripheral neuropathy in survivors of adult cancers: ASCO guideline update. *J Clin Oncol*. 2020;38(28):3325–3348.
7. Lotfollahi M, Wolf FA, Theis FJ. scGen predicts single-cell perturbation responses. *Nat Methods*. 2019;16(8):715–721.
8. Roohani Y, Huang K, Leskovec J. Predicting transcriptional outcomes of novel multigene perturbations with GEARS. *Nat Biotechnol*. 2024;42(6):927–935.
9. Ahlmann-Eltze C, Huber W, Anders S. Deep-learning-based gene perturbation effect prediction does not yet outperform simple linear baselines. *Nat Methods*. 2025;22(8):1657–1661. doi:10.1038/s41592-025-02772-6
10. Zhai M, Yang S, Lin S, Zhu H. Distinct gene expression patterns of ion channels and cytokines in rat primary sensory neurons during development of bone cancer and cancer pain. *Front Mol Neurosci*. 2021;14:665085. doi:10.3389/fnmol.2021.665085
11. Sun W, Hao Y, Li R, Ho IHT, et al. Comparative transcriptome of dorsal root ganglia reveals distinct etiologies of paclitaxel- and oxaliplatin-induced peripheral neuropathy in rats. *Neuroscience*. 2023;516:1–14. doi:10.1016/j.neuroscience.2023.02.007
12. Smyth GK. Linear models and empirical Bayes methods for assessing differential expression in microarray experiments. *Stat Appl Genet Mol Biol*. 2004;3:Article3.
13. Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *J R Stat Soc B*. 1995;57(1):289–300.
14. Luger NM, Mach DB, Sevcik MA, Mantyh PW. Bone cancer pain: from model to mechanism to therapy. *J Pain Symptom Manage*. 2005;29(5 Suppl):S32–S46.

---

## Tables

**Table 1. Module-level mean moderated t at day 7 (muscle-adjusted).**

| Module | Kind | n | mean t | mean logFC | q |
|---|---|---|---|---|---|
| CTRL interferon | control | 8 | +24.92 | +2.61 | 6.8 × 10⁻⁵ |
| Glial / immune | target | 5 | +11.41 | +2.01 | 0.0027 |
| Toll-like receptors | target | 3 | +6.36 | +0.87 | 0.0245 |
| Cytokine / chemokine | target | 4 | +5.88 | +1.38 | 0.0129 |
| Purinergic P2X | target | 3 | +0.69 | +0.04 | 0.272 |
| Prostanoid / COX | target | 4 | +1.07 | +0.29 | 0.315 |
| CTRL housekeeping | control | 11 | +0.65 | +0.04 | 0.524 |
| Injury / regeneration | target | 3 | −1.40 | −0.09 | 0.791 |
| Monoaminergic | target | 3 | −1.08 | −0.13 | 0.961 |
| Cannabinoid | target | 4 | −0.65 | −0.03 | 0.934 |
| TRP channels | target | 4 | −2.47 | −0.29 | 0.221 |
| Neuropeptides | target | 5 | −2.63 | −0.52 | 0.136 |
| Ca²⁺ channel α2δ | target | 4 | −3.51 | −0.24 | 0.136 |
| Opioid receptors | target | 4 | −3.76 | −0.40 | 0.077 |
| CTRL myelin / Schwann | control | 7 | −4.10 | −0.39 | 0.017 |
| Voltage-gated Na⁺ | target | 5 | −4.35 | −0.42 | 0.021 |
| Glutamate / NMDA | target | 3 | −4.35 | −0.38 | 0.060 |
| CTRL erythrocyte | control | 4 | −4.80 | −1.15 | 0.019 |
| CTRL pan-neuronal | control | 10 | −5.20 | −0.40 | 0.00031 |
| Neurotrophin | target | 3 | −5.29 | −0.52 | 0.033 |
| K⁺ / HCN channels | target | 5 | −5.64 | −0.48 | 0.0053 |
| GABAergic | target | 3 | −5.71 | −0.37 | 0.025 |

**Table 2. Cross-aetiology difference in canonical markers (exact permutation, 4,900 null pairs).**

| Gene | Bone (day 7) | Bone (day 14) | PTX | OXA | Largest Δ | q |
|---|---|---|---|---|---|---|
| Cd68 | +2.84 | +2.49 | −0.58 | +0.44 | +3.42 (D7 vs PTX) | 0.012 |
| Aif1 | +2.19 | +1.95 | −0.88 | −0.46 | +3.08 (D7 vs PTX) | 0.012 |
| Gal | −0.67 | +0.01 | +1.92 | −0.30 | −2.59 (D7 vs PTX) | 0.012 |
| Ccl2 | +0.90 | +0.84 | −0.23 | −0.49 | +1.39 (D7 vs OXA) | 0.017 |
| Il1b | +0.80 | +0.74 | −0.43 | −0.31 | +1.23 (D7 vs PTX) | 0.018 |
| Calca | −0.27 | −0.18 | +0.50 | +0.03 | −0.77 (D7 vs PTX) | 0.027 |
| Scn10a | −0.46 | −0.21 | +0.24 | −0.11 | −0.70 (D7 vs PTX) | 0.026 |
| Atf3 | +0.22 | +0.51 | +1.13 | +0.16 | — | ns |
| Gfap | +1.33 | +0.65 | +0.09 | +0.08 | — | ns |

**Table 3. Cross-aetiology difference in module mean logFC (exact permutation).**

| Module | Kind | logFC bone D7 | logFC bone D14 | logFC PTX | logFC OXA | Δ (D7 vs PTX) | q |
|---|---|---|---|---|---|---|---|
| Glial / immune | target | +2.05 | +1.30 | −0.64 | −0.13 | +2.70 | 0.0073 |
| Cytokine / chemokine | target | +1.40 | +1.46 | −0.31 | −0.35 | +1.70 | 0.0073 |
| Toll-like receptors | target | +0.85 | +0.68 | −0.42 | −0.14 | +1.28 | 0.0073 |
| Neurotrophin | target | −0.57 | −0.29 | +0.41 | — | −0.97 | 0.0073 |
| Ca²⁺ channel α2δ | target | −0.28 | −0.14 | +0.32 | — | −0.60 | 0.0073 |
| Neuropeptides | target | −0.52 | −0.18 | +0.83 | — | −1.35 | 0.0145 |
| CTRL pan-neuronal | control | −0.44 | −0.37 | +0.25 | — | −0.69 | 0.0145 |
| Glutamate / NMDA | target | −0.44 | — | +0.20 | — | −0.64 | 0.0162 |
| Voltage-gated Na⁺ | target | −0.49 | — | +0.16 | — | −0.65 | 0.0310 |
| GABAergic | target | −0.43 | — | +0.06 | — | −0.49 | 0.0310 |
| CTRL erythrocyte | control | −1.05 | — | — | +1.43 | −2.48 (D7 vs OXA) | 0.0073 |

**Table 4. Calibration of the in-silico perturbation layer (day 7, discovery).**

| Module | Kind | Reversal (knockout) | z vs degree-matched null | % of oracle | p | q |
|---|---|---|---|---|---|---|
| CTRL interferon | control | +0.228 | +4.61 | 75.4 | 0.0020 | 0.011 |
| Glial / immune | target | +0.155 | +3.27 | 61.5 | 0.0020 | 0.011 |
| CTRL pan-neuronal | control | −0.109 | −3.56 | 32.9 | 0.0020 | 0.011 |
| K⁺ / HCN channels | target | −0.092 | −2.87 | 36.5 | 0.0020 | 0.011 |
| Voltage-gated Na⁺ | target | −0.090 | −2.63 | 35.9 | 0.0080 | 0.020 |
| CTRL myelin / Schwann | control | −0.090 | −2.71 | 30.9 | 0.0080 | 0.020 |
| CTRL erythrocyte | control | +0.074 | +2.52 | 33.0 | 0.0040 | 0.016 |
| Neuropeptides | target | +0.067 | +2.56 | 26.9 | 0.0040 | 0.016 |

Spearman ρ between network and additive reversal = 0.991 (n = 132). Mean absolute reversal: network 0.057, additive 0.033, ridge 0.712. Mean percentage of oracle 20.5%.

**Table 5. Independent single-cell validation (GSE292522, B16F10 cancer-nerve-injury DRG).**

| Quantity | Value |
|---|---|
| Samples / cells | 4 tumour, 4 control; 55,967 cells |
| Pseudobulk DE (q < 0.05) | 1,517 up / 703 down |
| Top up-regulated genes | Sprr1a +1.55 log₂ (t = 54.5), Atf3 +0.84 — neuronal-injury programme |
| Interferon-stimulated genes | all flat (|logFC| < 0.11, q > 0.15) |
| Myeloid markers Cd68 / Aif1 | unchanged (Δ ≈ 0) |
| Myeloid cell fraction | 0.12% (control) vs 0.05% (B16F10); Mann–Whitney p = 0.97 |
| Concordance with Walker-256 | Spearman ρ = 0.098 (D7), 0.112 (D1), 0.030 (D14) |
