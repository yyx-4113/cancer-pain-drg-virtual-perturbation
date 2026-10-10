# A3 — Single-Cell Omics Methodology (Design) Review

**Manuscript:** `docs/manuscript/manuscript_final.md` (§2.8, §3.8, Table 5, Abstract, §5)
**Code reviewed:** `scripts/30_fetch_sc.py`, `scripts/31_process_sc.py`, `scripts/32_percelltype_sc.py`, `scripts/_limma.py`
**Artefacts recomputed from raw 10x matrices:** `sc_GSE292522_celltype_DE.csv`, `_celltype_scores.csv`, `_celltype_summary.json`, `sc_GSE292522_pseudobulk_DE.csv`, `sc_GSE292522_summary.json`, `sc_GSE292522_composition.csv`

Reviewer posture: first submission, no credit for effort. Every number below was
recomputed by me from `data/processed/GSE292522/*.mtx.gz` (55,335 features ×
55,967 barcodes across 8 samples) or derived from code I read line by line.

---

## Verdict up front

The per-compartment layer is built on a **statistically correct frame** (animal-level
pseudobulk, real limma-trend, no pseudoreplication) wrapped around an
**annotation whose bins do not mean what §3.8 says they mean**. Three of the
four headline compartment sentences in §3.8 fail on the paper's own artefacts:

1. The compartment where an interferon programme would most plausibly live —
   **myeloid** — was silently dropped by the ≥5-cells/animal rule (my recomputed
   myeloid cells per animal: 3, 10, 12, 5 in controls; 5, 3, 7, 1 in B16F10 — only
   5 of 8 animals clear the threshold). "No compartment shows ISG
   up-regulation" is therefore a claim about the compartments *that survived
   filtering*, and the immune compartment is not among them.
2. "**The only** significant ISG changes were Ifit2 in neurons and Isg15 in
   satellite glia" is contradicted by the deposited table itself: Ifit1 in
   SatelliteGlia, q = 0.0037 (my recomputation: 0.0036), is also significant and is
   not mentioned. The Isg15 result it does cite does not survive leave-one-out:
   refitting limma-trend per drop, q < 0.05 in **1 of 8** refits; the control
   animal Sample4 carries an ISG score of 0.9501 versus 0.42–0.45 in the other
   three controls.
3. "The neuronal-injury programme was **neuron-intrinsic**" is contradicted by the
   deposited table: Gal is up in Schwann (+0.453, q = 0.0062) and Endothelium
   (+0.671, q = 0.0001), both robust under leave-one-out (8/8 and 7/8 refits at
   q < 0.05). And the "Neuron" bin itself is a mixed bin: its median cell carries
   only 16–49 pan-neuronal counts (Snap25+Tubb3+Syp) and 3–6 myelin-gene counts,
   with Schwann as the runner-up label for 68–78% of its cells. Rebuilding the
   neuron pseudobulk from a strict, independently defined neuronal subset
   (476–759 cells per animal, all 8 animals) roughly halves the effects: Gal
   +0.781 (q = 0.353), Scn10a −0.135 (q = 0.59).

**What genuinely survives** — and it is not nothing: the animal-level pseudobulk
machinery is correct and reproduces the deposited numbers exactly; galanin
up-regulation is real at the per-cell level (mean Gal counts per neuron-bin cell
rise 2.3–3.0×, from 4.4–11.1 to 17.8–26.6, so it is not a within-compartment
normalisation artefact); Gal and Scn10a point estimates are leave-one-out stable
in the bins as constructed; and the narrow claim "no ISG is up-regulated at
q < 0.05 in any tested compartment" I verified exhaustively (29–38 ISGs × 8
compartments, zero hits). The injury programme is real. The compartment
*attribution* of it is not established by this analysis.

---

## Findings

### Finding 1 — CRITICAL. The myeloid compartment was silently excluded, so "absence in every compartment" is not supported

**【Problem】** The sentence "No compartment — neurons, neurofilament-positive or
peptidergic neuronal subtypes, satellite glia, Schwann cells, endothelium,
fibroblasts or erythrocytes — showed interferon-stimulated-gene up-regulation"
omits, without disclosure, the two compartments most relevant to an interferon
claim: myeloid cells and T cells, both of which were defined in the code and
both of which were dropped by the ≥5-cells/animal rule.

**【Evidence】** `32_percelltype_sc.py:131-138` builds a Myeloid and a TCell
pseudobulk whenever a sample has ≥5 assigned cells, and `32_percelltype_sc.py:149-150`
drops the compartment unless ≥6 samples qualify. I recomputed per-sample assigned
cell counts from the raw matrices (identical argmax rule, identical marker sets):

| sample | group | Myeloid | TCell | NP | TH | Erythrocyte |
|---|---|---|---|---|---|---|
| Sample1 Con rep1 | Con | 3 | 0 | 1 | 0 | 5 |
| Sample2 Con rep2 | Con | 10 | 0 | 0 | 0 | 4 |
| Sample3 Con rep3 | Con | 12 | 0 | 0 | 0 | 22 |
| Sample4 Con rep4 | Con | 5 | 0 | 0 | 0 | 3 |
| Sample5 B16F10 rep1 | B16F10 | 5 | 0 | 0 | 2 | 7 |
| Sample6 B16F10 rep2 | B16F10 | 3 | 0 | 1 | 2 | 5 |
| Sample7 B16F10 rep3 | B16F10 | 7 | 1 | 1 | 3 | 26 |
| Sample8 B16F10 rep4 | B16F10 | 1 | 0 | 0 | 1 | 15 |

Myeloid: 5 of 8 samples have ≥5 cells → compartment dropped. TCell: 1 cell in the
entire dataset → dropped. Neuron.NP: 3 cells total → dropped. Neuron.TH: 8 cells
total → dropped. Erythrocyte: 6 samples qualify, but the two dropped ones are both
controls → tested at **2 controls vs 4 tumours**. None of these drops appears in
§2.8 or §3.8; the whole-DRG myeloid fraction is 0.12% vs 0.05% per
`sc_GSE292522_summary.json`, i.e. myeloid cells exist but are just under the
per-animal threshold. The abstract compounds this: "this held in every DRG
compartment tested, confirming the interferon programme is model-specific."

**【Why it matters】** In the Walker-256 discovery data, the interferon programme
co-segregates with myeloid infiltration (§3.3: Cd68 +2.86, Aif1 +2.18). A reader
of §3.8 will reasonably conclude that the immune compartment was examined and
was ISG-negative. It was not examined. The claim "the interferon programme is
absent at the compartment level" (§3.8, last sentence of the second paragraph)
is currently supported only for non-immune compartments, and in the B16F10 model
the myeloid compartment is precisely where a residual IFN signal would hide. Note
this does not *weaken* the model-specificity thesis — a myeloid fraction of
0.05–0.12% is itself strong evidence against Walker-256-like infiltration — but
the ISG-absence claim and the infiltration-absence claim are different claims,
and only the second is actually supported at compartment level.

**【Specific fix】** Add to §3.8: "The myeloid and T-cell compartments and the NP and
TH neuronal subtypes could not be tested: they contained fewer than five assigned
cells in at least two animals (myeloid cells: 1–12 per animal; T cells: 1 cell in
total), so the per-compartment ISG claim covers the seven non-immune compartments
plus NF and PE neurons, not the immune compartment. The near-absence of assigned
myeloid cells in both groups (whole-DRG fraction 0.12% control vs 0.05% B16F10)
is itself the compartment-level finding: the Walker-256 myeloid infiltration has
no counterpart in this model." Change the abstract's "confirming" to "consistent
with". Additionally deposit a table of per-compartment per-animal assigned cell
counts (the analysis computes them; they are not in any artefact).

### Finding 2 — CRITICAL. "The only significant ISG changes" is false on the paper's own table, and the Isg15 result it names does not survive leave-one-out

**【Problem】** §3.8 states "the only significant ISG changes were down-regulation of
Ifit2 in neurons (q = 0.046) and Isg15 in satellite glia (q = 0.007)". Both halves
of the sentence are wrong: a third significant ISG change exists in the deposited
table and is unmentioned, and the Isg15 result is a one-animal effect.

**【Evidence】** (a) `results/tables/sc_GSE292522_celltype_DE.csv` row
`SatelliteGlia,Ifit1,-0.523,-4.73,0.0037`. I reproduced this independently
(my refit with the project's own `_limma.fit_limma`: logFC −0.523, t −4.74,
q 0.0036). Ifit1 is in the FOCUS gene list (`32_percelltype_sc.py:80-82`), so this
is not a hidden row — it is a row the text ignores while claiming exclusivity.
(b) Leave-one-out: I refit limma-trend with the project's `_limma` implementation
on a fixed 3,000-gene universe, dropping each animal in turn. SatelliteGlia/Isg15:
full-data q = 0.083 on that universe; LOO q range 0.0097–0.8039; **q < 0.05 in 1 of
8 refits**. Welch t on the same values: full p = 0.0385, max LOO p = 0.1397. The
driver is visible in the deposited `sc_GSE292522_celltype_scores.csv`: control ISG
scores are 0.4266, 0.4229, 0.4472, **0.9501** (Sample4) — one control animal is
2.1–2.2× the others. (c) SatelliteGlia/Ifit1: q < 0.05 in **0 of 8** LOO refits;
full-data Welch p = 0.1418. (d) The one counterweight: the direction is not
purely one-animal — all four controls exceed all four B16F10 samples (that is why
the MWU p = 0.0286, the two-sided floor at n=4+4, C(8,4)=70), and in my full
39-gene ISG screen the SatelliteGlia compartment is the only one with a clearly
negative median ISG logFC (−0.195), with Cxcl10 (−1.63), Ifit1, Ifit2, Isg15,
Rsad2, Gbp2, Ifi47 all down at q < 0.10. The composite direction is real; the
published gene-level q values are not.

**【Why it matters】** A reviewer who checks the deposited CSV against the sentence
finds an immediate factual error. A reviewer who checks stability finds that one
of the two named results is a single-animal effect. Either one costs credibility
for the entire section, including the parts that are correct.

**【Specific fix】** Replace the sentence with: "Interferon-stimulated genes were
never up-regulated at q < 0.05 in any tested compartment (38 ISGs × 8
compartments); the only directions of change were down-regulation — genome-wide
in neurons (Ifit2, q = 0.046) and broadly in the satellite-glia compartment, where
7 of 29 tested ISGs were down at q < 0.10 (Isg15, Ifit1, Ifit2, Cxcl10, Rsad2,
Gbp2, Ifi47), although the individual gene-level calls in that compartment are
sensitive to a single control animal (leave-one-out q < 0.05 in 1/8 refits for
Isg15) and the compartment's composition is itself uncertain (Finding 3)."

### Finding 3 — CRITICAL. The "SatelliteGlia" compartment contains essentially no satellite glia by its own defining marker

**【Problem】** The compartment that carries the paper's named ISG finding has a
Gfap content of zero, and its assignments are coin-flips against Schwann.

**【Evidence】** I recomputed, from the raw matrices, the raw counts of the four
glial markers inside the argmax-SatelliteGlia bin (25–108 cells per animal):

| sample | n_cells | Gfap | Slc1a3 | Aqp4 | Sox10 |
|---|---|---|---|---|---|
| Sample1 Con | 46 | **0** | 44 | 48 | 48 |
| Sample2 Con | 35 | **0** | 34 | 45 | 37 |
| Sample3 Con | 37 | **0** | 31 | 39 | 43 |
| Sample4 Con | 25 | **3** | 26 | 36 | 29 |
| Sample5 B16F10 | 36 | **0** | 32 | 41 | 39 |
| Sample6 B16F10 | 44 | **2** | 30 | 48 | 22 |
| Sample7 B16F10 | 108 | **2** | 69 | 86 | 91 |
| Sample8 B16F10 | 30 | **0** | 27 | 37 | 42 |

Total Gfap across the entire compartment, all eight animals: **7 counts**.
Satellite glia are the most Gfap-expressing cells in a healthy DRG; a bin with
7 Gfap counts in ~361 cells is not a satellite-glia bin. Consistently, the median
top1−top2 annotation margin in this bin is 0.041–0.078 (per sample), versus
1.0–1.4 for the "Neuron" bin — the score difference between the winning and
runner-up marker sets is a few percent of a log unit. The marker sets themselves
share Sox10 (`32_percelltype_sc.py:57-58`), so they cannot discriminate the two
glial bins by construction. This also explains why Gfap's row in the deposited
`celltype_DE.csv` for SatelliteGlia is logFC 0.000, t = 0.000: the gene is
constant-zero in both groups.

**【Why it matters】** §3.8 attributes an ISG signal to "satellite glia" (q = 0.007)
and §2.8 tells the reader cells were assigned to canonical DRG types. The bin so
labelled has no demonstrable satellite-glia identity; whatever cells it contains
are more plausibly Schwann-lineage or ambiguous. A compartment-level claim
inherited the label of a compartment that was never validated. Note the whole-DRG
Gfap signal exists (whole-DRG pseudobulk Gfap up, q = 0.006, §3.8), so Gfap+
cells exist in the data — they are simply not in this bin.

**【Specific fix】** Either validate the bin or rename it. Minimal validation to
report in §2.8: per-compartment median counts of the compartment's own defining
markers, per animal. If the Gfap content remains ~0, rename the bin (e.g.
"Sox10/Slc1a3-low glial") or merge SatelliteGlia with Schwann, re-run, and report
the merged result. A paste-ready sentence: "Because the satellite-glia and
Schwann marker sets share Sox10 and the satellite-glia-assigned cells carried
almost no Gfap transcript (7 counts in 361 cells), we treat the two glial bins as
a single Schwann-lineage compartment and do not interpret them separately."

### Finding 4 — MAJOR. The "Neuron" compartment is a mixed bin; on a strict neuronal definition the compartment effects halve and lose significance

**【Problem】** The bin that carries "galanin was up +1.47 in neurons (q < 0.001)" and
"Scn10a (Nav1.8) was down in neurons (q = 0.012)" is not a validated neuronal
compartment, and the headline effect sizes are inflated roughly twofold relative
to a defensible neuronal subset.

**【Evidence】** I profiled the argmax-Neuron bin (4,226–5,840 cells per sample):
median pan-neuronal counts (Snap25+Tubb3+Syp) per assigned cell **16–49**; only
**47.6–62.5%** of bin cells have ≥20 such counts; median myelin-gene counts
(Mpz+Pmp22+Plp1+Mag) **3–6 per cell** — i.e. most cells in the "neuron" bin carry
myelin transcript; the runner-up label is **Schwann for 67.9–77.9%** of bin cells
(and SatelliteGlia for a further 14–24%). Neurons are the largest cells in the
DRG and are known to be under-captured in whole-tissue 10x, so a neuronal bin
comprising 72–80% of all captured cells is implausible on biological grounds as
well. I then rebuilt the neuron pseudobulk from a strict, independent criterion —
(Snap25+Tubb3+Syp ≥ 20 counts) AND (Mpz+Pmp22+Plp1+Mag < 5) — which yields
476–759 cells in **every one of the 8 animals** (no animal dropped), and ran the
identical limma-trend: Gal **+0.781, t = 3.31, q = 0.353** (deposited bin: +1.472,
t = 10.17, q = 0.0001); Scn10a **−0.135, q = 0.59** (deposited: −0.247, q = 0.012);
Ifit2 −0.016 (null); Snap25 −0.143, Nefl −0.139 (null); Cd68/Aif1 null; 73 of
27,228 genes at q < 0.05. Direction is preserved for every checked gene; the
magnitudes and significance are not. Relatedly, regressing each compartment's
gene-wise logFC on the whole-DRG pseudobulk logFC gives slope **0.856,
Spearman ρ = 0.904** for Neuron and slope **1.011, ρ = 0.827** for Neuron.NF —
the "neuronal" compartment largely re-measures the whole-DRG contrast, as expected
for a bin holding ~75% of all cells, so it cannot by itself localise anything.

**【Why it matters】** The two compartment-level neuronal claims in §3.8 (+1.47 Gal
"neuron-intrinsic"; Scn10a down q = 0.012) are the quantitative basis for "the
injury programme localises to neurons". With the strict definition the effects
remain directionally present but are no longer significant, and the claimed
localisation is contradicted within the same table by Gal up in Schwann and
Endothelium (Finding 5). The honest statement is "the injury programme is present
in neuronal cells (per-cell Gal counts rise 2.3–3.0×) and is not separable from a
whole-tissue shift at this sample size", not "neuron-intrinsic".

**【Specific fix】** (i) Deposit the per-bin composition diagnostics (median
pan-neuronal and myelin counts, runner-up label distribution) as a supplementary
table. (ii) Re-state the neuronal findings from a validated neuronal subset. A
paste-ready sentence: "Restricting neurons to cells with ≥20 pan-neuronal counts
and fewer than 5 myelin-gene counts (476–759 cells per animal, all eight animals
retained), galanin remained directionally up (+0.78 log₂, nominal t = 3.31, not
significant after within-compartment correction) and Scn10a directionally down
(−0.14, not significant), so we report the neuronal localisation of the injury
programme as directionally supported but not significant at n = 4 + 4."
(iii) Do not use the phrase "neuron-intrinsic" anywhere; it is contradicted by
Finding 5.

### Finding 5 — MAJOR. "Neuron-intrinsic" is contradicted by the paper's own table, and the neuronal bins carry measurable cross-contamination signal

**【Problem】** Galanin is significantly up in the Schwann (+0.453, q = 0.0062) and
Endothelial (+0.671, q = 0.0001) compartments — both robust under leave-one-out
(Schwann 8/8 refits at q < 0.05; Endothelial 7/8, worst q = 0.053) — so the injury
programme is not neuron-restricted even on the paper's own bins; and trace
myeloid/glial genes move inside the neuronal bins with large t statistics,
indicating sensitivity to a handful of contaminated cells.

**【Evidence】** My refits reproduce the deposited Schwann/Gal and Endothelial/Gal
rows exactly (6.01/0.0062 and 6.29/0.0001). In the deposited table, the
Neuron.NF bin also shows Cx3cr1 up (+0.019, t = 5.56, q = 0.006), Gfap up (+0.001,
t = 3.74, q = 0.047) and Il1b up (+0.002, t = 3.89, q = 0.040) — genes of
myeloid and glial identity "up-regulated" inside a neuronal bin with logFC ≤
0.019. In the Endothelial bin, Lyz2 is "down" (−0.814, q = 0.0008), i.e. the bin's
small myeloid contamination differs between groups. These are the signature of
doublets and ambient RNA moving with group composition, not of regulation in the
labelled cell type.

**【Why it matters】** "The neuronal-injury programme was neuron-intrinsic" (§3.8)
and Table 5's "Neuronal-injury localisation" row are the section's positive
claims. A compartment analysis whose bins demonstrably exchange signal across
labels cannot support a localisation claim of this strength, and the
contradicting rows (Gal in Schwann/Endothelium) sit in the paper's own deposited
artefact.

**【Specific fix】** Replace "Conversely the neuronal-injury programme was
neuron-intrinsic: galanin was up +1.47 log₂ in neurons (q < 0.001) and +2.54 in
neurofilament-positive neurons (q < 0.001)" with: "Conversely the injury programme
was strongest in neuronal bins: galanin rose +1.47 log₂ in the neuronal bin and
+2.54 in neurofilament-positive neurons (both q < 0.001), with smaller rises in
Schwann (+0.45, q = 0.006) and endothelial (+0.67, q < 0.001) bins that are
consistent with partial cross-contamination of bins or with a broader injury
response; we therefore report the neuronal enrichment of the injury programme as
graded rather than exclusive."

### Finding 6 — MAJOR. The erythrocyte "compartment" is not a cell type, was tested at 2 vs 4 animals, and neither fact is disclosed

**【Problem】** §3.8 lists erythrocytes among the compartments that showed no ISG
up-regulation. The bin's cells have a median of 513–1,163 detected genes
(incompatible with enucleated erythrocytes), and the compartment lost both of its
weakest control samples to the ≥5-cell rule, leaving a 2-control vs 4-tumour test.

**【Evidence】** My per-compartment QC on the argmax-Erythrocyte bin: median
detected genes per cell 513–1,163 and median UMIs 1,487–6,281 across the eight
samples (per-animal cell counts in Finding 1's table; Sample2 with 4 cells and
Sample4 with 3 cells are dropped as controls). The limma-trend test therefore ran
on n = 2 + 4 (df_res = 4) over 14,601 genes, finding 0 genes at q < 0.05. The
per-animal cell counts are also wildly uneven (3–26 cells), so the pseudobulk
varies by orders of magnitude more than for other bins.

**【Why it matters】** A compartment of 3–26 RBC-derived, RNA-contaminated
barcodes contributes no biological information; including it in the "no ISG
up-regulation" list manufactures the appearance of a more complete compartment
survey than occurred, and the unbalanced 2-vs-4 design is never mentioned. This
is the kind of detail that, once noticed, makes a reader re-price every other
undisclosed choice in the section.

**【Specific fix】** Drop the erythrocyte compartment from the per-compartment claim,
or report it with its design: "the erythrocyte bin (RBC-derived barcodes with
ambient RNA, 3–26 cells per animal, testable in only 2 of 4 controls) is reported
for completeness and carries no inferential weight." A paste-ready replacement
for §2.8: "Compartments present in fewer than six animals, or with fewer than five
cells in an animal, were dropped; this removed the myeloid, T-cell, NP and TH
compartments entirely and left the erythrocyte compartment with two of four
control animals."

### Finding 7 — MODERATE. §2.8 misdescribes the implementation in four specific, checkable ways

**【Problem】** The Methods text for the single-cell layer does not match the code
that produced the results.

**【Evidence】**
(a) **"converted to log₂(CPM + 1)"** — both scripts compute counts per **10
thousand**, not per million: `31_process_sc.py:126` and `32_percelltype_sc.py:154`
(`div(...) * 1e4`). All log₂ fold-changes in §3.8 are CP10K units, ~1.7% of CPM
scale for the +1 offset; the logFC values themselves are unaffected because the
scale factor is constant across samples, but the Methods text is wrong.
(b) **"55,967 cells after quality filtering"** (§2.8) — no QC of any kind exists in
`30/31/32`; 55,967 is exactly the sum of raw barcodes in the eight deposited
matrices (6,514+7,322+6,019+5,817+6,175+7,480+9,826+6,814 — I recomputed each).
The sentence attributes filtering to this analysis that no code performs. (If the
GEO deposit happens to be CellRanger's filtered matrix, the filtering was the
submitters', and should be described as such.)
(c) **"assigned to one of eight compartments — neurons, and within them NF, NP,
PE and TH subtypes, plus satellite glia, Schwann cells, endothelium, fibroblasts
and erythrocytes"** — that list has ten labels for "eight"; the code defines 12
compartments (8 broad + 4 subtypes) and tests 8, after silently dropping Myeloid,
TCell, NP and TH (Finding 1).
(d) **Annotation by "the marker set with the highest mean log-expression"** is
accurate as far as it goes but omits the two consequential details: scores are
computed on raw counts without any per-cell size normalisation
(`32_percelltype_sc.py:97-102`), and there is no minimum-score floor, so a
zero-marker cell defaults to whichever set indexes first — Neuron
(`32_percelltype_sc.py:120`). In this dataset the failure mode is negligible (I
counted 0, 0, 0, 0, 0, 0, 1, 1 all-zero cells per sample), but the absence of a
floor is still an unstated design choice. Relatedly, three listed marker genes do
not exist in the features file at all and are silently skipped by the scoring
(`marker_score` drops missing genes): **PgP9.5** (mis-capitalised in 32; "Pgp9.5"
in 31 — the correct mouse symbol Uchl1 is absent from both), **Oct6** (Pou3f1),
**Blbp** (Fabp7). I verified all three are absent from the features file.

**【Why it matters】** Each is small; together they mean a reader cannot reconstruct
the analysis from the Methods, which is the one thing a Methods section is for.
(b) in particular is a reproducibility claim the artefacts do not support.

**【Specific fix】** (a) Replace "CPM" with "counts per 10,000 (CP10K)" in §2.8 and
Table 5 context. (b) Replace "55,967 cells after quality filtering" with
"55,967 cell barcodes as deposited (GEO-processed matrices; no additional
quality filtering was applied)". (c) Replace the compartment list with the
truthful one, including the drops (see Finding 6's fix). (d) Append to §2.8:
"Marker genes absent from the reference annotation (Uchl1/Pgp9.5, Pou3f1/Oct6,
Fabp7/Blbp) were dropped from their sets; scores used raw counts without
size-factor normalisation, and cells with no marker counts (≤2 per sample) were
assigned to the first marker set."

### Finding 8 — MODERATE. BH scope and gene-universe inconsistency between the two single-cell scripts

**【Problem】** Multiplicity correction is applied within each compartment over all
tested genes — defensible, but unstated — and the per-compartment pipeline uses a
different gene universe from the whole-DRG pipeline it is said to mirror.

**【Evidence】** `32_percelltype_sc.py:160-161`: `fit["q"] = bh(fit["P"].values)`
where `fit` covers every gene with nonzero sum in that compartment (my
recomputation: 16,077–34,076 genes depending on compartment). So q values are
BH-within-compartment. Meanwhile `31_process_sc.py:99-100` collapses duplicate
gene symbols by summation before fitting (55,335 feature rows → 34,603 unique
genes), but `32_percelltype_sc.py` does not deduplicate at all — I verified the
deposited features file contains **105 duplicated symbols** (all with distinct
Ensembl IDs, e.g. Ptp4a1, Arhgef4, Sept2). None of the FOCUS or ISG genes is
duplicated (verified), so no deposited focus-gene q is affected, but the
compartment-wide DE counts (e.g. the 1,890 q<0.05 genes in Neuron.NF printed by
the script; my deduplicated refit gives 1,886) use a different universe and
double-count 105 genes, and the Methods describe neither the BH scope nor the
universe.

**【Why it matters】** "1,890 DE genes" style statements are not reproducible from
the deposited artefacts (the per-compartment full fits are not deposited, only 27
FOCUS genes per compartment), and the same BH letter (q) means different things
in `sc_GSE292522_pseudobulk_DE.csv` (over 34,603 genes) and
`sc_GSE292522_celltype_DE.csv` (per-compartment, ~14.6k–34k genes).

**【Specific fix】** Add to §2.8: "Benjamini–Hochberg correction was applied within
each compartment over all genes detected in that compartment's pseudobulk
(14,600–34,100 genes); duplicate gene symbols in the 10x reference were retained
in the per-compartment fit (105 duplicated symbols), so compartment-wide gene
counts are not directly comparable to the deduplicated whole-DRG fit. Full
per-compartment fits are deposited." And deduplicate in 32 for consistency.

### Finding 9 — the claim that does survive, stated precisely

**【Problem】** "No compartment showed ISG up-regulation" is the section's central
negative claim and needs an exact threshold and an exact scope.

**【Evidence】** I refit all eight surviving compartments with the project's
`_limma.fit_limma` (trend = True, BH within compartment) and extracted every
member of the paper's 39-gene ISG set (29–38 testable per compartment): **zero
ISGs are up-regulated at q < 0.05 in any compartment**. At the looser q < 0.10,
Neuron.NF shows Mx1 (+0.004), Mx2 (+0.012) and Ifi35 (+0.159) up, and the median
ISG logFC is ≈ 0 in the neuronal bins (−0.0017 to +0.0009), negative only in
SatelliteGlia (−0.195) and weakly negative in Endothelium (−0.068) and Schwann
(−0.009). The composite ISG score (all 39 genes) is lower in B16F10 than control
in 6 of 7 compartments with an MWU-testable result (`sc_GSE292522_celltype_summary.json`;
p = 0.0286 for SatelliteGlia, which equals the 4+4 floor 2/C(8,4)). The direction
of the claim is therefore correct in the data as analysed.

**【Why it matters】** This is the one compartment-level sentence that I could not
break — provided it is stated at q < 0.05 and scoped to the tested compartments.
Its two vulnerabilities are scope (Finding 1) and the bins' composition (Findings
3–4); neither flips its sign, and the ISG median-logFC pattern (flat-to-negative
everywhere, most negative in the glial bin) is a genuine, compartment-resolved
observation.

**【Specific fix】** State it with the threshold and scope: "At q < 0.05 no
interferon-stimulated gene was up-regulated in any of the eight tested
compartments (38 ISGs; the closest approach is Ifi35, +0.16 log₂ at q < 0.10 in
neurofilament-positive neurons); composite ISG scores were flat to lower in
tumour animals in every compartment with a testable composite (n ≥ 3 per group;
the erythrocyte bin retained only two control animals and carries no
inferential weight)."

---

## § Stands up

Things I suspected on first read of the code and then verified to be correct —
each with the evidence that convinced me:

1. **Pseudobulk construction is correct.** Orientation is genes × cells as
   loaded (`mmread(...).tocsc()`), per-compartment sums are taken over cells
   within each animal (`mat[:, mask].sum(axis=1)`), CP10K normalisation is
   per-sample, and the design matrix is intercept + group over animals. I
   rebuilt all compartment vectors independently from the raw matrices for two
   samples spanning both groups and compared against the analysis cache:
   **maximum absolute difference 0.0 counts** for every compartment in both
   samples. The inference unit is genuinely the animal; there is no cell-level
   test anywhere in 31 or 32.
2. **Pseudoreplication is genuinely avoided.** All tests in the per-compartment
   layer operate on 4+4 (or fewer) animal-level pseudobulks; the MWU
   discrimination tests are on per-animal scores (`32_percelltype_sc.py:189-198`),
   and I confirmed the reported p = 0.0286 for SatelliteGlia equals the exact
   two-sided floor 2/70 = 0.0286 for n = 4+4 — i.e. the test is at its resolution
   limit, not over-resolved.
3. **The deposited compartment numbers reproduce exactly.** My independent refits
   reproduce the deposited `celltype_DE.csv` values to rounding: Neuron/Gal
   +1.472, t 10.17, q 0.0001; Neuron/Scn10a −0.247, t −5.15, q 0.012; NF/Gal
   +2.535, t 18.14, q ≈ 0; SatelliteGlia/Isg15 −0.629, t −4.54, q 0.0071;
   SatelliteGlia/Ifit1 −0.523, t −4.74, q 0.0036; Endothelial/Lyz2 −0.813,
   q 0.0008. The brief's "Neuron.NF had 1,890 DE genes" reproduced at 1,886–1,888
   depending on symbol deduplication.
4. **Gal up-regulation is not a within-compartment normalisation artefact.** This
   was my prior going in (mandate item 5): a compartment normalised to its own
   total could show a spurious "up" if other transcripts fell. It is not the case
   here: mean raw Gal counts per neuron-bin cell rise from 4.4–11.1 (controls) to
   17.8–26.6 (B16F10), 2.3–3.0×, with the fraction of Gal-positive neurons rising
   from 0.29–0.40 to 0.42–0.66. The signal exists at the count level, before any
   CP10K division.
5. **The ≥5-cell / ≥6-animal filter is implemented exactly as §2.8 states.**
   `32_percelltype_sc.py:133` (`mask.sum() >= 5`) and `:149-150`
   (`len(present) < 6`). The Methods sentence is accurate about the rule; what is
   missing is disclosure of its consequences (Findings 1 and 6).
6. **Leave-one-out stability of the two strongest compartment claims.** NF/Gal:
   LOO logFC range +2.420 to +2.669, zero sign flips, Welch p < 0.0016 in every
   refit, limma q < 0.05 in 8/8. Neuron/Gal: LOO +1.307 to +1.600, 8/8 q < 0.05.
   Scn10a in both bins: zero sign flips, 8/8 q < 0.05. The injury-programme
   numbers the paper leads with are not one-animal effects.
7. **The downloader (30) is sound.** Range-request probe, refusal to stream
   without byte-range support, size verification, and full tar-member CRC
   verification — the integrity discipline in this script is above the standard
   for this kind of paper.

---

## § Questions for the authors

1. What, exactly, does "55,967 cells after quality filtering" (§2.8) refer to?
   Were the GEO-deposited matrices CellRanger-filtered output, and is any QC
   claimed by this paper performed anywhere in the code? If not, the sentence
   should be corrected rather than the code.
2. What is the basis for labelling the SatelliteGlia bin as satellite glia when
   its Gfap content is 7 counts in 361 cells? Do the authors have cluster-level
   or reference-mapping evidence (not in the repository) that these 25–108
   cells/animal are SGCs rather than Schwann-lineage cells?
3. Was the myeloid compartment's exclusion by the ≥5-cell rule a deliberate
   design choice or an oversight? If deliberate, why is it not reported; if
   oversight, would the authors re-run with a ≥3-cell threshold (5/8 to 7/8
   animals per group, at the cost of noisier pseudobulks) before reasserting the
   compartment-level ISG claim?
4. The "Neuron" bin's median cell carries 3–6 myelin-gene counts and its
   runner-up label is Schwann for ~70–78% of cells. Do the authors interpret
   these as neurons with ambient myelin RNA, neuron–glia doublets, or something
   else — and was any doublet/ambient assessment considered?
5. Why does `32_percelltype_sc.py` not deduplicate gene symbols when
   `31_process_sc.py` does, given both are described as using "exactly the same"
   pipeline logic?
6. The manuscript reports Ifit2-in-neurons and Isg15-in-satellite-glia as "the
   only significant ISG changes" while the deposited table contains
   Ifit1-in-satellite-glia at q = 0.0037. Was this row seen and excluded for a
   reason, or missed?
7. Per-compartment full DE fits (beyond the 27 FOCUS genes) are not deposited,
   yet compartment-wide statements ("no compartment…") rest on them. Will the
   full per-compartment fits and the per-compartment per-animal cell-count
   table be deposited?

---

## § What I actually checked

**Files read in full:** `scripts/30_fetch_sc.py`, `31_process_sc.py`,
`32_percelltype_sc.py`, `_limma.py`; `docs/manuscript/manuscript_final.md`
(Abstract, §2.8, §3.8, §4, §5, Table 5); `results/tables/sc_GSE292522_celltype_DE.csv`
(all 198 rows), `_celltype_scores.csv` (all 189 rows), `_celltype_summary.json`,
`sc_GSE292522_summary.json`, `sc_GSE292522_composition.csv`,
`negative_result_quantification.json`; `data/processed/GSE292522/*_features.tsv.gz`
(symbol columns, duplicates, marker-gene presence).

**Recomputed from the eight raw 10x matrices** (independent re-implementation of
the annotation, marker scoring, pseudobulk aggregation and limma-trend; limma
statistics verified against the project's own `_limma` implementation, which I
read line by line):

- Per-sample per-compartment assigned cell counts for all 12 defined
  compartments (Finding 1 table) — reproduces the analysis' own assignment
  counts exactly (e.g. Neuron 5,159/5,840/4,606/4,589/4,226/5,366/5,489/4,890).
- Per-compartment QC: median UMIs, median detected genes, quartiles,
  top1−top2 annotation margins, runner-up label distributions (Findings 3, 4, 6).
- Strict-neuron pseudobulk DE (Finding 4): Gal +0.781/q 0.353, Scn10a
  −0.135/q 0.59, 73/27,228 genes at q < 0.05.
- Full per-compartment DE refits (real `_limma`): reproduces every deposited
  focus-gene value I checked to ≤ 0.002 in logFC and ≤ 0.0002 in q; compartment
  q<0.05 counts: Neuron 1,558; Neuron.NF 1,886; Neuron.PE 349; Schwann 228;
  SatelliteGlia 24; Endothelial 21; Fibroblast 0; Erythrocyte 0.
- Exhaustive 39-gene ISG screen per compartment (Finding 9): zero ISG up at
  q < 0.05 anywhere; Ifi35/Mx1/Mx2 up at q < 0.10 in Neuron.NF; SatelliteGlia is
  the only compartment with a negative median ISG logFC (−0.195).
- Leave-one-animal-out for 18 compartment×gene headline pairs, two ways
  (Welch t; limma-trend refit on a fixed 3,000-gene universe with the project's
  `_limma`): Gal and Scn10a stable (8/8); SatelliteGlia Isg15 1/8, Ifit1 0/8;
  Ifit2-in-neurons sign-stable but Welch p rises to 0.128; Schwann/Endothelial
  Gal robust (8/8, 7/8).
- Per-cell Gal counts within the neuron bin per animal (normalisation-artefact
  check): Con 7.9/4.4/11.1/7.5 vs B16F10 26.6/17.8/21.6/25.2 mean counts/cell.
- Regression of compartment logFC on whole-DRG logFC: Neuron slope 0.856
  (ρ = 0.904), Neuron.NF slope 1.011 (ρ = 0.827), Schwann slope 0.428 (ρ = 0.400).
- Features-file forensics: 55,335 rows, 105 duplicated symbols (none in FOCUS or
  ISG sets), PgP9.5/Pgp9.5/Oct6/Blbp absent, all eight feature files identical
  and identically ordered.
- BH scope: within-compartment over all detected genes, from
  `32_percelltype_sc.py:160-161`; whole-DRG BH over 34,603 deduplicated genes
  from `31_process_sc.py:99-134`.
- "55,967 cells after quality filtering": equals the raw barcode sum; no QC code
  exists anywhere in 30/31/32.

**Discrepancies vs the manuscript found and stated above:** (i) "the only
significant ISG changes" omits Ifit1-SatelliteGlia (deposited q = 0.0037);
(ii) "neuron-intrinsic" contradicted by Gal in Schwann and Endothelium (both
q ≤ 0.006, LOO-robust); (iii) "CPM" is CP10K; (iv) "after quality filtering"
unsupported by code; (v) the compartment list in §2.8 miscounts and omits the
dropped myeloid/T-cell/NP/TH compartments; (vi) the erythrocyte compartment was
tested at 2 vs 4 animals, undisclosed; (vii) compartment q values are
BH-within-compartment, unstated.

**Not checked (stated for completeness):** identity of the argmax-Neuron bin's
non-neuronal majority beyond the pan-neuronal/myelin criteria used here (full
clustering annotation was out of scope for review); leave-one-out of the
compartment-wide q<0.05 gene counts (only headline genes were LOO-tested);
whether the GEO deposit's own metadata documents upstream QC.
