#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
31_process_sc.py  --  Independent single-cell validation of the GSE149648
(bulk, Walker-256 tibia) findings using GSE292522 (10x scRNA-seq of mouse
DRG: 4 DRG_Con vs 4 DRG_B16F10 melanoma-inoculated).

Why this dataset
----------------
GSE149648 (our discovery) is a *bulk* DRG profile of Walker-256 tibial cancer.
Its dominant day-7 signal was a type-I interferon programme with myeloid
infiltration and satellite gliosis. GSE292522 is an *independent* bone-cancer
pain model (B16F10 melanoma in the femur/tibia) profiled at single-cell
resolution. If the same interferon / myeloid programme appears here, the bulk
finding is not a Walker-256- or batch-specific artefact. The single-cell layer
additionally lets us ask whether "myeloid infiltration" is real cellular
enrichment rather than a compositional read-out.

Three deliverables
------------------
1. Pseudobulk DE (B16F10 n=4 vs Con n=4) via limma-trend -- exactly the same
   model as the discovery pipeline -- on the summed per-sample counts.
2. Convergent validation: Spearman(B16F10-vs-Con logFC, Walker-256 D7 logFC)
   over shared genes. High agreement = the bulk programme reproduces.
3. Cell-type composition from canonical marker scoring (no clustering), to test
   whether the myeloid marker up-regulation reflects an increased myeloid
   fraction.
"""
from __future__ import annotations
import os, sys, json, gzip
import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy import sparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _limma import fit_limma, bh, muscle_score

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SC = os.path.join(ROOT, "data", "processed", "GSE292522")
TAB = os.path.join(ROOT, "results", "tables")
os.makedirs(TAB, exist_ok=True)

# ---- sample map (parsed from GSE292522_series_matrix.txt / filenames) ----
SAMPLES = {
    "Sample1_DRG_Con_rep1":    ("Con",      "GSM8860632_Sample1_DRG_Con_rep1"),
    "Sample2_DRG_Con_rep2":    ("Con",      "GSM8860633_Sample2_DRG_Con_rep2"),
    "Sample3_DRG_Con_rep3":    ("Con",      "GSM8860634_Sample3_DRG_Con_rep3"),
    "Sample4_DRG_Con_rep4":    ("Con",      "GSM8860635_Sample4_DRG_Con_rep4"),
    "Sample5_DRG_B16F10_rep1": ("B16F10",   "GSM8860636_Sample5_DRG_B16F10_rep1"),
    "Sample6_DRG_B16F10_rep2": ("B16F10",   "GSM8860637_Sample6_DRG_B16F10_rep2"),
    "Sample7_DRG_B16F10_rep3": ("B16F10",   "GSM8860638_Sample7_DRG_B16F10_rep3"),
    "Sample8_DRG_B16F10_rep4": ("B16F10",   "GSM8860639_Sample8_DRG_B16F10_rep4"),
}

# curated mouse DRG cell-type marker sets (gene symbols)
MARKERS = {
    "Neuron":        ["Snap25","Syp","Tubb3","Nefl","Nefh","Th","Scn10a","Scn9a",
                     "Trpv1","Calca","Tac1","Npy","Gal","Pgp9.5","Prph"],
    "SatelliteGlia": ["Gfap","Sox10","Slc1a3","Aqp4","Kcnj10","Glul","Blbp","Fgfr3"],
    "Schwann":       ["Mpz","Pmp22","Oct6","Mag","Plp1","Sox10","Dhh","Ncmap"],
    "Myeloid":       ["Cd68","Aif1","Itgam","Lyz2","Cx3cr1","Fcgr1","Cd14","Itgax",
                     "Ctss","Ly6c2","Ly6g","Spi1"],
    "TCell":         ["Cd3e","Cd3d","Cd4","Cd8a","Trbc1","Trac"],
    "Endothelial":   ["Pecam1","Cdh5","Cldn5","Vwf","Emcn"],
    "Fibroblast":    ["Col1a1","Col1a2","Pdgfra","Dcn","Lum","Col3a1"],
    "Erythrocyte":   ["Hbb-bs","Hbb-bt","Hba-a1","Hba-a2","Hbb-y"],
}
# interferon-stimulated gene (ISG) signature (mouse)
ISG = ["Ifit1","Ifit2","Ifit3","Isg15","Isg20","Mx1","Mx2","Oas1a","Oas1b","Oas2",
       "Oas3","Cxcl10","Ccl5","Stat1","Stat2","Irf7","Irf9","Ifih1","Ddx58","Rsad2",
       "Trim21","Bst2","Gbp2","Gbp3","Gbp4","Gbp5","Ifi44","Ifi47","Parp9","Usp18",
       "Slfn1","Slfn2","Slfn4","Slfn5","Sp100","Ifi35","Ifi204","Samhd1"]


def gz_read_lines(path):
    with gzip.open(path, "rt") as fh:
        return [ln.rstrip("\n").split("\t") for ln in fh]


def load_sample(prefix, gsm):
    feat = gz_read_lines(os.path.join(SC, gsm + "_features.tsv.gz"))
    sym = [r[1] if len(r) > 1 else r[0] for r in feat]
    mat = mmread(os.path.join(SC, gsm + "_matrix.mtx.gz")).tocsc()
    return mat, sym


def main():
    print("[31] loading 8 samples and building pseudobulk ...", flush=True)
    pbulk = {}          # sample -> pandas Series(gene symbol -> summed counts)
    compo = {}          # sample -> dict(celltype -> fraction of cells)
    ncells = {}
    # marker symbols present
    all_markers = sorted({g for s in MARKERS.values() for g in s})
    for prefix, (grp, gsm) in SAMPLES.items():
        mat, sym = load_sample(prefix, gsm)
        ng, nc = mat.shape
        ncells[prefix] = nc
        # pseudobulk
        s = pd.Series(np.asarray(mat.sum(axis=1)).ravel(), index=sym, dtype=float)
        s = s.groupby(s.index, axis=0).sum()
        pbulk[prefix] = s
        # cell-type composition on marker genes only (log1p per cell)
        mset = [g for g in all_markers if g in set(sym)]
        idx = [sym.index(g) for g in mset]
        sub = mat[idx, :].toarray().astype(np.float32)   # markers x cells
        sub = np.log1p(sub)
        # score per cell = mean log1p over each marker set
        scores = {}
        for ct, genes in MARKERS.items():
            gi = [mset.index(g) for g in genes if g in mset]
            scores[ct] = sub[gi, :].mean(axis=0) if gi else np.zeros(nc)
        S = np.vstack([scores[ct] for ct in MARKERS])     # ntype x cells
        assign = np.argmax(S, axis=0)
        vals = S.max(axis=0)
        frac = {}
        for i, ct in enumerate(MARKERS):
            frac[ct] = float((assign == i).mean())
        compo[prefix] = frac
        print(f"  {prefix} [{grp}] cells={nc} myeloidfrac={frac['Myeloid']:.3f} "
              f"neuronfrac={frac['Neuron']:.3f}", flush=True)

    # ---- pseudobulk matrix ----
    PB = pd.DataFrame(pbulk)            # gene x sample
    PB = PB.loc[PB.sum(axis=1) > 0]
    # log2(counts per 10,000 + 1)
    cp10k = PB.div(PB.sum(axis=0), axis=1) * 1e4 + 1.0
    L = np.log2(cp10k)
    groups = [SAMPLES[c][0] for c in PB.columns]
    design = np.column_stack([np.ones(8),
                              np.array([1.0 if g == "B16F10" else 0.0 for g in groups])])
    fit = fit_limma(L, design, coef_idx=1, trend=True)
    fit["q"] = bh(fit["P"].values)
    fit = fit.sort_values("P")
    fit.to_csv(os.path.join(TAB, "sc_GSE292522_pseudobulk_DE.csv"))

    # ---- focused: interferon + myeloid markers ----
    focus = ISG + ["Cd68","Aif1","Itgam","Lyz2","Cx3cr1","Gfap","Sox10","Mpz",
                  "Snap25","Scn9a","Calca","Tac1","Gal","Il1b","Tnf","Ccl2"]
    focus = [g for g in focus if g in fit.index]
    foc = fit.loc[focus, ["logFC","t","P","q"]].copy()
    foc = foc.round(3)
    foc.to_csv(os.path.join(TAB, "sc_GSE292522_focus_markers.csv"))

    # ---- convergent validation vs bulk Walker-256 D7 ----
    bulk_path = os.path.join(TAB, "DE_D7_vs_Sham_adj.csv")
    conv = None
    if os.path.exists(bulk_path):
        bk = pd.read_csv(bulk_path)
        gcol = "gene" if "gene" in bk.columns else bk.columns[0]
        bk = bk.set_index(gcol)
        shared = [g for g in L.index if g in bk.index]
        # proper Spearman between B16F10-vs-Con pseudobulk logFC and bulk D7 logFC
        lfc_sc = L.loc[shared, [c for c in PB.columns if SAMPLES[c][0]=="B16F10"]].mean(axis=1) \
                 - L.loc[shared, [c for c in PB.columns if SAMPLES[c][0]=="Con"]].mean(axis=1)
        bklfc = bk.loc[shared, "logFC"] if "logFC" in bk.columns else None
        if bklfc is not None:
            from scipy.stats import spearmanr
            r, p = spearmanr(lfc_sc.values, bklfc.values)
            conv = {"n_shared_genes": len(shared),
                    "spearman_rho": float(r), "spearman_p": float(p),
                    "shared_with_bulk_D7": len(shared)}
    # ---- composition summary ----
    comp_df = pd.DataFrame(compo).T
    comp_df.index.name = "sample"
    comp_df["group"] = [SAMPLES[i][0] for i in comp_df.index]
    comp_df.to_csv(os.path.join(TAB, "sc_GSE292522_composition.csv"))

    # myeloid enrichment test across 4+4 samples (Wilcoxon on per-sample fraction)
    from scipy.stats import mannwhitneyu
    my = comp_df.groupby("group")["Myeloid"].mean()
    # two-sided is the primary report; "less" is the directional test that
    # matches the biology (B16F10 myeloid fraction is lower). "greater" would
    # test the opposite direction and is not reported.
    mw_two = mannwhitneyu(comp_df.loc[comp_df.group=="B16F10","Myeloid"],
                          comp_df.loc[comp_df.group=="Con","Myeloid"], alternative="two-sided")
    mw_less = mannwhitneyu(comp_df.loc[comp_df.group=="B16F10","Myeloid"],
                           comp_df.loc[comp_df.group=="Con","Myeloid"], alternative="less")

    summary = {
        "dataset": "GSE292522 (10x scRNA-seq mouse DRG, B16F10 melanoma model)",
        "n_samples": 8, "groups": {"Con": 4, "B16F10": 4},
        "total_cells": int(sum(ncells.values())),
        "cells_per_sample": ncells,
        "pseudobulk_DE": {
            "n_genes": int(len(fit)),
            "n_up_B16F10_q05": int(((fit.logFC > 0) & (fit.q < 0.05)).sum()),
            "n_down_B16F10_q05": int(((fit.logFC < 0) & (fit.q < 0.05)).sum()),
        },
        "ISG_up_in_B16F10": [g for g in ISG if g in fit.index and fit.loc[g,"logFC"]>0 and fit.loc[g,"q"]<0.05],
        "myeloid_fraction_mean": {"Con": float(my["Con"]), "B16F10": float(my["B16F10"])},
        "myeloid_MWU_p_twosided": float(mw_two.pvalue),
        "myeloid_MWU_p_less": float(mw_less.pvalue),
        "convergent_validation_vs_bulk_D7": conv,
    }
    with open(os.path.join(TAB, "sc_GSE292522_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2, default=str)
    print(json.dumps(summary, indent=2, default=str))
    print("DONE")


if __name__ == "__main__":
    main()
