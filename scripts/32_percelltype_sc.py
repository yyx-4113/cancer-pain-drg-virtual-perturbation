#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
32_percelltype_sc.py -- Stratified single-cell re-analysis of GSE292522.

Why this script exists
----------------------
The original 31_process_sc.py collapsed every cell into a single pseudobulk and
estimated composition by assigning each cell to the marker set with the single
highest mean log-expression. That answers "is there a bulk interferon programme?"
but it cannot say *where* a programme lives, and a whole-DRG average dilutes any
compartment-specific signal.

This script asks a sharper question: does the interferon / myeloid / neuronal-
injury programme appear in *any* DRG compartment, or is it absent everywhere?
It (1) annotates every cell into a broad type and, within neurons, into a
sensory subtype, using curated DRG marker sets (Usoskin 2015, ref [16]); and
(2) for each compartment builds a *per-sample* pseudobulk (sum of counts of the
cells of that type in that animal) and runs the SAME limma-trend model across
the 4 tumour vs 4 control animals -- so the unit of inference is the animal,
not the cell, which avoids single-cell pseudoreplication.

Deliverables
------------
- sc_GSE292522_celltype_DE.csv   : per cell-type limma-trend of focus genes
- sc_GSE292522_celltype_scores.csv : per-sample composite ISG / injury / glial scores
- sc_GSE292522_celltype_summary.json : fractions, discrimination tests
"""
from __future__ import annotations
import os, sys, json, gzip
import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.stats import mannwhitneyu
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _limma import fit_limma, bh

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SC = os.path.join(ROOT, "data", "processed", "GSE292522")
TAB = os.path.join(ROOT, "results", "tables")
os.makedirs(TAB, exist_ok=True)

SAMPLES = {
    "Sample1_DRG_Con_rep1":    ("Con",    "GSM8860632_Sample1_DRG_Con_rep1"),
    "Sample2_DRG_Con_rep2":    ("Con",    "GSM8860633_Sample2_DRG_Con_rep2"),
    "Sample3_DRG_Con_rep3":    ("Con",    "GSM8860634_Sample3_DRG_Con_rep3"),
    "Sample4_DRG_Con_rep4":    ("Con",    "GSM8860635_Sample4_DRG_Con_rep4"),
    "Sample5_DRG_B16F10_rep1": ("B16F10", "GSM8860636_Sample5_DRG_B16F10_rep1"),
    "Sample6_DRG_B16F10_rep2": ("B16F10", "GSM8860637_Sample6_DRG_B16F10_rep2"),
    "Sample7_DRG_B16F10_rep3": ("B16F10", "GSM8860638_Sample7_DRG_B16F10_rep3"),
    "Sample8_DRG_B16F10_rep4": ("B16F10", "GSM8860639_Sample8_DRG_B16F10_rep4"),
}

BROAD = {
    "Neuron":        ["Snap25","Syp","Tubb3","Nefl","Nefh","Th","Scn10a","Scn9a",
                     "Trpv1","Calca","Tac1","Npy","Gal","PgP9.5","Prph"],
    "SatelliteGlia": ["Gfap","Sox10","Slc1a3","Aqp4","Kcnj10","Glul","Blbp","Fgfr3"],
    "Schwann":       ["Mpz","Pmp22","Oct6","Mag","Plp1","Sox10","Dhh","Ncmap"],
    "Myeloid":       ["Cd68","Aif1","Itgam","Lyz2","Cx3cr1","Fcgr1","Cd14","Itgax",
                     "Ctss","Ly6c2","Ly6g","Spi1"],
    "TCell":         ["Cd3e","Cd3d","Cd4","Cd8a","Trbc1","Trac"],
    "Endothelial":   ["Pecam1","Cdh5","Cldn5","Vwf","Emcn"],
    "Fibroblast":    ["Col1a1","Col1a2","Pdgfra","Dcn","Lum","Col3a1"],
    "Erythrocyte":   ["Hbb-bs","Hbb-bt","Hba-a1","Hba-a2","Hbb-y"],
}
# neuronal subtypes (within Neuron only)
NEURON_SUB = {
    "NF":    ["Nefh","Nefl","Prph","Map2","Tubb3"],
    "NP":    ["P2rx3","Ret","Sst","Npy","Runx1"],          # non-peptidergic / nociceptor-ish
    "PE":    ["Calca","Tac1","Trpv1","Scn10a","Scn9a"],    # peptidergic
    "TH":    ["Th","Dbh","Ddc"],                            # tyrosine-hydroxylase
}
ISG = ["Ifit1","Ifit2","Ifit3","Isg15","Isg20","Mx1","Mx2","Oas1a","Oas1b","Oas2",
       "Oas3","Cxcl10","Ccl5","Stat1","Stat2","Irf7","Irf9","Ifih1","Ddx58","Rsad2",
       "Trim21","Bst2","Gbp2","Gbp3","Gbp4","Gbp5","Ifi44","Ifi47","Parp9","Usp18",
       "Slfn1","Slfn2","Slfn4","Slfn5","Sp100","Ifi35","Ifi204","Samhd1"]
INJURY = ["Sprr1a","Atf3","Gadd45a","Ecel1","Serpinb1a","Jun","Fos","Socs3"]
GLIAL = ["Gfap","Sox10","Slc1a3","Aqp4"]

FOCUS = ["Cd68","Aif1","Itgam","Lyz2","Cx3cr1","Gfap","Sox10","Mpz",
         "Snap25","Scn9a","Scn10a","Calca","Tac1","Gal","Il1b","Tnf","Ccl2",
         "Ifit1","Ifit2","Ifit3","Isg15","Oas1a","Oas1b","Bst2","Mx1","Stat1","Irf7"]


def gz_lines(path):
    with gzip.open(path, "rt") as fh:
        return [ln.rstrip("\n").split("\t") for ln in fh]


def load(prefix, gsm):
    feat = gz_lines(os.path.join(SC, gsm + "_features.tsv.gz"))
    sym = [r[1] if len(r) > 1 else r[0] for r in feat]
    mat = mmread(os.path.join(SC, gsm + "_matrix.mtx.gz")).tocsc()
    return mat, sym


def marker_score(mat_sym, idx_map, geneset):
    gi = [idx_map[g] for g in geneset if g in idx_map]
    if not gi:
        return None
    sub = mat_sym[gi, :].toarray().astype(np.float32)
    return np.log1p(sub).mean(axis=0)


def main():
    print("[32] loading + annotating 8 samples (aggregating per compartment) ...",
          flush=True)
    # aggregate per-compartment summed counts during load to avoid holding raw matrices
    COMPARTMENTS = list(BROAD) + [f"Neuron.{s}" for s in NEURON_SUB]
    agg = {c: {} for c in COMPARTMENTS}     # compartment -> prefix -> vector
    celltype_frac = {}
    canon_sym = None
    for prefix, (grp, gsm) in SAMPLES.items():
        mat, sym = load(prefix, gsm)
        if canon_sym is None:
            canon_sym = sym
        ng, nc = mat.shape
        S_broad = np.vstack([marker_score(mat, {g:i for i,g in enumerate(sym)}, gset)
                             for gset in BROAD.values()])      # 8 x cells
        assign_broad = np.array(list(BROAD))[S_broad.argmax(axis=0)]
        broad_frac = {ct: float((assign_broad == ct).mean()) for ct in BROAD}
        celltype_frac[prefix] = broad_frac
        neu_mask = assign_broad == "Neuron"
        S_sub = np.vstack([marker_score(mat, {g:i for i,g in enumerate(sym)}, gset)
                           for gset in NEURON_SUB.values()])   # 4 x cells
        assign_sub = np.array(["none"]*nc, dtype=object)
        if neu_mask.any():
            ss = S_sub[:, neu_mask]
            assign_sub[neu_mask] = np.array(list(NEURON_SUB))[ss.argmax(axis=0)]
        # sum counts per compartment (dense over gene symbols)
        for ctype in BROAD:
            mask = (assign_broad == ctype)
            if mask.sum() >= 5:
                agg[ctype][prefix] = np.asarray(mat[:, mask].sum(axis=1)).ravel()
        for s in NEURON_SUB:
            mask = (assign_broad == "Neuron") & (assign_sub == s)
            if mask.sum() >= 5:
                agg[f"Neuron.{s}"][prefix] = np.asarray(mat[:, mask].sum(axis=1)).ravel()
        del mat
        print(f"  {prefix} [{grp}] n={nc} neuron={broad_frac['Neuron']:.2f} "
              f"myeloid={broad_frac['Myeloid']:.4f} glia={broad_frac['SatelliteGlia']:.3f}",
              flush=True)

    # ---- per-sample x cell-type pseudobulk ----
    print("[32] per-compartment limma-trend across 4+4 animals ...", flush=True)
    de_rows = []
    score_records = []
    for ctype in COMPARTMENTS:
        present = [p for p in SAMPLES if p in agg[ctype]]
        if len(present) < 6:
            continue
        M = pd.DataFrame({p: agg[ctype][p] for p in present}, index=canon_sym)
        M = M.loc[M.sum(axis=1) > 0]
        cp10k = M.div(M.sum(axis=0), axis=1) * 1e4 + 1.0
        L = np.log2(cp10k)
        groups = [SAMPLES[p][0] for p in L.columns]
        design = np.column_stack([np.ones(len(groups)),
                  np.array([1.0 if g == "B16F10" else 0.0 for g in groups])])
        fit = fit_limma(L, design, coef_idx=1, trend=True)
        fit["q"] = bh(fit["P"].values)
        fit = fit.sort_values("P")
        for g in FOCUS:
            if g in fit.index:
                r = fit.loc[g]
                de_rows.append({"cell_type": ctype, "gene": g,
                                "logFC": round(float(r.logFC),3),
                                "t": round(float(r.t),2),
                                "q": round(float(r.q),4)})
        for gset_name, gset in [("ISG", ISG), ("injury", INJURY), ("glial", GLIAL)]:
            gs = [g for g in gset if g in L.index]
            if not gs:
                continue
            means = L.loc[gs].mean(axis=0)
            for p in L.columns:
                score_records.append({"cell_type": ctype, "score": gset_name,
                                       "sample": p, "group": SAMPLES[p][0],
                                       "value": round(float(means[p]),4)})
        print(f"  [{ctype}] n_samples={len(present)} DE genes(q<.05)={int((fit.q<0.05).sum())}",
              flush=True)

    de_df = pd.DataFrame(de_rows)
    de_df.to_csv(os.path.join(TAB, "sc_GSE292522_celltype_DE.csv"), index=False)

    sc_df = pd.DataFrame(score_records)
    sc_df.to_csv(os.path.join(TAB, "sc_GSE292522_celltype_scores.csv"), index=False)

    # ---- discrimination: does ANY compartment show ISG up? ----
    # per cell-type, MWU of ISG score (B16F10 vs Con)
    isg_disc = {}
    for ctype in sc_df.cell_type.unique():
        sub = sc_df[(sc_df.cell_type==ctype) & (sc_df.score=="ISG")]
        con = sub[sub.group=="Con"].value.values
        b16 = sub[sub.group=="B16F10"].value.values
        if len(con)>=3 and len(b16)>=3:
            mw = mannwhitneyu(b16, con, alternative="two-sided")
            isg_disc[ctype] = {"mean_Con": round(float(con.mean()),4),
                               "mean_B16F10": round(float(b16.mean()),4),
                               "MWU_p": round(float(mw.pvalue),4)}

    summary = {
        "dataset": "GSE292522 (10x scRNA-seq mouse DRG, B16F10 model)",
        "annotation": "marker-score assignment (Usoskin 2015 DRG markers, ref [15])",
        "inference_unit": "animal (per-sample pseudobulk per compartment), n=4+4",
        "celltype_ISG_discrimination": isg_disc,
        "note": "per-compartment limma-trend; see sc_GSE292522_celltype_DE.csv",
    }
    with open(os.path.join(TAB, "sc_GSE292522_celltype_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2, default=str)
    print(json.dumps(summary, indent=2, default=str)[:2000])
    print("DONE 32")

if __name__ == "__main__":
    main()
