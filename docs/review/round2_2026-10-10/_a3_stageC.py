# -*- coding: utf-8 -*-
"""A3 stage-C: decisive test -- does the neuronal claim survive a clean
neuron-only pseudobulk built from a STRICT, independent criterion?
strict neuron = (Snap25+Tubb3+Syp >= 20 counts) AND (Mpz+Pmp22+Plp1+Mag < 5).
Also profiles the argmax-Neuron bin (runner-up sets, pan-neuronal content).
"""
import os, sys, json, gzip
import numpy as np, pandas as pd
from scipy.io import mmread

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))
SC = os.path.join(ROOT, "data", "processed", "GSE292522")
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from _limma import fit_limma, bh

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

def gz_lines(p):
    with gzip.open(p, "rt") as fh:
        return [ln.rstrip("\n").split("\t") for ln in fh]

def mscore(mat, im, gs):
    gi = [im[g] for g in gs if g in im]
    return np.log1p(mat[gi, :].toarray().astype(np.float32)).mean(axis=0)

def main():
    names = list(BROAD)
    strict_agg, canon = {}, None
    bin_profile = {}
    for prefix, (grp, gsm) in SAMPLES.items():
        sym = [r[1] for r in gz_lines(os.path.join(SC, gsm+"_features.tsv.gz"))]
        if canon is None:
            canon = sym
        mat = mmread(os.path.join(SC, gsm+"_matrix.mtx.gz")).tocsc()
        im = {g: i for i, g in enumerate(sym)}
        nc = mat.shape[1]
        S = np.vstack([mscore(mat, im, gs) for gs in BROAD.values()])
        assign = np.array(names)[S.argmax(axis=0)]
        pan = np.asarray(mat[[im[g] for g in ("Snap25","Tubb3","Syp")], :].sum(axis=0)).ravel()
        mye = np.asarray(mat[[im[g] for g in ("Mpz","Pmp22","Plp1","Mag") if g in im], :].sum(axis=0)).ravel()
        strict = (pan >= 20) & (mye < 5)
        # profile the argmax-Neuron bin
        nb = assign == "Neuron"
        order = np.argsort(-S, axis=0)
        runner = np.array(names)[order[1]]
        rn = pd.Series(runner[nb]).value_counts(normalize=True)
        bin_profile[prefix] = {
            "n_bin": int(nb.sum()),
            "frac_with_pan_ge20": float((pan[nb] >= 20).mean()),
            "median_pan_counts": float(np.median(pan[nb])),
            "median_mye_counts": float(np.median(mye[nb])),
            "runner_up_frac": {k: round(float(v), 4) for k, v in rn.head(3).items()}}
        strict_agg[prefix] = np.asarray(mat[:, strict].sum(axis=1)).ravel()
        del mat
        print(f"done {prefix} strict_n={int(strict.sum())}", flush=True)

    df = pd.DataFrame(strict_agg, index=canon)
    df = df.loc[df.sum(axis=1) > 0]
    L = np.log2(df.div(df.sum(axis=0), axis=1)*1e4 + 1.0)
    g = [SAMPLES[c][0] for c in L.columns]
    des = np.column_stack([np.ones(len(g)), np.array([1.0 if x == "B16F10" else 0.0 for x in g])])
    f = fit_limma(L, des, coef_idx=1, trend=True)
    f["q"] = bh(f["P"].values)
    res = {"n_samples": 8, "n_genes": int(len(f)), "n_q05": int((f.q < 0.05).sum()),
           "cells_per_sample": {p: int((pd.read_csv(os.devnull) if False else 0)) for p in []}}
    for gene in ["Gal","Scn10a","Snap25","Nefl","Ifit2","Ifit1","Isg15","Cx3cr1","Gfap","Il1b","Cd68","Aif1"]:
        if gene in f.index:
            r = f.loc[gene]
            res[gene] = {"logFC": round(float(r.logFC),3), "t": round(float(r.t),2), "q": round(float(r.q),4)}
    out = {"strict_neuron_DE": res, "argmax_neuron_bin_profile": bin_profile}
    print(json.dumps(out, indent=1))

if __name__ == "__main__":
    main()
