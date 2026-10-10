# -*- coding: utf-8 -*-
"""A3 independent recompute for the per-compartment (32_percelltype) analysis.

Independent of scripts/_limma.py and 32: re-implements the annotation +
pseudobulk + moderated-t pipeline from scratch (same maths, own code) so that
agreement is informative. Outputs JSON to stdout.
"""
import os, sys, json, gzip, math
import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy import stats
from scipy.special import digamma, polygamma

ROOT = ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))
SC = os.path.join(ROOT, "data", "processed", "GSE292522")

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
NEURON_SUB = {
    "NF": ["Nefh","Nefl","Prph","Map2","Tubb3"],
    "NP": ["P2rx3","Ret","Sst","Npy","Runx1"],
    "PE": ["Calca","Tac1","Trpv1","Scn10a","Scn9a"],
    "TH": ["Th","Dbh","Ddc"],
}
ISG = ["Ifit1","Ifit2","Ifit3","Isg15","Isg20","Mx1","Mx2","Oas1a","Oas1b","Oas2",
       "Oas3","Cxcl10","Ccl5","Stat1","Stat2","Irf7","Irf9","Ifih1","Ddx58","Rsad2",
       "Trim21","Bst2","Gbp2","Gbp3","Gbp4","Gbp5","Ifi44","Ifi47","Parp9","Usp18",
       "Slfn1","Slfn2","Slfn4","Slfn5","Sp100","Ifi35","Ifi204","Samhd1"]

def gz_lines(p):
    with gzip.open(p, "rt") as fh:
        return [ln.rstrip("\n").split("\t") for ln in fh]

def lowess_trend_s2(Y, design, frac=0.4):
    """limma-trend moderated statistics, own implementation (sm lowess)."""
    import statsmodels.api as sm
    n, p = design.shape
    df_res = n - p
    beta, *_ = np.linalg.lstsq(design, Y.T, rcond=None)
    resid = Y.T - design @ beta
    s2 = (resid ** 2).sum(axis=0) / df_res
    XtX_inv = np.linalg.pinv(design.T @ design)
    uns = math.sqrt(XtX_inv[1, 1])
    lfc = beta[1]
    Amean = Y.mean(axis=1)
    ok = np.isfinite(s2) & (s2 > 0)
    base = np.log(s2[ok])
    fitv = sm.nonparametric.lowess(base, Amean[ok], frac=frac, return_sorted=False)
    s2_scaled = np.exp(base - fitv)
    z = np.log(s2_scaled)
    e = z - digamma(df_res/2.0) + math.log(df_res/2.0)
    ebar = e.mean()
    var_e = ((e-ebar)**2).sum()/max(len(e)-1,1) - polygamma(1, df_res/2.0)
    if var_e <= 0:
        s20, df_prior = math.exp(ebar), 1e6
    else:
        lo, hi = 1e-4, 1e5
        f = lambda x: polygamma(1, x/2.0) - var_e
        for _ in range(200):
            mid = math.sqrt(lo*hi)
            if f(mid) > 0: lo = mid
            else: hi = mid
        df_prior = math.sqrt(lo*hi)
        s20 = math.exp(ebar + digamma(df_prior/2.0) - math.log(df_prior/2.0))
    s2_post = np.full_like(s2, np.nan)
    s2_post[ok] = (df_prior*s20 + df_res*s2[ok])/(df_prior+df_res)
    se = np.sqrt(np.maximum(s2_post,1e-12))*uns
    t = lfc/np.maximum(se,1e-12)
    pv = 2*stats.t.sf(np.abs(t), df_res+df_prior)
    return lfc, t, pv

def bh(p):
    p = np.asarray(p, float); n = len(p)
    order = np.argsort(p)
    q = p[order]*n/(np.arange(n)+1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n); out[order] = np.clip(q,0,1)
    return out

def fit(L, drop_idx=None):
    cols = list(L.columns)
    if drop_idx is not None:
        keep = [i for i in range(len(cols)) if i != drop_idx]
        L = L.iloc[:, keep]
    grp = np.array([1.0 if SAMPLES[c][0]=="B16F10" else 0.0 for c in L.columns])
    design = np.column_stack([np.ones(len(grp)), grp])
    Y = L.values.astype(float)
    lfc, t, pv = lowess_trend_s2(Y, design)
    q = bh(pv)
    return pd.DataFrame({"logFC": lfc, "t": t, "P": pv, "q": q}, index=L.index)

def main():
    out = {"marker_gene_presence": {}, "per_sample_compartment_cells": {},
           "zero_score_cells_assigned": {}, "feat_consistency": {}, "dup_symbols": None}
    comps = list(BROAD) + [f"Neuron.{s}" for s in NEURON_SUB]
    agg = {c: {} for c in comps}
    canon = None
    for prefix, (grp, gsm) in SAMPLES.items():
        feat = gz_lines(os.path.join(SC, gsm+"_features.tsv.gz"))
        sym = [r[1] for r in feat]
        if canon is None:
            canon = sym
            from collections import Counter
            dups = {k: v for k, v in Counter(sym).items() if v > 1}
            out["dup_symbols"] = {"n_rows": len(sym), "n_dup_symbols": len(dups),
                                  "examples": dict(list(dups.items())[:5])}
        out["feat_consistency"][prefix] = (sym == canon)
        mat = mmread(os.path.join(SC, gsm+"_matrix.mtx.gz")).tocsc()
        ng, nc = mat.shape
        idx_map = {g: i for i, g in enumerate(sym)}
        # marker presence
        if prefix == "Sample1_DRG_Con_rep1":
            for ct, gs in list(BROAD.items()) + list(NEURON_SUB.items()):
                out["marker_gene_presence"][ct] = {
                    "missing": [g for g in gs if g not in idx_map],
                    "n_used": sum(1 for g in gs if g in idx_map)}
        def mscore(gset):
            gi = [idx_map[g] for g in gset if g in idx_map]
            if not gi:
                return np.zeros(nc)
            return np.log1p(mat[gi, :].toarray().astype(np.float32)).mean(axis=0)
        Sb = np.vstack([mscore(g) for g in BROAD.values()])
        ab = np.array(list(BROAD))[Sb.argmax(axis=0)]
        maxb = Sb.max(axis=0)
        out["zero_score_cells_assigned"][prefix] = {
            "n_cells": int(nc),
            "n_allzero_max_score": int((maxb <= 0).sum()),
            "allzero_assigned_to": {ct: int(((maxb <= 0) & (ab == ct)).sum()) for ct in BROAD}}
        # per-compartment cell counts BEFORE the >=5 filter
        cc = {}
        for ct in BROAD:
            cc[ct] = int((ab == ct).sum())
        neu = ab == "Neuron"
        Ss = np.vstack([mscore(g) for g in NEURON_SUB.values()])
        asu = np.array(["none"]*nc, dtype=object)
        if neu.any():
            asu[neu] = np.array(list(NEURON_SUB))[Ss[:, neu].argmax(axis=0)]
        for s in NEURON_SUB:
            cc[f"Neuron.{s}"] = int((neu & (asu == s)).sum())
        out["per_sample_compartment_cells"][prefix] = cc
        # store pseudobulk vectors
        for ct in BROAD:
            m = ab == ct
            if m.sum() >= 5:
                agg[ct][prefix] = np.asarray(mat[:, m].sum(axis=1)).ravel()
        for s in NEURON_SUB:
            m = neu & (asu == s)
            if m.sum() >= 5:
                agg[f"Neuron.{s}"][prefix] = np.asarray(mat[:, m].sum(axis=1)).ravel()
        del mat
        print(f"loaded {prefix}", flush=True)

    # ---- per-compartment DE, full re-run ----
    res = {}
    for ctype in comps:
        present = [p for p in SAMPLES if p in agg[ctype]]
        if len(present) < 6:
            res[ctype] = {"n_samples": len(present), "dropped": True}
            continue
        M = pd.DataFrame({p: agg[ctype][p] for p in present}, index=canon)
        M = M.loc[M.sum(axis=1) > 0]
        cp = M.div(M.sum(axis=0), axis=1)*1e4 + 1.0
        L = np.log2(cp)
        fitr = fit(L)
        nq = int((fitr.q < 0.05).sum())
        entry = {"n_samples": len(present), "groups": {SAMPLES[p][0] for p in present},
                 "n_genes": int(len(fitr)), "n_q05": nq}
        # full-ISG check (all 39, not just FOCUS)
        isg = fitr.loc[[g for g in ISG if g in fitr.index]]
        entry["ISG_up_q05"] = [g for g in isg.index if isg.loc[g, "logFC"] > 0 and isg.loc[g, "q"] < 0.05]
        entry["ISG_dn_q05"] = [g for g in isg.index if isg.loc[g, "logFC"] < 0 and isg.loc[g, "q"] < 0.05]
        for g in ["Gal", "Scn10a", "Ifit1", "Ifit2", "Isg15", "Gfap", "Cd68"]:
            if g in fitr.index:
                r = fitr.loc[g]
                entry[g] = {"logFC": round(float(r.logFC), 3), "t": round(float(r.t), 2),
                            "q": round(float(r.q), 4)}
        entry_main = dict(entry)
        res[ctype] = entry_main
        out["compartment_DE"] = res
        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                               "_a3_recompute_partial.json"), "w") as fh:
            json.dump(out, fh, indent=1, default=str)
        # leave-one-out on key genes
        loo = {}
        for g in ["Gal", "Ifit1", "Isg15", "Scn10a"]:
            if g not in fitr.index:
                continue
            loo[g] = []
            for d in range(len(L.columns)):
                fl = fit(L, drop_idx=d)
                if g in fl.index:
                    r = fl.loc[g]
                    loo[g].append({"dropped": L.columns[d],
                                   "logFC": round(float(r.logFC), 3),
                                   "q": round(float(r.q), 4)})
        entry["LOO"] = loo
        res[ctype] = entry
        out["compartment_DE"] = res
        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                               "_a3_recompute_partial.json"), "w") as fh:
            json.dump(out, fh, indent=1, default=str)
        print(f"done {ctype}", flush=True)

    out["compartment_DE"] = res
    print(json.dumps(out, indent=1, default=str))

if __name__ == "__main__":
    main()
