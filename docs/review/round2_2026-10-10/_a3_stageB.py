# -*- coding: utf-8 -*-
"""A3 stage-B verification: annotation validity + compositional artefacts.

All numbers cited in the A3 review come from THIS script or _a3_recompute.py.
Outputs one JSON to stdout.
"""
import os, sys, json, gzip, math
import numpy as np
import pandas as pd
from scipy.io import mmread

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))
SC = os.path.join(ROOT, "data", "processed", "GSE292522")
TAB = os.path.join(ROOT, "results", "tables")

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
NEURON_SUB = {"NF": ["Nefh","Nefl","Prph","Map2","Tubb3"],
              "NP": ["P2rx3","Ret","Sst","Npy","Runx1"],
              "PE": ["Calca","Tac1","Trpv1","Scn10a","Scn9a"],
              "TH": ["Th","Dbh","Ddc"]}

def gz_lines(p):
    with gzip.open(p, "rt") as fh:
        return [ln.rstrip("\n").split("\t") for ln in fh]

def mscore(mat, im, gs):
    gi = [im[g] for g in gs if g in im]
    return np.log1p(mat[gi, :].toarray().astype(np.float32)).mean(axis=0)

def main():
    out = {"qc_per_compartment": {}, "satglia_markers": {}, "neuron_purity": {},
           "margin_summary": {}, "gal_per_cell": {}, "global_resid": {}}
    # whole-DRG logFC from deposited pseudobulk DE (for slope comparison)
    pbde = pd.read_csv(os.path.join(TAB, "sc_GSE292522_pseudobulk_DE.csv"), index_col=0)
    lfc_bulk = pbde["logFC"]
    comp_lfc = {c: {} for c in ["Neuron", "Neuron.NF", "Schwann"]}
    agg = {c: {} for c in list(BROAD) + [f"Neuron.{s}" for s in NEURON_SUB]}
    canon = None
    for prefix, (grp, gsm) in SAMPLES.items():
        sym = [r[1] for r in gz_lines(os.path.join(SC, gsm+"_features.tsv.gz"))]
        if canon is None:
            canon = sym
        mat = mmread(os.path.join(SC, gsm+"_matrix.mtx.gz")).tocsc()
        im = {g: i for i, g in enumerate(sym)}
        nc = mat.shape[1]
        S = np.vstack([mscore(mat, im, gs) for gs in BROAD.values()])
        order = np.argsort(-S, axis=0)
        top1, top2 = order[0], order[1]
        assign = np.array(list(BROAD))[top1]
        marg = S[top1, np.arange(nc)] - S[top2, np.arange(nc)]
        umi = np.asarray(mat.sum(axis=0)).ravel()
        det = np.asarray((mat > 0).sum(axis=0)).ravel()
        # per-compartment QC
        q = {}
        for ct in list(BROAD):
            m = assign == ct
            if m.sum() == 0:
                q[ct] = {"n": 0}
                continue
            q[ct] = {"n": int(m.sum()),
                     "median_umi": float(np.median(umi[m])),
                     "median_genes": float(np.median(det[m])),
                     "p25_genes": float(np.percentile(det[m], 25)),
                     "p75_genes": float(np.percentile(det[m], 75)),
                     "median_margin": float(np.median(marg[m]))}
        neu = assign == "Neuron"
        Ssub = np.vstack([mscore(mat, im, gs) for gs in NEURON_SUB.values()])
        asu = np.array(["none"]*nc, dtype=object)
        asu[neu] = np.array(list(NEURON_SUB))[Ssub[:, neu].argmax(axis=0)]
        for s in NEURON_SUB:
            m = neu & (asu == s)
            q[f"Neuron.{s}"] = {"n": int(m.sum()),
                                "median_umi": float(np.median(umi[m])) if m.sum() else None,
                                "median_genes": float(np.median(det[m])) if m.sum() else None,
                                "median_margin": float(np.median(marg[m])) if m.sum() else None}
        out["qc_per_compartment"][prefix] = q
        # neuron purity proxy: pan-neuronal high & myelin-low cells
        pan = np.asarray(mat[[im[g] for g in ("Snap25","Tubb3","Syp") if g in im], :].sum(axis=0)).ravel()
        mye = np.asarray(mat[[im[g] for g in ("Mpz","Pmp22","Plp1","Mag") if g in im], :].sum(axis=0)).ravel()
        true_neu = (pan >= 20) & (mye < 5)
        land = pd.Series(assign[true_neu]).value_counts().to_dict()
        out["neuron_purity"][prefix] = {
            "n_true_neuron_cells": int(true_neu.sum()),
            "n_argmax_Neuron": int(neu.sum()),
            "true_neu_in_argmax_Neuron": int((assign[true_neu] == "Neuron").sum()),
            "true_neu_land": {k: int(v) for k, v in land.items()},
            "argmax_Neuron_median_genes": q["Neuron"]["median_genes"]}
        # SatGlia marker counts inside the compartment
        m_sg = assign == "SatelliteGlia"
        gf = {}
        for g in ("Gfap","Slc1a3","Aqp4","Sox10"):
            if g in im:
                gf[g] = int(mat[im[g], :][:, m_sg].sum())
        out["satglia_markers"][prefix] = {"n_cells": int(m_sg.sum()), "counts_in_compartment": gf}
        # margin summary
        out["margin_summary"][prefix] = {
            "frac_margin_lt_0p1": float((marg < 0.1).mean()),
            "frac_margin_lt_0p3": float((marg < 0.3).mean()),
            "median_margin_by_type": {ct: float(np.median(marg[assign == ct]))
                                      for ct in BROAD if (assign == ct).sum() > 0}}
        # Gal per-cell mean counts in argmax-Neuron cells
        if "Gal" in im:
            gal = np.asarray(mat[im["Gal"], :][:, neu].todense()).ravel()
            out["gal_per_cell"][prefix] = {
                "group": grp, "n_neuron_cells": int(neu.sum()),
                "mean_gal_counts_per_neuron_cell": float(gal.mean()),
                "median_gal": float(np.median(gal)),
                "frac_gal_positive_neurons": float((gal > 0).mean())}
        # pseudobulk vectors for slope comparison
        for ct in ["Neuron", "Neuron.NF", "Schwann"]:
            m = (assign == ct) if ct in BROAD else (neu & (asu == ct.split(".")[1]))
            if m.sum() >= 5:
                agg[ct][prefix] = np.asarray(mat[:, m].sum(axis=1)).ravel()
        del mat
        print(f"done {prefix}", flush=True)

    # slope of compartment logFC vs whole-DRG logFC
    for ct, d in agg.items():
        pres = [p for p in SAMPLES if p in d]
        if len(pres) < 6:
            continue
        M = pd.DataFrame({p: d[p] for p in pres}, index=canon)
        M = M.groupby(M.index).sum()
        M = M.loc[M.sum(axis=1) > 0]
        L = np.log2(M.div(M.sum(axis=0), axis=1)*1e4 + 1.0)
        lfc = L.loc[:, [p for p in pres if SAMPLES[p][0] == "B16F10"]].mean(axis=1) - \
              L.loc[:, [p for p in pres if SAMPLES[p][0] == "Con"]].mean(axis=1)
        shared = [g for g in L.index if g in lfc_bulk.index]
        x = lfc_bulk.loc[shared].values
        y = lfc.loc[shared].values
        slope = float(np.polyfit(x, y, 1)[0])
        rho = float(pd.Series(x).corr(pd.Series(y), method="spearman"))
        out["global_resid"][ct] = {"n_genes": len(shared), "slope_vs_bulk_lfc": round(slope, 4),
                                   "spearman_rho_vs_bulk": round(rho, 4)}

    # ---- npz cross-validation: recomputed vectors vs cached ones ----
    npz_path = os.path.join(ROOT, "results", "review_a3", "compartment_pseudobulk.npz")
    if os.path.exists(npz_path):
        z = np.load(npz_path, allow_pickle=True)
        nsym = pd.Series(canon)
        checks = {}
        for prefix in ["Sample1_DRG_Con_rep1", "Sample7_DRG_B16F10_rep3"]:
            sym1 = [r[1] for r in gz_lines(os.path.join(SC, SAMPLES[prefix][1]+"_features.tsv.gz"))]
            mat = mmread(os.path.join(SC, SAMPLES[prefix][1]+"_matrix.mtx.gz")).tocsc()
            im = {g: i for i, g in enumerate(sym1)}
            nc = mat.shape[1]
            S = np.vstack([mscore(mat, im, gs) for gs in BROAD.values()])
            assign = np.array(list(BROAD))[S.argmax(axis=0)]
            neu = assign == "Neuron"
            Ssub = np.vstack([mscore(mat, im, gs) for gs in NEURON_SUB.values()])
            asu = np.array(["none"]*nc, dtype=object)
            asu[neu] = np.array(list(NEURON_SUB))[Ssub[:, neu].argmax(axis=0)]
            vec = {}
            for ct in BROAD:
                m = assign == ct
                if m.sum() >= 5:
                    vec[ct] = np.asarray(mat[:, m].sum(axis=1)).ravel()
            for s in NEURON_SUB:
                m = neu & (asu == s)
                if m.sum() >= 5:
                    vec[f"Neuron.{s}"] = np.asarray(mat[:, m].sum(axis=1)).ravel()
            diffs = {}
            for ct, v in vec.items():
                key = f"{ct}||{prefix}"
                if key in z.files:
                    cached = z[key]
                    if len(cached) == len(v):
                        d = float(np.abs(cached - v).max())
                    else:
                        d = f"len mismatch {len(cached)} vs {len(v)}"
                else:
                    d = "missing_in_npz"
                diffs[ct] = d
            checks[prefix] = diffs
            del mat
        out["npz_crosscheck_maxabsdiff"] = checks

    # ---- full DE reproduction per compartment with the real _limma ----
    sys.path.insert(0, os.path.join(ROOT, "scripts"))
    from _limma import fit_limma, bh
    z2 = np.load(npz_path, allow_pickle=True) if os.path.exists(npz_path) else None
    sym_all = pd.Series(canon)
    def vecs_from_npz(ct):
        pres = [p for p in SAMPLES if z2 is not None and f"{ct}||{p}" in z2.files]
        return pres, {p: pd.Series(z2[f"{ct}||{p}"]).groupby(sym_all.values).sum().values
                      for p in pres} if pres else (None, None)
    de_repro = {}
    for ct in ["Neuron", "SatelliteGlia", "Schwann", "Endothelial", "Fibroblast",
               "Erythrocyte", "Neuron.NF", "Neuron.PE"]:
        if z2 is None:
            break
        pres = [p for p in SAMPLES if f"{ct}||{p}" in z2.files]
        if len(pres) < 6:
            de_repro[ct] = {"n_samples": len(pres), "dropped": True}
            continue
        df = pd.DataFrame({p: pd.Series(z2[f"{ct}||{p}"]).groupby(sym_all.values).sum()
                           for p in pres}, index=canon)
        df = df.loc[df.sum(axis=1) > 0]
        L = np.log2(df.div(df.sum(axis=0), axis=1)*1e4 + 1.0)
        g = [SAMPLES[p][0] for p in pres]
        des = np.column_stack([np.ones(len(g)), np.array([1.0 if x == "B16F10" else 0.0 for x in g])])
        f = fit_limma(L, des, coef_idx=1, trend=True)
        f["q"] = bh(f["P"].values)
        entry = {"n_samples": len(pres), "groups_present": sorted(set(g)),
                 "n_genes": int(len(f)), "n_q05": int((f.q < 0.05).sum())}
        for gene in ["Gal", "Scn10a", "Ifit1", "Ifit2", "Isg15", "Cx3cr1", "Gfap", "Il1b", "Lyz2", "Snap25"]:
            if gene in f.index:
                r = f.loc[gene]
                entry[gene] = {"logFC": round(float(r.logFC), 3), "t": round(float(r.t), 2),
                               "q": round(float(r.q), 4)}
        de_repro[ct] = entry
        print(f"DE repro done {ct}: n_q05={entry.get('n_q05')}", flush=True)
    out["de_reproduction"] = de_repro
    print(json.dumps(out, indent=1))

if __name__ == "__main__":
    main()
