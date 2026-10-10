#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
_A3_panel_verify_sc.py -- R3 reviewer A3 (single-cell / pseudobulk methodology)
independent verification of the GSE292522 claims in manuscript section 3.8.
Replicates the annotation logic of scripts/31 and scripts/32 FROM RAW DATA
and recomputes every compartment-level claim. Reads no other reviewer file.
"""
import gzip, os, json
import numpy as np
import pandas as pd
from scipy.io import mmread

ROOT = r"D:/2026.10/不同部位癌痛患者对药物敏感性生信分析/cancer_pain_virtual_perturbation"
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

# exactly as in scripts/32 (note: "PgP9.5" casing, "Blbp", "Oct6" aliases)
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
MARKERS31 = {
    "Neuron": ["Pgp9.5"],   # script 31 casing, to compare
}

def gz_lines(p):
    with gzip.open(p, "rt") as fh:
        return [ln.rstrip("\n").split("\t") for ln in fh]

rows_counts, rows_gfap, purity_rows, marker_presence = [], [], [], None
total = 0
for prefix, (grp, gsm) in SAMPLES.items():
    feat = gz_lines(os.path.join(SC, gsm + "_features.tsv.gz"))
    sym = [r[1] if len(r) > 1 else r[0] for r in feat]
    mat = mmread(os.path.join(SC, gsm + "_matrix.mtx.gz")).tocsc()
    ng, nc = mat.shape
    total += nc
    idx = {g: i for i, g in enumerate(sym)}
    if marker_presence is None:
        allm = sorted({g for s in list(BROAD.values()) + list(MARKERS31.values()) for g in s})
        marker_presence = {g: (g in idx) for g in allm}
    # marker scores exactly as script 32: mean log1p over present markers
    S = np.vstack([
        (lambda gi: np.log1p(mat[gi, :].toarray().astype(np.float32)).mean(axis=0)
         if gi else np.zeros(nc))([idx[g] for g in gs if g in idx])
        for gs in BROAD.values()])
    assign = np.array(list(BROAD))[S.argmax(axis=0)]
    counts = {ct: int((assign == ct).sum()) for ct in BROAD}
    rows_counts.append({"sample": prefix, "group": grp, "n_cells": nc, **counts})
    # SatelliteGlia bin: absolute Gfap expression
    sg_mask = assign == "SatelliteGlia"
    n_sg = int(sg_mask.sum())
    if n_sg:
        lib = np.asarray(mat[:, sg_mask].sum(axis=1)).ravel().sum()
        gfap_cnt = float(mat[idx["Gfap"], :][:, sg_mask].sum())
        gfap_det = int((mat[idx["Gfap"], :][:, sg_mask].toarray().ravel() > 0).sum())
        gfap_cpm = gfap_cnt / lib * 1e4 if lib else np.nan
        # mean log1p Gfap per cell in bin
        gfap_meanlog = float(np.log1p(mat[idx["Gfap"], :][:, sg_mask].toarray().ravel()).mean())
    else:
        gfap_cpm = gfap_meanlog = np.nan; gfap_det = 0
    rows_gfap.append({"sample": prefix, "group": grp, "n_SatGlia": n_sg,
                      "Gfap_counts": gfap_cnt, "Gfap_cells_detected": gfap_det,
                      "Gfap_CPM_in_bin": round(gfap_cpm, 4),
                      "Gfap_meanlog1p_per_cell": round(gfap_meanlog, 5)})
    # neuron-bin purity under the manuscript criterion
    neu_mask = assign == "Neuron"
    if neu_mask.any():
        pan = ["Snap25", "Tubb3", "Syp"]
        myelin = ["Mpz", "Pmp22", "Mag", "Plp1", "Mbp"]
        pan_pos = np.zeros(nc, bool)
        for g in pan:
            if g in idx:
                pan_pos |= (mat[idx[g], :].toarray().ravel() > 0)
        myl_neg = np.ones(nc, bool)
        for g in myelin:
            if g in idx:
                myl_neg &= ~(mat[idx[g], :].toarray().ravel() > 0)
        sel = neu_mask
        pur = float((pan_pos[sel] & myl_neg[sel]).mean())
        pur_any = float((pan_pos[sel]).mean())
        purity_rows.append({"sample": prefix, "n_neuron_bin": int(sel.sum()),
                            "purity_strict": round(pur, 4),
                            "pan_positive_only": round(pur_any, 4)})

cnt = pd.DataFrame(rows_counts)
print("=== total cells:", total)
cnt["Myeloid_frac"] = cnt.Myeloid / cnt.n_cells
print(cnt[["sample","group","n_cells","Neuron","SatelliteGlia","Schwann","Myeloid",
           "TCell","Endothelial","Fibroblast","Erythrocyte","Myeloid_frac"]].to_string(index=False))
print("\nMyeloid total across 8 animals:", int(cnt.Myeloid.sum()))
print("TCell total:", int(cnt.TCell.sum()))
print("Myeloid frac range: %.5f%% - %.5f%%" % (100*cnt.Myeloid_frac.min(), 100*cnt.Myeloid_frac.max()))
print("Animals with >=5 myeloid cells:", int((cnt.Myeloid >= 5).sum()))
print("Animals with >=5 T cells:", int((cnt.TCell >= 5).sum()))
print("Animals with >=5 erythrocyte cells:", int((cnt.Erythrocyte >= 5).sum()))
print("Animals with >=5 NP-candidate? n/a")

# cross-check against deposited composition.csv
comp = pd.read_csv(os.path.join(TAB, "sc_GSE292522_composition.csv")).set_index("sample")
maxd = 0.0
for _, r in cnt.iterrows():
    for ct in BROAD:
        d = abs(comp.loc[r["sample"], ct] - r[ct] / r["n_cells"])
        maxd = max(maxd, d)
print("max |deposited frac - recomputed frac|:", maxd)

# myeloid fraction statistics (per-animal)
from scipy.stats import mannwhitneyu
from statsmodels.stats.power import TTestIndPower
con = cnt.loc[cnt.group == "Con", "Myeloid_frac"].values
b16 = cnt.loc[cnt.group == "B16F10", "Myeloid_frac"].values
d_obs = (con.mean() - b16.mean()) / np.sqrt(((len(con)-1)*con.var(ddof=1)+(len(b16)-1)*b16.var(ddof=1))/(len(con)+len(b16)-2))
mw = mannwhitneyu(b16, con, alternative="two-sided")
mde = TTestIndPower().solve_power(effect_size=None, nobs1=4, alpha=0.05, power=0.8, ratio=1.0, alternative="two-sided")
print("\nmyeloid Cohen d = %.4f | MWU U1 = %s  p2 = %.4f | MDE d = %.4f" % (d_obs, mw.statistic, mw.pvalue, mde))
mw1 = mannwhitneyu(b16, con, alternative="less")
print("MWU one-sided less p = %.4f ; U stats: U1=%s U2=%s" % (mw1.pvalue, mw.statistic, mw1.statistic if hasattr(mw1,'statistic') else ''))

gf = pd.DataFrame(rows_gfap)
print("\n=== Gfap in SatelliteGlia bin ===")
print(gf.to_string(index=False))

pu = pd.DataFrame(purity_rows)
print("\n=== Neuron-bin purity (criterion variants) ===")
print(pu.to_string(index=False))
print("pooled strict purity:", round(pd.concat([pu]*1).purity_strict.mean(), 4))

print("\n=== marker genes absent from features ===")
print({g: v for g, v in marker_presence.items() if not v})
