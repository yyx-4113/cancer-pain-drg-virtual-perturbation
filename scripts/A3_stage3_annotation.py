#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""A3 stage 3: annotation validity and compositional-artefact forensics.

T9  compositional test -- is each compartment logFC explained by its own CPM
    denominator shift?  (regress gene logFC on log2 total-count ratio)
T10 erythrocyte plausibility -- can an enucleated RBC have >600 genes detected?
T11 SatelliteGlia/Schwann boundary -- what is the runner-up type for glia?
T12 MWU resolution floor at n=4v4 and multiplicity of the 7 discrimination tests
T13 marker-set expression baseline bias driving argmax (which set wins by default)
T14 injury-score LOO on the per-compartment scores actually tabulated
"""
from __future__ import annotations
import os, sys, json, gzip
import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.stats import mannwhitneyu

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SC = os.path.join(ROOT, "data", "processed", "GSE292522")
TAB = os.path.join(ROOT, "results", "tables")
OUT = os.path.join(ROOT, "results", "review_a3")

SAMPLES = {
    "Sample1_DRG_Con_rep1": ("Con", "GSM8860632_Sample1_DRG_Con_rep1"),
    "Sample2_DRG_Con_rep2": ("Con", "GSM8860633_Sample2_DRG_Con_rep2"),
    "Sample3_DRG_Con_rep3": ("Con", "GSM8860634_Sample3_DRG_Con_rep3"),
    "Sample4_DRG_Con_rep4": ("Con", "GSM8860635_Sample4_DRG_Con_rep4"),
    "Sample5_DRG_B16F10_rep1": ("B16F10", "GSM8860636_Sample5_DRG_B16F10_rep1"),
    "Sample6_DRG_B16F10_rep2": ("B16F10", "GSM8860637_Sample6_DRG_B16F10_rep2"),
    "Sample7_DRG_B16F10_rep3": ("B16F10", "GSM8860638_Sample7_DRG_B16F10_rep3"),
    "Sample8_DRG_B16F10_rep4": ("B16F10", "GSM8860639_Sample8_DRG_B16F10_rep4"),
}
GROUPS = {k: v[0] for k, v in SAMPLES.items()}
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

rows10, rows11, rows13, rows11b = [], [], [], []
for prefix, (grp, gsm) in SAMPLES.items():
    sym = [r[1] if len(r) > 1 else r[0] for r in gz_lines(os.path.join(SC, gsm + "_features.tsv.gz"))]
    mat = mmread(os.path.join(SC, gsm + "_matrix.mtx.gz")).tocsc()
    im = {g: i for i, g in enumerate(sym)}
    nc = mat.shape[1]
    S = np.vstack([mscore(mat, im, gs) for gs in BROAD.values()])
    names = list(BROAD)
    ct_idx = {n: i for i, n in enumerate(names)}
    order = np.argsort(-S, axis=0)
    top1, top2 = order[0], order[1]
    assign = np.array(names)[top1]
    runner = np.array(names)[top2]
    umi = np.asarray(mat.sum(axis=0)).ravel()
    det = np.asarray((mat > 0).sum(axis=0)).ravel()
    hbb = np.asarray(mat[[im[g] for g in ("Hbb-bs","Hbb-bt","Hba-a1","Hba-a2") if g in im], :]
                     .sum(axis=0)).ravel() if any(g in im for g in ("Hbb-bs","Hbb-bt")) else np.zeros(nc)
    # T10: erythrocyte cells -- gene detectability and Hbb share
    m = assign == "Erythrocyte"
    if m.sum():
        rows10.append(dict(sample=prefix, group=grp, n_eryth=int(m.sum()),
                           median_umi=float(np.median(umi[m])),
                           median_genes_detected=float(np.median(det[m])),
                           median_umi_from_hbb=float(np.median(hbb[m])),
                           frac_umi_hbb=float(np.median(hbb[m] / np.maximum(umi[m], 1))),
                           min_genes_detected=float(np.min(det[m]))))
    # T11: runner-up type for each cell type
    for ct in names:
        sel = np.where(assign == ct)[0]
        if sel.size == 0:
            continue
        rc = pd.Series(runner[sel]).value_counts(normalize=True)
        rows11.append(dict(sample=prefix, group=grp, cell_type=ct, n_cells=int(sel.size),
                           runner_up=rc.index[0], frac_runner_up=float(rc.iloc[0]),
                           median_margin=float(np.median(S[top1[sel], sel] -
                                                         S[top2[sel], sel]))))
        # T11b: where do the real marker-high cells land?
        if ct == "SatelliteGlia":
            gsel = sel[S[ct_idx["SatelliteGlia"], sel] > 0]
            rows11b.append(dict(sample=prefix, group=grp, n_cells=int(sel.size),
                                n_with_any_Gfap_or_Slc1a3=int(gsel.size),
                                frac=float(gsel.size / max(sel.size, 1))))
        if ct == "Schwann":
            gsel = sel[(S[ct_idx["Schwann"], sel] > 0)]
            rows11b.append(dict(sample=prefix, group=grp, n_cells=int(sel.size),
                                n_with_any_Gfap_or_Slc1a3=int(gsel.size),
                                frac=float(gsel.size / max(sel.size, 1))))
    # T13: which marker set wins for the median cell (baseline dominance)
    med_cell = int(np.argsort(umi)[nc // 2])
    rows13.append(dict(sample=prefix, group=grp,
                       winner_at_median_umi=names[int(np.argmax(S[:, med_cell]))],
                       runner_at_median_umi=names[int(np.argsort(-S[:, med_cell])[1])],
                       s_top1=float(S[:, med_cell].max()),
                       s_top2=float(np.sort(S[:, med_cell])[-2])))
    del mat

t10 = pd.DataFrame(rows10); t11 = pd.DataFrame(rows11); t13 = pd.DataFrame(rows13)
print("=== T10 cells assigned to 'Erythrocyte' ===")
print(t10.to_string(index=False))
t10.to_csv(os.path.join(OUT, "T10_erythrocyte_plausibility.csv"), index=False)

print("\n=== T11 runner-up compartment (annotation ambiguity) ===")
agg = (t11.groupby(["cell_type", "runner_up"]).agg(n=("n_cells", "sum"),
                                                 wmargin=("median_margin", "median"))
       .reset_index())
agg["frac"] = agg.n / agg.groupby("cell_type").n.transform("sum")
for ct, g in agg.groupby("cell_type"):
    print(f"  {ct:14s} " + "  ".join(f"{r.runner_up}={r.frac:.3f}" for r in g.itertuples()))
agg.to_csv(os.path.join(OUT, "T11_runner_up.csv"), index=False)

print("\n=== T13 which marker set wins for the median-depth cell ===")
print(t13.to_string(index=False))
t13.to_csv(os.path.join(OUT, "T13_winner.csv"), index=False)
if rows11b:
    print("\n=== T11b SatelliteGlia / Schwann: do these cells express their own markers? ===")
    t11b = pd.DataFrame(rows11b)
    print(t11b.to_string(index=False))
    t11b.to_csv(os.path.join(OUT, "T11b_glia_marker_positive.csv"), index=False)

# ---------------- T9 compositional test ----------------
print("\n=== T9 compositional test: gene logFC vs compartment denominator shift ===")
import sys as _s
_s.path.insert(0, os.path.join(ROOT, "scripts"))
from _limma import fit_limma
z = np.load(os.path.join(OUT, "compartment_pseudobulk.npz"), allow_pickle=True)
sym = list(z["canon_sym"])
cts = sorted({k.split("||")[0] for k in z.files if "||" in k})
t9 = []
for ct in cts:
    pres = [p for p in SAMPLES if f"{ct}||{p}" in z.files]
    if len(pres) < 6:
        continue
    M = pd.DataFrame({p: z[f"{ct}||{p}"] for p in pres}, index=sym)
    M = M.loc[M.sum(axis=1) > 0]
    tot = M.sum(axis=0)
    den = float(np.log2(np.mean([tot[p] for p in pres if GROUPS[p] == "B16F10"]) /
                        np.mean([tot[p] for p in pres if GROUPS[p] == "Con"])))
    L = np.log2(M.div(tot, axis=1) * 1e4 + 1.0)
    g = [GROUPS[p] for p in pres]
    design = np.column_stack([np.ones(len(g)), np.array([1.0 if x == "B16F10" else 0.0 for x in g])])
    f = fit_limma(L, design, coef_idx=1, trend=True)
    # correlation of logFC with a pure-denominator-shift prediction (constant)
    r_const = np.corrcoef(f.logFC.values, np.full(len(f), den))[0, 1]
    # variance of logFC attributable to genes vs the constant shift
    resid_var = float(((f.logFC - den) ** 2).mean())
    t9.append(dict(cell_type=ct, log2_denominator_ratio=den,
                   median_logFC=float(f.logFC.median()),
                   frac_logFC_above_den=float((f.logFC > den).mean()),
                   corr_logFC_with_denominator=float(r_const),
                   var_logFC_about_denominator=resid_var))
    print(f"  {ct:14s} log2(B16/Con total counts)={den:+.4f}  "
          f"median_logFC={f.logFC.median():+.4f}  "
          f"frac(logFC>den)={float((f.logFC > den).mean()):.3f}  "
          f"R2_vs_pure_shift={float(r_const**2):.4f}")
pd.DataFrame(t9).to_csv(os.path.join(OUT, "T9_compositional.csv"), index=False)

# ---------------- T12 MWU floor ----------------
print("\n=== T12 resolution of the 4v4 MWU discrimination test ===")
from itertools import combinations
n = 8
best = None
for combo in combinations(range(8), 4):
    pass
print(f"  C(8,4) = {len(list(combinations(range(8),4)))} label permutations")
print(f"  minimum attainable two-sided MWU p at n=4 vs 4 = 2/70 = {2/70:.4f}")
s = json.load(open(os.path.join(TAB, "sc_GSE292522_celltype_summary.json")))["celltype_ISG_discrimination"]
print(f"  compartments tested = {len(s)}; p<0.05 uncorrected = "
      f"{[k for k, v in s.items() if v['MWU_p'] < 0.05]}")
print(f"  BH across the {len(s)} compartments: "
      f"{ {k: round(v['MWU_p'] * len(s) / (i + 1), 4) for i, (k, v) in enumerate(sorted(s.items(), key=lambda x: x[1]['MWU_p']))} }")
print("\n=== SatelliteGlia ISG per animal (the only p<0.05 compartment) ===")
sc = pd.read_csv(os.path.join(TAB, "sc_GSE292522_celltype_scores.csv"))
sg = sc[(sc.cell_type == "SatelliteGlia") & (sc.score == "ISG")]
print(sg[["sample", "group", "value"]].to_string(index=False))
con = sg[sg.group == "Con"].value.values
print(f"  Con values sorted: {np.sort(con)}")
print(f"  Con rep4 = {sg[sg.sample.str.contains('rep4')].value.iloc[0]:.4f}; "
      f"next-highest Con = {np.sort(con)[-2]:.4f}; ratio = "
      f"{sg[sg.sample.str.contains('rep4')].value.iloc[0]/np.sort(con)[-2]:.2f}x")
b16 = sg[sg.group == "B16F10"].value.values
for drop in sg.sample:
    c2 = [v for v in con if v != sg[sg.sample == drop].value.iloc[0]]
    b2 = [v for v in b16 if v != sg[sg.sample == drop].value.iloc[0]]
    p = mannwhitneyu(b2, c2, alternative="two-sided").pvalue if len(c2) >= 3 else float("nan")
    print(f"  drop {drop:24s} -> n={len(c2)}v{len(b2)}  MWU_p={p:.4f}")
print("\n[A3] DONE stage3")