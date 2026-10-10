#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""A3 stage 4: gene-level leave-one-animal-out for the headline compartment genes."""
from __future__ import annotations
import os, sys
import numpy as np, pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from _limma import fit_limma, bh

OUT = os.path.join(ROOT, "results", "review_a3")
z = np.load(os.path.join(OUT, "compartment_pseudobulk.npz"), allow_pickle=True)
sym = pd.Series(list(z["canon_sym"]))
G = {"Sample1_DRG_Con_rep1": "Con", "Sample2_DRG_Con_rep2": "Con",
     "Sample3_DRG_Con_rep3": "Con", "Sample4_DRG_Con_rep4": "Con",
     "Sample5_DRG_B16F10_rep1": "B16F10", "Sample6_DRG_B16F10_rep2": "B16F10",
     "Sample7_DRG_B16F10_rep3": "B16F10", "Sample8_DRG_B16F10_rep4": "B16F10"}
keys = list(G)

def fit(ct, samples):
    df = pd.DataFrame({p: pd.Series(z[f"{ct}||{p}"]).groupby(sym.values).sum()
                       for p in samples})
    df = df.loc[df.sum(axis=1) > 0]
    L = np.log2(df.div(df.sum(axis=0), axis=1) * 1e4 + 1.0)
    g = [G[p] for p in samples]
    des = np.column_stack([np.ones(len(g)),
                           np.array([1.0 if x == "B16F10" else 0.0 for x in g])])
    f = fit_limma(L, des, coef_idx=1, trend=True)
    f["q"] = bh(f["P"].values)
    return f

HEAD = [("Neuron", "Gal"), ("Neuron", "Scn10a"), ("Neuron", "Ifit2"), ("Neuron", "Cx3cr1"),
        ("Neuron.NF", "Gal"), ("Neuron.NF", "Scn10a"), ("Neuron.NF", "Il1b"),
        ("Neuron.NF", "Cx3cr1"), ("SatelliteGlia", "Isg15"), ("SatelliteGlia", "Ifit1"),
        ("SatelliteGlia", "Calca"), ("Schwann", "Gal"), ("Schwann", "Ccl2"),
        ("Endothelial", "Lyz2"), ("Endothelial", "Gal"), ("Neuron.PE", "Gal")]
print(f"{'compartment':14s} {'gene':7s} {'fullLFC':>8s} {'LOOmin':>8s} {'LOOmax':>8s} "
      f"{'flip':>5s} {'fullq':>7s} {'q<0.05':>8s}", flush=True)
rows = []
for ct, g in HEAD:
    pres = [p for p in keys if f"{ct}||{p}" in z.files]
    if len(pres) < 6:
        continue
    f0 = fit(ct, pres)
    if g not in f0.index:
        continue
    l0, q0 = float(f0.loc[g].logFC), float(f0.loc[g].q)
    lf, qs = [], []
    for d in pres:
        sub = [p for p in pres if p != d]
        if len(sub) < 6:
            continue
        f = fit(ct, sub)
        if g in f.index:
            lf.append(float(f.loc[g].logFC)); qs.append(float(f.loc[g].q))
    flip = sum(1 for v in lf if (v > 0) != (l0 > 0))
    nq = sum(1 for v in qs if v < 0.05)
    print(f"{ct:14s} {g:7s} {l0:+8.3f} {min(lf):+8.3f} {max(lf):+8.3f} {flip:5d} "
          f"{q0:7.4f} {nq:5d}/{len(qs)}", flush=True)
    rows.append(dict(cell_type=ct, gene=g, full_logFC=l0, full_q=q0, loo_min=min(lf),
                     loo_max=max(lf), sign_flips=flip, n_q05=nq, n_loo=len(qs),
                     n_cells_min=min(int(np.load(os.path.join(OUT, "compartment_cellcounts.json"),
                                                 allow_pickle=True) if False else [0])) if False else 0))
pd.DataFrame(rows).to_csv(os.path.join(OUT, "A3_LOO_gene_level.csv"), index=False)
print("DONE stage4", flush=True)