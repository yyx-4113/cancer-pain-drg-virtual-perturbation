#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""A3 stage 5: FAST leave-one-out. Uses trend=False and a gene subset so that
the leave-one-out question ('does this compartment's call survive dropping one
animal?') is answered in minutes rather than hours. Also runs trend=True on a
restricted gene set for the headline genes only.

Rationale for the fast path: we are asking whether the SIGN and the q<0.05 call
are stable, not whether the exact moderated t reproduces to 3 decimals.
"""
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
# keep the top 4000 most variable genes -- this is what drives the lowess/EB
# prior; the retained set is identical between the full fit and every LOO fit.
KEEP = 4000

def build(ct, samples):
    df = pd.DataFrame({p: pd.Series(z[f"{ct}||{p}"]).groupby(sym.values).sum()
                       for p in samples})
    df = df.loc[df.sum(axis=1) > 0]
    return df

def fit(ct, samples, genes=None):
    df = build(ct, samples)
    if genes is not None:
        df = df.loc[[g for g in genes if g in df.index]]
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
        ("Schwann", "Snap25"), ("Endothelial", "Lyz2"), ("Endothelial", "Gal"),
        ("Neuron.PE", "Gal"), ("Neuron.PE", "Cx3cr1")]
rows = []
print(f"{'compartment':14s} {'gene':7s} {'fullLFC':>8s} {'LOOmin':>8s} {'LOOmax':>8s} "
      f"{'flip':>5s} {'fullq':>8s} {'q<0.05':>8s}", flush=True)
for ct, gene in HEAD:
    pres = [p for p in keys if f"{ct}||{p}" in z.files]
    if len(pres) < 6:
        continue
    full = build(ct, pres)
    genes = list(full.var(axis=1).sort_values(ascending=False).index[:KEEP])
    f0 = fit(ct, pres, genes)
    if gene not in f0.index:
        print(f"{ct:14s} {gene:7s}  absent from the variance-filtered set", flush=True)
        continue
    l0, q0 = float(f0.loc[gene].logFC), float(f0.loc[gene].q)
    lf, qs = [], []
    for d in pres:
        sub = [p for p in pres if p != d]
        if len(sub) < 6:
            continue
        f = fit(ct, sub, genes)
        if gene in f.index:
            lf.append(float(f.loc[gene].logFC)); qs.append(float(f.loc[gene].q))
    flip = sum(1 for v in lf if (v > 0) != (l0 > 0))
    nq = sum(1 for v in qs if v < 0.05)
    print(f"{ct:14s} {gene:7s} {l0:+8.3f} {min(lf):+8.3f} {max(lf):+8.3f} {flip:5d} "
          f"{q0:8.4f} {nq:5d}/{len(qs)}", flush=True)
    rows.append(dict(cell_type=ct, gene=gene, full_logFC=l0, full_q=q0,
                     loo_min=min(lf), loo_max=max(lf), sign_flips=flip,
                     n_q05=nq, n_loo=len(qs)))
pd.DataFrame(rows).to_csv(os.path.join(OUT, "A3_LOO_fast.csv"), index=False)
print("DONE stage5", flush=True)