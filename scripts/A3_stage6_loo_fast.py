#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""A3 stage 6: fast leave-one-animal-out at gene level.

The question is 'does the compartment-level call survive dropping any single
animal?'. We answer it on the log2(CPM+1) scale with plain group means and a
Welch t-test, which needs no lowess refit. We ALSO recompute the moderated-t
version once per (compartment, gene) on a fixed 3000-gene universe, so the
significance call can be checked under the model the paper actually used.
"""
from __future__ import annotations
import os, sys, json
import numpy as np, pandas as pd
from scipy import stats

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

def logcpm(ct, samples):
    df = pd.DataFrame({p: pd.Series(z[f"{ct}||{p}"]).groupby(sym.values).sum()
                       for p in samples})
    df = df.loc[df.sum(axis=1) > 0]
    return np.log2(df.div(df.sum(axis=0), axis=1) * 1e4 + 1.0)

HEAD = [("Neuron", "Gal"), ("Neuron", "Scn10a"), ("Neuron", "Ifit2"), ("Neuron", "Cx3cr1"),
        ("Neuron.NF", "Gal"), ("Neuron.NF", "Scn10a"), ("Neuron.NF", "Il1b"),
        ("Neuron.NF", "Cx3cr1"), ("Neuron.NF", "Gfap"),
        ("SatelliteGlia", "Isg15"), ("SatelliteGlia", "Ifit1"), ("SatelliteGlia", "Calca"),
        ("Schwann", "Gal"), ("Schwann", "Ccl2"), ("Schwann", "Snap25"),
        ("Endothelial", "Lyz2"), ("Endothelial", "Gal"), ("Neuron.PE", "Gal")]

rows = []
print(f"{'compartment':14s} {'gene':7s} {'fullLFC':>8s} {'LOOmin':>8s} {'LOOmax':>8s} "
      f"{'flip':>5s} {'Welch p full':>12s} {'maxLOO p':>9s}", flush=True)
for ct, gene in HEAD:
    pres = [p for p in keys if f"{ct}||{p}" in z.files]
    if len(pres) < 6:
        continue
    L = logcpm(ct, pres)
    if gene not in L.index:
        continue
    c = [p for p in pres if G[p] == "Con"]
    b = [p for p in pres if G[p] == "B16F10"]
    l0 = float(L.loc[gene, b].mean() - L.loc[gene, c].mean())
    _, p0 = stats.ttest_ind(L.loc[gene, b], L.loc[gene, c], equal_var=False)
    lf, ps = [], []
    for d in pres:
        bb = [p for p in b if p != d]
        cc = [p for p in c if p != d]
        if len(bb) < 3 or len(cc) < 3:
            continue
        lv = float(L.loc[gene, bb].mean() - L.loc[gene, cc].mean())
        _, pv = stats.ttest_ind(L.loc[gene, bb], L.loc[gene, cc], equal_var=False)
        lf.append(lv); ps.append(float(pv))
    flip = sum(1 for v in lf if (v > 0) != (l0 > 0))
    print(f"{ct:14s} {gene:7s} {l0:+8.3f} {min(lf):+8.3f} {max(lf):+8.3f} {flip:5d} "
          f"{p0:12.4f} {max(ps):9.4f}", flush=True)
    rows.append(dict(cell_type=ct, gene=gene, full_logFC=l0, welch_p=p0,
                     loo_min=min(lf), loo_max=max(lf), sign_flips=flip,
                     max_loo_p=max(ps), n_loo=len(lf)))

# --- moderated-t version on a fixed universe, for the significance call only ---
print()
print("=== moderated limma-trend q, refit per LOO on a fixed 3000-gene universe ===", flush=True)
KEEP = 3000
print(f"{'compartment':14s} {'gene':7s} {'fullq':>8s} {'minLOOq':>8s} {'maxLOOq':>8s} "
      f"{'n q<.05':>9s}", flush=True)
rows2 = []
for ct, gene in HEAD:
    pres = [p for p in keys if f"{ct}||{p}" in z.files]
    if len(pres) < 6:
        continue
    L0 = logcpm(ct, pres)
    genes = list(L0.var(axis=1).sort_values(ascending=False).index[:KEEP])
    def q_of(samples):
        df = pd.DataFrame({p: pd.Series(z[f"{ct}||{p}"]).groupby(sym.values).sum()
                           for p in samples})
        df = df.loc[df.sum(axis=1) > 0]
        df = df.loc[[g for g in genes if g in df.index]]
        L = np.log2(df.div(df.sum(axis=0), axis=1) * 1e4 + 1.0)
        g = [G[p] for p in samples]
        des = np.column_stack([np.ones(len(g)),
                               np.array([1.0 if x == "B16F10" else 0.0 for x in g])])
        f = fit_limma(L, des, coef_idx=1, trend=True)
        f["q"] = bh(f["P"].values)
        return (float(f.loc[gene].q) if gene in f.index else np.nan)
    q0 = q_of(pres)
    qs = []
    for d in pres:
        sub = [p for p in pres if p != d]
        if len(sub) >= 6:
            qs.append(q_of(sub))
    nq = sum(1 for v in qs if v == v and v < 0.05)
    print(f"{ct:14s} {gene:7s} {q0:8.4f} {min(qs):8.4f} {max(qs):8.4f} {nq:5d}/{len(qs)}",
          flush=True)
    rows2.append(dict(cell_type=ct, gene=gene, full_q=q0, loo_q_min=min(qs),
                      loo_q_max=max(qs), n_q05=nq, n_loo=len(qs)))

pd.DataFrame(rows).to_csv(os.path.join(OUT, "A3_LOO_welch.csv"), index=False)
pd.DataFrame(rows2).to_csv(os.path.join(OUT, "A3_LOO_limma_q.csv"), index=False)
print("DONE stage6", flush=True)