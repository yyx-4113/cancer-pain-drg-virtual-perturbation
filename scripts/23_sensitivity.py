#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
23_sensitivity.py  (robustness of every analysis choice)
========================================================
The manuscript reports module rankings produced by a pipeline with several
free choices: diffusion restart probability alpha, k of the k-NN graph, the
number of high-variance genes, and the lowess span of the limma-trend prior.
None of these has a canonical value. A result that only holds at one point of
that grid is a result about the grid, not about biology.

This script sweeps the grid and reports, for each setting, the Spearman
correlation of the module ranking with the reference configuration, plus how
many of the reference-significant modules stay significant. It also reports
Spearman(network model, additive model), which is the head-to-head the 2025
Nature Methods benchmark insists on.

Outputs -> results/tables/sensitivity_grid.csv, sensitivity_summary.json
"""
from __future__ import annotations

import os
import sys
import json

import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)

from _limma import (build_network, diffusion_operator, cos_batch, bh,      # noqa: E402
                    fit_limma, muscle_score, TARGETS, CONTROL_MODULES)

RAW = os.path.join(ROOT, "data", "raw")
TAB = os.path.join(ROOT, "results", "tables")
os.makedirs(TAB, exist_ok=True)

SEED = 20261009
REF = {"alpha": 0.55, "topk": 12, "hv": 3000, "frac": 0.4}
ALPHAS = [0.30, 0.45, 0.55, 0.70, 0.85]
TOPKS = [5, 8, 12, 20]
HVS = [1000, 2000, 3000]
FRACS = [0.4, 0.2, 0.3, 0.5, 0.6]   # reference first so it is available

ALL_MODS = dict(TARGETS)
ALL_MODS.update(CONTROL_MODULES)

sys.path.insert(0, HERE)
from _limma import MUSCLE_GENES                                            # noqa: E402


def log(m):
    print(f"[23] {m}", flush=True)


def load_study(dirname: str):
    """Reuse 22_etiology_comparison's reader so the two scripts cannot drift
    in how they parse the (non-uniform) GEO member files."""
    import importlib
    m = importlib.import_module("22_etiology_comparison")
    return m


def main():
    m22 = load_study("GSE149648")
    Lb_all = None
    # ---- load both studies exactly as script 22 does ----
    bone_fpkm, fb = m22.read_fpkm("GSE149648")
    chemo_fpkm, fc = m22.read_fpkm("GSE160543")
    gb = {g: m22.group_from_filename(fb[g], "bone") for g in bone_fpkm.columns}
    gc = {g: m22.group_from_filename(fc[g], "chemo") for g in chemo_fpkm.columns}
    keep_b = ((bone_fpkm >= 1.0).sum(axis=1) >= 4) & (bone_fpkm.mean(axis=1) >= 0.5)
    keep_c = ((chemo_fpkm >= 1.0).sum(axis=1) >= 4) & (chemo_fpkm.mean(axis=1) >= 0.5)
    Lb = m22.log2fpkm(bone_fpkm).loc[keep_b]
    Lc = m22.log2fpkm(chemo_fpkm).loc[keep_c]
    shared = sorted(set(Lb.index) & set(Lc.index))
    Lb, Lc = Lb.loc[shared], Lc.loc[shared]
    log(f"shared genes: {len(shared)}")

    ca14 = [s for s in Lb.columns if gb.get(s) == "D14"]
    ca7 = [s for s in Lb.columns if gb.get(s) == "D7"]
    sh = [s for s in Lb.columns if gb.get(s) == "Sham"]
    ptx = [s for s in Lc.columns if gc.get(s) == "PTX"]
    veh = [s for s in Lc.columns if gc.get(s) == "Vehicle"]

    def axis_of(L, ca, co, genes):
        return (L[co].mean(axis=1) - L[ca].mean(axis=1)).reindex(genes).values

    MUST = sorted({g for gs in ALL_MODS.values() for g in gs} & set(shared))

    def rank_of(L, ca, co, alpha, topk, hv, sign=-1.0):
        genes, W = build_network(L[list(ca) + list(co)], hv=hv, topk=topk,
                                 must_include=MUST)
        P = diffusion_operator(W, alpha=alpha, exact=True)
        ax = axis_of(L, ca, co, genes)
        out = {}
        for mod, gs in ALL_MODS.items():
            present = [g for g in gs if g in genes]
            if len(present) < 3:
                continue
            idx = [genes.index(g) for g in present]
            seed = np.zeros(len(genes))
            seed[idx] = sign
            out[mod] = float(cos_batch((P @ seed).reshape(-1, 1), ax)[0])
        return pd.Series(out)

    # ---------------- 1. network / diffusion grid ----------------
    rows = []
    ref_rank = None
    log("=== grid: alpha x topk x hv (GSE149648 D14 knockout) ===")
    for hv in HVS:
        for topk in TOPKS:
            for alpha in ALPHAS:
                r = rank_of(Lb, ca14, sh, alpha, topk, hv, -1.0)
                if ref_rank is None and (alpha, topk, hv) == \
                        (REF["alpha"], REF["topk"], REF["hv"]):
                    pass
                rows.append({"setting": "discovery_D14_KO", "alpha": alpha,
                             "topk": topk, "hv": hv, "frac": REF["frac"],
                             "n_modules": len(r),
                             **{f"rev_{k}": v for k, v in r.items()}})
                log(f"  alpha={alpha} topk={topk} hv={hv} -> {len(r)} modules")
    grid = pd.DataFrame(rows)
    ref = grid[(grid.alpha == REF["alpha"]) & (grid.topk == REF["topk"])
               & (grid.hv == REF["hv"])].iloc[0]
    refcols = [c for c in grid.columns if c.startswith("rev_")]
    rho = []
    for _, row in grid.iterrows():
        rho.append(stats.spearmanr(row[refcols].astype(float).values,
                                   ref[refcols].astype(float).values)[0])
    grid["spearman_vs_reference"] = rho
    # top-5 overlap
    def top5(row):
        s = row[refcols].astype(float)
        return list(s.sort_values(ascending=False).index[:5])
    ref_top5 = set(top5(ref))
    grid["top5_overlap"] = [len(set(top5(r)) & ref_top5) for _, r in grid.iterrows()]

    # ---------------- 2. additive head-to-head ----------------
    genes, W = build_network(Lb[ca14 + sh], hv=REF["hv"], topk=REF["topk"])
    P = diffusion_operator(W, alpha=REF["alpha"], exact=True)
    ax = axis_of(Lb, ca14, sh, genes)
    net, add = {}, {}
    for mod, gs in ALL_MODS.items():
        present = [g for g in gs if g in genes]
        if len(present) < 3:
            continue
        idx = [genes.index(g) for g in present]
        seed = np.zeros(len(genes))
        seed[idx] = -1.0
        net[mod] = float(cos_batch((P @ seed).reshape(-1, 1), ax)[0])
        add[mod] = float(cos_batch(seed.reshape(-1, 1), ax)[0])
    rho_na = stats.spearmanr(pd.Series(net), pd.Series(add))[0]
    log(f"\nSpearman(network, additive) at reference = {rho_na:.4f}")

    # ---------------- 3. limma-trend lowess span ----------------
    log("=== grid: lowess span (limma-trend) ===")
    de_rows = []
    ref_t = None
    for frac in FRACS + ["none"]:
        Xc = Lb[ca14 + sh]
        gvec = np.array([1.0] * len(ca14) + [0.0] * len(sh))
        msc = muscle_score(Xc).values.astype(float)
        msc -= msc.mean()
        design = np.column_stack([np.ones(len(ca14) + len(sh)), gvec, msc])
        r = fit_limma(Xc, design, coef_idx=1, trend=(frac != "none"),
                      frac=(0.4 if frac == "none" else frac))
        r["q"] = bh(r["P"].values)
        if frac == REF["frac"]:
            ref_t = r["t"]
        rho_t = np.nan if (ref_t is None or frac == "none") else \
            stats.spearmanr(r["t"].reindex(ref_t.index), ref_t)[0]
        de_rows.append({"setting": "limma_trend_span", "frac": frac,
                        "n_q05": int((r["q"] < 0.05).sum()),
                        "n_q10": int((r["q"] < 0.10).sum()),
                        "n_P01": int((r["P"] < 0.01).sum()),
                        "spearman_t_vs_reference": (float(rho_t)
                                                    if np.isfinite(rho_t) else np.nan)})
        log(f"  frac={frac}: q<0.05={int((r['q'] < 0.05).sum())} "
            f"spearman(t)={rho_t:.4f}" if np.isfinite(rho_t)
            else f"  frac={frac}: q<0.05={int((r['q'] < 0.05).sum())} "
                 f"spearman(t)=NA")
    de_grid = pd.DataFrame(de_rows)

    # ---------------- 4. cross-aetiology Delta grid ----------------
    log("=== grid: alpha x topk, cross-aetiology Delta (Bone_D14 vs PTX) ===")
    d_rows = []
    # the two networks must share one node set, otherwise the intersection of
    # two independent top-variance selections removes the target genes
    v2 = (Lb[ca14 + sh].var(axis=1).add(Lc[ptx + veh].var(axis=1),
                                        fill_value=0.0) / 2.0).loc[shared] \
        .sort_values(ascending=False)
    NET_GENES = list(v2.index[:REF["hv"]]) + \
        [g for g in MUST if g not in set(v2.index[:REF["hv"]])]
    for alpha in ALPHAS:
        for topk in [5, 12, 20]:
            _, Wa = build_network(Lb[ca14 + sh], topk=topk, genes=NET_GENES)
            _, Wb = build_network(Lc[ptx + veh], topk=topk, genes=NET_GENES)
            common = NET_GENES
            Pa = diffusion_operator(Wa, alpha=alpha, exact=True)
            Pb = diffusion_operator(Wb, alpha=alpha, exact=True)
            axa = axis_of(Lb, ca14, sh, common)
            axb = axis_of(Lc, ptx, veh, common)
            d = {}
            for mod, gs in ALL_MODS.items():
                present = [g for g in gs if g in common]
                if len(present) < 3:
                    continue
                seed = np.zeros(len(common))
                seed[[common.index(g) for g in present]] = -1.0
                d[mod] = float(cos_batch((Pa @ seed).reshape(-1, 1), axa)[0]
                               - cos_batch((Pb @ seed).reshape(-1, 1), axb)[0])
            d_rows.append({"setting": "delta_BoneD14_vs_PTX", "alpha": alpha,
                           "topk": topk, "hv": REF["hv"],
                           **{f"d_{k}": v for k, v in d.items()}})
            log(f"  alpha={alpha} topk={topk} -> {len(d)} modules")
    dgrid = pd.DataFrame(d_rows)
    dref = dgrid[(dgrid.alpha == REF["alpha"]) & (dgrid.topk == REF["topk"])].iloc[0]
    dcols = [c for c in dgrid.columns if c.startswith("d_")]
    dgrid["spearman_vs_reference"] = [
        stats.spearmanr(r[dcols].astype(float).values,
                        dref[dcols].astype(float).values)[0]
        for _, r in dgrid.iterrows()]

    out = pd.concat([grid, dgrid], ignore_index=True, sort=False)
    out.to_csv(os.path.join(TAB, "sensitivity_grid.csv"), index=False)
    de_grid.to_csv(os.path.join(TAB, "sensitivity_limma_trend.csv"), index=False)

    summ = {
        "reference": REF,
        "spearman_network_vs_additive": float(rho_na),
        "grid_discovery_min_spearman": float(grid["spearman_vs_reference"].min()),
        "grid_discovery_median_spearman": float(grid["spearman_vs_reference"].median()),
        "grid_discovery_min_top5_overlap": int(grid["top5_overlap"].min()),
        "grid_delta_min_spearman": float(dgrid["spearman_vs_reference"].min()),
        "grid_delta_median_spearman": float(dgrid["spearman_vs_reference"].median()),
        "de_q05_by_frac": {str(r["frac"]): int(r["n_q05"])
                           for _, r in de_grid.iterrows()},
        "n_grid_points": int(len(grid) + len(dgrid)),
    }
    with open(os.path.join(TAB, "sensitivity_summary.json"), "w") as fh:
        json.dump(summ, fh, indent=2)
    log(json.dumps(summ, indent=2))
    log("DONE")


if __name__ == "__main__":
    main()
