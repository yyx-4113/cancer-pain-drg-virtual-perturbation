#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
21_validate_pipeline.py
=======================
Pipeline validation on the SECOND, independent dataset (GSE249443).

Purpose
-------
1. Verify that the functions used in 20_multistage_pipeline.py run correctly
   end-to-end on a real matrix before trusting the discovery results
   (code audit, not just output audit).
2. Produce an independent replication cohort for the manuscript:
   GSE249443 = rat DRG, tibial MRMT-1 bone cancer pain, BCP (n=3) vs
   sham (n=3), FPKM. Same species (rat), same tissue (DRG), same pain
   modality (bone cancer pain) as the discovery set -> a clean replication,
   not a confounded cross-tissue comparison.

Key precedent from earlier runs: in GSE249443, the top "DEGs" were skeletal
muscle / sarcomere genes (Tnni2, Acta1, Ckm, Myl1; log2FC -4 to -5.4), i.e.
unequal muscle contamination between groups. After muscle-score adjustment
the q<0.05 count fell to 0. That is a *finding about this dataset*, and it is
reported as such rather than being silently dropped.
"""
from __future__ import annotations

import os
import sys
import json
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _limma import fit_limma, bh, muscle_score, build_network, \
    diffusion_operator, cos_batch, ridge_baseline, TARGETS  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
PROC = os.path.join(ROOT, "data", "processed")
TAB = os.path.join(ROOT, "results", "tables")
os.makedirs(TAB, exist_ok=True)
SEED = 20261009
rng = np.random.default_rng(SEED)


def log(m): print(f"[21] {m}", flush=True)


def main():
    p = os.path.join(PROC, "expr_log2_bone.csv")
    if not os.path.exists(p):
        log(f"missing {p}; run 10_bone_end_to_end.py first")
        sys.exit(2)
    L = pd.read_csv(p, index_col=0)
    log(f"loaded GSE249443: {L.shape[0]} genes x {L.shape[1]} samples; cols={list(L.columns)}")

    case = [c for c in L.columns if "sBCP" not in c]
    ctrl = [c for c in L.columns if "sBCP" in c]
    log(f"BCP={case}  sham={ctrl}")

    # ---- confounding diagnostic -------------------------------------------
    ms = muscle_score(L)
    log(f"muscle score: BCP mean={ms[case].mean():+.3f}  sham mean={ms[ctrl].mean():+.3f}  "
        f"Welch p={stats.ttest_ind(ms[case], ms[ctrl], equal_var=False).pvalue:.4f}")

    # ---- DE with and without muscle adjustment -----------------------------
    out = {}
    for tag, use_m in [("unadj", False), ("adj", True)]:
        cols = case + ctrl
        X = L[cols]
        g = np.array([1.0] * len(case) + [0.0] * len(ctrl))
        if use_m:
            mc = ms[cols].values.astype(float); mc -= mc.mean()
            design = np.column_stack([np.ones(len(cols)), g, mc])
        else:
            design = np.column_stack([np.ones(len(cols)), g])
        r = fit_limma(X, design, coef_idx=1, trend=True)
        r["q"] = bh(r["P"].values)
        r = r.sort_values("P")
        r.to_csv(os.path.join(TAB, f"VALID_GSE249443_DE_{tag}.csv"))
        out[tag] = r
        log(f"  {tag}: q<0.05={int((r['q']<0.05).sum())} q<0.10={int((r['q']<0.10).sum())}")

    # ---- what are the top genes, and are they muscle? ----------------------
    for tag in ["unadj", "adj"]:
        top = out[tag].head(15)
        log(f"  top-15 ({tag}): {', '.join(top.index)}")
        mus = MUSCLE_SET.intersection(set(top.index))
        log(f"    of which skeletal-muscle genes: {sorted(mus) if mus else 'none'}")

    # ---- virtual perturbation (replication of the discovery direction) -----
    genes, W = build_network(L, hv=2500, topk=10)
    P = diffusion_operator(W, alpha=0.55, n_iter=80)
    pain_axis = (L[ctrl].mean(axis=1) - L[case].mean(axis=1)).reindex(genes).values
    rows = []
    for mod, tlist in TARGETS.items():
        present = [g for g in tlist if g in genes]
        if not present:
            continue
        idx = np.array([genes.index(g) for g in present])
        for ptype, sign in [("knockout", -1.0), ("overexpression", 1.0)]:
            seed = np.zeros(len(genes)); seed[idx] = sign
            s_net = float(cos_batch((P @ seed).reshape(-1, 1), pain_axis)[0])
            s_add = float(cos_batch(seed.reshape(-1, 1), pain_axis)[0])
            s_ridge = ridge_baseline(L, genes, present, sign, pain_axis)
            Ind = np.zeros((len(genes), 500))
            for k in range(500):
                Ind[rng.choice(len(genes), len(present), replace=False), k] = 1.0
            sp = cos_batch(P @ (Ind * sign), pain_axis)
            pe = (np.sum(sp <= s_net) + 1) / 501 if sign < 0 else (np.sum(sp >= s_net) + 1) / 501
            rows.append({"module": mod, "n": len(present), "perturbation": ptype,
                         "reversal_network": s_net, "reversal_additive": s_add,
                         "reversal_ridge": s_ridge, "p_empirical": pe,
                         "z_vs_perm": (s_net - sp.mean()) / (sp.std() + 1e-12)})
    pert = pd.DataFrame(rows)
    pert["q_empirical"] = bh(pert["p_empirical"].values)
    pert = pert.sort_values("p_empirical")
    pert.to_csv(os.path.join(TAB, "VALID_GSE249443_perturbation.csv"), index=False)
    log(pert.head(24).to_string(index=False))

    summary = {"dataset": "GSE249443", "n_case": len(case), "n_ctrl": len(ctrl),
               "n_genes": int(L.shape[0]),
               "muscle_score_BCP": float(ms[case].mean()),
               "muscle_score_sham": float(ms[ctrl].mean()),
               "n_sig_q05_unadj": int((out["unadj"]["q"] < 0.05).sum()),
               "n_sig_q05_adj": int((out["adj"]["q"] < 0.05).sum())}
    with open(os.path.join(TAB, "VALID_GSE249443_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    log(json.dumps(summary, indent=2))
    log("DONE")


MUSCLE_SET = {"Acta1", "Ckm", "Myl1", "Myh1", "Myh2", "Tnni2", "Tnnc2", "Tnnt3",
              "Des", "Mybph", "Actn3", "Mybpc1", "Myoz1", "Tnnc1", "Mylpf",
              "Myh4", "Atp2a1", "Casq1", "Mb", "Actc1"}

if __name__ == "__main__":
    main()
