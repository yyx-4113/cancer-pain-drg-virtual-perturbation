#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
34_concordance_per_day.py
-------------------------
Per-day cross-model concordance between the B16F10 single-cell pseudobulk
log2FC (GSE292522, tumour vs control) and the Walker-256 bulk log2FC
(day 1 / day 7 / day 14 vs sham, muscle-adjusted limma).

Source of the manuscript Section 3.8 sentence: "highest at day 1
(rho = 0.112) and vanishing by day 14 (rho = 0.030)".
Regenerates results/tables/sc_GSE292522_concordance_per_day.json from
results/tables/ artefacts only (no raw data needed).
"""
from __future__ import annotations
import os, json
import pandas as pd
from scipy.stats import spearmanr

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAB = os.path.join(ROOT, "results", "tables")


def main():
    pb = pd.read_csv(os.path.join(TAB, "sc_GSE292522_pseudobulk_DE.csv"))
    pb = pb.set_index("Unnamed: 0")
    lfc_sc = pb["logFC"]
    out = {"method": ("Spearman of B16F10-vs-Con pseudobulk log2FC against "
                      "Walker-256 bulk log2FC (muscle-adjusted limma) over "
                      "shared genes; sc gene universe as in script 31")}
    for day in ["D1", "D7", "D14"]:
        bk = pd.read_csv(os.path.join(TAB, f"DE_{day}_vs_Sham_adj.csv"))
        gcol = "gene" if "gene" in bk.columns else bk.columns[0]
        bk = bk.set_index(gcol)
        shared = [g for g in lfc_sc.index if g in bk.index]
        r, p = spearmanr(lfc_sc.loc[shared].values,
                         bk.loc[shared, "logFC"].values)
        out[day] = {"n_shared_genes": len(shared), "spearman_rho": round(float(r), 4),
                    "spearman_p": float(p)}
    with open(os.path.join(TAB, "sc_GSE292522_concordance_per_day.json"), "w") as fh:
        json.dump(out, fh, indent=2)
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
