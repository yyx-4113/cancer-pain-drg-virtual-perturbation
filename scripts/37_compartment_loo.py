"""Lean leave-one-animal-out sensitivity for the per-compartment ISG programme
section 3.8 relies on, using the deposited per-sample composite ISG scores
(results/tables/sc_GSE292522_celltype_scores.csv). No matrix reload."""
import os, json
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind

ROOT = "D:/2026.10/不同部位癌痛患者对药物敏感性生信分析/cancer_pain_virtual_perturbation"
d = pd.read_csv(os.path.join(ROOT, "results", "tables", "sc_GSE292522_celltype_scores.csv"))
out = {"note": "leave-one-animal-out of the per-compartment composite ISG score",
       "method": "Welch t-test, tumour vs control, ISG score", "tests": {}}

for ct in ["Neuron", "SatelliteGlia", "Schwann", "Neuron.NF", "Neuron.PE"]:
    s = d[(d.cell_type == ct) & (d.score == "ISG")]
    if len(s) < 8:
        continue
    con = s[s.group == "Con"].set_index("sample")["value"]
    tum = s[s.group == "B16F10"].set_index("sample")["value"]
    full_lfc = float(tum.mean() - con.mean())
    full_p = float(ttest_ind(tum.values, con.values, equal_var=False).pvalue)
    loo = []
    for drop in list(con.index) + list(tum.index):
        c = con.drop(drop) if drop in con.index else con
        t = tum.drop(drop) if drop in tum.index else tum
        lfc = float(t.mean() - c.mean())
        p = float(ttest_ind(t.values, c.values, equal_var=False).pvalue)
        loo.append((drop, lfc, p))
    lfcs = np.array([l for _, l, _ in loo])
    ps = np.array([p for _, _, p in loo])
    out["tests"][ct] = {
        "full_mean_diff_tumour_minus_con": round(full_lfc, 4),
        "full_p": round(full_p, 4),
        "loo_diff_range": [round(float(lfcs.min()), 4), round(float(lfcs.max()), 4)],
        "loo_p_max": round(float(ps.max()), 4),
        "loo_runs_sign_consistent": int((np.sign(lfcs) == np.sign(full_lfc)).sum()),
        "loo_runs_p_lt_0.05": int((ps < 0.05).sum()),
        "loo_runs_total": len(loo),
    }
    print("%-14s full diff=%+.3f p=%.4f | LOO diff[%+.3f,%+.3f] p_max=%.4f "
          "sign-consistent %d/%d, p<0.05 in %d/%d"
          % (ct, full_lfc, full_p, lfcs.min(), lfcs.max(), ps.max(),
             out["tests"][ct]["loo_runs_sign_consistent"], len(loo),
             out["tests"][ct]["loo_runs_p_lt_0.05"], len(loo)))

with open(os.path.join(ROOT, "results", "tables", "sc_GSE292522_compartment_loo.json"), "w") as fh:
    json.dump(out, fh, indent=2)
print("written sc_GSE292522_compartment_loo.json")
