#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
22b_patch_artifact.py --  Recompute the BH-adjusted q-values of
etiology_bootstrap_differences.csv with the CORRECT independent-test count.

Background
----------
Script 22 computes, for every (module, aetiology-pair, perturbation-direction),
a reversal-score Delta and two p-values (exact label permutation; degree-matched
random-gene-set). Knockout and overexpression are exact mirror images
(reversal_KO = -reversal_OE, and therefore identical exact-permutation p,
verified max|diff| = 0). They are therefore ONE independent hypothesis, not two.

The pipeline's heavy compute (networks, permutation nulls, bootstrap, random
sets) already finished and wrote the authoritative p-values to
etiology_bootstrap_differences.csv. Only the final BH step needed correcting, so
this script re-derives q from those p-values without re-running the 16-minute
pipeline. The dedup logic mirrors the (now fixed) block in 22_etiology_comparison.py.
"""
import os, json
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAB = os.path.join(ROOT, "results", "tables")
boot = pd.read_csv(os.path.join(TAB, "etiology_bootstrap_differences.csv"))


def bh(p):
    p = np.asarray(p, dtype=float)
    m = len(p)
    o = np.argsort(p)
    q = p[o] * m / (np.arange(m) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(m)
    out[o] = np.clip(q, 0, 1)
    return out


ded = boot[boot.perturbation == "knockout"].copy()
for sc in boot["network_scope"].unique():
    m = ded["network_scope"] == sc
    ded.loc[m, "q_permutation"] = bh(ded.loc[m, "p_permutation_exact"].values)
    ded.loc[m, "q_randset"] = bh(ded.loc[m, "p_randset"].values)

qmap = ded[["network_scope", "module", "aetiology_A", "aetiology_B",
            "q_permutation", "q_randset"]]
boot = boot.merge(qmap, on=["network_scope", "module", "aetiology_A",
                            "aetiology_B"], how="left", suffixes=("", "_ded"))
for base in ["q_permutation", "q_randset"]:
    dc = base + "_ded"
    if dc in boot.columns:
        boot[base] = boot.pop(dc)
boot = boot.sort_values(["p_permutation_exact", "module"])
boot.to_csv(os.path.join(TAB, "etiology_bootstrap_differences.csv"), index=False)

# ---- update summary json ----
sj = os.path.join(TAB, "etiology_summary.json")
summ = json.load(open(sj)) if os.path.exists(sj) else {}
sel = boot[boot.network_scope == "per_contrast"]
ded_all = boot[boot.perturbation == "knockout"]
ded_pc = ded_all[ded_all.network_scope == "per_contrast"]
summ.update({
    "n_rows_raw": int(len(boot)),
    "n_independent_tests_per_scope": int(len(ded_pc)),
    "n_independent_tests_total": int(len(ded_all)),
    "n_sig_perm_q05_percontrast_indep": int((ded_pc.q_permutation < 0.05).sum()),
    "n_sig_randset_q05_percontrast_indep": int((ded_pc.q_randset < 0.05).sum()),
    "min_q_permutation_per_contrast": float(ded_pc.q_permutation.min()),
    "min_q_randset_per_contrast": float(ded_pc.q_randset.min()),
    "randset_sig_breakdown_percontrast": {
        "target": int(((ded_pc.q_randset < 0.05) & (ded_pc.kind == "target")).sum()),
        "control": int(((ded_pc.q_randset < 0.05) & (ded_pc.kind == "control")).sum()),
    },
    "bh_note": ("Knockout and overexpression are exact mirror images (identical "
                "exact-permutation p); BH is applied over the deduplicated "
                "knockout set = 132 independent tests per scope."),
})
json.dump(summ, open(sj, "w"), indent=2, default=str)
print("patched. per_contrast min q_perm=%.4f (n_sig=%d); randset min q=%.4f (n_sig=%d: %d target / %d control)" % (
    ded_pc.q_permutation.min(), (ded_pc.q_permutation < 0.05).sum(),
    ded_pc.q_randset.min(), (ded_pc.q_randset < 0.05).sum(),
    ((ded_pc.q_randset < 0.05) & (ded_pc.kind == "target")).sum(),
    ((ded_pc.q_randset < 0.05) & (ded_pc.kind == "control")).sum()))
