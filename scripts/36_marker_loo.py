"""Leave-one-animal-out sensitivity for the headline cross-aetiology markers
(the exhaustive-permutation endpoint, which had no deposited LOO refit).
Reuses the loader in scripts/22_etiology_comparison.py."""
import os, sys, json
import numpy as np
import pandas as pd

ROOT = "D:/2026.10/不同部位癌痛患者对药物敏感性生信分析/cancer_pain_virtual_perturbation"
sys.path.insert(0, os.path.join(ROOT, "scripts"))
import importlib.util
spec = importlib.util.spec_from_file_location("et22", os.path.join(ROOT, "scripts", "22_etiology_comparison.py"))
et = importlib.util.module_from_spec(spec)
spec.loader.exec_module(et)           # __main__ guard means main() does not run

bone_fpkm, fb = et.read_fpkm("GSE149648")
chemo_fpkm, fc = et.read_fpkm("GSE160543")
gb = {g: et.group_from_filename(fb[g], "bone") for g in bone_fpkm.columns}
gc = {g: et.group_from_filename(fc[g], "chemo") for g in chemo_fpkm.columns}
keep_b = ((bone_fpkm >= 1.0).sum(axis=1) >= 4) & (bone_fpkm.mean(axis=1) >= 0.5)
keep_c = ((chemo_fpkm >= 1.0).sum(axis=1) >= 4) & (chemo_fpkm.mean(axis=1) >= 0.5)
Lb = et.log2fpkm(bone_fpkm).loc[keep_b]
Lc = et.log2fpkm(chemo_fpkm).loc[keep_c]
shared = sorted(set(Lb.index) & set(Lc.index))
Lb, Lc = Lb.loc[shared], Lc.loc[shared]


def logfc(L, cases, controls, gene):
    return float(L.loc[gene, cases].mean() - L.loc[gene, controls].mean())


def loo_delta(L, cases, controls, genes):
    base = {g: logfc(L, cases, controls, g) for g in genes}
    out = {}
    for drop in list(cases) + list(controls):
        c2 = [c for c in cases if c != drop]
        k2 = [c for c in controls if c != drop]
        out[drop] = {g: logfc(L, c2, k2, g) for g in genes}
    return base, out


BONE_CASES = [s for s in Lb.columns if gb.get(s) == "D7"]
BONE_CTRL = [s for s in Lb.columns if gb.get(s) == "Sham"]
PTX_CASES = [s for s in Lc.columns if gc.get(s) == "PTX"]
PTX_CTRL = [s for s in Lc.columns if gc.get(s) == "Vehicle"]
GENES = ["Cd68", "Aif1", "Gal", "Ccl2", "Il1b", "Scn10a"]

res = {"pair": "Bone_D7 vs PTX", "genes": {}, "n_animals_each": 4}
for g in GENES:
    if g not in Lb.index or g not in Lc.index:
        continue
    ab, a_loo = loo_delta(Lb, BONE_CASES, BONE_CTRL, [g])
    bb, b_loo = loo_delta(Lc, PTX_CASES, PTX_CTRL, [g])
    d0 = ab[g] - bb[g]
    # LOO delta combining one drop from either study
    deltas = []
    for drop_a in a_loo:
        deltas.append(a_loo[drop_a][g] - bb[g])
    for drop_b in b_loo:
        deltas.append(ab[g] - b_loo[drop_b][g])
    deltas = np.array(deltas)
    res["genes"][g] = {"delta_full": round(d0, 4),
                       "loo_min": round(float(deltas.min()), 4),
                       "loo_max": round(float(deltas.max()), 4),
                       "sign_stable": bool(np.all(np.sign(deltas) == np.sign(d0)))}
    print("%-7s D=%+.3f  LOO [%+.3f, %+.3f]  sign stable=%s"
          % (g, d0, deltas.min(), deltas.max(), res["genes"][g]["sign_stable"]))

with open(os.path.join(ROOT, "results", "tables", "etiology_marker_loo.json"), "w") as fh:
    json.dump(res, fh, indent=2)
print("written etiology_marker_loo.json")
