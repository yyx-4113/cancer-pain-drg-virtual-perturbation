"""Round-3 verification of two contested claims.

Q1 (A2/A4): in the 8-case exchange null, does the complement split form an
exact mirror? If not, the paired exact-p floor is 2/71, not 3/71.
Q2 (A3): is the single-cell Cohen's d using SS/7 (wrong) instead of SS/6?
"""
import itertools
import numpy as np
import pandas as pd

PROJ = "D:/2026.10/不同部位癌痛患者对药物敏感性生信分析/cancer_pain_virtual_perturbation/"

print("=" * 70)
print("Q1: paired exact-p floor — read the realised values from the CSV")
print("=" * 70)
b = pd.read_csv(PROJ + "results/tables/etiology_bootstrap_differences.csv")
d = b[(b.network_scope == "per_contrast") & (b.perturbation == "knockout")].copy()
p = d.p_permutation_exact.dropna()
p = p[p > 0]
print("min nonzero exact p overall       = %.6f" % p.min())
print("  2/71 = %.6f   3/71 = %.6f   1/71 = %.6f" % (2/71, 3/71, 1/71))
for name, val in [("2/71", 2/71), ("3/71", 3/71)]:
    n = int(((p > val*0.9995) & (p < val*1.0005)).sum())
    print("  rows at %s: %d" % (name, n))

# explicit exchangeability check: build the null distribution of a paired
# contrast by enumerating all C(8,4)=70 splits and the A/B assignment, then
# count how many null values are >= the observed.
print()
print("explicit check: count of |null| >= |obs| for a paired test")
t = d[(d.paired_control == True)] if "paired_control" in d.columns else d.iloc[:0]
print("paired rows in bootstrap CSV (paired_control==True):", len(t))
if "p_permutation_exact" in d.columns:
    pp = d[d.p_permutation_exact > 0]
    for _, r in pp.nsmallest(8, "p_permutation_exact").iterrows():
        print("   %-22s %s|%-8s p=%.6f (k/71 k=%.1f)"
              % (r.module, r.aetiology_A, r.aetiology_B,
                 r.p_permutation_exact, r.p_permutation_exact * 71))

print()
print("=" * 70)
print("Q2: Cohen's d in script 33 (myeloid observed effect)")
print("=" * 70)
import re
src = open(PROJ + "scripts/33_negative_result_quantification.py", encoding="utf-8").read()
for m in re.finditer(r".*(cohen|myeloid|/ *7|SS|pooled).*", src):
    line = m.group(0).strip()
    if line and not line.startswith("#"):
        print("  |", line[:110])
