"""Decisive breakdown of the exact-p floor: paired vs unpaired rows."""
import numpy as np
import pandas as pd

PROJ = "D:/2026.10/不同部位癌痛患者对药物敏感性生信分析/cancer_pain_virtual_perturbation/"
b = pd.read_csv(PROJ + "results/tables/etiology_bootstrap_differences.csv")
d = b[(b.network_scope == "per_contrast") & (b.perturbation == "knockout")].copy()
print("columns:", [c for c in d.columns])
print("contrast pairs present:",
      sorted(set(zip(d.aetiology_A, d.aetiology_B))))

for flag in [True, False]:
    sub = d[d.paired_control == flag]
    p = sub.p_permutation_exact.dropna()
    p = p[p > 0]
    print("\npaired_control == %s : %d rows" % (flag, len(sub)))
    print("  min p = %.6f" % p.min())
    for val, lab in [(1/71, "1/71"), (2/71, "2/71"), (3/71, "3/71"),
                     (3/4901, "3/4901"), (2/4901, "2/4901")]:
        n1 = int((abs(p - val) < 1e-9).sum())
        n2 = int(((p > val*0.9995) & (p < val*1.0005)).sum())
        print("    ==%s: %d exact ; ~%s: %d" % (lab, n1, lab, n2))
    print("  smallest 6:")
    for _, r in sub.nsmallest(6, "p_permutation_exact").iterrows():
        print("    %-22s %s|%-8s p=%.6f x71=%.2f x4901=%.2f"
              % (r.module, r.aetiology_A, r.aetiology_B,
                 r.p_permutation_exact, r.p_permutation_exact*71,
                 r.p_permutation_exact*4901))

# marker level
print("\n" + "="*60)
m = pd.read_csv(PROJ + "results/tables/etiology_marker_cross.csv")
print("marker paired min p =", m[m.paired_control == True].p_permutation_exact.min())
print("marker unpaired min p =", m[m.paired_control == False].p_permutation_exact.min())
