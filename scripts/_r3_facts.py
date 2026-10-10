"""Round-3 fact gathering for the manuscript fixes."""
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.power import TTestIndPower
import json, os

PV = "D:/2026.10/不同部位癌痛患者对药物敏感性生信分析/cancer_pain_virtual_perturbation/results/tables/"

print("=" * 72)
print("A. paired/unpaired exact-p floor for the module-logFC endpoint")
print("=" * 72)
mod = pd.read_csv(PV + "etiology_module_logfc_cross.csv")
print("module_logfc cols:", list(mod.columns))
for flag in [True, False]:
    sub = mod[mod.paired_control == flag]
    p = sub.p_permutation_exact.dropna(); p = p[p > 0]
    print("paired=%s n=%d min p=%.6f  x71=%.2f x4901=%.2f"
          % (flag, len(sub), p.min(), p.min()*71, p.min()*4901))

print()
print("=" * 72)
print("B. §3.2 median |logFC|: all-gene vs q<0.05 set")
print("=" * 72)
# per-day DE with muscle adjustment
for day, f in [("D1", "DE_D1_vs_Sham_adj.csv"),
               ("D7", "DE_D7_vs_Sham_adj.csv"),
               ("D14", "DE_D14_vs_Sham_adj.csv")]:
    try:
        d = pd.read_csv(PV + f)
        c = [x for x in d.columns if x.lower() in ("logfc", "log2fc")]
        q = [x for x in d.columns if x.lower() in ("q", "qvalue", "adj.p.val", "adj_p_val")]
        lc = d[c[0]] if c else d.iloc[:, 1]
        qv = d[q[0]] if q else None
        print("%s: all-gene median|logFC|=%.4f" % (day, lc.abs().median()))
        if qv is not None:
            sig = lc[qv < 0.05]
            print("      sig-set (q<0.05) n=%d median|logFC|=%.4f"
                  % (len(sig), sig.abs().median()))
    except Exception as e:
        print(day, "ERR", e)

print()
print("=" * 72)
print("C. single-cell composition shift")
print("=" * 72)
comp = pd.read_csv(PV + "sc_GSE292522_composition.csv")
print("cols:", list(comp.columns))
print(comp.to_string(index=False))

print()
print("=" * 72)
print("D. T3 recompute (Cohen's d) + T5 MDE conventions")
print("=" * 72)
tp = TTestIndPower()
d_mde4 = tp.solve_power(nobs1=4, alpha=0.05, power=0.8, ratio=1.0, alternative="two-sided")
d_mde2 = tp.solve_power(nobs1=2, alpha=0.05, power=0.8, ratio=1.0, alternative="two-sided")
print("solve_power(n=4/group) d=%.4f ; solve_power(n=2/group) d=%.4f" % (d_mde4, d_mde2))
con = comp[comp.group == "Con"].Myeloid.values
b16 = comp[comp.group == "B16F10"].Myeloid.values
ss = ((con-con.mean())**2).sum() + ((b16-b16.mean())**2).sum()
sd7 = np.sqrt(((con-con.mean())**2).sum()/7 + ((b16-b16.mean())**2).sum()/7)
sd6 = np.sqrt(ss/(len(con)+len(b16)-2))
print("pooled_sd /7 = %.5f  (old)  d=%.4f" % (sd7, (con.mean()-b16.mean())/sd7))
print("pooled_sd /6 = %.5f  (correct) d=%.4f" % (sd6, (con.mean()-b16.mean())/sd6))

bt = pd.read_csv(PV + "etiology_bootstrap_differences.csv")
pc = bt[(bt.network_scope == "per_contrast") & (bt.perturbation == "knockout")].copy()
se = (pc.ci95_hi - pc.ci95_lo)/(2*1.96)
med_se = float(se.median())
mde_norm = (1.959964 + 0.841621) * med_se
mde_t = d_mde2 * med_se
tg = pc[pc.kind == "target"]
print("median bootstrap SE of Delta = %.5f" % med_se)
print("MDE normal-approx (2.802*SE) = %.4f" % mde_norm)
print("MDE t-based (%.3f*SE)        = %.4f" % (d_mde2, mde_t))
print("max target |Delta| = %.4f" % tg.delta_observed.abs().max())
print("median target |Delta| = %.4f ; ratio to MDE_norm = %.2f"
      % (tg.delta_observed.abs().median(), tg.delta_observed.abs().median()/mde_norm))
for lab, mde in [("0.1046(old)", 0.1046), ("MDE_norm %.4f" % mde_norm, mde_norm)]:
    print("  targets exceeding %s: %d/%d" % (lab, (tg.delta_observed.abs() > mde).sum(), len(tg)))
print("n exact p<0.05 = %d/%d" % ((pc.p_permutation_exact < 0.05).sum(), len(pc)))

print()
print("=" * 72)
print("E. anti-selectivity excluding the interferon control module")
print("=" * 72)
pc["absz"] = pc.z_vs_randset.abs()
for excl in [None, "CTRL interferon"]:
    sub = pc if excl is None else pc[pc.module != excl]
    t = sub[sub.kind == "target"]; c = sub[sub.kind == "control"]
    u = mannwhitneyu(c.absz, t.absz, alternative="two-sided")
    print("excl=%-16s n_ctl=%d median|z| ctl=%.3f tgt=%.3f p=%.5f"
          % (str(excl), len(c), c.absz.median(), t.absz.median(), u.pvalue))
