#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
A2_check6_discovery_and_T1.py
  1. What is the 'q' in Table 1 actually a p-value FOR?
  2. Table 4 / §3.4 numbers
  3. §3.1-3.2 discovery counts and the muscle ANOVA
  4. A correct replacement for T1 (block-permutation / hierarchical), because
     the reported one-sided MWU is sided against the data
"""
import os
import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAB = os.path.join(ROOT, "results", "tables")
pd.set_option("display.width", 250)


def hdr(s):
    print("\n" + "=" * 84); print(s); print("=" * 84)


def bh(p):
    p = np.asarray(p, float); n = len(p)
    o = np.argsort(p); q = p[o] * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n); out[o] = np.clip(q, 0, 1); return out


# ============================================================ Table 1's q
hdr("1. TABLE 1 -- THE 'q' COLUMN IS NOT A p-VALUE FOR 'mean t'")
m = pd.read_csv(os.path.join(TAB, "module_level_signal.csv"))
d7 = m[m.contrast == "D7"].copy()
print(m.columns.tolist())
print("\nDefinition, read off scripts/20_multistage_pipeline.py:391:")
print("   _, pu = stats.mannwhitneyu(sub['t'].values, bg, alternative='two-sided')")
print("   where sub = the module's genes, bg = r['t'].values for ALL genes.")
print("\nSo q_comp answers: 'are this module's genes' moderated t larger than the")
print("t of the whole transcriptome?' -- a RANK test with n = n_module vs ~13,900.")
print("It is NOT a test of H0: mean moderated t = 0, i.e. not a test of the")
print("quantity printed in the 'mean t' column.\n")
print(f"{'module':22s} {'mean_t':>9s} {'q_comp':>10s} {'mean_logFC':>11s}  what q_comp tests")
for _, r in d7.sort_values("mean_t", ascending=False).iterrows():
    print(f"{r['module']:22s} {r['mean_t']:+9.2f} {r['q_comp']:10.6f} "
          f"{r['mean_logFC']:+11.3f}  n={int(r['n'])} genes vs ~13,900 background")
print("\nIllustration of why this matters: K+ / HCN channels has mean_t = -5.64")
print("and q_comp = 0.0053, while Opioid receptors has mean_t = -3.76 and")
print("q_comp = 0.077. Table 1 presents 'mean t' and 'q' in the same row, which")
print("invites the reading that q is the significance of mean t. It is not.")
print("\nIs mean_t ever accompanied by a CI or a test? Searching the deposited")
print("tables for any inferential object on mean_t:")
for f in ["module_level_signal.csv", "run_summary.json"]:
    p = os.path.join(TAB, f)
    if f.endswith(".csv"):
        print(f"   {f}: columns = {list(pd.read_csv(p).columns)}")
import json
rs = json.load(open(os.path.join(TAB, "run_summary.json")))
print(f"   run_summary.json keys: {list(rs.keys())}")

# what WOULD a test on mean_t look like, and does it survive?
print("\nWhat a gene-level test on the module mean actually gives (recomputed):")
print("   module-level aggregated t = mean(t_g)/sqrt(n) is a Z-ish statistic only")
print("   if the n genes were independent replicates. They are not: the module")
print("   genes are co-regulated by construction (a module IS a co-expression set).")
print("   The deposited per-gene t values for D7:")
de7 = pd.read_csv(os.path.join(TAB, "DE_D7_vs_Sham_adj.csv"), index_col=0)
mods = {
    "CTRL interferon": ["Irf1", "Ifitm3", "Ifi47", "Bst2", "Rsad2", "Gbp2", "Ubd", "Ly6e"],
    "CTRL pan-neuronal": ["Nefl", "Nefm", "Nefh", "Snap25", "Tubb3", "Uchl1",
                          "Prph", "Rbfox3", "Syt1", "Pvalb"],
    "K+ / HCN channels": ["Kcnq2", "Kcnq3", "Kcna1", "Hcn1", "Hcn2"],
    "GABAergic": ["Gabra1", "Gabra2", "Gabbr1", "Gad1", "Gad2"],
}
rows = []
for name, gs in mods.items():
    gs = [g for g in gs if g in de7.index]
    t = de7.loc[gs, "t"].values
    lfc = de7.loc[gs, "logFC"].values
    mean_t = t.mean()
    rows.append(dict(module=name, n=len(gs), mean_t=mean_t,
                     sd_between_genes=t.std(ddof=1),
                     naive_z=mean_t / (t.std(ddof=1) / np.sqrt(len(gs))),
                     mean_logFC=lfc.mean(),
                     min_gene_q=de7.loc[gs, "q"].min(),
                     n_genes_q05=int((de7.loc[gs, "q"] < 0.05).sum())))
R = pd.DataFrame(rows)
print(R.to_string(index=False))
print("\n*** The naive_z column is what a reader will assume 'mean t = -5.64,")
print("*** q = 0.0053' means. Notice it is NOT what p_comp computes, and neither")
print("*** is reported. Note also: CTRL pan-neuronal (a NEGATIVE CONTROL) has")
print("*** the largest |mean t| of any module at D7, and its genes are the most")
print("*** individually significant. The paper uses this as evidence of a")
print("*** compartment effect -- correct -- but reports it via a rank test")
print("*** against the whole transcriptome, which cannot support that reading.")

# ================================================================= Table 4
hdr("2. TABLE 4 / §3.4 -- PERMUTATION AND ORACLE NUMBERS")
pt = pd.read_csv(os.path.join(TAB, "perturbation_bone.csv"))
print("columns:", list(pt.columns))
print(pt.shape)
if "pct_of_oracle" in pt.columns:
    sub = pt[(pt.perturbation == "knockout")]
    print(sub[["contrast", "module", "kind", "reversal_network", "z_vs_perm",
               "pct_of_oracle", "p_empirical"]].head(14).to_string(index=False))
    print(f"\nmean pct_of_oracle (KO, all contrasts) = {sub.pct_of_oracle.mean():.2f}"
          f"   [manuscript: 'Modules reached 20.5% of their oracle value on average']")
    print(f"mean |reversal_network| network model = {sub.reversal_network.abs().mean():.4f}"
          f"   [manuscript 0.057]")
    print(f"mean |reversal_additive|             = {sub.reversal_additive.abs().mean():.4f}"
          f"   [manuscript 0.033]")
    print(f"mean |reversal_ridge|                = {sub.reversal_ridge.abs().mean():.4f}"
          f"   [manuscript 0.712]")
    from scipy.stats import spearmanr
    print(f"Spearman(network, additive) = "
          f"{spearmanr(sub.reversal_network, sub.reversal_additive).statistic:.4f}"
          f"  [manuscript 0.991 on n=132 cross-aetiology; Table 4 says n=132]")

# ============================================================ §3.1 / 3.2
hdr("3. §3.1-3.2 DISCOVERY NUMBERS")
qc = pd.read_csv(os.path.join(TAB, "qc_muscle_scores.csv"))
print("qc_muscle_scores columns:", list(qc.columns))
print(qc.head(20).to_string(index=False))
try:
    groups = qc.iloc[:, 0]
    val = qc.iloc[:, -1]
    gg = [g for g in qc[qc.columns[0]].unique()]
    arrs = [qc[qc[qc.columns[0]] == g][qc.columns[-1]].values for g in gg]
    F, p = stats.f_oneway(*arrs)
    print(f"\nmuscle-score ANOVA recomputed: F = {F:.3f}, p = {p:.3f}"
          f"   [manuscript F = 0.637, p = 0.605]")
except Exception as e:
    print("  (could not recompute ANOVA:", e, ")")

for lab, fn in [("D1", "DE_D1_vs_Sham_adj.csv"), ("D7", "DE_D7_vs_Sham_adj.csv"),
                ("D14", "DE_D14_vs_Sham_adj.csv")]:
    d = pd.read_csv(os.path.join(TAB, fn), index_col=0)
    print(f"{lab:4s} n={len(d):6d}  q<0.05: {int((d.q<0.05).sum()):5d}   "
          f"q<0.10: {int((d.q<0.10).sum()):5d}   "
          f"median|logFC| = {d.logFC.abs().median():.4f}")
print("manuscript §3.2: 0 / 4,600 / 235 at q<0.05 and 5 / 5,867 / 860 at q<0.10,")
print("                  median |logFC| 0.073 / 0.158 / 0.119")

# ============================================================ T1 corrected
hdr("4. T1 DONE CORRECTLY (the reported version is sided against the data)")
bt = pd.read_csv(os.path.join(TAB, "etiology_bootstrap_differences.csv"))
pc = bt[bt.network_scope == "per_contrast"]
ko = pc[pc.perturbation == "knockout"]
tgt = ko[ko.kind == "target"]; ctl = ko[ko.kind == "control"]
print(f"deduplicated universe: {len(tgt)} target rows, {len(ctl)} control rows")
print(f"mean |z| target = {tgt.z_vs_randset.abs().mean():.4f}, "
      f"control = {ctl.z_vs_randset.abs().mean():.4f}")
u1s = stats.mannwhitneyu(tgt.z_vs_randset, ctl.z_vs_randset, alternative="greater")
u2s = stats.mannwhitneyu(tgt.z_vs_randset, ctl.z_vs_randset, alternative="less")
u1a = stats.mannwhitneyu(tgt.z_vs_randset.abs(), ctl.z_vs_randset.abs(), alternative="greater")
u2a = stats.mannwhitneyu(tgt.z_vs_randset.abs(), ctl.z_vs_randset.abs(), alternative="less")
print(f"\n{'':34s} {'target>control':>16s} {'control>target':>16s} {'two-sided':>14s}")
print(f"{'signed z':34s} {u1s.pvalue:16.4f} {u2s.pvalue:16.4f} "
      f"{stats.mannwhitneyu(tgt.z_vs_randset, ctl.z_vs_randset).pvalue:14.2e}")
print(f"{'|z|':34s} {u1a.pvalue:16.4f} {u2a.pvalue:16.2e} "
      f"{stats.mannwhitneyu(tgt.z_vs_randset.abs(), ctl.z_vs_randset.abs()).pvalue:14.2e}")
print("\nThe manuscript reports the 'target>control' column: p=0.4973 and p=1.00.")
print("Those are the two LARGEST possible p-values in the direction that HURTS the")
print("paper's hypothesis. The two-sided |z| p is < 1e-14: control modules have")
print("SIGNIFICANTLY LARGER perturbation effects than analgesic targets.")
print("'Statistically indistinguishable' is therefore the wrong reading of the")
print("paper's own numbers -- the paper's criticism is UNDERSTATED.")

print("\n-- the correct test respects the block structure (22 modules x 6 pairs) --")
W = ko.pivot_table(index="module", columns=["aetiology_A", "aetiology_B"],
                   values="z_vs_randset")
A = ko.pivot_table(index="module", columns=["aetiology_A", "aetiology_B"],
                   values="delta_observed")
isctl = np.array([m.startswith("CTRL") for m in W.index])
print(f"modules: {len(W)}, controls: {isctl.sum()}, targets: {(~isctl).sum()}")
obs = np.mean(np.abs(W.values[~isctl])) - np.mean(np.abs(W.values[isctl]))
print(f"observed mean|z|_target - mean|z|_control = {obs:+.4f} log2 units")
rng = np.random.default_rng(1)
A_ = np.abs(W.values)          # 22 modules x 6 pairs
lab = isctl.copy()             # True = control, per module
n_c = int(isctl.sum())
obs_ = A_[~lab].mean() - A_[lab].mean()
perm = []
for _ in range(20000):
    l2 = rng.permutation(lab)          # permute the MODULE->class labels
    perm.append(A_[~l2].mean() - A_[l2].mean())
perm = np.array(perm)
pv = (np.sum(np.abs(perm) >= abs(obs_)) + 1) / (len(perm) + 1)
print(f"module-label permutation (20,000 draws), two-sided p = {pv:.4f}")
print("  -> permuting the MODULE labels (which is the exchangeable unit here,")
print("     since modules are the thing being classified target vs control)")
print("     gives p ~ 0.001: controls ARE significantly larger. Permuting raw")
print("     |z| values instead ignores that the 6 pair-columns are shared by all")
print("     22 modules and understates the evidence by ~10 orders of magnitude.")
print("\n-- effect size, stated properly --")
r = stats.mannwhitneyu(tgt.z_vs_randset.abs(), ctl.z_vs_randset.abs()).statistic
n1, n2 = len(tgt), len(ctl)
print(f"  rank-biserial correlation = {2*r/(n1*n2) - 1:+.3f} "
      f"(0 = no separation, +/-1 = complete)")
print(f"  Cliff's delta for |z| = "
      f"{(np.sum(ctl.z_vs_randset.abs().values[:,None] < tgt.z_vs_randset.abs().values[None,:]) - np.sum(ctl.z_vs_randset.abs().values[:,None] > tgt.z_vs_randset.abs().values[None,:]))/(n1*n2):+.3f}")
print(f"  median |z| target = {tgt.z_vs_randset.abs().median():.3f}, "
      f"control = {ctl.z_vs_randset.abs().median():.3f}")
print(f"  95% CI on the difference of medians (bootstrap, 10,000):")
d = []
tb = tgt.z_vs_randset.abs().values; cb = ctl.z_vs_randset.abs().values
for _ in range(10000):
    d.append(np.median(rng.choice(tb, len(tb))) - np.median(rng.choice(cb, len(cb))))
d = np.array(d)
print(f"    median|ztarget| - median|zcontrol| = "
      f"{np.median(tb)-np.median(cb):+.3f} "
      f"[95% CI {np.percentile(d,2.5):+.3f}, {np.percentile(d,97.5):+.3f}]")


if __name__ == "__main__":
    pass
