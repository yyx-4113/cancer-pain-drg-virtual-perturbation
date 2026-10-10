#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""A2 (biostatistician) independent recompute — round 3 panel.
All numbers recomputed from results/tables artefacts; nothing taken from
prior reviews. Run from project root."""
import os, json
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, binomtest
from statsmodels.stats.power import TTestIndPower

TAB = os.path.join("results", "tables")
bt = pd.read_csv(os.path.join(TAB, "etiology_bootstrap_differences.csv"))
nq = json.load(open(os.path.join(TAB, "negative_result_quantification.json")))

print("=== 0. universe sizes ===")
print("total rows:", len(bt), "| scopes:", bt.network_scope.unique().tolist())
pc_all = bt[bt.network_scope == "per_contrast"]
print("per_contrast rows:", len(pc_all))
pc = pc_all[pc_all.perturbation == "knockout"].copy()
print("dedup KO per_contrast:", len(pc),
      "| targets:", (pc.kind == "target").sum(),
      "| controls:", (pc.kind == "control").sum())
print("modules:", pc.module.nunique(),
      "| pairs:", pc.groupby(["aetiology_A", "aetiology_B"]).ngroups)
print("paired rows (KO):", pc.paired_control.sum(),
      "| n_null unpaired:", pc.loc[~pc.paired_control, "n_null"].unique(),
      "| n_null paired:", pc.loc[pc.paired_control, "n_null"].unique())

print("\n=== 1. KO/OE mirror check ===")
koi = pc.set_index(["module", "aetiology_A", "aetiology_B"])
oei = pc_all[pc_all.perturbation == "overexpression"].set_index(
    ["module", "aetiology_A", "aetiology_B"])
d_ko, d_oe = koi.delta_observed, oei.delta_observed
r = np.corrcoef(d_ko, d_oe)[0, 1]
print("corr(delta_KO, delta_OE) = %.6f" % r)
print("max |delta_KO + delta_OE| = %.3e" % float((d_ko + d_oe).abs().max()))
print("max |p_perm KO - p_perm OE| = %.3e" % float((koi.p_permutation_exact
                                                   - oei.p_permutation_exact).abs().max()))

print("\n=== 2. T1 Mann-Whitney |z| ===")
pc["absz"] = pc.z_vs_randset.abs()
tgt, ctl = pc[pc.kind == "target"], pc[pc.kind == "control"]
u = mannwhitneyu(ctl.absz, tgt.absz, alternative="two-sided")
print("median |z| target  = %.4f (json %.4f)" % (tgt.absz.median(),
      nq["T1_target_vs_control"]["median_abs_z_target"]))
print("median |z| control = %.4f (json %.4f)" % (ctl.absz.median(),
      nq["T1_target_vs_control"]["median_abs_z_control"]))
print("mean   |z| target  = %.4f | control = %.4f" % (tgt.absz.mean(), ctl.absz.mean()))
print("MWU two-sided p = %.6f (json %.6f)" % (u.pvalue,
      nq["T1_target_vs_control"]["mwu_abs_z_two_sided_p"]))
print("erythrocyte |z| max row:")
zmx = pc.loc[pc.z_vs_randset.abs().idxmax()]
print("  ", zmx.module, zmx.aetiology_A, "vs", zmx.aetiology_B,
      "z=%.4f p_randset=%.6g" % (zmx.z_vs_randset, zmx.p_randset))
erow = pc[(pc.module == "CTRL erythrocyte") & (pc.aetiology_A == "PTX")
          & (pc.aetiology_B == "OXA")]
print("  erythrocyte PTX vs OXA z = %.4f" % erow.z_vs_randset.iloc[0])
prow = pc[(pc.module == "CTRL pan-neuronal") & (pc.aetiology_A == "Bone_D14")
          & (pc.aetiology_B == "Bone_D7")]
print("  pan-neuronal Bone_D14 vs Bone_D7 z = %.4f" % prow.z_vs_randset.iloc[0])

print("\n=== 3. T2 binomial ===")
sig = pc[pc.q_randset < 0.05]
k_ctl = int((sig.kind == "control").sum())
base = len(ctl) / len(pc)
bp = binomtest(k_ctl, len(sig), base, alternative="greater").pvalue
print("randset sig total=%d targets=%d controls=%d (json: %d/%d/%d)" % (
    len(sig), (sig.kind == "target").sum(), k_ctl,
    nq["T2_randset_null_selectivity"]["randset_significant_total"],
    nq["T2_randset_null_selectivity"]["randset_significant_targets"],
    nq["T2_randset_null_selectivity"]["randset_significant_controls"]))
print("base rate = %d/%d = %.4f" % (len(ctl), len(pc), base))
print("binomial one-sided p = %.5f (json %.5f)" % (bp,
      nq["T2_randset_null_selectivity"]["binomial_overrepresentation_p"]))
excl = int((((sig.ci95_lo > 0) & (sig.ci95_hi > 0)) |
            ((sig.ci95_lo < 0) & (sig.ci95_hi < 0))).sum())
print("of the %d randset-sig, CI excluding zero: %d" % (len(sig), excl))
print("min p_randset in dedup universe = %.6g  (2/501 = %.6g)" % (
    pc.p_randset.min(), 2 / 501))
print("n tests at p_randset floor 2/501:", int((pc.p_randset <= 2/501 + 1e-12).sum()))

print("\n=== 4. T4 percentile ===")
allz = pc.copy()
allz["rank_abs"] = allz.absz.rank(ascending=False, method="average")
allz["pct_abs"] = (allz.rank_abs - 1) / (len(allz) - 1) * 100
tg = allz[allz.kind == "target"]
print("median pct_absz targets = %.1f (json %.1f)" % (
    tg.pct_abs.median(), nq["T4_target_rank_percentile"]["median_pct_absz"]))

print("\n=== 5. T5 power ===")
tp = TTestIndPower()
d22 = tp.solve_power(nobs1=2, alpha=0.05, power=0.8, ratio=1.0,
                     alternative="two-sided")
d44 = tp.solve_power(nobs1=4, alpha=0.05, power=0.8, ratio=1.0,
                     alternative="two-sided")
print("solve_power n=2+2 d = %.4f (json %.3f) | n=4+4 d = %.4f" % (
    d22, nq["T5_reversal_score_power"]["mde_cohens_d_unpaired_n2n2"], d44))
se = (pc.ci95_hi - pc.ci95_lo) / (2 * 1.96)
med_se = se.median()
print("median bootstrap SE (from CI width/3.92) = %.5f (json %.5f)" % (
    med_se, nq["T5_reversal_score_power"]["median_bootstrap_se_delta"]))
mde = d22 * med_se
print("MDE = %.4f (json %.4f)" % (mde, nq["T5_reversal_score_power"]["mde_delta_units"]))
print("normal-approx alternative MDE = 2.802 x SE = %.4f" % (2.8019 * med_se))
mt = pc[pc.kind == "target"].copy()
mt["absd"] = mt.delta_observed.abs()
mx = mt.loc[mt.absd.idxmax()]
print("max |target delta| = %.4f (%s, %s vs %s)" % (
    mx.absd, mx.module, mx.aetiology_A, mx.aetiology_B))
print("ratio max/MDE = %.3f (json %.2f)" % (mx.absd / mde,
      nq["T5_reversal_score_power"]["max_delta_over_mde"]))
nsig = int((pc.p_permutation_exact < 0.05).sum())
print("exact p<0.05 count = %d (json %d)" % (
    nsig, nq["T5_reversal_score_power"]["n_exact_permutation_p_lt_0.05"]))
print("NOTE json interpretation string says 26/132 (stale hardcode in script 33)")

print("\n=== 6. glial/immune row + q floor ===")
gi = pc[(pc.module == "Glial / immune") & (pc.aetiology_A == "Bone_D7")
        & (pc.aetiology_B == "PTX")].iloc[0]
print("delta=%.4f CI=[%.4f, %.4f] p=%.7f q=%.4f" % (
    gi.delta_observed, gi.ci95_lo, gi.ci95_hi,
    gi.p_permutation_exact, gi.q_permutation))
print("min q_permutation (132) = %.4f" % pc.q_permutation.min())
print("n q_permutation<0.05:", int((pc.q_permutation < 0.05).sum()))
# implied ge from p: p=(ge+1)/(n_null+1)
for row in [gi]:
    ge = row.p_permutation_exact * (row.n_null + 1) - 1
    print("implied ge for glial/immune = %.2f (n_null=%d)" % (ge, row.n_null))

print("\n=== 7. permutation p floors in the CSV ===")
un = pc[~pc.paired_control]; pr = pc[pc.paired_control]
print("unpaired: min p = %.7f ; 3/4901 = %.7f ; n at floor = %d" % (
    un.p_permutation_exact.min(), 3 / 4901,
    int((un.p_permutation_exact <= 3/4901 + 1e-12).sum())))
print("paired:   min p = %.7f ; 2/71 = %.7f ; 3/71 = %.7f ; n at 2/71 = %d ; n at 3/71 = %d" % (
    pr.p_permutation_exact.min(), 2/71, 3/71,
    int((pr.p_permutation_exact <= 2/71 + 1e-12).sum()),
    int(((pr.p_permutation_exact - 3/71).abs() < 1e-12).sum())))
print("paired p value counts:")
print(pr.p_permutation_exact.value_counts().head(8).to_string())

print("\n=== 8. synthetic check: 8-case exchange-null floor ===")
# replicate scripts/22 paired-null code path on synthetic data
rng = np.random.default_rng(7)
ge_counts = {}
for trial in range(300):
    ng = 60
    dA = rng.normal(size=ng); dB = rng.normal(size=ng)
    c0a = rng.normal(size=ng); c0b = rng.normal(size=ng)
    A = rng.normal(size=(ng, 4)); B = rng.normal(size=(ng, 4))
    C8 = np.hstack([A, B])
    from itertools import combinations
    combs = list(combinations(range(8), 4))
    null = np.empty(len(combs))
    for r_i, c in enumerate(combs):
        s1 = np.zeros(8, bool); s1[list(c)] = True
        ra = np.corrcoef(dA, c0a - C8[:, s1].mean(axis=1))[0, 1]
        rb = np.corrcoef(dB, c0b - C8[:, ~s1].mean(axis=1))[0, 1]
        null[r_i] = ra - rb
    s_true = np.zeros(8, bool); s_true[:4] = True
    obs_i = combs.index((0, 1, 2, 3))
    obs = null[obs_i]
    ge = int((np.abs(null) >= abs(obs) - 1e-12).sum())
    ge_counts[ge] = ge_counts.get(ge, 0) + 1
print("distribution of ge (count of |null|>=|obs|) over 300 synthetic trials:")
print("  ", dict(sorted(ge_counts.items())))
print("  -> if ge=1 occurs, attainable floor is 2/71=%.5f, NOT the manuscript's 3/71" % (2/71))
# also: complement split value vs -obs
combs = list(combinations(range(8), 4))
ng = 60
dA = rng.normal(size=ng); dB = rng.normal(size=ng)
c0a = rng.normal(size=ng); c0b = rng.normal(size=ng)
A = rng.normal(size=(ng, 4)); B = rng.normal(size=(ng, 4))
C8 = np.hstack([A, B])
s1 = np.zeros(8, bool); s1[:4] = True
ra = np.corrcoef(dA, c0a - C8[:, s1].mean(axis=1))[0, 1]
rb = np.corrcoef(dB, c0b - C8[:, ~s1].mean(axis=1))[0, 1]
obs = ra - rb
s2 = ~s1
ra2 = np.corrcoef(dA, c0a - C8[:, s2].mean(axis=1))[0, 1]
rb2 = np.corrcoef(dB, c0b - C8[:, s2].mean(axis=1))[0, 1]
print("complement-split null = %.5f vs -obs = %.5f (differ -> no mirror in 8-case null)" % (
    ra2 - rb2, -obs))

print("\n=== 9. unpaired mirror-in-null check (synthetic) ===")
from itertools import combinations as comb2
n_ok = 0
for trial in range(50):
    ng = 40
    dA = rng.normal(size=ng); dB = rng.normal(size=ng)
    La = rng.normal(size=(ng, 8)); Lb = rng.normal(size=(ng, 8))
    cbs = list(comb2(range(8), 4))
    sa = np.array([np.corrcoef(dA, La[:, [j for j in range(8) if j not in c]].mean(axis=1)
                               - La[:, list(c)].mean(axis=1))[0, 1] for c in cbs])
    sb = np.array([np.corrcoef(dB, Lb[:, [j for j in range(8) if j not in c]].mean(axis=1)
                               - Lb[:, list(c)].mean(axis=1))[0, 1] for c in cbs])
    obs_i = cbs.index((0, 1, 2, 3))
    null = (sa[:, None] - sb[None, :]).ravel()
    obs = sa[obs_i] - sb[obs_i]
    ge = int((np.abs(null) >= abs(obs) - 1e-12).sum())
    if ge >= 2:
        n_ok += 1
print("ge>=2 in %d/50 trials (mirror pair always present -> floor 3/4901)" % n_ok)

print("\n=== 10. bootstrap CI checks ===")
print("max n_degenerate_boot =", int(bt.n_degenerate_boot.max()))
cover = ((pc.ci95_lo <= pc.delta_observed) & (pc.delta_observed <= pc.ci95_hi)).mean()
print("CI covers observed (132):", cover)
n05 = int((pc.p_permutation_exact < 0.05).sum())
n_ci = int((((pc.ci95_lo > 0) & (pc.ci95_hi > 0)) |
            ((pc.ci95_lo < 0) & (pc.ci95_hi < 0)))[pc.p_permutation_exact < 0.05].sum())
print("among %d exact-p<0.05 tests, CI excludes zero: %d" % (n05, n_ci))
# asymmetry of percentile CI (glial/immune)
print("glial/immune CI asymmetry: down=%.4f up=%.4f (point %.4f)" % (
    gi.delta_observed - gi.ci95_lo, gi.ci95_hi - gi.delta_observed, gi.delta_observed))

print("\n=== 11. network vs additive rho universes ===")
h2h = json.load(open(os.path.join(TAB, "head_to_head_baselines.json")))
print("discovery head_to_head spearman_network_vs_additive =",
      h2h.get("spearman_network_vs_additive"), "| keys:", list(h2h.keys())[:12])

print("\n=== 12. bootstrap mean vs delta_observed coherence ===")
bias = (pc.delta_boot_mean - pc.delta_observed).abs()
print("max |boot_mean - obs| = %.4f ; median = %.4f" % (bias.max(), bias.median()))
