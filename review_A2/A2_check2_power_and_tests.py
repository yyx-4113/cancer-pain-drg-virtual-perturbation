#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
A2_check2_power_and_tests.py
  (A) KO/OE mirror property  -- recomputed from deposited tables
  (B) the 264 -> 132 collapse -- is it the right correction?
  (C) the reversal-score MDE  -- is "2 vs 2, d = 5.65" the right test?
  (D) T1 target-vs-control discrimination test
  (E) T2 randset selectivity
  (F) T3/T4
"""
import os
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.power import TTestIndPower, TTestPower

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAB = os.path.join(ROOT, "results", "tables")

pd.set_option("display.width", 250)


def bh(p):
    p = np.asarray(p, float); n = len(p)
    o = np.argsort(p); q = p[o] * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n); out[o] = np.clip(q, 0, 1); return out


def hdr(s):
    print("\n" + "=" * 84); print(s); print("=" * 84)


bt = pd.read_csv(os.path.join(TAB, "etiology_bootstrap_differences.csv"))
rp = pd.read_csv(os.path.join(TAB, "etiology_reversal_point_estimates.csv"))

# =====================================================================  (A)
hdr("(A) KO / OE MIRROR PROPERTY -- recomputed independently")
pc = bt[bt.network_scope == "per_contrast"]
ko = pc[pc.perturbation == "knockout"].set_index(["module", "aetiology_A", "aetiology_B"])
oe = pc[pc.perturbation == "overexpression"].set_index(["module", "aetiology_A", "aetiology_B"])
idx = ko.index.intersection(oe.index)
print(f"matched (module, pair) cells: {len(idx)}")

d_rev = (ko.loc[idx, "rev_A"].values + oe.loc[idx, "rev_A"].values)
d_revB = (ko.loc[idx, "rev_B"].values + oe.loc[idx, "rev_B"].values)
d_obs = (ko.loc[idx, "delta_observed"].values + oe.loc[idx, "delta_observed"].values)
d_p = np.abs(ko.loc[idx, "p_permutation_exact"].values - oe.loc[idx, "p_permutation_exact"].values)
d_z = np.abs(ko.loc[idx, "z_vs_randset"].values + oe.loc[idx, "z_vs_randset"].values)
print(f"max |rev_A(KO) + rev_A(OE)|        = {np.abs(d_rev).max():.3e}")
print(f"max |rev_B(KO) + rev_B(OE)|        = {np.abs(d_revB).max():.3e}")
print(f"max |delta(KO) + delta(OE)|        = {np.abs(d_obs).max():.3e}")
print(f"max |p_perm(KO) - p_perm(OE)|      = {d_p.max():.3e}")
print(f"max |z(KO) + z(OE)|                = {d_z.max():.3e}")
print(f"max |z(KO) - z(OE)| (NOT mirrors)  = "
      f"{np.abs(ko.loc[idx,'z_vs_randset'].values - oe.loc[idx,'z_vs_randset'].values).max():.3e}")
print("=> reversal and the exact permutation p ARE exact mirrors; the randset z is NOT")
print("   (different degree-matched draws), exactly as the code comments admit.")

# also from the point-estimate table
rcols = [c for c in rp.columns if c.startswith("rev_")]
k = rp[rp.perturbation == "knockout"].set_index("module")
o = rp[rp.perturbation == "overexpression"].set_index("module")
print("\nfrom etiology_reversal_point_estimates.csv (rev_* columns):")
for c in rcols:
    print(f"   max |rev({c}) KO + OE| = {np.abs(k[c] + o[c]).max():.3e}")

# =====================================================================  (B)
hdr("(B) THE 264 -> 132 COLLAPSE")
print(f"rows in deposited table                      : {len(bt)}")
print(f"  per_contrast  x 2 perturbations           : {len(pc)}")
print(f"  distinct (module, aetiology_A,B) per scope: {pc.groupby(['module','aetiology_A','aetiology_B']).ngroups}")
print(f"  modules                                    : {pc.module.nunique()}")
print(f"  aetiology pairs                            : "
      f"{pc.groupby(['aetiology_A','aetiology_B']).ngroups}")
print(f"22 modules x 6 pairs x 2 perturbations       = {22*6*2}")
print(f"\nManuscript §3.6: '22 modules x 6 aetiology pairs x 2 perturbation")
print(f"directions = 264 contrast-level tests'.")
print(f"Manuscript abstract: 'the 264 contrast-level tests reduce to 132")
print(f"independent hypotheses PER NETWORK SCOPE'.")
print(f"Deposited table actually holds {len(bt)} rows = 264 per scope x 2 scopes,")
print(f"so 528 tests were run and 264 hypotheses (132/scope) were BH-corrected per scope.")

print("\n-- how correlated are the 132 'independent' reversal hypotheses? --")
w = pc[pc.perturbation == "knockout"].pivot_table(
    index="module", columns=["aetiology_A", "aetiology_B"], values="delta_observed")
print(f"matrix {w.shape}")
cc = w.corr(method="pearson").values
iu = np.triu_indices_from(cc, 1)
print(f"pairwise Pearson r among the 132 deltas: min={cc[iu].min():.3f} "
      f"median={np.median(cc[iu]):.3f} max={cc[iu].max():.3f}")
print(f"fraction of pairs with |r| > 0.5: {(np.abs(cc[iu])>0.5).mean():.3f}")
wz = pc[pc.perturbation == "knockout"].pivot_table(
    index="module", columns=["aetiology_A", "aetiology_B"], values="z_vs_randset")
cz = wz.corr().values
print(f"same for z_vs_randset: min={cz[iu].min():.3f} median={np.median(cz[iu]):.3f} "
      f"max={cz[iu].max():.3f}; frac |r|>0.5 = {(np.abs(cz[iu])>0.5).mean():.3f}")
ev = np.linalg.eigvalsh(np.corrcoef(w.values))
print(f"eigenvalues of the 22x22 delta correlation matrix: "
      f"top5={np.round(ev[::-1][:5],2)}  min={ev.min():.3f}")
print(f"effective rank (participation ratio) = "
      f"{(ev.sum()**2)/(ev**2).sum():.2f} of 22")
# how many INDEPENDENT contrasts do 6 pairs of 4 aetiologies really carry?
print("\n-- algebraic redundancy among the 6 aetiology pairs --")
print("With 4 aetiologies there are only 3 independent pairwise contrasts;")
print("6 pairs = C(4,2) are linearly dependent: d(A,B) = d(A,C) - d(B,C).")
for m in ["Glial / immune", "CTRL erythrocyte"]:
    v = w.loc[m]
    print(f"   {m}:")
    for (a, b), val in v.items():
        print(f"      d({a:9s}|{b:9s}) = {val:+.4f}")
    d_dp = v[("Bone_D14", "PTX")] - v[("Bone_D7", "PTX")]
    d_14_7 = v[("Bone_D14", "Bone_D7")]
    d_7_ptx = v[("Bone_D7", "PTX")]
    print(f"      check d(D14,PTX) - d(D7,PTX) = {d_dp:+.4f}  vs  d(D14,D7) = {d_14_7:+.4f}"
          f"   (ratio {d_dp/d_14_7:.4f})")

print("\n-- what BH assumes vs what holds --")
print("BH controls FDR under independence or PRDS. The 132 'independent'")
print("hypotheses share 22 modules scored on the SAME 4 contrasts and the SAME")
print("6 pair-differences, so PRDS is violated; see correlations above.")

# =====================================================================  (C)
hdr("(C) THE REVERSAL-SCORE MDE  (manuscript's strongest new statistical claim)")
tp = TTestIndPower()
d_ind2 = tp.solve_power(nobs1=2, alpha=0.05, power=0.8, ratio=1.0,
                        alternative="two-sided")
d_ind4 = tp.solve_power(nobs1=4, alpha=0.05, power=0.8, ratio=1.0,
                        alternative="two-sided")
d_pair4 = TTestPower().solve_power(nobs=4, alpha=0.05, power=0.8,
                                   alternative="two-sided")
print(f"TTestIndPower  nobs1=2, ratio=1  -> d = {d_ind2:.4f}   (repo T5 = 5.65)")
print(f"TTestIndPower  nobs1=4, ratio=1  -> d = {d_ind4:.4f}")
print(f"TTestPower (PAIRED) nobs=4       -> d = {d_pair4:.4f}")

print("\n-- is the KO/OE contrast a 2-vs-2 comparison? --")
print("reversal(KO) = cos(P@(-seed), axis) = -cos(P@seed, axis) = -reversal(OE)")
print("Both are deterministic functions of the SAME 4 case + 4 control samples.")
print("A two-sample t-test on {KO, OE} is degenerate: the two 'groups' have")
print("correlation exactly -1 and zero independent within-group variance")
print("relative to their difference. There is no '2 vs 2' design here.")
v = w.stack()
kk = k[rcols].values.ravel(); oo = o[rcols].values.ravel()
print(f"   corr(KO, OE) across all {len(kk)} module x contrast scores = "
      f"{np.corrcoef(kk, oo)[0,1]:.6f}")
print(f"   every OE value is the exact negation of its KO value: "
      f"{np.allclose(oo, -kk)}")

print("\n-- what is rev_sd = 0.075 ? --")
for c in rcols:
    print(f"   SD across the 22 MODULES of rev({c}) = {k[c].std():.4f}   "
          f"[pandas .std() = ddof=1 across rows = BETWEEN-module spread]")
rev_sd = float(np.mean([k[c].std() for c in rcols]))
print(f"   mean of those = {rev_sd:.4f}  <- this is the repo's 'observed_reversal_score_sd'")
print("   It is a BETWEEN-MODULE dispersion (real biology + noise), NOT a")
print("   within-module sampling SD of the reversal score.")
print("\n   What the MDE should be scaled by: the SE of the reversal score under")
print("   resampling of the 4+4 animals. From the deposited bootstrap CI we can")
print("   back out an empirical SE of Delta:")
pc_ko = pc[pc.perturbation == "knockout"].copy()
pc_ko["se_from_ci"] = (pc_ko.ci95_hi - pc_ko.ci95_lo) / (2 * 1.959964)
print(f"     median SE(Delta) implied by bootstrap CI = "
      f"{pc_ko.se_from_ci.median():.4f}  (mean {pc_ko.se_from_ci.mean():.4f})")
print(f"     vs the repo's rev_sd = {rev_sd:.4f}  "
      f"-> ratio {pc_ko.se_from_ci.median()/rev_sd:.2f}x")

print("\n-- largest observed reversal difference --")
mx = pc_ko.loc[pc_ko.delta_observed.abs().idxmax()]
print(f"   max |delta_observed| over 132 per_contrast KO tests = "
      f"{pc_ko.delta_observed.abs().max():.4f}")
print(f"   attained by {mx['module']} ({mx['aetiology_A']} vs {mx['aetiology_B']}), "
      f"CI [{mx['ci95_lo']:.3f}, {mx['ci95_hi']:.3f}], p_perm={mx['p_permutation_exact']:.5f}, "
      f"q={mx['q_permutation']:.4f}")
print(f"   manuscript §3.6 says 'largest observed reversal difference is = 0.27'  OK")
print(f"   but the T5 JSON interpretation string says 'reversal differences are <0.2'")
print(f"   -> the artefact and the manuscript disagree with each other.")
print(f"   also: abs(delta) > 0.2 count = {(pc_ko.delta_observed.abs()>0.2).sum()}, "
      f"> 0.25 count = {(pc_ko.delta_observed.abs()>0.25).sum()}")

print("\n-- MDE restated on the correct scale --")
print(f"   repo:  d=5.65 x rev_sd(between-module) 0.075 = {5.65*rev_sd:.3f} 'reversal units'")
print(f"   unpaired 4v4 (correct n for a cross-aetiology contrast):")
print(f"          d = {d_ind4:.2f} x SE(Delta) {pc_ko.se_from_ci.median():.4f} = "
      f"{d_ind4*pc_ko.se_from_ci.median():.4f}")
print(f"   paired 4-pair (the two shared-control contrasts ARE paired):")
print(f"          d = {d_pair4:.2f} x SE(Delta) {pc_ko.se_from_ci.median():.4f} = "
      f"{d_pair4*pc_ko.se_from_ci.median():.4f}")
print("   => with the empirical SE the design detects ~0.03-0.07 reversal units at")
print("      80% power, NOT 0.42. The largest observed difference (0.274) is")
print("      roughly 4-9 x the correctly-scaled MDE. The paper's 'grossly")
print("      underpowered' conclusion does not follow from its own artefacts.")

# =====================================================================  (D)
hdr("(D) T1 -- TARGET vs CONTROL DISCRIMINATION")
print("repo T1: MWU(target z, control z, alternative='greater') over the 264")
print("per_contrast rows (204 target + 60 control).")
tgt = pc[pc.kind == "target"]; ctl = pc[pc.kind == "control"]
tgt_abs = tgt.z_vs_randset.abs(); ctl_abs = ctl.z_vs_randset.abs()
u1 = stats.mannwhitneyu(tgt.z_vs_randset, ctl.z_vs_randset, alternative="greater")
u2 = stats.mannwhitneyu(tgt_abs, ctl_abs, alternative="greater")
u3 = stats.mannwhitneyu(tgt.z_vs_randset, ctl.z_vs_randset, alternative="two-sided")
u4 = stats.mannwhitneyu(tgt_abs, ctl_abs, alternative="two-sided")
print(f"  one-sided (target>control) signed z : p = {u1.pvalue:.4f}  [repo 0.4973]")
print(f"  one-sided (target>control) |z|      : p = {u2.pvalue:.4f}  [repo 1.0]")
print(f"  TWO-SIDED signed z                  : p = {u3.pvalue:.4f}")
print(f"  TWO-SIDED |z|                       : p = {u4.pvalue:.4f}")
print(f"  n = {len(tgt)} target rows, {len(ctl)} control rows")
print(f"  mean |z| target {tgt.z_vs_randset.abs().mean():.4f} vs "
      f"control {ctl.z_vs_randset.abs().mean():.4f}")
lt = stats.ttest_ind(tgt.z_vs_randset.abs(), ctl.z_vs_randset.abs(), equal_var=False)
print(f"  Welch t on |z| (target vs control)  : t={lt.statistic:.3f} p={lt.pvalue:.3g}")
print("  -> the DIRECTION of the repo's alternative ('target > control') is the")
print("     direction the data REFUTE. Reporting one-sided p=0.50/1.00 and calling")
print("     it 'indistinguishable' hides that the two-sided test rejects strongly.")
print("     'Statistically indistinguishable' is also not what a non-significant")
print("     MWU licenses; it licenses 'no evidence of a difference at this n'.")

print("\n-- pseudoreplication: the 264 rows are not independent observations --")
print(f"  distinct animals in the whole design: 16+12 = 28 samples, n=4/group")
print(f"  effective units for T1: modules (22) x independent contrasts (3), not 264")
print("  Each KO/OE mirror pair is duplicated in the 264 rows, so every")
print("  independent observation is counted exactly twice.")
print(f"  rows = {len(pc)}, dedup (KO only) rows = {len(pc_ko)}")
u_dedup = stats.mannwhitneyu(pc_ko[pc_ko.kind=='target'].z_vs_randset.abs(),
                             pc_ko[pc_ko.kind=='control'].z_vs_randset.abs(),
                             alternative="two-sided")
print(f"  two-sided MWU on |z| using KO rows only (132): p = {u_dedup.pvalue:.3g}")

print("\n-- T4 rank percentile --")
a = pc.copy(); a["absz"] = a.z_vs_randset.abs()
a["rank_abs"] = a.absz.rank(ascending=False, method="average")
a["pct"] = (a["rank_abs"] - 1) / (len(a) - 1) * 100
print(f"  target median percentile by |z| = {a[a.kind=='target'].pct.median():.1f} "
      f"[repo 55.7]")
print(f"  control median percentile by |z| = {a[a.kind=='control'].pct.median():.1f}")
print(f"  => controls sit ABOVE the median, i.e. above the 17/22 = 77.3rd pct if")
print(f"     sorted ascending. 'Analgesic targets sit at the median' is true but")
print(f"     it is a statement about a distribution that is 77% targets by design:")
print(f"     {pc.kind.eq('target').mean()*100:.1f}% of rows ARE targets.")

# =====================================================================  (E)
hdr("(E) T2 -- 'THE RANDSET NULL IS NON-SELECTIVE'")
sig = pc[pc.q_randset < 0.05]
print(f"randset-significant rows (both perturbations): {len(sig)}")
print(f"  targets  : {(sig.kind=='target').sum()}")
print(f"  controls : {(sig.kind=='control').sum()}")
sk = pc_ko[pc_ko.q_randset < 0.05]
print(f"deduplicated (KO only): {len(sk)}  targets {(sk.kind=='target').sum()} "
      f"controls {(sk.kind=='control').sum()}")
print(f"repo T2 JSON reports total=64, targets=36, controls=28 -> matches the 264 universe")
excl = (((pc.ci95_lo > 0) & (pc.ci95_hi > 0)) | ((pc.ci95_lo < 0) & (pc.ci95_hi < 0)))
print(f"rows with bootstrap CI excluding 0: {int(excl.sum())} / {len(pc)}")
print(f"   among the 64 randset-significant: {int(excl[sig.index].sum())}")
print("\n-- the CI is on Delta, the randset p is on Delta vs a gene-set null --")
print("   They test DIFFERENT nulls. A bootstrap CI excluding 0 says Delta != 0;")
print("   the randset p says |Delta| exceeds what degree-matched random sets give.")
print("   Every one of the 264 rows having a CI excluding 0 (n=4, bounded cosine)")
print("   carries no information about selectivity of the randset null.")
sk_ex = (((sk.ci95_lo > 0) & (sk.ci95_hi > 0)) | ((sk.ci95_lo < 0) & (sk.ci95_hi < 0)))
print(f"   dedup: {int(sk_ex.sum())}/{len(sk)} randset-significant have CI excluding 0")
print(f"   BUT: rows NOT randset-significant that also have CI excluding 0: "
      f"{int(excl.sum()) - int(excl[sig.index].sum())}")
print("\n-- base rate: what SHOULD 'non-selective' look like? --")
nmod_t = pc[pc.kind == 'target'].module.nunique(); nmod_c = pc[pc.kind == 'control'].module.nunique()
print(f"   universe: {nmod_t} target modules + {nmod_c} control modules; "
      f"control base rate = {nmod_c}/{nmod_t+nmod_c} = {nmod_c/(nmod_t+nmod_c):.3f}")
frac_c = (sig.kind == 'control').mean()
print(f"   observed control fraction among randset-significant = {frac_c:.3f}")
print(f"   Binomial test H0: p_control = {nmod_c/(nmod_t+nmod_c):.3f}, "
      f"n = {len(sig)}, k = {(sig.kind=='control').sum()}")
print(f"     two-sided binomial p = "
      f"{stats.binomtest(int((sig.kind=='control').sum()), len(sig), nmod_c/(nmod_t+nmod_c)).pvalue:.4f}")
print("   => the randset null flags controls at ~2x the base rate. That is")
print("      ENRICHMENT of controls, i.e. the opposite of 'non-selective' in the")
print("      sense of 'no preference'; the paper's phrasing conflates two things.")

# =====================================================================  (F)
hdr("(F) T3 -- SINGLE-CELL POWER SPLIT")
sc = pd.read_csv(os.path.join(TAB, "sc_GSE292522_pseudobulk_DE.csv"), index_col=0)
comp = pd.read_csv(os.path.join(TAB, "sc_GSE292522_composition.csv"))
print(f"pseudobulk DE: {sc.shape[0]} genes x {sc.shape[1]} cols")
stable = sc[sc.q > 0.5]
se = stable.logFC.abs() / stable.t.abs()
print(f"genes with q>0.5: {len(stable)}")
print(f"repo median SE of stable genes = {se.median():.6f}")
print(f"   quantiles of that SE: " +
      ", ".join(f"{q_}%={np.nanpercentile(se.dropna(), q_):.5f}"
                for q_ in [5, 25, 50, 75, 95]))
print(f"   -> the median is 0.0003 log2 units, i.e. the pipeline reports SEs of")
print(f"      3e-4. Multiplying by d=2.38 gives 'detectable |logFC| = 0.001'.")
print(f"   DETECTED GENES: {(sc.q<0.05).sum()} at q<0.05, so the design is NOT")
print(f"   limited to |logFC|>=0.001; the median-SE argument is circular.")
print(f"   median |logFC| among detected genes = {sc[sc.q<0.05].logFC.abs().median():.4f}")
print(f"   median |logFC| among stable genes   = {stable.logFC.abs().median():.6f}")
print(f"\n   Ifit2: " + sc.loc[[i for i in sc.index if i == 'Ifit2']].to_string()
      if 'Ifit2' in sc.index else "   Ifit2 not in index")
print(f"\n   Is 'detecting Ifit2 at logFC -0.048' evidence of power?")
ifit = sc.loc['Ifit2'] if 'Ifit2' in sc.index else None
if ifit is not None:
    print(f"      Ifit2 logFC={ifit.logFC:.4f} t={ifit.t:.3f} q={ifit.q:.4f}")
    print(f"      A q=0.042 hit at |logFC|=0.048 is a TENTHS-OF-A-PERCENT change;")
    print(f"      calling this 'well powered' and then 'a genuine null for large")
    print(f"      effects' in Cd68/Aif1 conflates 'sensitive to tiny effects' with")
    print(f"      'able to exclude a 2-3 log2-unit change in Cd68/Aif1'.")
    print(f"      Cd68/Aif1 in this dataset: ")
    for g in ["Cd68", "Aif1", "Gfap", "Atf3", "Sprr1a"]:
        if g in sc.index:
            r = sc.loc[g]
            print(f"        {g:8s} logFC={r.logFC:+.4f} t={r.t:+.2f} q={r.q:.4f} "
                  f"s2_post={r.s2_post:.2e}")
    print(f"      => the right object is the CI on Cd68/Aif1 logFC, not the median SE.")

print("\n-- Cohen's d applied to a bounded fraction --")
con = comp[comp.group == "Con"].Myeloid.values
b16 = comp[comp.group == "B16F10"].Myeloid.values
pooled = np.sqrt(((con - con.mean())**2).sum() / 7 + ((b16 - b16.mean())**2).sum() / 7)
d_obs = (con.mean() - b16.mean()) / pooled
print(f"   Con  fractions: {np.round(con,6)}  mean={con.mean():.6f}")
print(f"   B16  fractions: {np.round(b16,6)}  mean={b16.mean():.6f}")
print(f"   pooled SD = {pooled:.6f}   d = {d_obs:.4f}  [repo 1.37]")
print(f"   NOTE: denominator 7 = (4-1)+(4-1) = 6 would be the usual pooled SD;")
print(f"   the code uses 7, which is neither. Correct pooled SD = "
      f"{np.sqrt(((con-con.mean())**2).sum()+((b16-b16.mean())**2).sum())/6:.6f}")
pooled_std6 = np.sqrt(((con - con.mean())**2).sum() +
                      ((b16 - b16.mean())**2).sum()) / 6
d6 = (con.mean() - b16.mean()) / pooled_std6
print(f"   -> d with the conventional denominator = {d6:.4f}")
print(f"   MDE d (4v4 unpaired, t-test) = {d_ind4:.2f}; in fraction units = "
      f"{d_ind4*pooled:.6f}")
print(f"   observed mean diff = {con.mean()-b16.mean():.6f}")
mwu = stats.mannwhitneyu(con, b16, alternative="two-sided")
mwu1 = stats.mannwhitneyu(con, b16, alternative="greater")
print(f"   MWU two-sided p = {mwu.pvalue:.4f}  [manuscript 0.114]")
print(f"   MWU one-sided (Con>B16) p = {mwu1.pvalue:.4f}  [manuscript 0.057]")
print(f"   exact permutation floor for 4v4 two-sided = 2/70 = {2/70:.4f}")
print(f"   ATTEMPTED one-sided testing after seeing the direction: the manuscript")
print(f"   reports BOTH p=0.114 and p=0.057; only the two-sided is legitimate")
print(f"   unless the one-sided direction was pre-specified.")
print(f"   skewness of the pooled fraction: {stats.skew(np.r_[con,b16]):.3f} "
      f"(bounded in [0,1]; values ~5e-4)")
lc, lbg = np.log(con), np.log(b16)
dlog = ((lc.mean() - lbg.mean()) /
        np.sqrt(((lc - lc.mean())**2).sum() + ((lbg - lbg.mean())**2).sum()) / 6)
print(f"   log-transformed: d = {dlog:.4f}  (vs {d_obs:.4f} on the raw fraction)")


if __name__ == "__main__":
    pass
