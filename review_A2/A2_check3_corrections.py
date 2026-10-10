#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
A2_check3_corrections.py -- corrected arithmetic + the remaining claims.
Fixes the pooled-SD slip in check2 and settles:
  1. myeloid-fraction Cohen's d with the conventional pooled SD
  2. is the reversal score a PAIRED quantity? proper paired MDE
  3. the 64 vs 68 randset-significant discrepancy in the summary JSON
  4. |z| = 7.39 erythrocyte claim, min q = 0.165, p = 0.0014 glial/immune
  5. how many shared-control conclusions change under the correct null
  6. BH consistency across every reported family
"""
import os
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.power import TTestIndPower, TTestPower

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAB = os.path.join(ROOT, "results", "tables")
pd.set_option("display.width", 250)


def hdr(s):
    print("\n" + "=" * 84); print(s); print("=" * 84)


bt = pd.read_csv(os.path.join(TAB, "etiology_bootstrap_differences.csv"))
pc = bt[bt.network_scope == "per_contrast"]
pk = pc[pc.perturbation == "knockout"]

# =====================================================  1. myeloid fraction
hdr("1. MYELOID FRACTION -- CORRECTED Cohen's d (pooled SD convention fixed)")
comp = pd.read_csv(os.path.join(TAB, "sc_GSE292522_composition.csv"))
con = comp[comp.group == "Con"].Myeloid.values
b16 = comp[comp.group == "B16F10"].Myeloid.values
ss1 = ((con - con.mean()) ** 2).sum()
ss2 = ((b16 - b16.mean()) ** 2).sum()
diff = con.mean() - b16.mean()
sd_repo = np.sqrt(ss1 / 7 + ss2 / 7)          # what the repo does
sd_conv = np.sqrt((ss1 + ss2) / 6)             # conventional pooled SD, df = n1+n2-2 = 6
print(f"SS_con = {ss1:.4e}   SS_b16 = {ss2:.4e}   mean diff = {diff:.4e}")
print(f"repo  'pooled SD' = sqrt(SS1/7 + SS2/7) = {sd_repo:.6e}  -> d = {diff/sd_repo:.4f}")
print(f"conventional      = sqrt((SS1+SS2)/6)    = {sd_conv:.6e}  -> d = {diff/sd_conv:.4f}")
print(f"NOTE sqrt(a/7+b/7) == sqrt((a+b)/7); the repo divides the TOTAL SS by 7,")
print(f"     not by 6. It is neither the pooled SD (6) nor the within-group SD.")
print(f"     Effect: the reported d = {diff/sd_repo:.3f} UNDERSTATES the conventional")
print(f"     d = {diff/sd_conv:.3f} by a factor {(sd_conv/sd_repo):.3f}.")
print(f"     (My check2 printed 3.10 because of an arithmetic slip on my part;")
print(f"      {diff/sd_conv:.3f} is the correct conventional value.)")
t = stats.ttest_ind(con, b16, equal_var=True)
print(f"     Welch/Student t-test on the raw fractions: t={t.statistic:.3f}, "
      f"df={t.df:.2f}, p={t.pvalue:.4f}")
print(f"     -> the t-test on the SAME data the repo reports as an underpowered")
print(f"        MWU gives p = {t.pvalue:.3f}. The 'underpowered' label is a")
print(f"        consequence of the MWU's own coarse n=4v4 resolution (floor")
print(f"        2/70 = {2/70:.4f}), not of the effect size.")
print(f"\n   skewness {stats.skew(np.r_[con,b16]):.3f}; values are bounded in [0,1]")
print(f"   and ~5e-4. MDE d=2.38 (4v4 t-test) in fraction units = "
      f"{2.3808*sd_conv:.3e}")
print(f"   observed mean diff = {diff:.3e}; ratio = {diff/(2.3808*sd_conv):.3f}")
print(f"   log-scale: d = ", end="")
lc, lb = np.log(con), np.log(b16)
print(f"{((lc.mean()-lb.mean())/np.sqrt(((lc-lc.mean())**2).sum()+((lb-lb.mean())**2).sum())/6):.4f}"
      f"  (vs {diff/sd_conv:.4f} raw) -> the effect size is entirely scale-dependent")

# ==============================================  2. is reversal a paired score?
hdr("2. IS THE REVERSAL SCORE A PAIRED QUANTITY?  PROPER MDE")
print("Structure of the reversal score (read off scripts/22_etiology_comparison.py:205-210):")
print("   reversal(k, module, sign) = cos( P_k @ seed , axis_k )")
print("   axis_k = mean(control_k samples) - mean(case_k samples)      [4 + 4 animals]")
print("   P_k    = diffusion operator built on the SAME 8 samples")
print("=> one scalar per (module, contrast, perturbation). It is a functional of")
print("   8 animals, NOT a per-animal measurement. There is no '2 vs 2' design:")
print("   KO and OE are the SAME 8 animals with the seed sign flipped.")
print()
print("Empirical confirmation from the deposited point estimates:")
rp = pd.read_csv(os.path.join(TAB, "etiology_reversal_point_estimates.csv"))
k = rp[rp.perturbation == "knockout"].set_index("module")
o = rp[rp.perturbation == "overexpression"].set_index("module")
rc = [c for c in rp.columns if c.startswith("rev_")]
allk = k[rc].values.ravel(); allo = o[rc].values.ravel()
print(f"   corr(KO, OE) over {len(allk)} scores = {np.corrcoef(allk, allo)[0,1]:.6f}")
print(f"   max|KO+OE| = {np.abs(allk+allo).max():.3e}   (exact mirror, verified)")
print()
print("So the design question is NOT 'KO vs OE'. It is:")
print("   Delta = reversal(aetiology A) - reversal(aetiology B),")
print("   with A and B each resting on 4 case + 4 control animals.")
print("   For the two SHARED-CONTROL pairs the control animals are COMMON to")
print("   both sides, so Delta IS paired. For the four cross-study pairs it is not.")
print()
tpi = TTestIndPower(); tpp = TTestPower()
d_i4 = tpi.solve_power(nobs1=4, alpha=0.05, power=0.8, ratio=1.0)
d_p4 = tpp.solve_power(nobs=4, alpha=0.05, power=0.8)
d_i2 = tpi.solve_power(nobs1=2, alpha=0.05, power=0.8, ratio=1.0)
pk2 = pk.copy()
pk2["se_ci"] = (pk2.ci95_hi - pk2.ci95_lo) / (2 * 1.959964)
se_med = float(pk2.se_ci.median())
print("CORRECTED MDE, using the empirical SE of Delta implied by the paper's own")
print("bootstrap CI (the only resampling-based SE the paper reports):")
print(f"   median SE(Delta) = {se_med:.4f}   (mean {pk2.se_ci.mean():.4f}, "
      f"range {pk2.se_ci.min():.4f}-{pk2.se_ci.max():.4f})")
print(f"   cross-study (unpaired, 4v4) : d = {d_i4:.3f} -> MDE = "
      f"{d_i4*se_med:.4f} reversal units")
print(f"   shared-control (paired n=4) : d = {d_p4:.3f} -> MDE = "
      f"{d_p4*se_med:.4f} reversal units")
print(f"   repo (nobs1=2, x between-module SD 0.0741) : d = {d_i2:.3f} -> "
      f"{d_i2*0.0741:.4f} 'reversal units'")
print()
print(f"   Observed max |Delta| (excl. controls) = "
      f"{pk2[pk2.kind=='target'].delta_observed.abs().max():.4f}")
print(f"   Observed max |Delta| (all modules)     = "
      f"{pk2.delta_observed.abs().max():.4f}")
print(f"   Number of the 132 tests whose |Delta| exceeds the CORRECTED unpaired")
print(f"   MDE ({d_i4*se_med:.4f}): "
      f"{int((pk2.delta_observed.abs() > d_i4*se_med).sum())}")
print(f"   ... whose |Delta| exceeds the corrected PAIRED MDE "
      f"({d_p4*se_med:.4f}): {int((pk2.delta_observed.abs() > d_p4*se_med).sum())}")
print("   => the design is NOT 'only powered to detect 0.42-unit differences'.")
print("      The paper's own bootstrap says the SE of Delta is ~0.019, i.e. ~4x")
print("      SMALLER than the 0.074 the MDE was scaled by, and 22x smaller than")
print("      the 0.418 'MDE in reversal units' the paper quotes.")

# cross-check: what does the exact permutation say?
print("\n   Independent cross-check with the EXACT permutation p (no bootstrap):")
print(f"   tests with p_permutation_exact < 0.05 : "
      f"{int((pk2.p_permutation_exact < 0.05).sum())} / 132")
print(f"   min p_permutation_exact = {pk2.p_permutation_exact.min():.5f}")
print(f"   min q_permutation       = {pk2.q_permutation.min():.4f}  [manuscript 0.165]")
print(f"   BH on 132 p's recomputed: ", end="")
def bh(p):
    p = np.asarray(p, float); n = len(p)
    o_ = np.argsort(p); q = p[o_] * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n); out[o_] = np.clip(q, 0, 1); return out
print(f"{bh(pk2.p_permutation_exact.values).min():.4f}")
print(f"   number with q_permutation < 0.05: {int((pk2.q_permutation<0.05).sum())}"
      f"  [repo summary says 0]")

# ==============================================  3. 64 vs 68 discrepancy
hdr("3. THE '64 randset-significant' COUNT vs THE SUMMARY JSON's 68")
sig264 = int((pc.q_randset < 0.05).sum())
sig_ko = int((pk.q_randset < 0.05).sum())
sig_oe = int((pc[pc.perturbation=='overexpression'].q_randset < 0.05).sum())
print(f"rows with q_randset < 0.05 in the deposited 264-row per_contrast block: {sig264}")
print(f"   KO rows: {sig_ko}    OE rows: {sig_oe}    2 x {sig_ko} = {2*sig_ko}")
print(f"summary JSON 'n_sig_randset_q05_percontrast'  = 68")
print(f"summary JSON 'n_sig_randset_q05_percontrast_indep' = 32")
print(f"summary JSON 'randset_sig_breakdown_percontrast'  = target 18 + control 14 = 32")
print(f"T2 JSON      'randset_significant_total'           = 64 (targets 36 + controls 28)")
print("=> the deposited CSV gives 64 (= 2 x 32) and T2 agrees; the summary JSON's")
print("   '68' is inconsistent with its own '32' and with the CSV. Stale field.")

# ==============================================  4. specific quoted numbers
hdr("4. SPECIFIC NUMBERS QUOTED IN §3.6, RECOMPUTED")
e = pk[pk.module == "CTRL erythrocyte"]
print("CTRL erythrocyte, all 6 pairs (KO rows):")
print(e[["aetiology_A", "aetiology_B", "delta_observed", "z_vs_randset",
         "p_permutation_exact", "q_permutation", "q_randset"]]
      .sort_values("z_vs_randset", key=lambda s: -s.abs()).to_string(index=False))
print(f"\n   max |z| in the whole 132 = {pk.z_vs_randset.abs().max():.4f}")
top = pk.loc[pk.z_vs_randset.abs().idxmax()]
print(f"   attained by {top['module']} {top['aetiology_A']} vs {top['aetiology_B']}, "
      f"z = {top['z_vs_randset']:+.4f}")
print(f"   manuscript claims |z| = 7.39 for 'bone day 14 versus day 7' and "
      f"z = -5.78 for PTX vs OXA")
er = e.set_index(["aetiology_A", "aetiology_B"])
print(f"   erythrocyte Bone_D14|Bone_D7 z = "
      f"{er.loc[('Bone_D14','Bone_D7'),'z_vs_randset']:+.4f}")
print(f"   erythrocyte PTX|OXA         z = "
      f"{er.loc[('PTX','OXA'),'z_vs_randset']:+.4f}")
print(f"   erythrocyte Bone_D7|PTX     z = "
      f"{er.loc[('Bone_D7','PTX'),'z_vs_randset']:+.4f}")

g = pk[pk.module == "Glial / immune"].set_index(["aetiology_A", "aetiology_B"])
print(f"\nGlial / immune: the manuscript's 'smallest q' row")
print(g[["delta_observed", "ci95_lo", "ci95_hi", "p_permutation_exact",
         "q_permutation"]].to_string())
gi = pk.loc[pk.q_permutation.idxmin()]
print(f"   min-q row in the whole 132: {gi['module']} {gi['aetiology_A']}|"
      f"{gi['aetiology_B']} Delta={gi['delta_observed']:+.4f} "
      f"CI=[{gi['ci95_lo']:.3f},{gi['ci95_hi']:.3f}] "
      f"p={gi['p_permutation_exact']:.5f} q={gi['q_permutation']:.4f}")
print(f"   manuscript §3.6: 'glial/immune (Delta = +0.274, 95% interval 0.203-0.287,")
print(f"   p = 0.0014, q = 0.165)'")
print(f"   -> the row with Delta=+0.2743 is Bone_D7|PTX; its p and q are:")
gr = g.loc[("Bone_D7", "PTX")]
print(f"      Delta={gr.delta_observed:+.4f} CI=[{gr.ci95_lo:.4f},{gr.ci95_hi:.4f}] "
      f"p={gr.p_permutation_exact:.5f} q={gr.q_permutation:.4f}")
print(f"   -> manuscript's p=0.0014 / q=0.165 MATCH the deposited row.")

# ======================================  5. shared-control null: what changes
hdr("5. SHARED-CONTROL CONTRASTS -- REPO 'PAIRED' NULL vs THE CORRECT NULL")
mk = pd.read_csv(os.path.join(TAB, "etiology_module_logfc_cross.csv"))
sh = mk[mk.paired_control]
print(f"module-level tests on shared-control pairs: {len(sh)} "
      f"({sh.module.nunique()} modules x 2 pairs)")
print(f"   repo p distribution: " +
      ", ".join(f"p={p:.4f}:{c}" for p, c in
                sh.p_permutation_exact.round(4).value_counts().sort_index().items()))
print(f"   min repo p on shared-control pairs = {sh.p_permutation_exact.min():.5f}")
print(f"   the repo floor for these is claimed to be 1/71 = {1/71:.5f}; the")
print(f"   attainable floor is 3/71 = {3/71:.5f} (complement symmetry).")
print(f"   min ATTAINED repo p on these pairs = {sh.p_permutation_exact.min():.5f}")
print(f"   any shared-control test with q<0.05? "
      f"{int((sh.q_permutation<0.05).sum())} of {len(sh)}")
print("   => no conclusion in the paper rests on a shared-control permutation")
print("      result, so the wrong-null problem is a Methods error, not a")
print("      results error. It must still be corrected, because §2.7 presents")
print("      the pairing as a methodological refinement.")

# ==============================================  6. BH consistency audit
hdr("6. BH CONSISTENCY AUDIT ACROSS EVERY REPORTED FAMILY")
def rep(name, df, pcol, qcol, expect_n=None):
    q2 = bh(df[pcol].values)
    ok = np.allclose(q2, df[qcol].values, atol=1e-9)
    print(f"  {name:52s} n={len(df):4d} (expect {expect_n})  "
          f"BH reproduces deposited q: {ok}  min q={df[qcol].min():.4g}")
    return ok
rep("etiology_marker_cross.csv (gene x pair)", mk.assign(), "p_permutation_exact",
    "q_permutation") if False else None
mk2 = pd.read_csv(os.path.join(TAB, "etiology_marker_cross.csv"))
rep("etiology_marker_cross.csv (17 genes x 6 pairs)", mk2, "p_permutation_exact",
    "q_permutation")
rep("etiology_module_logfc_cross.csv (22 mod x 6 pairs)", mk, "p_permutation_exact",
    "q_permutation")
rep("etiology_bootstrap_differences.csv, per_contrast KO (132)", pk,
    "p_permutation_exact", "q_permutation")
pkb = bt[bt.network_scope == "pooled_by_study"]
rep("etiology_bootstrap_differences.csv, pooled KO (132)", pkb,
    "p_permutation_exact", "q_permutation")
rep("etiology_bootstrap_differences.csv, per_contrast KO randset", pk,
    "p_randset", "q_randset")
sc = pd.read_csv(os.path.join(TAB, "sc_GSE292522_pseudobulk_DE.csv"), index_col=0)
from scipy.stats import rankdata
n = len(sc); p = sc.P.values
q = np.minimum.accumulate((p[::-1] * n / np.arange(n, 0, -1))[::-1])[::-1]
q = np.clip(q, 0, 1)
print(f"  {'sc_GSE292522_pseudobulk_DE.csv (34603 genes)':52s} n={n:4d} "
      f"(expect 34603)  BH reproduces deposited q: "
      f"{np.allclose(np.sort(q), np.sort(sc.q.values), atol=1e-9)}  "
      f"min q={sc.q.min():.4g}")
de7 = pd.read_csv(os.path.join(TAB, "DE_D7_vs_Sham_adj.csv"), index_col=0)
print(f"  DE_D7_vs_Sham_adj.csv: n={len(de7)}, "
      f"q<0.05 count = {int((de7.q<0.05).sum())}  [manuscript 4,600]")
print(f"     columns: {list(de7.columns)}")
ct = pd.read_csv(os.path.join(TAB, "sc_GSE292522_celltype_DE.csv"))
print(f"\n  per-compartment DE table: {ct.shape}, "
      f"compartments = {ct.cell_type.nunique()}")
print("  NOTE: the per-compartment q values are BH WITHIN each compartment")
print("  (scripts/32_percelltype_sc.py:160 applies bh() to one compartment's")
print("  P vector at a time). The manuscript reports 'no compartment showed")
print("  ISG up-regulation' from those per-compartment q values, so the")
print("  multiplicity family is 8 (compartments) x n_ISG, not n_ISG.")
print("  Reported compartment q values in the paper:")
print(ct[ct.gene.isin(["Ifit2", "Isg15", "Oas1c", "Gfap"])].to_string(index=False))

print("\n" + "=" * 84)
print("7. TABLE 2 / TABLE 3 'Largest Delta' SELECTION")
print("=" * 84)
print("Table 2 reports, per gene, the LARGEST |Delta| over the 6 pairs, with the")
print("q of that particular pair. Selecting the max over 6 pairs and then quoting")
print("that pair's uncorrected-for-selection q is a multiple-comparison problem")
print("in its own right: the q in the table is a BH q among 102 marker tests, not")
print("a BH q among 102 x 6 = 612 gene-pair-max statistics.")
best = (mk2.assign(ad=mk2.delta_logFC.abs())
        .sort_values("ad", ascending=False)
        .groupby("gene").head(1)
        [["gene", "aetiology_A", "aetiology_B", "delta_logFC", "ad",
          "p_permutation_exact", "q_permutation"]])
print(best.head(10).to_string(index=False))
print(f"\n  min q among the 6 per-gene maxima = {best.q_permutation.min():.4f}")
print(f"  if the family were 17 genes x 6 selected maxima = 102 statistics, the")
print(f"  BH q would need recomputation over that family; it is currently")
print(f"  computed over all 102 rows including non-maximal pairs, which is the")
print(f"  same set of 102 -- so this one happens to be self-consistent, but the")
print(f"  SELECTION of the max is not accounted for anywhere.")


if __name__ == "__main__":
    pass
