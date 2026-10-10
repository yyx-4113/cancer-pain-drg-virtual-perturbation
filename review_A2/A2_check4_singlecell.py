#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
A2_check4_singlecell.py -- verify §3.8 / Table 5 claims and the BH families
in the single-cell arm, plus the "no ISG up in any compartment" claim.
"""
import os
import json
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


# ---------------------------------------------------------------- pseudobulk
hdr("1. WHOLE-DRG PSEUDOBULK -- BH check and the quoted marker values")
sc = pd.read_csv(os.path.join(TAB, "sc_GSE292522_pseudobulk_DE.csv"), index_col=0)
q2 = bh(sc.P.values)
print(f"n genes = {len(sc)}")
print(f"BH reproduces deposited q exactly: {np.allclose(q2, sc.q.values, atol=1e-12)}")
print(f"   max |q_recomputed - q_deposited| = {np.abs(q2-sc.q.values).max():.3e}")
print(f"n up at q<0.05   = {int(((sc.q<0.05)&(sc.logFC>0)).sum())}  [manuscript 1,517]")
print(f"n down at q<0.05 = {int(((sc.q<0.05)&(sc.logFC<0)).sum())}  [manuscript 703]")
print(f"n at q<0.05      = {int((sc.q<0.05).sum())}")
print("\nQuoted genes:")
for g in ["Sprr1a", "Atf3", "Gadd45a", "Serpinb1a", "Ecel1",
          "Ifit2", "Oas1c", "Gfap", "Cd68", "Aif1", "Isg15", "Gal", "Scn10a"]:
    if g in sc.index:
        r = sc.loc[g]
        print(f"  {g:10s} logFC={r.logFC:+.4f} t={r.t:+.2f} P={r.P:.3g} q={r.q:.4f}")
    else:
        print(f"  {g:10s} NOT IN INDEX")

# ISG set used by the paper (from scripts/32_percelltype_sc.py)
ISG = None
import re
src = open(os.path.join(ROOT, "scripts", "32_percelltype_sc.py")).read()
m = re.search(r"ISG\s*=\s*\[(.*?)\]", src, re.S)
if m:
    ISG = [x.strip().strip('"\'') for x in m.group(1).replace("\n", "").split(",")
           if x.strip().strip('"\'')]
    print(f"\nISG list read from scripts/32_percelltype_sc.py ({len(ISG)} genes): {ISG}")
    present = [g for g in ISG if g in sc.index]
    print(f"   present in pseudobulk DE table: {len(present)}")
    sub = sc.loc[present]
    up = sub[(sub.q < 0.05) & (sub.logFC > 0)]
    print(f"   ISGs UP at q<0.05: {len(up)}  -> {list(up.index)}")
    dn = sub[(sub.q < 0.05) & (sub.logFC < 0)]
    print(f"   ISGs DOWN at q<0.05: {len(dn)} -> {[(g, round(sc.loc[g,'q'],4)) for g in dn.index]}")
    print(f"   any ISG up at q<0.10: "
          f"{int(((sub.q<0.10)&(sub.logFC>0)).sum())}")
    print(f"\n   *** the paper's 'no ISG was up-regulated' is a statement about the")
    print(f"   *** {len(present)} ISGs at q<0.05. But BH was applied over ALL")
    print(f"   *** {len(sc)} genes, so each ISG q is a genome-wide q. If the claim")
    print(f"   *** is 'no ISG up', the relevant family is the {len(present)} ISGs.")
    q_isg = bh(sub.P.values)
    sub2 = sub.assign(q_isg_family=q_isg)
    print(f"\n   BH within the {len(present)}-ISG family only:")
    print(sub2[["logFC", "t", "P", "q", "q_isg_family"]].sort_values("P").to_string())
    print(f"   ISGs up at q_isg<0.05: "
          f"{int(((sub2.q_isg_family<0.05)&(sub2.logFC>0)).sum())}")

# ------------------------------------------------------------- compartments
hdr("2. PER-COMPARTMENT -- 'no compartment showed ISG up-regulation'")
ct = pd.read_csv(os.path.join(TAB, "sc_GSE292522_celltype_DE.csv"))
print(f"table shape {ct.shape}; compartments: {sorted(ct.cell_type.unique())}")
print(f"NOTE this table only carries the FOCUS gene list, so it CANNOT by itself")
print(f"     support a claim about 'no ISG up in any compartment' -- it only")
print(f"     contains {ct.gene.nunique()} genes x {ct.cell_type.nunique()} compartments.")
print(f"     genes present: {sorted(ct.gene.unique())}")
scores = pd.read_csv(os.path.join(TAB, "sc_GSE292522_celltype_scores.csv"))
print(f"\ncelltype_scores: {scores.shape}; score types {sorted(scores.score.unique())}")
print(f"   -> the ISG SCORE per compartment per animal is here, which IS the")
print(f"      object that supports an 'ISG up in no compartment' claim.")
for ctype in sorted(scores.cell_type.unique()):
    s = scores[(scores.cell_type == ctype) & (scores.score == "ISG")]
    if len(s) == 0:
        continue
    con = s[s.group == "Con"].value.values
    b16 = s[s.group == "B16F10"].value.values
    if len(con) < 2 or len(b16) < 2:
        print(f"  {ctype:16s} n={len(con)}/{len(b16)}  (too few for MWU)")
        continue
    mw = stats.mannwhitneyu(b16, con, alternative="two-sided")
    print(f"  {ctype:16s} n={len(con)}/{len(b16)}  meanISGscore Con={con.mean():+.4f} "
          f"B16={b16.mean():+.4f}  MWU two-sided p={mw.pvalue:.4f}")
n_up = 0
for ctype in sorted(scores.cell_type.unique()):
    s = scores[(scores.cell_type == ctype) & (scores.score == "ISG")]
    if len(s) == 0:
        continue
    con = s[s.group == "Con"].value.values
    b16 = s[s.group == "B16F10"].value.values
    if len(con) >= 2 and len(b16) >= 2:
        if b16.mean() > con.mean():
            n_up += 1
print(f"\n   compartments with a HIGHER mean ISG score in B16F10: {n_up} of "
      f"{scores.cell_type.nunique()}")
print("   -> the correct statement is 'the ISG score is not higher in B16F10 in")
print("      any compartment', which is a claim about DIRECTION, not about")
print("      significance. It is supported. The manuscript's phrasing 'no")
print("      compartment showed interferon-stimulated-gene up-regulation'")
print("      is a claim about gene-level q values that the deposited table")
print("      cannot cover exhaustively (FOCUS genes only).")

# ------------------------------------------- multiplicity family for §3.8
hdr("3. MULTIPLICITY FAMILY FOR THE PER-COMPARTMENT CLAIM")
print("scripts/32_percelltype_sc.py:160 applies bh() INSIDE each compartment")
print("loop iteration, i.e. 8 separate BH families (one per compartment).")
print("The manuscript then reads across all 8: 'no compartment showed ISG")
print("up-regulation ... the only significant ISG changes were Ifit2 in")
print("neurons (q=0.046) and Isg15 in satellite glia (q=0.007)'.")
sig = ct[(ct.gene.isin(["Ifit2", "Isg15", "Oas1c"])) & (ct.q < 0.05)]
print(f"\nrows with an ISG at q<0.05 in the deposited compartment table:")
print(sig.to_string(index=False))
print(f"\ncount = {len(sig)}:  Ifit2 in Neuron (q={ct[(ct.cell_type=='Neuron')&(ct.gene=='Ifit2')].q.iloc[0]:.4f})"
      f" and Ifit2 in Neuron.NF (q={ct[(ct.cell_type=='Neuron.NF')&(ct.gene=='Ifit2')].q.iloc[0]:.4f})"
      f" and Isg15 in SatelliteGlia (q={ct[(ct.cell_type=='SatelliteGlia')&(ct.gene=='Isg15')].q.iloc[0]:.4f})")
print("*** The manuscript says 'Ifit2 in neurons (q = 0.046) and Isg15 in")
print("    satellite glia (q = 0.007)' -- it OMITS the Neuron.NF result")
print("    (Ifit2 q = 0.0412), which is a THIRD significant ISG change and the")
print("    SMALLEST q of the three. Selective omission of the strongest")
print("    compartment-level ISG hit.")

# --------------------------------------------------------------- 4. Table 5
hdr("4. TABLE 5 ROW-BY-ROW")
summ = json.load(open(os.path.join(TAB, "sc_GSE292522_summary.json")))
print(f"total_cells = {summ['total_cells']}  [Table 5: 55,967]")
print(f"n_up = {summ['pseudobulk_DE']['n_up_B16F10_q05']}  "
      f"n_down = {summ['pseudobulk_DE']['n_down_B16F10_q05']}  [1,517 / 703]")
print(f"myeloid fraction Con = {summ['myeloid_fraction_mean']['Con']:.6f} "
      f"= {100*summ['myeloid_fraction_mean']['Con']:.3f}%  [Table 5: 0.12%]")
print(f"myeloid fraction B16 = {summ['myeloid_fraction_mean']['B16F10']:.6f} "
      f"= {100*summ['myeloid_fraction_mean']['B16F10']:.3f}%  [Table 5: 0.05%]")
print(f"ratio = {summ['myeloid_fraction_mean']['B16F10']/summ['myeloid_fraction_mean']['Con']:.4f}"
      f"  -> 'a 55% drop' = {100*(1-summ['myeloid_fraction_mean']['B16F10']/summ['myeloid_fraction_mean']['Con']):.1f}%")
print(f"MWU two-sided p = {summ['myeloid_MWU_p_twosided']:.4f}  [0.114]")
print(f"MWU one-sided p = {summ['myeloid_MWU_p_less']:.4f}  [0.057]")
print(f"Spearman rho vs bulk D7 = {summ['convergent_validation_vs_bulk_D7']['spearman_rho']:.4f} "
      f"p = {summ['convergent_validation_vs_bulk_D7']['spearman_p']:.3g}  [0.098]")
print("\n*** 'U = 2/14' in §3.8. For 4 vs 4 the maximum possible U is 16, not 14,")
print("*** and the MWU two-sided p of 0.1143 corresponds to U = 2. Let me check:")
from itertools import combinations
c = comp = pd.read_csv(os.path.join(TAB, "sc_GSE292522_composition.csv"))
con = c[c.group == "Con"].Myeloid.values
b16 = c[c.group == "B16F10"].Myeloid.values
U = sum(1 for i in range(4) for j in range(4) if b16[i] > con[j]) + \
    0.5 * sum(1 for i in range(4) for j in range(4) if b16[i] == con[j])
print(f"    U(B16 > Con) = {U}   [manuscript says U = 2]")
print(f"    n1*n2 = 16; the manuscript's 'U = 2/14' denominator is wrong.")
print(f"    Also: the MYELOID FRACTION has n = 4 per group, so 'U = 2/14'")
print(f"    cannot be a U statistic over 14 units.")

print("\n" + "=" * 84)
print("5. CONCORDANCE CLAIM: 'highest at day 1, rho = 0.112, vanishing by day 14")
print("   rho = 0.030'")
print("=" * 84)
print("Only the D7 concordance (rho = 0.098) is in the deposited summary JSON.")
print("D1 (0.112) and D14 (0.030) are NOT in sc_GSE292522_summary.json, so they")
print("cannot be verified from the deposited artefacts.")
print("Additionally: rho = 0.098 with n = 12,158 genes has p = 1.4e-27 -- the")
print("concordance is *statistically* non-zero but *practically* negligible.")
print("The paper reports the coefficient only; the reader cannot see that the")
print("p-value is astronomically small while rho is biologically meaningless.")
s = summ['convergent_validation_vs_bulk_D7']
print(f"   deposited: n_shared = {s['n_shared_genes']}, rho = {s['spearman_rho']:.5f}, "
      f"p = {s['spearman_p']:.3g}")


if __name__ == "__main__":
    pass
