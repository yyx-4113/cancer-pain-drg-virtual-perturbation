#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""A3 stage 2: reproduce published compartment DE, then stress it.

Tests:
  T1  exact reproduction of sc_GSE292522_celltype_DE.csv from cached pseudobulk
  T2  BH scope: within-compartment vs global (count q<0.05 both ways)
  T3  global shift test: median logFC per compartment == denominator artefact?
  T4  leave-one-animal-out for every headline gene
  T5  circularity: FOCUS genes that are ALSO compartment marker genes
  T6  duplicated symbols among ISG/FOCUS
  T7  near-zero-variance genes (Gfap in SatelliteGlia)
  T8  depth confound: per-compartment total counts and per-cell UMI by group
"""
from __future__ import annotations
import os, sys, json
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from _limma import fit_limma, bh
from scipy.stats import mannwhitneyu, ttest_ind

TAB = os.path.join(ROOT, "results", "tables")
OUT = os.path.join(ROOT, "results", "review_a3")
os.makedirs(OUT, exist_ok=True)

SAMPLES = {
    "Sample1_DRG_Con_rep1": "Con", "Sample2_DRG_Con_rep2": "Con",
    "Sample3_DRG_Con_rep3": "Con", "Sample4_DRG_Con_rep4": "Con",
    "Sample5_DRG_B16F10_rep1": "B16F10", "Sample6_DRG_B16F10_rep2": "B16F10",
    "Sample7_DRG_B16F10_rep3": "B16F10", "Sample8_DRG_B16F10_rep4": "B16F10",
}
BROAD = {"Neuron","SatelliteGlia","Schwann","Myeloid","TCell",
         "Endothelial","Fibroblast","Erythrocyte"}
NEURON_SUB = {"NF","NP","PE","TH"}
FOCUS = ["Cd68","Aif1","Itgam","Lyz2","Cx3cr1","Gfap","Sox10","Mpz",
         "Snap25","Scn9a","Scn10a","Calca","Tac1","Gal","Il1b","Tnf","Ccl2",
         "Ifit1","Ifit2","Ifit3","Isg15","Oas1a","Oas1b","Bst2","Mx1","Stat1","Irf7"]
BROAD_MARKERS = {
    "Neuron": ["Snap25","Syp","Tubb3","Nefl","Nefh","Th","Scn10a","Scn9a",
               "Trpv1","Calca","Tac1","Npy","Gal","PgP9.5","Prph"],
    "SatelliteGlia": ["Gfap","Sox10","Slc1a3","Aqp4","Kcnj10","Glul","Blbp","Fgfr3"],
    "Schwann": ["Mpz","Pmp22","Oct6","Mag","Plp1","Sox10","Dhh","Ncmap"],
    "Myeloid": ["Cd68","Aif1","Itgam","Lyz2","Cx3cr1","Fcgr1","Cd14","Itgax",
                "Ctss","Ly6c2","Ly6g","Spi1"],
    "TCell": ["Cd3e","Cd3d","Cd4","Cd8a","Trbc1","Trac"],
    "Endothelial": ["Pecam1","Cdh5","Cldn5","Vwf","Emcn"],
    "Fibroblast": ["Col1a1","Col1a2","Pdgfra","Dcn","Lum","Col3a1"],
    "Erythrocyte": ["Hbb-bs","Hbb-bt","Hba-a1","Hba-a2","Hbb-y"],
}
NEURON_SUB_MARK = {
    "NF": ["Nefh","Nefl","Prph","Map2","Tubb3"],
    "NP": ["P2rx3","Ret","Sst","Npy","Runx1"],
    "PE": ["Calca","Tac1","Trpv1","Scn10a","Scn9a"],
    "TH": ["Th","Dbh","Ddc"],
}

z = np.load(os.path.join(OUT, "compartment_pseudobulk.npz"), allow_pickle=True)
sym = list(z["canon_sym"])
COMPARTMENTS = [c for c in list(BROAD) + [f"Neuron.{s}" for s in NEURON_SUB]
                if f"{c}||Sample1_DRG_Con_rep1" in z or any(
                    k.startswith(f"{c}||") for k in z.files)]
COMPARTMENTS = sorted({k.split("||")[0] for k in z.files if "||" in k})
print("[T0] compartments recovered:", COMPARTMENTS)

celln = json.load(open(os.path.join(OUT, "compartment_cellcounts.json")))
report = {}

# ---------------------------------------------------------------- T1 reproduce
def build(ct, samples):
    M = pd.DataFrame({p: z[f"{ct}||{p}"] for p in samples}, index=sym)
    M = M.loc[M.sum(axis=1) > 0]
    cp10k = M.div(M.sum(axis=0), axis=1) * 1e4 + 1.0
    return M, np.log2(cp10k)

def fit_of(L, samples):
    g = [SAMPLES[p] for p in samples]
    design = np.column_stack([np.ones(len(g)),
                              np.array([1.0 if x == "B16F10" else 0.0 for x in g])])
    f = fit_limma(L, design, coef_idx=1, trend=True)
    f["q"] = bh(f["P"].values)
    return f, g

pub = pd.read_csv(os.path.join(TAB, "sc_GSE292522_celltype_DE.csv"))
repro, de_all = [], {}
print("\n[T1] reproduce published compartment DE")
for ct in COMPARTMENTS:
    pres = [p for p in SAMPLES if f"{ct}||{p}" in z.files]
    if len(pres) < 6:
        print(f"  {ct}: only {len(pres)} samples -> dropped by script")
        continue
    M, L = build(ct, pres)
    f, g = fit_of(L, pres)
    de_all[ct] = (f, L, pres, M)
    n_con, n_b16 = g.count("Con"), g.count("B16F10")
    nq = int((f.q < 0.05).sum())
    print(f"  {ct}: n={len(pres)} ({n_con}Con/{n_b16}B16) genes={len(f)} "
          f"DE(q<.05)={nq}  median_logFC={f.logFC.median():+.4f}  "
          f"mean_logFC={f.logFC.mean():+.4f}")
    for gene in FOCUS:
        if gene in f.index:
            r = f.loc[gene]
            repro.append(dict(cell_type=ct, gene=gene, logFC=round(float(r.logFC), 3),
                              t=round(float(r.t), 2), q=round(float(r.q), 4),
                              n_con=n_con, n_b16=n_b16, nDE_q05=nq))
rep = pd.DataFrame(repro).merge(pub, on=["cell_type", "gene"], how="outer",
                                suffixes=("_mine", "_pub"))
rep["d_logFC"] = (rep.logFC_mine - rep.logFC_pub).abs()
rep["d_q"] = (rep.q_mine - rep.q_pub).abs()
bad = rep[(rep.d_logFC > 1e-6) | (rep.d_q > 1e-6)]
print(f"  rows mine={len(rep)} published={len(pub)}  mismatches={len(bad)}")
if len(bad):
    print(bad[["cell_type","gene","logFC_mine","logFC_pub","q_mine","q_pub"]].to_string())
report["T1_reproduced"] = bool(len(bad) == 0)
report["T1_rows"] = int(len(rep))

# ---------------------------------------------------------------- T2 BH scope
print("\n[T2] BH scope: within-compartment vs across-all-compartments")
allp = []
for ct, (f, L, pres, M) in de_all.items():
    d = f[["P"]].copy(); d["ct"] = ct
    allp.append(d.reset_index().rename(columns={"index": "gene"}))
allp = pd.concat(allp, ignore_index=True)
allp["q_global"] = bh(allp["P"].values)
n_within = sum(int((f.q < 0.05).sum()) for f, _, _, _ in de_all.values())
n_global = int((allp.q_global < 0.05).sum())
print(f"  total gene tests across compartments = {len(allp)}")
print(f"  DE at q<0.05 with BH WITHIN compartment   = {n_within}")
print(f"  DE at q<0.05 with BH ACROSS compartments  = {n_global}")
rep["q_global"] = [float(allp[(allp.ct == r.cell_type) & (allp.gene == r.gene)].q_global.iloc[0])
                   if (allp.ct == r.cell_type).any() and
                      ((allp.ct == r.cell_type) & (allp.gene == r.gene)).any() else np.nan
                   for r in rep.itertuples()]
rep["sig_within"] = rep.q_pub < 0.05
rep["sig_global"] = rep.q_global < 0.05
print("  FOCUS-gene calls that change family under global BH:")
ch = rep[rep.sig_within != rep.sig_global]
print(ch[["cell_type","gene","q_pub","q_global","sig_within","sig_global"]].to_string()
      if len(ch) else "   (none change)")
report["T2_n_within"] = n_within; report["T2_n_global"] = n_global
report["T2_focus_flips"] = int(len(ch))
rep.to_csv(os.path.join(OUT, "T2_reproduction_and_BH.csv"), index=False)

# ---------------------------------------------------------------- T3 global shift
print("\n[T3] global (denominator) shift per compartment")
print("   if median logFC != 0 the compartment CPM denominator differs by group")
rows = []
for ct, (f, L, pres, M) in de_all.items():
    tot = M.sum(axis=0)
    con = [tot[p] for p in pres if SAMPLES[p] == "Con"]
    b16 = [tot[p] for p in pres if SAMPLES[p] == "B16F10"]
    ncell = [celln[ct][p] for p in pres]
    ncell_c = [celln[ct][p] for p in pres if SAMPLES[p] == "Con"]
    ncell_b = [celln[ct][p] for p in pres if SAMPLES[p] == "B16F10"]
    lfc_ratio = np.log2(np.mean(b16) / np.mean(con))
    rows.append(dict(cell_type=ct, n_con=len(con), n_b16=len(b16),
                     median_logFC=float(f.logFC.median()),
                     frac_pos=float((f.logFC > 0).mean()),
                     log2_total_b16_over_con=float(lfc_ratio),
                     tot_con_mean=float(np.mean(con)), tot_b16_mean=float(np.mean(b16)),
                     cells_con=float(np.mean(ncell_c)),
                     cells_b16=float(np.mean(ncell_b)),
                     counts_per_cell_con=float(np.mean(con) / max(np.mean(ncell_c), 1)),
                     counts_per_cell_b16=float(np.mean(b16) / max(np.mean(ncell_b), 1))))
t3 = pd.DataFrame(rows)
print(t3.to_string(index=False))
t3.to_csv(os.path.join(OUT, "T3_global_shift.csv"), index=False)

# ---------------------------------------------------------------- T4 LOO
print("\n[T4] leave-one-animal-out on headline genes")
HEAD = [("Neuron", "Gal"), ("Neuron.NF", "Gal"), ("Neuron", "Scn10a"),
        ("Neuron", "Ifit2"), ("Neuron", "Cx3cr1"), ("SatelliteGlia", "Isg15"),
        ("SatelliteGlia", "Ifit1"), ("Schwann", "Gal"), ("Schwann", "Ccl2"),
        ("Endothelial", "Lyz2"), ("Endothelial", "Gal"), ("Neuron.NF", "Il1b"),
        ("Neuron.NF", "Scn10a"), ("Fibroblast", "Gal")]
loo = []
for ct, gene in HEAD:
    if ct not in de_all:
        continue
    f0, L0, pres, M0 = de_all[ct]
    full = (float(f0.loc[gene].logFC), float(f0.loc[gene].q)) if gene in f0.index else (np.nan, np.nan)
    for drop in pres:
        sub = [p for p in pres if p != drop]
        if len(sub) < 6:
            continue
        M, L = build(ct, sub)
        f, g = fit_of(L, sub)
        if gene in f.index:
            loo.append(dict(cell_type=ct, gene=gene, dropped=drop,
                            logFC=float(f.loc[gene].logFC), t=float(f.loc[gene].t),
                            q_within=float(f.loc[gene].q),
                            logFC_full=full[0], q_full=full[1]))
loo = pd.DataFrame(loo)
loo.to_csv(os.path.join(OUT, "T4_leave_one_out.csv"), index=False)
for (ct, gene), grp in loo.groupby(["cell_type", "gene"], sort=False):
    d = grp.logFC - grp.logFC_full
    flip = int(((d > 0) != (grp.logFC_full > 0)).sum())
    qx = int((grp.q_within < 0.05).sum())
    print(f"  {ct:14s} {gene:8s} full logFC={grp.logFC_full:+.3f} q={grp.q_full:.4f} | "
          f"LOO logFC [{grp.logFC.min():+.3f},{grp.logFC.max():+.3f}] "
          f"sign-flips={flip}/8  still q<.05 in {qx}/8 drops")
report["T4_loops"] = int(len(loo))

# ---------------------------------------------------------------- T5 circularity
print("\n[T5] circularity: FOCUS genes that are also the compartment's own markers")
circ = []
for ct, (f, L, pres, M) in de_all.items():
    mk = set(BROAD_MARKERS.get(ct, []))
    if ct.startswith("Neuron."):
        mk = set(NEURON_SUB_MARK[ct.split(".")[1]]) | set(BROAD_MARKERS["Neuron"])
    ov = sorted(mk & set(FOCUS))
    if ov:
        sub = f.loc[[g for g in ov if g in f.index]]
        sig = sub[(sub.q < 0.05)]
        circ.append(dict(cell_type=ct, overlap_genes=",".join(ov),
                         n_overlap=len(ov), n_sig=sig.shape[0],
                         sig_genes=",".join(f"{g}(q={q:.3f})" for g, q in sig.q.items())))
        print(f"  {ct:14s} overlap={ov}")
        print(f"      of which q<0.05: {list(sig.index)}")
t5 = pd.DataFrame(circ)
t5.to_csv(os.path.join(OUT, "T5_circularity.csv"), index=False)

# ---------------------------------------------------------------- T6 duplicates
print("\n[T6] duplicated gene symbols")
dup = pd.Series(sym)[pd.Series(sym).duplicated()].unique().tolist()
print(f"  n duplicated symbols in features.tsv = {len(dup)}")
print(f"  ISG/FOCUS genes that are duplicated   = "
      f"{sorted(set(dup) & (set(FOCUS) | set(['Ifit1','Ifit3','Isg15','Oas1a','Bst2','Mx1','Stat1','Irf7'])))}")
report["T6_n_dup"] = len(dup)

# ---------------------------------------------------------------- T7 zero-variance
print("\n[T7] near-zero-variance genes retained in the DE table")
for ct, (f, L, pres, M) in de_all.items():
    v = L.var(axis=1)
    lo = v[(v < 1e-8)]
    if len(lo):
        print(f"  {ct}: {len(lo)} genes with variance < 1e-8 (all-samples identical)")
gf = []
for ct, (f, L, pres, M) in de_all.items():
    for gene in ("Gfap", "Sox10", "Aif1", "Cd68"):
        if gene in L.index:
            gf.append(dict(cell_type=ct, gene=gene,
                           values=";".join(f"{v:.3f}" for v in L.loc[gene].values),
                           var=float(L.loc[gene].var())))
print(pd.DataFrame(gf).to_string(index=False))
pd.DataFrame(gf).to_csv(os.path.join(OUT, "T7_marker_values.csv"), index=False)

# ---------------------------------------------------------------- T8 depth
print("\n[T8] per-compartment sequencing depth by group (Mann-Whitney on totals)")
for ct, (f, L, pres, M) in de_all.items():
    tot = M.sum(axis=0)
    c = [tot[p] for p in pres if SAMPLES[p] == "Con"]
    b = [tot[p] for p in pres if SAMPLES[p] == "B16F10"]
    u = mannwhitneyu(b, c, alternative="two-sided")
    cc = [celln[ct][p] for p in pres if SAMPLES[p] == "Con"]
    cb = [celln[ct][p] for p in pres if SAMPLES[p] == "B16F10"]
    print(f"  {ct:14s} total counts Con={np.mean(c):.3g} B16={np.mean(b):.3g} "
          f"MWU_p={u.pvalue:.4f} | cells Con={np.mean(cc):.0f} B16={np.mean(cb):.0f}")
print("\n[A3] DONE stage2")
json.dump(report, open(os.path.join(OUT, "A3_stage2_report.json"), "w"), indent=2)