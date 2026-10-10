#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Regenerate Tables 2, 3, 4 of the manuscript directly from the deposited
artefacts, so the table cells are guaranteed to trace to the CSVs (the 98-checker
only validates the CSVs, not the table cells, so hand-transcribed tables can
drift). Output is printed and written to results/tables/_generated_tables.md.
"""
import os, pandas as pd
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
T=os.path.join(ROOT,"results","tables")
def r2(x): return f"{x:+.2f}"
def r3(x): return f"{x:+.3f}"

mk=pd.read_csv(os.path.join(T,"etiology_marker_de.csv")).set_index(["gene","contrast"])
mc=pd.read_csv(os.path.join(T,"etiology_marker_cross.csv"))
mx=pd.read_csv(os.path.join(T,"etiology_module_logfc_cross.csv"))
vp=pd.read_csv(os.path.join(T,"virtual_perturbation_scores.csv"))

# ---- Table 2: canonical markers ----
genes=["Cd68","Aif1","Gal","Ccl2","Il1b","Calca","Scn10a","Atf3","Gfap"]
studies=["Bone_D7","Bone_D14","PTX","OXA"]
study_lbl={"Bone_D7":"Bone (day 7)","Bone_D14":"Bone (day 14)","PTX":"PTX","OXA":"OXA"}
def marker_logfc(g,c):
    # per-study simple-mean logFC extracted from the cross-aetiology CSV, the
    # same artefact that supplies Delta and q, so every Table-2 cell traces to
    # one source instead of mixing de.csv (muscle-adjusted) with cross.csv.
    rowA=mc[(mc.gene==g)&(mc.aetiology_A==c)]
    if len(rowA): return float(rowA.iloc[0]["logFC_A"])
    rowB=mc[(mc.gene==g)&(mc.aetiology_B==c)]
    return float(rowB.iloc[0]["logFC_B"])
def largest_delta(g):
    # restrict to bone-vs-CIPN contrasts (Bone_D7/D14 vs PTX/OXA): the only
    # contrasts the manuscript reports in this column. Report only when it
    # survives q<0.05; otherwise mark ns (matches the manuscript convention).
    sub=mc[(mc.gene==g)&(mc.aetiology_A.isin(["Bone_D7","Bone_D14"]))&(mc.aetiology_B.isin(["PTX","OXA"]))]
    i=sub["delta_logFC"].abs().idxmax()
    r=sub.loc[i]
    if r["q_permutation"]<0.05:
        return r["delta_logFC"], r["aetiology_A"], r["aetiology_B"], r["q_permutation"]
    return None
L=[]
L.append("| Gene | Bone (day 7) | Bone (day 14) | PTX | OXA | Largest Δ | q |")
L.append("|---|---|---|---|---|---|---|")
for g in genes:
    cells=[g]+[r2(marker_logfc(g,c)) for c in studies]
    res=largest_delta(g)
    if res is None:
        cells.append("—"); cells.append("ns")
    else:
        d,a,b,q=res
        tag=f"{a.replace('Bone_','')} vs {b}" if a.startswith('Bone') else f"{a} vs {b}"
        cells.append(f"{r2(d)} ({tag})"); cells.append(f"{q:.3f}")
    L.append("| "+" | ".join(cells)+" |")
table2="\n".join(L)

# ---- Table 3: module mean logFC ----
mods=["Glial / immune","Cytokine / chemokine","Toll-like receptors","Neurotrophin",
      "Ca2+ channel a2d","Neuropeptides","CTRL pan-neuronal","Glutamate / NMDA",
      "Voltage-gated Na+","GABAergic","CTRL erythrocyte"]
def mod_logfc(m,study):
    # value of `study` from cross csv (as A or B)
    rA=mx[(mx.module==m)&(mx.aetiology_A==study)]
    rB=mx[(mx.module==m)&(mx.aetiology_B==study)]
    vals=[]
    if len(rA): vals.append(float(rA.iloc[0]["logFC_A"]))
    if len(rB): vals.append(float(rB.iloc[0]["logFC_B"]))
    assert vals and max(vals)-min(vals)<1e-6, (m,study,vals)
    return vals[0]
def mod_delta(m):
    r=mx[(mx.module==m)&(mx.aetiology_A=="Bone_D7")&(mx.aetiology_B=="PTX")].iloc[0]
    return r["delta_logFC"], r["q_permutation"]
# kind lookup
kind_map=mx.drop_duplicates("module").set_index("module")["kind"].to_dict()
L=[]
L.append("| Module | Kind | logFC bone D7 | logFC bone D14 | logFC PTX | logFC OXA | Δ (D7 vs PTX) | q |")
L.append("|---|---|---|---|---|---|---|---|")
for m in mods:
    k=kind_map.get(m,"target")
    cells=[m,k]+[r2(mod_logfc(m,s)) for s in studies]
    d,q=mod_delta(m); cells.append(r2(d)); cells.append(f"{q:.4f}")
    L.append("| "+" | ".join(cells)+" |")
table3="\n".join(L)

# ---- Table 4: perturbation D7 (knockout) ----
v7=vp[(vp.contrast=="D7")&(vp.perturbation=="knockout")].set_index("module")
order=["CTRL interferon","Glial / immune","CTRL pan-neuronal","K+ / HCN channels",
       "Voltage-gated Na+","CTRL myelin / Schwann","CTRL erythrocyte","Neuropeptides"]
L=[]
L.append("| Module | Kind | Reversal (knockout) | z vs degree-matched null | % of oracle | p | q |")
L.append("|---|---|---|---|---|---|---|")
for m in order:
    r=v7.loc[m]
    cells=[m,r["kind"],r2(float(r["reversal_network"])),f"{float(r['z_vs_perm']):+.2f}",
           f"{float(r['pct_of_oracle']):.1f}",f"{float(r['p_empirical']):.4f}",
           f"{float(r['q_empirical']):.3f}"]
    L.append("| "+" | ".join(cells)+" |")
table4="\n".join(L)

out=f"## Table 2\n\n{table2}\n\n## Table 3\n\n{table3}\n\n## Table 4\n\n{table4}\n"
open(os.path.join(T,"_generated_tables.md"),"w").write(out)
print(out)
PY = None
