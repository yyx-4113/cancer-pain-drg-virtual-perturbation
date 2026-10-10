"""Verify Gfap expression (UMI) in the satellite-glia bin per animal,
and the neuron-bin marker-positive fraction (purity) with an explicit
criterion. Reuses script 31's loader and marker sets."""
import os, sys, gzip, numpy as np, pandas as pd
from scipy.io import mmread

ROOT = "D:/2026.10/不同部位癌痛患者对药物敏感性生信分析/cancer_pain_virtual_perturbation"
SC = os.path.join(ROOT, "data", "processed", "GSE292522")
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from importlib import import_module
m31 = import_module("31_process_sc") if False else None

# replicate the SAMPLES + MARKERS from script 31
SAMPLES = {
    "Sample1_DRG_Con_rep1": ("Con", "GSM8860632_Sample1_DRG_Con_rep1"),
    "Sample2_DRG_Con_rep2": ("Con", "GSM8860633_Sample2_DRG_Con_rep2"),
    "Sample3_DRG_Con_rep3": ("Con", "GSM8860634_Sample3_DRG_Con_rep3"),
    "Sample4_DRG_Con_rep4": ("Con", "GSM8860635_Sample4_DRG_Con_rep4"),
    "Sample5_DRG_B16F10_rep1": ("B16F10", "GSM8860636_Sample5_DRG_B16F10_rep1"),
    "Sample6_DRG_B16F10_rep2": ("B16F10", "GSM8860637_Sample6_DRG_B16F10_rep2"),
    "Sample7_DRG_B16F10_rep3": ("B16F10", "GSM8860638_Sample7_DRG_B16F10_rep3"),
    "Sample8_DRG_B16F10_rep4": ("B16F10", "GSM8860639_Sample8_DRG_B16F10_rep4"),
}
MARKERS = {
    "Neuron": ["Snap25","Syp","Tubb3","Nefl","Nefh","Th","Scn10a","Scn9a",
               "Trpv1","Calca","Tac1","Npy","Gal","Pgp9.5","Prph"],
    "SatelliteGlia": ["Gfap","Sox10","Slc1a3","Aqp4","Kcnj10","Glul","Blbp","Fgfr3"],
    "Schwann": ["Mpz","Pmp22","Oct6","Mag","Plp1","Sox10","Dhh","Ncmap"],
    "Myeloid": ["Cd68","Aif1","Itgam","Lyz2","Cx3cr1","Fcgr1","Cd14","Itgax",
                "Ctss","Ly6c2","Ly6g","Spi1"],
    "TCell": ["Cd3e","Cd3d","Cd4","Cd8a","Trbc1","Trac"],
    "Endothelial": ["Pecam1","Cdh5","Cldn5","Vwf","Emcn"],
    "Fibroblast": ["Col1a1","Col1a2","Pdgfra","Dcn","Lum","Col3a1"],
    "Erythrocyte": ["Hbb-bs","Hbb-bt","Hba-a1","Hba-a2","Hbb-y"],
}
CT_ORDER = list(MARKERS)

def gz_lines(p):
    with gzip.open(p, "rt") as fh:
        return [ln.rstrip("\n").split("\t") for ln in fh]

rows = []
for prefix, (grp, gsm) in SAMPLES.items():
    feat = gz_lines(os.path.join(SC, gsm + "_features.tsv.gz"))
    sym = [r[1] if len(r) > 1 else r[0] for r in feat]
    mat = mmread(os.path.join(SC, gsm + "_matrix.mtx.gz")).tocsc()
    all_markers = sorted({g for s in MARKERS.values() for g in s})
    mset = [g for g in all_markers if g in set(sym)]
    idx = [sym.index(g) for g in mset]
    sub = np.log1p(mat[idx, :].toarray().astype(np.float32))
    scores = {}
    for ct, genes in MARKERS.items():
        gi = [mset.index(g) for g in genes if g in mset]
        scores[ct] = sub[gi, :].mean(axis=0) if gi else np.zeros(mat.shape[1])
    S = np.vstack([scores[ct] for ct in CT_ORDER])
    assign = np.argmax(S, axis=0)
    sg_i = CT_ORDER.index("SatelliteGlia")
    ne_i = CT_ORDER.index("Neuron")
    sg_mask = assign == sg_i
    ne_mask = assign == ne_i
    # Gfap raw UMI in the SG bin
    gfap_linked = False
    if "Gfap" in sym:
        g_row = sym.index("Gfap")
        gfap_umi = np.asarray(mat[g_row, :].todense()).ravel()
        gfap_sg = float(gfap_umi[sg_mask].sum()); gfap_sg_cells_pos = int((gfap_umi[sg_mask] > 0).sum())
        gfap_all = float(gfap_umi.sum())
    else:
        gfap_sg = gfap_sg_cells_pos = gfap_all = -1
    # purity of neuron bin under a naive independent-ish criterion:
    # Snap25>0 & Syp>0 & Tubb3>0 & Plp1==0 & Mpz==0
    def umi(g):
        return np.asarray(mat[sym.index(g), :].todense()).ravel() if g in sym else np.zeros(mat.shape[1])
    crit = (umi("Snap25") > 0) & (umi("Syp") > 0) & (umi("Tubb3") > 0) & (umi("Plp1") == 0) & (umi("Mpz") == 0)
    pur_ne = float(crit[ne_mask].mean()) if ne_mask.sum() else float("nan")
    rows.append(dict(sample=prefix, group=grp, ncells=int(mat.shape[1]),
                     sg_cells=int(sg_mask.sum()), sg_gfap_umi=gfap_sg,
                     sg_gfap_pos_cells=gfap_sg_cells_pos, gfap_total_umi=gfap_all,
                     ne_cells=int(ne_mask.sum()), purity_ne=round(pur_ne, 4)))
    print(rows[-1], flush=True)

df = pd.DataFrame(rows)
df.to_csv(os.path.join(ROOT, "results/tables/sc_GSE292522_gfap_purity_check.csv"), index=False)
print("\nSG bin total Gfap UMI across 8 animals:", df.sg_gfap_umi.sum())
print("animals with zero Gfap UMI in SG bin:", int((df.sg_gfap_umi == 0).sum()))
print("mean neuron purity (naive criterion): %.3f  range %.3f-%.3f"
      % (df.purity_ne.mean(), df.purity_ne.min(), df.purity_ne.max()))
