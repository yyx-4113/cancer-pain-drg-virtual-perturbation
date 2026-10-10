#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""A3 independent recomputation of GSE292522 per-compartment quantities.

Replicates 32_percelltype_sc.py exactly (same marker sets, same score, same
argmax rule) but ALSO records the diagnostics the manuscript never reports:
  - per-compartment per-animal CELL COUNTS
  - top1-top2 annotation margin (assignment confidence)
  - cells with ZERO marker counts (which argmax silently sends to "Neuron")
  - n_genes / UMI / mito per cell, and per-compartment
  - doublet proxy: cells scoring high on >=2 mutually exclusive sets
  - per-compartment per-animal summed counts cached to .npz for LOO analysis
Also records which samples are PRESENT per compartment (>=5 cells rule).
"""
from __future__ import annotations
import os, sys, json, gzip, hashlib
import numpy as np
import pandas as pd
from scipy.io import mmread

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SC = os.path.join(ROOT, "data", "processed", "GSE292522")
OUT = os.path.join(ROOT, "results", "review_a3")
os.makedirs(OUT, exist_ok=True)

SAMPLES = {
    "Sample1_DRG_Con_rep1":    ("Con",    "GSM8860632_Sample1_DRG_Con_rep1"),
    "Sample2_DRG_Con_rep2":    ("Con",    "GSM8860633_Sample2_DRG_Con_rep2"),
    "Sample3_DRG_Con_rep3":    ("Con",    "GSM8860634_Sample3_DRG_Con_rep3"),
    "Sample4_DRG_Con_rep4":    ("Con",    "GSM8860635_Sample4_DRG_Con_rep4"),
    "Sample5_DRG_B16F10_rep1": ("B16F10", "GSM8860636_Sample5_DRG_B16F10_rep1"),
    "Sample6_DRG_B16F10_rep2": ("B16F10", "GSM8860637_Sample6_DRG_B16F10_rep2"),
    "Sample7_DRG_B16F10_rep3": ("B16F10", "GSM8860638_Sample7_DRG_B16F10_rep3"),
    "Sample8_DRG_B16F10_rep4": ("B16F10", "GSM8860639_Sample8_DRG_B16F10_rep4"),
}

BROAD = {
    "Neuron":        ["Snap25","Syp","Tubb3","Nefl","Nefh","Th","Scn10a","Scn9a",
                     "Trpv1","Calca","Tac1","Npy","Gal","PgP9.5","Prph"],
    "SatelliteGlia": ["Gfap","Sox10","Slc1a3","Aqp4","Kcnj10","Glul","Blbp","Fgfr3"],
    "Schwann":       ["Mpz","Pmp22","Oct6","Mag","Plp1","Sox10","Dhh","Ncmap"],
    "Myeloid":       ["Cd68","Aif1","Itgam","Lyz2","Cx3cr1","Fcgr1","Cd14","Itgax",
                     "Ctss","Ly6c2","Ly6g","Spi1"],
    "TCell":         ["Cd3e","Cd3d","Cd4","Cd8a","Trbc1","Trac"],
    "Endothelial":   ["Pecam1","Cdh5","Cldn5","Vwf","Emcn"],
    "Fibroblast":    ["Col1a1","Col1a2","Pdgfra","Dcn","Lum","Col3a1"],
    "Erythrocyte":   ["Hbb-bs","Hbb-bt","Hba-a1","Hba-a2","Hbb-y"],
}
NEURON_SUB = {
    "NF": ["Nefh","Nefl","Prph","Map2","Tubb3"],
    "NP": ["P2rx3","Ret","Sst","Npy","Runx1"],
    "PE": ["Calca","Tac1","Trpv1","Scn10a","Scn9a"],
    "TH": ["Th","Dbh","Ddc"],
}
ISG = ["Ifit1","Ifit2","Ifit3","Isg15","Isg20","Mx1","Mx2","Oas1a","Oas1b","Oas2",
       "Oas3","Cxcl10","Ccl5","Stat1","Stat2","Irf7","Irf9","Ifih1","Ddx58","Rsad2",
       "Trim21","Bst2","Gbp2","Gbp3","Gbp4","Gbp5","Ifi44","Ifi47","Parp9","Usp18",
       "Slfn1","Slfn2","Slfn4","Slfn5","Sp100","Ifi35","Ifi204","Samhd1"]


def gz_lines(path):
    with gzip.open(path, "rt") as fh:
        return [ln.rstrip("\n").split("\t") for ln in fh]


def marker_score(mat, idx_map, geneset):
    gi = [idx_map[g] for g in geneset if g in idx_map]
    if not gi:
        return None
    sub = mat[gi, :].toarray().astype(np.float32)
    return np.log1p(sub).mean(axis=0)


def main():
    # ---- feature-file integrity: are all 8 identical / same order? ----
    feat_hashes = {}
    for prefix, (grp, gsm) in SAMPLES.items():
        p = os.path.join(SC, gsm + "_features.tsv.gz")
        feat_hashes[prefix] = hashlib.md5(open(p, "rb").read()).hexdigest()
    print("[A3] feature-file md5:", json.dumps(feat_hashes, indent=2), flush=True)
    uniq = len(set(feat_hashes.values()))
    print(f"[A3] distinct feature files: {uniq}", flush=True)

    COMPARTMENTS = list(BROAD) + [f"Neuron.{s}" for s in NEURON_SUB]
    agg = {c: {} for c in COMPARTMENTS}
    celln = {c: {} for c in COMPARTMENTS}      # per-compartment CELL COUNTS
    qc_rows, margin_rows, zero_rows, doublet_rows = [], [], [], []
    canon_sym = None
    dup_sym = None

    for prefix, (grp, gsm) in SAMPLES.items():
        feat = gz_lines(os.path.join(SC, gsm + "_features.tsv.gz"))
        sym = [r[1] if len(r) > 1 else r[0] for r in feat]
        mat = mmread(os.path.join(SC, gsm + "_matrix.mtx.gz")).tocsc()
        if canon_sym is None:
            canon_sym = sym
            dup_sym = int(pd.Series(sym).duplicated().sum())
        ng, nc = mat.shape
        assert len(sym) == ng, f"symbol/row mismatch {len(sym)} vs {ng}"
        assert sym == canon_sym, f"feature order differs in {prefix}"

        idx_map = {g: i for i, g in enumerate(sym)}
        S = np.vstack([marker_score(mat, idx_map, gs) for gs in BROAD.values()])
        order = np.argsort(-S, axis=0)
        top1 = order[0, :]; top2 = order[1, :]
        assign_broad = np.array(list(BROAD))[top1]
        s1 = S[top1, np.arange(nc)]; s2 = S[top2, np.arange(nc)]
        marg = s1 - s2

        # cells with ZERO counts across the entire marker panel
        zero_marker = (S.max(axis=0) == 0)
        # per-cell QC
        umi = np.asarray(mat.sum(axis=0)).ravel()
        det = np.asarray((mat > 0).sum(axis=0)).ravel()
        mito = np.asarray(mat[[idx_map[g] for g in ("mt-Co1","mt-Co2","mt-Nd1",
                "mt-Nd2","mt-Cytb") if g in idx_map], :].sum(axis=0)).ravel() if any(
                g in idx_map for g in ("mt-Co1","mt-Co2")) else np.zeros(nc)

        for ct in BROAD:
            m = assign_broad == ct
            celln[ct][prefix] = int(m.sum())
            if m.sum() >= 5:
                agg[ct][prefix] = np.asarray(mat[:, m].sum(axis=1)).ravel()

        neu = assign_broad == "Neuron"
        assign_sub = np.array(["none"] * nc, dtype=object)
        if neu.any():
            Ssub = np.vstack([marker_score(mat, idx_map, gs) for gs in NEURON_SUB.values()])
            assign_sub[neu] = np.array(list(NEURON_SUB))[Ssub[:, neu].argmax(axis=0)]
        for s in NEURON_SUB:
            m = neu & (assign_sub == s)
            celln[f"Neuron.{s}"][prefix] = int(m.sum())
            if m.sum() >= 5:
                agg[f"Neuron.{s}"][prefix] = np.asarray(mat[:, m].sum(axis=1)).ravel()

        # ---- diagnostics ----
        for ct in BROAD:
            m = assign_broad == ct
            if m.sum() == 0:
                qc_rows.append(dict(sample=prefix, group=grp, cell_type=ct, n_cells=0))
                continue
            qc_rows.append(dict(sample=prefix, group=grp, cell_type=ct, n_cells=int(m.sum()),
                                median_umi=float(np.median(umi[m])),
                                median_genes=float(np.median(det[m])),
                                median_mito_pct=float(100*np.median(mito[m])/max(np.median(umi[m]),1)),
                                median_margin=float(np.median(marg[m]))))
        margin_rows.append(dict(sample=prefix, group=grp, n_cells=nc,
                                frac_margin_lt_0p1=float((marg < 0.1).mean()),
                                frac_margin_lt_0p3=float((marg < 0.3).mean()),
                                median_margin=float(np.median(marg)),
                                frac_zero_marker=int(zero_marker.mean()),
                                n_zero_marker=int(zero_marker.sum()),
                                median_umi=float(np.median(umi)),
                                median_genes=float(np.median(det))))
        # where do zero-marker cells go?
        zc = pd.Series(assign_broad[zero_marker]).value_counts().to_dict()
        zero_rows.append(dict(sample=prefix, group=grp, n_zero=int(zero_marker.sum()),
                              assigned_to=json.dumps({k: int(v) for k, v in zc.items()})))
        # doublet proxy: top2 score >= 0.5*top1 and both substantial
        dbl = (s2 >= 0.5 * s1) & (s1 > 0.3)
        doublet_rows.append(dict(sample=prefix, group=grp, n_cells=nc,
                                 n_ambiguous=int(dbl.sum()),
                                 frac_ambiguous=float(dbl.mean())))
        del mat
        print(f"  {prefix} [{grp}] ncells={nc} zeroMarker={int(zero_marker.sum())} "
              f"ambiguous={int(dbl.sum())} medMargin={np.median(marg):.3f}", flush=True)

    np.savez_compressed(os.path.join(OUT, "compartment_pseudobulk.npz"),
                        **{f"{c}||{p}": agg[c][p] for c in COMPARTMENTS for p in agg[c]},
                        canon_sym=np.array(canon_sym))
    json.dump({c: celln[c] for c in COMPARTMENTS},
              open(os.path.join(OUT, "compartment_cellcounts.json"), "w"), indent=2)
    pd.DataFrame(qc_rows).to_csv(os.path.join(OUT, "compartment_qc.csv"), index=False)
    pd.DataFrame(margin_rows).to_csv(os.path.join(OUT, "annotation_margin.csv"), index=False)
    pd.DataFrame(zero_rows).to_csv(os.path.join(OUT, "zero_marker_cells.csv"), index=False)
    pd.DataFrame(doublet_rows).to_csv(os.path.join(OUT, "ambiguous_cells.csv"), index=False)
    print(f"[A3] duplicate gene symbols in features.tsv: {dup_sym}", flush=True)
    print("[A3] DONE", flush=True)


if __name__ == "__main__":
    main()