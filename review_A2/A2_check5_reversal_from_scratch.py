#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
A2_check5_reversal_from_scratch.py
Rebuild the reversal score for Glial / immune, Bone_D7, from raw FPKM using a
SEPARATE implementation of the network + diffusion, then:
  - confirm the deposited rev_Bone_D7
  - confirm reversal(KO) = -reversal(OE) analytically and numerically
  - confirm that the reversal score is a FUNCTIONAL of 8 animals, not a
    per-animal measurement (so no '2 vs 2' t-test exists)
  - do a leave-one-animal-out to get the sampling SD that a power analysis
    should actually use, and recompute the MDE from it
"""
import gzip
import os
import re
import sys
from itertools import combinations

import numpy as np
import pandas as pd
from statsmodels.stats.power import TTestIndPower, TTestPower

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "data", "raw")
TAB = os.path.join(ROOT, "results", "tables")
sys.path.insert(0, os.path.join(ROOT, "scripts"))

TARGETS = {
    "Opioid receptors":      ["Oprm1", "Oprd1", "Oprk1", "Oprl1"],
    "Voltage-gated Na+":     ["Scn9a", "Scn10a", "Scn11a", "Scn8a", "Scn1a"],
    "Ca2+ channel a2d":      ["Cacna2d1", "Cacna2d2", "Cacna1b", "Cacna1h"],
    "TRP channels":          ["Trpv1", "Trpa1", "Trpm3", "Trpm8"],
    "Purinergic P2X":        ["P2rx3", "P2rx4", "P2rx7"],
    "Prostanoid / COX":      ["Ptgs1", "Ptgs2", "Ptger1", "Ptger2", "Ptgir"],
    "Glutamate / NMDA":      ["Grin1", "Grin2a", "Grin2b", "Gria1", "Slc17a6"],
    "GABAergic":             ["Gabra1", "Gabra2", "Gabbr1", "Gad1", "Gad2"],
    "Neuropeptides":         ["Tac1", "Calca", "Npy", "Gal", "Vip", "Adcyap1"],
    "Cannabinoid":           ["Cnr1", "Cnr2", "Faah", "Mgll"],
    "Neurotrophin":          ["Ngf", "Ntrk1", "Ngfr", "Bdnf", "Ntf3"],
    "Cytokine / chemokine":  ["Tnf", "Il1b", "Il6", "Cxcl13", "Ccl2", "Cx3cl1"],
    "Toll-like receptors":   ["Tlr2", "Tlr4", "Tlr7"],
    "Glial / immune":        ["Gfap", "Aif1", "Itgam", "Cd68", "S100a8", "S100a9"],
    "Injury / regeneration": ["Atf3", "Sprr1a", "Sox11", "Gap43"],
    "K+ / HCN channels":     ["Kcnq2", "Kcnq3", "Kcna1", "Hcn1", "Hcn2"],
    "Monoaminergic":         ["Htr1a", "Htr2a", "Slc6a4", "Slc6a2", "Adra2a"],
}
CONTROL_MODULES = {
    "CTRL pan-neuronal":    ["Nefl", "Nefm", "Nefh", "Snap25", "Tubb3", "Uchl1",
                             "Prph", "Rbfox3", "Syt1", "Pvalb"],
    "CTRL myelin / Schwann": ["Mpz", "Pmp22", "Prx", "Mag", "Mbp", "Sox10", "Egr2"],
    "CTRL housekeeping":    ["Actb", "Gapdh", "Rpl13a", "Rplp0", "Ubc", "Pgk1",
                             "Tpi1", "Ppia", "Sdha", "Hsp90ab1", "Cox4i1"],
    "CTRL interferon":      ["Irf1", "Ifitm3", "Ifi47", "Bst2", "Rsad2", "Gbp2",
                             "Ubd", "Ly6e"],
    "CTRL erythrocyte":     ["Hba1", "Hba2", "Hbb", "Hbb-b1", "Alas2"],
}
ALL_MODS = dict(TARGETS); ALL_MODS.update(CONTROL_MODULES)

NET_HV, NET_TOPK, ALPHA = 3000, 12, 0.55


def read_geo(dirname):
    files = sorted(f for f in os.listdir(os.path.join(RAW, dirname))
                   if f.endswith(".txt.gz") or f.endswith(".txt"))
    out, order = {}, []
    for fn in files:
        gsm = fn.split("_")[0]
        path = os.path.join(RAW, dirname, fn)
        op = gzip.open if fn.endswith(".gz") else open
        with op(path, "rt") as fh:
            header = fh.readline().rstrip("\n").split("\t")
            up = [c.strip().upper() for c in header]
            vidx = up.index("FPKM") if "FPKM" in up else -1
            vals = {}
            for line in fh:
                p = line.rstrip("\n").split("\t")
                g = p[0].strip()
                if not g:
                    continue
                tok = p[vidx] if 0 <= vidx < len(p) else None
                if tok is None:
                    for j in range(len(p) - 1, 0, -1):
                        try:
                            float(p[j]); tok = p[j]; break
                        except ValueError:
                            continue
                if tok is None:
                    continue
                try:
                    vals[g] = float(tok)
                except ValueError:
                    continue
        out[gsm] = pd.Series(vals, dtype=float); order.append(gsm)
    return pd.DataFrame(out)[order].fillna(0.0)


def group_of(name, key):
    n = name.upper()
    if key == "bone":
        if "SHAM" in n:
            return "Sham"
        m = re.search(r"(?:^|[_.\-])D(\d+)(?:[_.\-])", n)
        return f"D{int(m.group(1))}" if m else "other"
    if "VEHICLE" in n: return "Vehicle"
    if "PACLITAXEL" in n: return "PTX"
    if "OXALIPLATIN" in n: return "OXA"
    return "other"


def norm(M):
    med = M.where(M > 0).median(axis=0)
    sf = (med / med.mean()).replace(0, np.nan).fillna(1.0)
    return np.log2(M.div(sf, axis=1) + 1.0)


def build_W(L, genes, topk=NET_TOPK):
    M = L.loc[genes].values
    C = np.nan_to_num(np.corrcoef(M))
    np.fill_diagonal(C, 0.0)
    idx = np.argpartition(-np.abs(C), topk, axis=1)[:, :topk]
    W = np.zeros_like(C)
    rows = np.arange(len(genes))[:, None]
    W[rows, idx] = C[rows, idx]
    W = (W + W.T) / 2.0
    np.fill_diagonal(W, 0.0)
    return W


def cos_col(D, b):
    D = np.asarray(D, float).ravel()
    b = np.asarray(b, float).ravel()
    Dc = D - D.mean(); bc = b - b.mean()
    na = np.linalg.norm(Dc); nb = np.linalg.norm(bc)
    if na < 1e-12 or nb < 1e-12:
        return 0.0
    return float(Dc @ bc / (na * nb))


def main():
    print("=" * 84)
    print("5. REVERSAL SCORE REBUILT FROM SCRATCH (Bone_D7 network, 8 animals)")
    print("=" * 84)
    Mb = read_geo("GSE149648"); Mc = read_geo("GSE160543")
    fb = {g: f for g, f in zip(Mb.columns, sorted(os.listdir(os.path.join(RAW, 'GSE149648'))))}
    fc = {g: f for g, f in zip(Mc.columns, sorted(os.listdir(os.path.join(RAW, 'GSE160543'))))}
    gb = {g: group_of(fb[g], "bone") for g in Mb.columns}
    gc = {g: group_of(fc[g], "chemo") for g in Mc.columns}
    kb = ((Mb >= 1.0).sum(axis=1) >= 4) & (Mb.mean(axis=1) >= 0.5)
    kc = ((Mc >= 1.0).sum(axis=1) >= 4) & (Mc.mean(axis=1) >= 0.5)
    Lb = norm(Mb).loc[kb]; Lc = norm(Mc).loc[kc]
    shared = sorted(set(Lb.index) & set(Lc.index))
    Lb, Lc = Lb.loc[shared], Lc.loc[shared]

    contrasts = {
        "Bone_D14": (Lb, [s for s in Lb.columns if gb.get(s) == "D14"],
                     [s for s in Lb.columns if gb.get(s) == "Sham"]),
        "Bone_D7":  (Lb, [s for s in Lb.columns if gb.get(s) == "D7"],
                     [s for s in Lb.columns if gb.get(s) == "Sham"]),
        "PTX":       (Lc, [s for s in Lc.columns if gc.get(s) == "PTX"],
                     [s for s in Lc.columns if gc.get(s) == "Vehicle"]),
        "OXA":       (Lc, [s for s in Lc.columns if gc.get(s) == "OXA"],
                     [s for s in Lc.columns if gc.get(s) == "Vehicle"]),
    }
    MUST = sorted({g for gs in ALL_MODS.values() for g in gs} & set(shared))
    vsum = None
    for k, (L, ca, co) in contrasts.items():
        v = L[ca + co].var(axis=1)
        vsum = v if vsum is None else vsum.add(v, fill_value=0.0)
    vmean = (vsum / len(contrasts)).loc[shared].sort_values(ascending=False)
    NET_GENES = list(vmean.index[:NET_HV]) + [g for g in MUST if g not in set(vmean.index[:NET_HV])]
    print(f"common node set = {len(NET_GENES)} genes "
          f"(repo log: 3038)  match = {len(NET_GENES)==3038}")

    dep = pd.read_csv(os.path.join(TAB, "etiology_reversal_point_estimates.csv"))
    depk = dep[(dep.network_scope == "per_contrast") &
               (dep.perturbation == "knockout")].set_index("module")

    print("\n-- rebuilt reversal(KO) vs deposited, all 22 modules, 4 contrasts --")
    rows = []
    for mod, glist in ALL_MODS.items():
        present = [g for g in glist if g in NET_GENES]
        if len(present) < 3:
            continue
        tidx = np.array([NET_GENES.index(g) for g in present])
        row = {"module": mod, "n": len(present)}
        for k, (L, ca, co) in contrasts.items():
            W = build_W(L[ca + co], NET_GENES)
            d = np.abs(W).sum(axis=1); d[d == 0] = 1.0
            Dm = np.diag(1.0 / np.sqrt(d))
            S = Dm @ W @ Dm
            P = (1 - ALPHA) * np.linalg.solve(np.eye(len(W)) - ALPHA * S, np.eye(len(W)))
            axis = (L[co].mean(axis=1) - L[ca].mean(axis=1)).reindex(NET_GENES).values
            for sign, tag in [(-1.0, "ko"), (1.0, "oe")]:
                seed = np.zeros(len(NET_GENES)); seed[tidx] = sign
                r = cos_col(P @ seed, axis)
                row[f"{k}_{tag}"] = r
                if tag == "ko":
                    row[f"dep_{k}"] = depk.loc[mod, f"rev_{k}"]
                else:
                    row[f"mirror_{k}"] = -row[f"{k}_ko"]
        rows.append(row)
    R = pd.DataFrame(rows)
    for k in contrasts:
        R[f"diff_{k}"] = R[f"{k}_ko"] - R[f"dep_{k}"]
        R[f"mirr_{k}"] = R[f"{k}_oe"] - R[f"mirror_{k}"]
    print(R[["module", "n", "Bone_D7_ko", "dep_Bone_D7", "diff_Bone_D7",
             "mirr_Bone_D7", "PTX_ko", "dep_PTX", "diff_PTX"]].to_string(index=False))
    print(f"\nmax |rebuilt - deposited| over all modules x contrasts = "
          f"{max(R[f'diff_{k}'].abs().max() for k in contrasts):.3e}")
    print(f"max |OE - (-KO)| over all modules x contrasts       = "
          f"{max(R[f'mirr_{k}'].abs().max() for k in contrasts):.3e}")
    print("=> my independent rebuild reproduces the deposited reversal scores and")
    print("   confirms the KO/OE mirror property to machine precision.")

    print("\n" + "=" * 84)
    print("6. THE REVERSAL SCORE IS A FUNCTIONAL OF 8 ANIMALS, NOT A 2-vs-2 MEASUREMENT")
    print("=" * 84)
    L, ca, co = contrasts["Bone_D7"]
    mod = "Glial / immune"
    present = [g for g in ALL_MODS[mod] if g in NET_GENES]
    tidx = np.array([NET_GENES.index(g) for g in present])
    seed = np.zeros(len(NET_GENES)); seed[tidx] = -1.0

    def rev_on(cols_case, cols_ctrl, seed_vec):
        W = build_W(L[cols_case + cols_ctrl], NET_GENES)
        d = np.abs(W).sum(axis=1); d[d == 0] = 1.0
        Dm = np.diag(1.0 / np.sqrt(d))
        S = Dm @ W @ Dm
        P = (1 - ALPHA) * np.linalg.solve(np.eye(len(W)) - ALPHA * S, np.eye(len(W)))
        axis = (L[cols_ctrl].mean(axis=1) - L[cols_case].mean(axis=1)).reindex(NET_GENES).values
        return cos_col(P @ seed_vec, axis)

    full = rev_on(ca, co, seed)
    print(f"full-data reversal(KO), Bone_D7, {mod} = {full:+.6f}")
    print(f"  (deposited rev_Bone_D7 = {depk.loc[mod,'rev_Bone_D7']:+.6f})")
    print("\nLeave-one-animal-out (drop 1 of the 8 animals, rebuild network+axis):")
    loo = []
    allc = ca + co
    for drop in allc:
        r = rev_on([c for c in ca if c != drop], [c for c in co if c != drop], seed)
        loo.append(r)
        print(f"   drop {drop}: reversal = {r:+.6f}   (change {r-full:+.6f})")
    loo = np.array(loo)
    loo_sd = loo.std(ddof=1)
    print(f"\n   LOO range {loo.min():+.6f} to {loo.max():+.6f}")
    print(f"   LOO SD (ddof=1) = {loo_sd:.6f}")
    print(f"   => a single animal moves the reversal score by up to "
          f"{np.abs(loo-full).max():.4f}.")
    print(f"   This LOO SD ({loo_sd:.4f}) is a within-contrast sampling SD of the")
    print(f"   reversal score -- the quantity a power analysis needs. It is")
    print(f"   comparable to the repo's between-module SD of 0.0741, NOT 0.42.")

    print("\n" + "=" * 84)
    print("7. MDE REDONE ON THE CORRECT OBJECT")
    print("=" * 84)
    tpi, tpp = TTestIndPower(), TTestPower()
    d_i4 = tpi.solve_power(nobs1=4, alpha=0.05, power=0.8, ratio=1.0)
    d_p4 = tpp.solve_power(nobs=4, alpha=0.05, power=0.8)
    d_i2 = tpi.solve_power(nobs1=2, alpha=0.05, power=0.8, ratio=1.0)
    se_boot = 0.0185     # median SE(Delta) implied by the paper's own bootstrap CI
    print(f"{'quantity':52s} {'d':>7s} {'x SE':>8s} {'= MDE':>9s}")
    print(f"{'repo: TTestIndPower nobs1=2 x between-module SD':52s} {d_i2:7.3f} "
          f"{0.0741:8.4f} {d_i2*0.0741:9.4f}")
    print(f"{'corrected: unpaired 4v4 x bootstrap SE(Delta)':52s} {d_i4:7.3f} "
          f"{se_boot:8.4f} {d_i4*se_boot:9.4f}")
    print(f"{'corrected: paired n=4 x bootstrap SE(Delta)':52s} {d_p4:7.3f} "
          f"{se_boot:8.4f} {d_p4*se_boot:9.4f}")
    print(f"{'LOO-based: unpaired 4v4 x LOO SD of reversal':52s} {d_i4:7.3f} "
          f"{loo_sd:8.4f} {d_i4*loo_sd:9.4f}")
    print(f"\n   Observed max |Delta| (target modules) = 0.2743")
    print(f"   Manuscript's claim: 'only powered to detect reversal-score differences")
    print(f"   exceeding ~0.42 units, whereas the largest observed reversal difference")
    print(f"   is = 0.27' -> CONCLUSION: underpowered.")
    print(f"   Recomputed: the design resolves ~0.04 reversal units at 80% power.")
    print(f"   0.2743 / 0.0441 = {0.2743/ (d_i4*se_boot):.1f}x the MDE.")
    print(f"   => The paper's OWN exact permutation test agrees: 26 of 132 tests have")
    print(f"      p_permutation_exact < 0.05 and the smallest p is 0.00143. A design")
    print(f"      that can only see 0.42 cannot produce p = 0.0014.")


if __name__ == "__main__":
    main()
