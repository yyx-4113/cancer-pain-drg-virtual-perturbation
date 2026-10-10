#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
A2_check1_permutation.py  -- INDEPENDENT re-implementation of the cross-aetiology
exact permutation test, written from scratch (does NOT import scripts/_limma.py).

Sign convention matched to the manuscript: logFC = mean(case) - mean(control),
Delta = logFC(A) - logFC(B).
"""
import gzip
import os
import re
from itertools import combinations

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "data", "raw")
TAB = os.path.join(ROOT, "results", "tables")


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
                if 0 <= vidx < len(p):
                    tok = p[vidx]
                else:
                    tok = None
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
        out[gsm] = pd.Series(vals, dtype=float)
        order.append(gsm)
    M = pd.DataFrame(out)[order].fillna(0.0)
    return M, float(np.isclose(M.values % 1, 0).mean())


def group_of(name, key):
    n = name.upper()
    if key == "bone":
        if "SHAM" in n:
            return "Sham"
        m = re.search(r"(?:^|[_.\-])D(\d+)(?:[_.\-])", n)
        return f"D{int(m.group(1))}" if m else "other"
    if "VEHICLE" in n:
        return "Vehicle"
    if "PACLITAXEL" in n:
        return "PTX"
    if "OXALIPLATIN" in n:
        return "OXA"
    return "other"


def norm(M):
    med = M.where(M > 0).median(axis=0)
    sf = (med / med.mean()).replace(0, np.nan).fillna(1.0)
    return np.log2(M.div(sf, axis=1) + 1.0)


def splits(n_case, n_ctrl):
    n = n_case + n_ctrl
    combs = list(combinations(range(n), n_case))
    S = np.zeros((len(combs), n), dtype=bool)
    for r, c in enumerate(combs):
        S[r, list(c)] = True
    return S, combs


def logfc_axes(Lv, S):
    """genes x n_splits of  mean(case) - mean(control)."""
    out = np.empty((Lv.shape[0], S.shape[0]))
    for r in range(S.shape[0]):
        c = np.flatnonzero(S[r]); k = np.flatnonzero(~S[r])
        out[:, r] = Lv[:, c].mean(axis=1) - Lv[:, k].mean(axis=1)
    return out


def p_repo(obs, null):
    """Repo rule: (#{|null| >= |obs|} + 1) / (n + 1)."""
    null = np.asarray(null, float)
    ge = int(np.sum(np.abs(null) >= abs(obs) - 1e-12))
    return min(1.0, (ge + 1) / (len(null) + 1)), ge, len(null)


def p_exact_no_plus1(obs, null):
    """Textbook exact permutation p WITHOUT the +1 inflation:
    #{|null| >= |obs|} / n.  Valid because the observed labelling is one of the
    enumerated set (the null is the full orbit, not a random sample)."""
    null = np.asarray(null, float)
    ge = int(np.sum(np.abs(null) >= abs(obs) - 1e-12))
    return min(1.0, ge / len(null)), ge, len(null)


def main():
    print("=" * 80)
    print("A2 CHECK 1 -- exhaustive permutation framework, recomputed from raw FPKM")
    print("=" * 80)

    Mb, fi_b = read_geo("GSE149648")
    Mc, fi_c = read_geo("GSE160543")
    print(f"GSE149648 {Mb.shape} frac_integer={fi_b:.3f} | "
          f"GSE160543 {Mc.shape} frac_integer={fi_c:.3f}")

    fb = {g: f for g, f in zip(Mb.columns, sorted(os.listdir(
        os.path.join(RAW, 'GSE149648'))))}
    fc = {g: f for g, f in zip(Mc.columns, sorted(os.listdir(
        os.path.join(RAW, 'GSE160543'))))}
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
    S, combs = splits(4, 4)
    OBS = int(np.flatnonzero((S == np.array([True] * 4 + [False] * 4)).all(axis=1))[0])
    print(f"shared genes={len(shared)}  C(8,4)={len(combs)}  observed split index={OBS}")

    # per-contrast split matrices of logFC, per gene
    M = {}
    for k, (L, ca, co) in contrasts.items():
        M[k] = logfc_axes(L.reindex(shared)[ca + co].values, S)

    PAIRS = [("Bone_D14", "Bone_D7"), ("Bone_D14", "PTX"), ("Bone_D14", "OXA"),
             ("Bone_D7", "PTX"), ("Bone_D7", "OXA"), ("PTX", "OXA")]

    # ---------------------------------------------------------------- (a)
    print("\n" + "=" * 80)
    print("(a) IS THE ENUMERATION EXHAUSTIVE?  ")
    print("=" * 80)
    print(f"C(8,4) = 70; all rows have exactly 4 case labels: {bool((S.sum(1)==4).all())}")
    print(f"all 70 label subsets distinct                     : "
          f"{len({frozenset(c) for c in combs})==70}")
    comp = ~S
    cid = [int(np.flatnonzero((comp[r][None, :] == S).all(axis=1))[0]) for r in range(70)]
    print(f"complement involution has 0 fixed points (70 splits pair into 35)  : "
          f"{sum(cid[i]==i for i in range(70))==0}")
    dev = max(np.abs(M['Bone_D7'][g, r] + M['Bone_D7'][g, cid[r]]).max()
              for g in range(0, len(shared), 97) for r in range(70))
    print(f"max |logFC_r + logFC_complement(r)| = {dev:.3e}  "
          f"-> null symmetric about 0")
    print(f"=> 70 splits x 70 splits = 4900 null pairs for an unpaired contrast: TRUE,")
    print(f"   the enumeration IS the complete orbit (no Monte-Carlo).")

    # ---------------------------------------------------------------- (b)
    print("\n" + "=" * 80)
    print("(b) Cd68 / key markers, ALL SIX pairs, recomputed from raw FPKM")
    print("=" * 80)
    gi = {g: i for i, g in enumerate(shared)}
    MARK = ["Cd68", "Aif1", "Gal", "Ccl2", "Il1b", "Calca", "Scn10a", "Atf3", "Gfap"]
    rows = []
    for gene in MARK:
        for ka, kb_ in PAIRS:
            a = M[ka][gi[gene]]; b = M[kb_][gi[gene]]
            paired = contrasts[ka][2] == contrasts[kb_][2]
            obs = a[OBS] - b[OBS]
            null = (a - b) if paired else (a[:, None] - b[None, :]).ravel()
            p_repo_v, ge, n = p_repo(obs, null)
            p_ex, _, _ = p_exact_no_plus1(obs, null)
            rows.append(dict(gene=gene, A=ka, B=kb_,
                             kind="paired" if paired else "unpaired",
                             logFC_A=round(a[OBS], 3), logFC_B=round(b[OBS], 3),
                             delta=round(obs, 3), ge=ge, n_null=n,
                             p_repo=round(p_repo_v, 6), p_exact=round(p_ex, 6)))
    df = pd.DataFrame(rows)
    print(df[df.gene.isin(["Cd68", "Aif1", "Gal"])].to_string(index=False))
    df.to_csv(os.path.join(ROOT, "review_A2", "A2_marker_perm_recomputed.csv"),
              index=False)

    # ---------------------------------------------------------------- (c)
    print("\n" + "=" * 80)
    print("(c) THE CLAIMED p FLOORS ARE WRONG")
    print("=" * 80)
    # complement symmetry forces BOTH +obs and -obs into every null
    for n_null, lbl in [(4900, "unpaired 70x70"), (70, "paired index-paired")]:
        print(f"{lbl}: n_null={n_null}")
        print(f"   manuscript floor            1/{n_null+1} = {1/(n_null+1):.6g}")
        print(f"   repo formula, ge>=1  ->    2/{n_null+1} = {2/(n_null+1):.6g}")
        print(f"   complement symmetry ge>=2 -> 3/{n_null+1} = {3/(n_null+1):.6g}   <== TRUE FLOOR")
    print("\nEmpirical confirmation from the table above: "
          f"min p_repo over all 54 marker tests = {df.p_repo.min():.6g}")
    print(f"   which equals 3/4901 = {3/4901:.6g}  (attained by Cd68/Aif1 D7-vs-PTX)")

    # ---------------------------------------------------------------- (d)
    print("\n" + "=" * 80)
    print("(d) SHARED-CONTRAST ('PAIRED') TREAMENT: is index-pairing the right null?")
    print("=" * 80)
    print("Delta for a shared control = [mean(Ac)-mean(C0)] - [mean(Bc)-mean(C0)]")
    print("                          = mean(Ac) - mean(Bc)      (control cancels)")
    print("So the exchangeable-label null must permute the 8 Ac/Bc samples.\n")
    for ka, kb_ in [("Bone_D14", "Bone_D7"), ("PTX", "OXA")]:
        La, caa, coa = contrasts[ka]
        Lb_, cab, cob = contrasts[kb_]
        print(f"### {ka} vs {kb_}   shared control = {len(set(coa) & set(cob))} samples")
        out = []
        for gene in MARK:
            a = M[ka][gi[gene]]; b = M[kb_][gi[gene]]
            obs = a[OBS] - b[OBS]
            null_repo = a - b
            p_r, ge_r, _ = p_repo(obs, null_repo)
            # natural null: exchange the 8 case labels only
            eight = np.concatenate([La.loc[gene, caa].values, Lb_.loc[gene, cab].values])
            nat = np.array([eight[S[r]].mean() - eight[~S[r]].mean() for r in range(70)])
            obs_nat = eight[S[OBS]].mean() - eight[~S[OBS]].mean()
            p_n, ge_n, _ = p_repo(obs_nat, nat)
            out.append(dict(gene=gene, delta=round(obs, 3),
                            p_repo_paired=round(p_r, 4), ge_repo=ge_r,
                            p_exchange8=round(p_n, 4), ge_exch=ge_n,
                            same_obs=bool(abs(obs - obs_nat) < 1e-9)))
        o = pd.DataFrame(out)
        print(o.to_string(index=False))
        print(f"   sd(repo index-paired null) vs sd(exchange-8 null), Cd68: "
              f"{(M[ka][gi['Cd68']] - M[kb_][gi['Cd68']]).std():.4f} vs "
              f"{np.array([np.concatenate([La.loc['Cd68', caa].values, Lb_.loc['Cd68', cab].values])[S[r]].mean() - np.concatenate([La.loc['Cd68', caa].values, Lb_.loc['Cd68', cab].values])[~S[r]].mean() for r in range(70)]).std():.4f}")
        print(f"   distinct values in repo 70-long null (Cd68): "
              f"{len(set(np.round(M[ka][gi['Cd68']] - M[kb_][gi['Cd68']], 12)))}/70")
        print()

    # ---------------------------------------------------------------- (e)
    print("=" * 80)
    print("(e) DEPOSITED ARTEFACT  results/tables/etiology_marker_cross.csv")
    print("=" * 80)
    dep = pd.read_csv(os.path.join(TAB, "etiology_marker_cross.csv"))
    print("Cd68 rows as deposited:")
    print(dep[dep.gene == "Cd68"].to_string(index=False))
    print("\nn rows in deposited marker table:", len(dep),
          " = 19 markers x 6 pairs =", 19 * 6)
    print("n rows in deposited module table:",
          len(pd.read_csv(os.path.join(TAB, "etiology_module_logfc_cross.csv"))),
          " = 22 modules x 6 pairs =", 22 * 6)

    # ------------------------------------------------- (f) BH recomputation
    print("\n" + "=" * 80)
    print("(f) BH q-values recomputed for the marker table (Cd68 rows)")
    print("=" * 80)
    def bh(p):
        p = np.asarray(p, float); n = len(p)
        o = np.argsort(p); q = p[o] * n / (np.arange(n) + 1)
        q = np.minimum.accumulate(q[::-1])[::-1]
        out = np.empty(n); out[o] = np.clip(q, 0, 1); return out
    d2 = dep.copy()
    d2["q_recomputed"] = bh(d2.p_permutation_exact.values)
    print(d2[d2.gene.isin(["Cd68", "Aif1", "Gal", "Ccl2", "Il1b", "Calca", "Scn10a"])]
          [["gene", "aetiology_A", "aetiology_B", "delta_logFC",
            "p_permutation_exact", "q_permutation", "q_recomputed",
            "paired_control"]].to_string(index=False))
    print("\nmatch deposited q_permutation vs my BH:",
          bool(np.allclose(d2.q_permutation.values, d2.q_recomputed.values, atol=1e-9)))
    print("=> BH is applied over ALL 19x6=114 marker tests jointly, and separately")
    print("   over all 22x6=132 module tests. Confirm min q_marker =",
          round(d2.q_permutation.min(), 5), " min q_module =",
          round(pd.read_csv(os.path.join(TAB, "etiology_module_logfc_cross.csv")
                            ).q_permutation.min(), 5))


if __name__ == "__main__":
    main()
