#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
22_etiology_comparison.py  (core cross-aetiology result of the manuscript)
=========================================================================
Is the in-silico perturbation response of the DRG different between cancer
pain of different aetiology?

Contrast design
---------------
    GSE149648   rat DRG, Walker-256 tibial bone cancer (TCI),
                D7 / D14 vs Sham, GPL24782, n = 4 per group  -> "bone cancer pain"
    GSE160543   rat DRG, paclitaxel or oxaliplatin vs vehicle,
                GPL24782, n = 4 per group                    -> "CIPN"

Same species (Rattus norvegicus), same tissue (DRG), same sequencing platform
(GPL24782). What differs is the aetiology. That is what makes a cross-study
comparison interpretable at all; comparing rat DRG against human tumour would
not be.

Three changes relative to the previous revision, each made in response to a
specific reviewer criticism
--------------------------------------------------------------------------
(1) NETWORKS ARE BUILT PER CONTRAST, NOT POOLED.
    Previously PTX and OXA shared one network built on all 12 samples of
    GSE160543. Co-expression structure is itself treatment-dependent, so a
    shared network forces the two chemotherapies to be propagated through an
    identical graph and removes part of the very difference being tested.
    Each contrast now gets a network built on its own 8 samples. The pooled
    variant is still computed and reported as a sensitivity analysis.

(2) THE DIFFERENCE TEST IS AN EXACT LABEL PERMUTATION, NOT A BOOTSTRAP.
    At n = 4 vs 4 the entire label-assignment space is C(8,4) = 70 splits, so
    the null over the 70 x 70 = 4900 (A,B) split pairs can be enumerated
    exhaustively: no Monte-Carlo error, no seed dependence, and a hard floor
    on the p-value of 1/4901 rather than an artificial 1/2000. The bootstrap
    is retained ONLY for the confidence interval on Delta, which a
    permutation test does not provide.

(3) NEGATIVE-CONTROL MODULES ARE CARRIED THROUGH.
    If a pan-neuronal or erythrocyte gene set shows the same cross-aetiology
    Delta as the Na+ channel module, the Delta reflects compartment /
    composition change, not target-module specificity.

Outputs -> results/tables/etiology_*.csv
"""
from __future__ import annotations

import os
import re
import sys
import json
import gzip
from itertools import combinations

import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)

from _limma import (build_network, diffusion_operator, cos_batch,          # noqa: E402
                    ridge_baseline, muscle_score, bh, fit_limma,
                    TARGETS, CONTROL_MODULES,
                    degree_matched_sets, all_splits, perm_axes,
                    exact_two_sided_p, emp_p_two_sided)

RAW = os.path.join(ROOT, "data", "raw")
TAB = os.path.join(ROOT, "results", "tables")
os.makedirs(TAB, exist_ok=True)

SEED = 20261009
rng = np.random.default_rng(SEED)
N_BOOT = 2000
N_RANDSET = 500
NET_HV = 3000
NET_TOPK = 12
ALPHA = 0.55

ALL_MODS = dict(TARGETS)
ALL_MODS.update(CONTROL_MODULES)


def log(m):
    print(f"[22] {m}", flush=True)


# ---------------------------------------------------------------------------
def read_fpkm(dirname: str) -> tuple[pd.DataFrame, dict]:
    """Read a GEO RAW directory of per-sample files into a gene x sample FPKM
    matrix.

    The two studies ship different column layouts and GSE160543 is not even
    internally uniform (6 columns for most members, 13 for some, because of
    trailing empty tab-delimited fields). Reading a fixed column index picks
    up the Entrez id and yields a matrix of identifiers - a silent failure
    that produces zero fold change everywhere. The FPKM column is therefore
    located from the header, with a right-to-left numeric fallback, and the
    result is asserted to behave like expression.
    """
    files = sorted([f for f in os.listdir(os.path.join(RAW, dirname))
                    if f.endswith(".txt.gz") or f.endswith(".txt")])
    if not files:
        raise FileNotFoundError(f"no member files in {dirname}")
    series, order, fnames, skipped = {}, [], {}, []
    for fn in files:
        gsm = fn.split("_")[0]
        path = os.path.join(RAW, dirname, fn)
        op = gzip.open if fn.endswith(".gz") else open
        rows = {}
        try:
            with op(path, "rt") as fh:
                vidx = -1
                for line in fh:
                    p = line.rstrip("\n").split("\t")
                    if vidx < 0:
                        up = [c.strip().upper() for c in p]
                        if "FPKM" in up:
                            vidx = up.index("FPKM")
                            continue
                        vidx = -2
                    g = p[0].strip()
                    if not g:
                        continue
                    tok = None
                    if 0 <= vidx < len(p):
                        tok = p[vidx].strip()
                    else:
                        for j in range(len(p) - 1, 0, -1):
                            try:
                                float(p[j].strip())
                                tok = p[j].strip()
                                break
                            except ValueError:
                                continue
                    if tok is None:
                        continue
                    try:
                        rows[g] = float(tok)
                    except ValueError:
                        continue
        except (EOFError, OSError) as e:
            skipped.append((fn, type(e).__name__))
            log(f"  WARNING: skipped unreadable member {fn} ({e})")
            continue
        if not rows:
            raise RuntimeError(f"no numeric expression parsed from {fn}")
        series[gsm] = pd.Series(rows, dtype=float)
        fnames[gsm] = fn
        order.append(gsm)
    if skipped:
        log(f"  {len(skipped)} member(s) skipped -> n={len(order)}")
    M = pd.DataFrame(series)[order].fillna(0.0)
    frac_int = float(np.isclose(M.values % 1, 0).mean())
    log(f"  {dirname}: {M.shape[0]} genes x {M.shape[1]} samples; "
        f"frac integer-valued={frac_int:.3f}")
    if frac_int > 0.95:
        raise RuntimeError(f"{dirname}: matrix looks like identifiers, not FPKM")
    return M, fnames


def log2fpkm(fpkm: pd.DataFrame) -> pd.DataFrame:
    """FPKM is nominally length- and depth-normalised, but per-sample FPKM
    totals are not equal across groups in these series (a 6-9% global offset
    in the bone-cancer-pain set). Median-of-expressed-genes scaling is applied
    before the log transform - identical to the discovery pipeline - so the two
    studies are on the same footing.
    """
    med = fpkm.where(fpkm > 0).median(axis=0)
    sf = med / med.mean()
    sf = sf.replace(0, np.nan).fillna(1.0)
    return np.log2(fpkm.div(sf, axis=1) + 1.0)


def group_from_filename(name: str, key: str) -> str:
    """Group assignment from the member filename, which is unambiguous.

    A naive substring test such as `"tci_d1" in name` also matches `tci_d14`
    and collapses D14 into D1; the day token is therefore matched with
    delimiters.
    """
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


# ---------------------------------------------------------------------------
def run(scope: str, contrasts: dict, L_by_study: dict, common_genes: list,
        P_common: dict, W_ops: dict, est_rows: list, boot_rows: list):
    """Point estimates + exact-permutation / bootstrap / random-set tests for
    every module and every aetiology pair, under one network scope."""

    def reversal(k, tidx, sign):
        L, ca, co = contrasts[k]
        axis = (L[co].mean(axis=1) - L[ca].mean(axis=1)).reindex(common_genes).values
        seed = np.zeros(len(common_genes))
        seed[tidx] = sign
        return float(cos_batch((P_common[k] @ seed).reshape(-1, 1), axis)[0])

    def additive(k, tidx, sign):
        L, ca, co = contrasts[k]
        axis = (L[co].mean(axis=1) - L[ca].mean(axis=1)).reindex(common_genes).values
        seed = np.zeros(len(common_genes))
        seed[tidx] = sign
        return float(cos_batch(seed.reshape(-1, 1), axis)[0])

    # ---- exhaustive label-permutation axes, one matrix per contrast ----
    perm = {}
    for k, (L, ca, co) in contrasts.items():
        cols = list(ca) + list(co)
        Lv = L.reindex(common_genes)[cols].values
        sp = all_splits(len(ca), len(co))
        Ax = perm_axes(Lv, sp)                       # genes x n_splits
        true_mask = np.array([True] * len(ca) + [False] * len(co))
        obs_i = int(np.flatnonzero((sp == true_mask).all(axis=1))[0])
        perm[k] = {"axes": Ax, "obs_i": obs_i, "n_split": sp.shape[0],
                   "shares_control": co}
        log(f"  [{scope}] {k}: {sp.shape[0]} label splits, observed index {obs_i}")

    idx_of = {g: i for i, g in enumerate(common_genes)}
    gidx = np.array(sorted(idx_of.values()))

    pairs = []
    keys = list(contrasts)
    for i in range(len(keys)):
        for j in range(len(keys)):
            if i < j:
                pairs.append((keys[i], keys[j]))

    for ptype, sign in [("knockout", -1.0), ("overexpression", 1.0)]:
        for mod, tlist in ALL_MODS.items():
            present = [g for g in tlist if g in idx_of]
            if len(present) < 3:
                log(f"  skip {mod}: only {len(present)} of {len(tlist)} genes "
                    f"present on the common network")
                continue
            tidx = np.array([idx_of[g] for g in present])
            seed = np.zeros(len(common_genes))
            seed[tidx] = sign
            delta = {k: (P_common[k] @ seed).reshape(-1, 1) for k in contrasts}

            row = {"network_scope": scope, "module": mod,
                   "kind": "control" if mod.startswith("CTRL") else "target",
                   "n_targets": len(present), "genes_used": ";".join(present),
                   "perturbation": ptype}
            for k in contrasts:
                row[f"rev_{k}"] = reversal(k, tidx, sign)
                row[f"add_{k}"] = additive(k, tidx, sign)
                L, ca, co = contrasts[k]
                row[f"ridge_{k}"] = ridge_baseline(
                    L.reindex(common_genes)[list(ca) + list(co)],
                    common_genes, present, sign,
                    (L[co].mean(axis=1) - L[ca].mean(axis=1)).reindex(
                        common_genes).values, loso=True)
            est_rows.append(row)

            for ka, kb in pairs:
                sa = cos_batch(delta[ka], perm[ka]["axes"])[0]   # n_split
                sb = cos_batch(delta[kb], perm[kb]["axes"])[0]
                # Self-check: the observed split of the permutation matrix
                # must reproduce the point estimate exactly. Without this the
                # (historically real) possibility that the two code paths
                # normalise differently would pass unnoticed.
                for kk, ss in [(ka, sa), (kb, sb)]:
                    chk = reversal(kk, tidx, sign)
                    if abs(ss[perm[kk]["obs_i"]] - chk) > 1e-9:
                        raise RuntimeError(
                            f"{kk}/{mod}/{ptype}: permutation axis "
                            f"({ss[perm[kk]['obs_i']]:.6f}) disagrees with "
                            f"point estimate ({chk:.6f})")
                # Shared-control contrasts use the exchange-8-case null: the
                # control term of the pain axis is common to both sides, so
                # the only exchangeable labels are the eight case animals
                # (4 A-case + 4 B-case); enumerate all C(8,4) = 70 splits of
                # those eight with the shared control fixed. Pairing the
                # split index of two independent per-contrast enumerations
                # (the previous implementation) corresponds to no
                # randomisation of the eight case labels and is a different,
                # too-narrow null.
                paired = (list(perm[ka]["shares_control"])
                          == list(perm[kb]["shares_control"]))
                if paired:
                    La, caa, coa = contrasts[ka]
                    Lb, cab, cob = contrasts[kb]
                    A_case = La.reindex(common_genes)[list(caa)].values
                    B_case = Lb.reindex(common_genes)[list(cab)].values
                    c0a = La.reindex(common_genes)[list(coa)].mean(axis=1).values
                    c0b = Lb.reindex(common_genes)[list(cob)].mean(axis=1).values
                    C8 = np.hstack([A_case, B_case])   # genes x 8 case animals
                    sp8 = all_splits(4, 4)
                    null = np.empty(sp8.shape[0])
                    for r in range(sp8.shape[0]):
                        s1 = sp8[r]
                        ax_a = c0a - C8[:, s1].mean(axis=1)
                        ax_b = c0b - C8[:, ~s1].mean(axis=1)
                        ra = float(cos_batch(delta[ka], ax_a)[0])
                        rb = float(cos_batch(delta[kb], ax_b)[0])
                        null[r] = ra - rb
                    # self-check: the true-assignment split must reproduce
                    # the observed point estimates exactly
                    s_true = np.array([True] * 4 + [False] * 4)
                    chk_a = float(cos_batch(
                        delta[ka],
                        c0a - C8[:, s_true].mean(axis=1))[0])
                    chk_b = float(cos_batch(
                        delta[kb],
                        c0b - C8[:, ~s_true].mean(axis=1))[0])
                    if (abs(chk_a - sa[perm[ka]["obs_i"]]) > 1e-9
                            or abs(chk_b - sb[perm[kb]["obs_i"]]) > 1e-9):
                        raise RuntimeError(
                            f"{ka}/{kb}/{mod}/{ptype}: 8-case null observed "
                            f"split disagrees with point estimates")
                else:
                    null = (sa[:, None] - sb[None, :]).ravel()
                obs = sa[perm[ka]["obs_i"]] - sb[perm[kb]["obs_i"]]
                p_perm = exact_two_sided_p(obs, null)

                # ---- bootstrap CI for Delta (permutation gives no CI) ----
                La, caa, coa = contrasts[ka]
                Lb, cab, cob = contrasts[kb]
                Aa = La.reindex(common_genes).values
                Ab = Lb.reindex(common_genes).values
                ia_ca = [La.columns.get_loc(c) for c in caa]
                ia_co = [La.columns.get_loc(c) for c in coa]
                ib_ca = [Lb.columns.get_loc(c) for c in cab]
                ib_co = [Lb.columns.get_loc(c) for c in cob]
                dB = np.empty(N_BOOT)
                ndeg = 0
                for b in range(N_BOOT):
                    if paired:
                        pick = rng.choice(len(ia_co), len(ia_co), replace=True)
                        ctr_a = [ia_co[t] for t in pick]
                        ctr_b = [ib_co[t] for t in pick]
                    else:
                        ctr_a = rng.choice(ia_co, len(ia_co), replace=True)
                        ctr_b = rng.choice(ib_co, len(ib_co), replace=True)
                    a1 = (Aa[:, ctr_a].mean(1)
                          - Aa[:, rng.choice(ia_ca, len(ia_ca), replace=True)].mean(1))
                    a2 = (Ab[:, ctr_b].mean(1)
                          - Ab[:, rng.choice(ib_ca, len(ib_ca), replace=True)].mean(1))
                    if np.linalg.norm(a1 - a1.mean()) < 1e-9 or \
                       np.linalg.norm(a2 - a2.mean()) < 1e-9:
                        ndeg += 1                      # degenerate resample
                    dB[b] = cos_batch(delta[ka], a1)[0] - cos_batch(delta[kb], a2)[0]
                dB = dB[:N_BOOT - ndeg] if ndeg else dB

                # ---- random gene-set null: is THIS module special? ----
                # Batched: one (n x N_RANDSET) matmul per aetiology instead of
                # N_RANDSET separate matrix-vector products.
                msets = degree_matched_sets(
                    (np.abs(W_ops[ka]) + np.abs(W_ops[kb])) / 2.0,
                    tidx, N_RANDSET, rng, n_bins=10)
                Sr = np.zeros((len(common_genes), N_RANDSET))
                for r in range(N_RANDSET):
                    Sr[msets[r], r] = sign
                dR = (cos_batch(P_common[ka] @ Sr,
                                perm[ka]["axes"][:, perm[ka]["obs_i"]])
                      - cos_batch(P_common[kb] @ Sr,
                                  perm[kb]["axes"][:, perm[kb]["obs_i"]]))
                p_randset = emp_p_two_sided(obs, dR)

                boot_rows.append({
                    "network_scope": scope, "module": mod,
                    "kind": "control" if mod.startswith("CTRL") else "target",
                    "perturbation": ptype, "aetiology_A": ka, "aetiology_B": kb,
                    "paired_control": paired,
                    "rev_A": float(sa[perm[ka]["obs_i"]]),
                    "rev_B": float(sb[perm[kb]["obs_i"]]),
                    "delta_observed": float(obs),
                    "delta_boot_mean": float(dB.mean()),
                    "ci95_lo": float(np.percentile(dB, 2.5)),
                    "ci95_hi": float(np.percentile(dB, 97.5)),
                    "p_permutation_exact": p_perm,
                    "n_null": int(len(null)),
                    "p_bootstrap": float(min(1.0, max(
                        2.0 * min((dB <= 0).mean(), (dB >= 0).mean()),
                        1.0 / max(len(dB), 1)))),
                    "randset_mean": float(dR.mean()),
                    "randset_sd": float(dR.std()),
                    "z_vs_randset": float((obs - dR.mean()) / (dR.std() + 1e-12)),
                    "p_randset": p_randset,
                    "n_degenerate_boot": int(ndeg),
                })
    return pairs


# ---------------------------------------------------------------------------
def de_positive_controls(contrasts, L_by_study):
    """Differential expression of canonical injury / CIPN markers, so the
    reader can see that the two aetiologies really do differ at the level of
    individual genes before any network machinery is applied."""
    MARKERS = ["Atf3", "Sprr1a", "Gal", "Npy", "Calca", "Tac1", "Scn9a",
               "Scn10a", "Oprm1", "Cacna2d1", "Trpv1", "Gfap", "Aif1",
               "Cd68", "Ccl2", "Il1b", "Tnf", "Mpz", "Nefl"]
    rows = []
    for k, (L, ca, co) in contrasts.items():
        cols = list(ca) + list(co)
        X = L[cols]
        gvec = np.array([1.0] * len(ca) + [0.0] * len(co))
        ms = muscle_score(X).values.astype(float)
        mc = ms - ms.mean()
        design = np.column_stack([np.ones(len(cols)), gvec, mc])
        r = fit_limma(X, design, coef_idx=1, trend=True)
        r["q"] = bh(r["P"].values)
        for g in MARKERS:
            if g in r.index:
                rows.append({"contrast": k, "gene": g,
                             "logFC": round(float(r.loc[g, "logFC"]), 3),
                             "t": round(float(r.loc[g, "t"]), 3),
                             "P": float(f"{r.loc[g, 'P']:.3g}"),
                             "q": round(float(r.loc[g, "q"]), 4)})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
def module_logfc_splits(L, ca, co, gene_list, splits):
    """Mean log2 fold change (case - control) of every module, for every label
    split. Returns dict module -> (n_splits,) array.

    This is the better-powered cross-aetiology endpoint. A reversal cosine is
    a network-dependent, scale-free quantity whose magnitude is hard to
    interpret; the module's mean fold change is not. And because logFC does
    not involve the empirical-Bayes prior, it can be recomputed for all
    C(8,4) = 70 splits cheaply, which makes the permutation test exact.
    """
    cols = list(ca) + list(co)
    Lv = L.reindex(gene_list)[cols].values
    n = splits.shape[0]
    out = np.empty((Lv.shape[0], n))
    for r in range(n):
        case = np.flatnonzero(splits[r])
        ctrl = np.flatnonzero(~splits[r])
        out[:, r] = Lv[:, case].mean(axis=1) - Lv[:, ctrl].mean(axis=1)
    return out


def cross_module_logfc(contrasts, common_genes, markers=None):
    """Delta of module mean logFC between every pair of aetiologies, with an
    exact label-permutation p."""
    splits = {k: all_splits(len(ca), len(co)) for k, (_, ca, co) in contrasts.items()}
    obs_i = {k: int(np.flatnonzero(
        (splits[k] == np.array([True] * len(ca) + [False] * len(co))).all(axis=1))[0])
        for k, (_, ca, co) in contrasts.items()}
    mats = {}
    for k, (L, ca, co) in contrasts.items():
        mats[k] = module_logfc_splits(L, ca, co, common_genes, splits[k])
    mods = {m: [g for g in gs if g in common_genes] for m, gs in ALL_MODS.items()}
    mods = {m: g for m, g in mods.items() if len(g) >= 3}
    pos = {g: i for i, g in enumerate(common_genes)}
    keys = list(contrasts)
    rows = []
    for i in range(len(keys)):
        for j in range(len(keys)):
            if i >= j:
                continue
            ka, kb = keys[i], keys[j]
            La, caa, coa = contrasts[ka]
            Lb, cab, cob = contrasts[kb]
            paired = list(coa) == list(cob)
            if paired:
                # Exchange-8-case null. For a shared control the control term
                # cancels identically in the logFC difference
                #   Δ = [mean(A_case) - mean(C)] - [mean(B_case) - mean(C)]
                #     = mean(A_case) - mean(B_case),
                # so the only exchangeable labels are the eight case animals
                # (4 A-case + 4 B-case): enumerate all C(8,4) = 70 splits of
                # those eight with the shared control fixed. Pairing the split
                # INDEX of two independent 70-split enumerations (the previous
                # implementation) corresponds to no randomisation of the eight
                # case labels and is a different, too-narrow null.
                A_case = La.reindex(common_genes)[list(caa)].values   # genes x 4
                B_case = Lb.reindex(common_genes)[list(cab)].values   # genes x 4
                sp8 = all_splits(4, 4)                                # 70 x 8 bool
            for m, gs in mods.items():
                ri = np.array([pos[g] for g in gs])
                a = mats[ka][ri].mean(axis=0)
                b = mats[kb][ri].mean(axis=0)
                if paired:
                    C8 = np.hstack([A_case[ri], B_case[ri]])   # genes x 8 case animals
                    null = np.empty(sp8.shape[0])
                    for r in range(sp8.shape[0]):
                        s1 = sp8[r]
                        null[r] = (C8[:, s1].mean(axis=1).mean()
                                   - C8[:, ~s1].mean(axis=1).mean())
                else:
                    null = (a[:, None] - b[None, :]).ravel()
                obs = a[obs_i[ka]] - b[obs_i[kb]]
                rows.append({
                    "module": m,
                    "kind": "control" if m.startswith("CTRL") else "target",
                    "aetiology_A": ka, "aetiology_B": kb,
                    "logFC_A": float(a[obs_i[ka]]), "logFC_B": float(b[obs_i[kb]]),
                    "delta_logFC": float(obs),
                    "null_lo": float(np.percentile(null, 2.5)),
                    "null_hi": float(np.percentile(null, 97.5)),
                    "p_permutation_exact": exact_two_sided_p(obs, null),
                    "n_null": int(len(null)), "paired_control": paired,
                })
    df = pd.DataFrame(rows)
    df["q_permutation"] = bh(df["p_permutation_exact"].values)
    return df.sort_values("p_permutation_exact")


def cross_marker_logfc(contrasts, markers):
    """Same test at the level of individual canonical marker genes."""
    gene_list = [g for g in markers if all(g in L.index for L, _, _ in contrasts.values())]
    splits = {k: all_splits(len(ca), len(co)) for k, (_, ca, co) in contrasts.items()}
    obs_i = {k: int(np.flatnonzero(
        (splits[k] == np.array([True] * len(ca) + [False] * len(co))).all(axis=1))[0])
        for k, (_, ca, co) in contrasts.items()}
    mats = {k: module_logfc_splits(L, ca, co, gene_list, splits[k])
            for k, (L, ca, co) in contrasts.items()}
    keys = list(contrasts)
    rows = []
    for i in range(len(keys)):
        for j in range(len(keys)):
            if i >= j:
                continue
            ka, kb = keys[i], keys[j]
            La, caa, coa = contrasts[ka]
            Lb, cab, cob = contrasts[kb]
            paired = list(coa) == list(cob)
            if paired:
                # Exchange-8-case null (see cross_module_logfc): the shared
                # control cancels, so the eight case animals are the only
                # exchangeable labels. C(8,4) = 70 splits, control fixed.
                A_case = La.reindex(gene_list)[list(caa)].values   # genes x 4
                B_case = Lb.reindex(gene_list)[list(cab)].values
                sp8 = all_splits(4, 4)
            for gi, g in enumerate(gene_list):
                a = mats[ka][gi]
                b = mats[kb][gi]
                if paired:
                    c8 = np.concatenate([A_case[gi], B_case[gi]])   # 8 case animals
                    null = np.empty(sp8.shape[0])
                    for r in range(sp8.shape[0]):
                        s1 = sp8[r]
                        null[r] = c8[s1].mean() - c8[~s1].mean()
                else:
                    null = (a[:, None] - b[None, :]).ravel()
                obs = a[obs_i[ka]] - b[obs_i[kb]]
                rows.append({"gene": g, "aetiology_A": ka, "aetiology_B": kb,
                             "logFC_A": float(a[obs_i[ka]]),
                             "logFC_B": float(b[obs_i[kb]]),
                             "delta_logFC": float(obs),
                             "p_permutation_exact": exact_two_sided_p(obs, null),
                             "paired_control": paired})
    df = pd.DataFrame(rows)
    df["q_permutation"] = bh(df["p_permutation_exact"].values)
    return df.sort_values("p_permutation_exact")


# ---------------------------------------------------------------------------
def main():
    log("=== load both studies ===")
    bone_fpkm, fb = read_fpkm("GSE149648")
    chemo_fpkm, fc = read_fpkm("GSE160543")
    gb = {g: group_from_filename(fb[g], "bone") for g in bone_fpkm.columns}
    gc = {g: group_from_filename(fc[g], "chemo") for g in chemo_fpkm.columns}
    log(f"  GSE149648 groups={pd.Series(list(gb.values())).value_counts().to_dict()}")
    log(f"  GSE160543 groups={pd.Series(list(gc.values())).value_counts().to_dict()}")

    keep_b = ((bone_fpkm >= 1.0).sum(axis=1) >= 4) & (bone_fpkm.mean(axis=1) >= 0.5)
    keep_c = ((chemo_fpkm >= 1.0).sum(axis=1) >= 4) & (chemo_fpkm.mean(axis=1) >= 0.5)
    Lb = log2fpkm(bone_fpkm).loc[keep_b]
    Lc = log2fpkm(chemo_fpkm).loc[keep_c]
    shared = sorted(set(Lb.index) & set(Lc.index))
    Lb, Lc = Lb.loc[shared], Lc.loc[shared]
    log(f"  shared genes after filtering: {len(shared)}")

    for name, L, gmap in [("GSE149648", Lb, gb), ("GSE160543", Lc, gc)]:
        ms = muscle_score(L)
        df = pd.DataFrame({"g": [gmap[s] for s in L.columns], "ms": ms.values})
        log(f"  {name} muscle score by group: " +
            ", ".join(f"{k}={v.mean():+.2f}" for k, v in df.groupby("g")["ms"]))

    contrasts = {
        # D14 is the established bone-cancer-pain state. D1 shows essentially
        # no transcriptional response, so D7 and D14 are the informative
        # reference points and both are carried.
        "Bone_D14": (Lb, [s for s in Lb.columns if gb.get(s) == "D14"],
                     [s for s in Lb.columns if gb.get(s) == "Sham"]),
        "Bone_D7": (Lb, [s for s in Lb.columns if gb.get(s) == "D7"],
                    [s for s in Lb.columns if gb.get(s) == "Sham"]),
        "PTX": (Lc, [s for s in Lc.columns if gc.get(s) == "PTX"],
                [s for s in Lc.columns if gc.get(s) == "Vehicle"]),
        "OXA": (Lc, [s for s in Lc.columns if gc.get(s) == "OXA"],
                [s for s in Lc.columns if gc.get(s) == "Vehicle"]),
    }
    for k, (_, ca, co) in contrasts.items():
        if len(ca) == 0 or len(co) == 0:
            raise RuntimeError(f"{k}: empty case/control set - group labels are "
                               f"{sorted(set(gb.values()) | set(gc.values()))}")
        log(f"  {k}: case n={len(ca)} control n={len(co)}")

    # ---- cross-aetiology test on module mean logFC (primary endpoint) ----
    log("=== cross-aetiology Delta of module mean logFC (exact permutation) ===")
    modx = cross_module_logfc(contrasts, shared)
    modx.to_csv(os.path.join(TAB, "etiology_module_logfc_cross.csv"), index=False)
    log(modx.head(20)[["module", "kind", "aetiology_A", "aetiology_B",
                       "logFC_A", "logFC_B", "delta_logFC",
                       "p_permutation_exact", "q_permutation"]].to_string(index=False))

    log("=== cross-aetiology Delta of canonical markers (exact permutation) ===")
    MARKERS = ["Atf3", "Sprr1a", "Gal", "Npy", "Calca", "Tac1", "Scn9a",
               "Scn10a", "Oprm1", "Cacna2d1", "Trpv1", "Gfap", "Aif1",
               "Cd68", "Ccl2", "Il1b", "Tnf", "Mpz", "Nefl"]
    mkx = cross_marker_logfc(contrasts, MARKERS)
    mkx.to_csv(os.path.join(TAB, "etiology_marker_cross.csv"), index=False)
    log(mkx.head(24).to_string(index=False))

    pc = de_positive_controls(contrasts, {"b": Lb, "c": Lc})
    pc.to_csv(os.path.join(TAB, "etiology_marker_de.csv"), index=False)
    log("\n=== marker-level DE (muscle-adjusted limma-trend) ===")
    piv = pc.pivot(index="gene", columns="contrast", values="logFC")
    log(piv.round(2).to_string())

    est_rows, boot_rows = [], []
    summary = {"seed": SEED, "n_boot": N_BOOT, "n_randset": N_RANDSET,
               "alpha": ALPHA, "hv": NET_HV, "topk": NET_TOPK}

    # One node set for every contrast. Selecting high-variance genes
    # independently per contrast gives four different 3000-gene sets whose
    # intersection is only a few hundred genes, which silently deletes most
    # analgesic target genes from the graph. Ranking is by the MEAN variance
    # across the four contrasts so no single contrast dominates the choice,
    # and every module gene is force-included.
    MUST = sorted({g for gs in ALL_MODS.values() for g in gs} & set(shared))
    vsum = None
    for k, (L, ca, co) in contrasts.items():
        v = L[list(ca) + list(co)].var(axis=1)
        vsum = v if vsum is None else vsum.add(v, fill_value=0.0)
    vmean = (vsum / len(contrasts)).loc[shared].sort_values(ascending=False)
    NET_GENES = list(vmean.index[:NET_HV]) + \
        [g for g in MUST if g not in set(vmean.index[:NET_HV])]
    log(f"  common node set: {len(NET_GENES)} genes "
        f"(top-{NET_HV} by mean variance + {len(MUST)} module genes)")
    log("  module genes retained: " +
        ", ".join(f"{m}={len([g for g in gs if g in NET_GENES])}/{len(gs)}"
                  for m, gs in ALL_MODS.items()))

    for scope in ["per_contrast", "pooled_by_study"]:
        log(f"=== networks: scope={scope} ===")
        P_ops, W_ops = {}, {}
        for k, (L, ca, co) in contrasts.items():
            Ln = L[list(ca) + list(co)] if scope == "per_contrast" else (
                Lb if k.startswith("Bone") else Lc)
            genes, W = build_network(Ln, topk=NET_TOPK, genes=NET_GENES)
            W_ops[k] = W
            P_ops[k] = diffusion_operator(W, alpha=ALPHA, exact=True)
        common_genes = NET_GENES
        log(f"  nodes: {len(common_genes)} (identical for all {len(P_ops)} "
            f"contrasts by construction)")
        P_common = P_ops
        pairs = run(scope, contrasts, {"b": Lb, "c": Lc}, common_genes,
                    P_common, W_ops, est_rows, boot_rows)
        summary[f"common_network_genes_{scope}"] = len(common_genes)

    est = pd.DataFrame(est_rows)
    est.to_csv(os.path.join(TAB, "etiology_reversal_point_estimates.csv"), index=False)

    boot = pd.DataFrame(boot_rows)
    # --- correct multiplicity correction -----------------------------------
    # Knockout and overexpression are exact mirror images: reversal(KO) =
    # -reversal(OE) for every (module, aetiology pair), so the EXACT
    # permutation p is identical (verified max|diff| = 0 across all 132
    # module-pair tests). They are therefore ONE independent test, not two.
    # Applying BH over the raw 2 x N rows double-counts each hypothesis and
    # inflates the denominator. We instead apply BH over the deduplicated
    # knockout set *within each network scope* (132 independent tests per
    # scope, 264 across both scopes) and merge the q back to the OE twin.
    # The randset p is not exactly identical between KO and OE (different
    # degree-matched draws), but the question is the same, so the same
    # deduplication is used; this is the conservative choice.
    ded = boot[boot.perturbation == "knockout"].copy()
    for sc in boot["network_scope"].unique():
        m = ded["network_scope"] == sc
        ded.loc[m, "q_permutation"] = bh(ded.loc[m, "p_permutation_exact"].values)
        ded.loc[m, "q_randset"] = bh(ded.loc[m, "p_randset"].values)
    qmap = ded[["network_scope", "module", "aetiology_A", "aetiology_B",
                "q_permutation", "q_randset"]]
    boot = boot.merge(qmap, on=["network_scope", "module", "aetiology_A",
                                "aetiology_B"], how="left", suffixes=("", "_ded"))
    # qmap columns land as q_permutation / q_randset when boot lacks them, or
    # as q_permutation_ded / q_randset_ded when a clash exists; normalise both.
    for base in ["q_permutation", "q_randset"]:
        ded_col = base + "_ded"
        if ded_col in boot.columns:
            boot[base] = boot.pop(ded_col)
    boot = boot.sort_values(["p_permutation_exact", "module"])
    boot.to_csv(os.path.join(TAB, "etiology_bootstrap_differences.csv"), index=False)

    log("\n=== Delta ranked by exact permutation p (per_contrast) ===")
    sel = boot[boot.network_scope == "per_contrast"]
    log(sel[["module", "kind", "perturbation", "aetiology_A", "aetiology_B",
             "rev_A", "rev_B", "delta_observed", "ci95_lo", "ci95_hi",
             "p_permutation_exact", "q_permutation", "z_vs_randset",
             "p_randset"]].head(30).to_string(index=False))

    log("\n=== negative-control modules: do non-target sets move as much? ===")
    nc = sel[sel.kind == "control"][["module", "aetiology_A", "aetiology_B",
                                     "delta_observed", "p_permutation_exact",
                                     "q_permutation"]]
    log(nc.to_string(index=False))

    sel = boot[boot.network_scope == "per_contrast"]
    ded_all = boot[boot.perturbation == "knockout"]
    ded_pc = ded_all[ded_all.network_scope == "per_contrast"]
    summary.update({
        "studies": {"GSE149648": "bone cancer pain (Walker 256 tibia), GPL24782, rat DRG",
                    "GSE160543": "paclitaxel/oxaliplatin CIPN, GPL24782, rat DRG"},
        "shared_genes": len(shared),
        "contrasts": {k: {"case": len(v[1]), "control": len(v[2])}
                      for k, v in contrasts.items()},
        "pairs_tested": [f"{a}|{b}" for a, b in pairs],
        "n_rows_raw": int(len(boot)),
        "n_independent_tests_per_scope": int(len(ded_pc)),
        "n_independent_tests_total": int(len(ded_all)),
        "n_sig_perm_q05_indep": int((ded_all.q_permutation < 0.05).sum()),
        "n_sig_perm_q05_percontrast_indep": int((ded_pc.q_permutation < 0.05).sum()),
        "n_sig_randset_q05_percontrast_indep": int((ded_pc.q_randset < 0.05).sum()),
        "min_q_permutation_per_contrast": float(ded_pc.q_permutation.min()),
        "min_q_randset_per_contrast": float(ded_pc.q_randset.min()),
        "max_degenerate_boot": int(boot.n_degenerate_boot.max()),
    })
    with open(os.path.join(TAB, "etiology_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2, default=str)
    log(json.dumps(summary, indent=2, default=str))
    log("DONE")


if __name__ == "__main__":
    main()
