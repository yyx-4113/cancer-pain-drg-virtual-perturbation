#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
20_multistage_pipeline.py  (discovery cohort)
=============================================
Cancer-pain DRG multi-stage transcriptome: confounding-controlled differential
expression + in-silico perturbation (virtual knockout / overexpression) of
analgesic target modules, with permutation testing and baseline head-to-head.

All statistical / perturbation primitives live in _limma.py so that this
discovery script and 21_validate_pipeline.py execute identical code (guards
against version drift between manuscript numbers and analysis code).

Discovery dataset: GSE149648
    Walker 256 mammary carcinoma cells injected into rat tibia (TCI).
    DRG harvested Day 1 / Day 7 / Day 14 post-TCI; Sham (PBS) at Day 14.
    n = 4 per group, 16 samples, raw counts, platform GPL24782.

Design decisions (each is a claim the manuscript must defend)
------------------------------------------------------------
D1  Small-n RNA-seq: moderated t with empirical-Bayes variance shrinkage
    (Smyth 2004) plus a limma-trend intensity-dependent prior. A naive Welch
    test at n=4 has ~2-6 residual df and is badly calibrated at both tails.

D2  DRG dissection is systematically contaminated by paraspinal skeletal
    muscle. If contamination differs between groups, the apparent "pain
    signature" is a "muscle signature". We fit BOTH an unadjusted and a
    muscle-adjusted model and treat genes whose signal vanishes after
    adjustment as confounded rather than biological.

D3  Virtual perturbation: beyond differential expression, we diffuse a seed
    perturbation over a data-derived co-expression graph and score how far the
    perturbed cancer-pain profile moves back toward sham (reversal score =
    cosine with the pain->sham axis). Size-matched random gene sets build the
    empirical null.

D4  Baseline head-to-head: 2025 Nature Methods benchmarking showed perturbation
    models often fail to beat simple linear/additive baselines. We co-report
    (i) additive (no propagation) and (ii) ridge-linear alongside the network.

Performance note
----------------
The diffusion map is linear: delta_inf = (1-a)(I - a S)^-1 seed. We build the
operator once per dataset, after which every real and permuted perturbation is
a single matrix-vector product. This reduces an otherwise intractable
permutation test (102 conditions x 1000 permutations x 80 iterations of a
3000x3000 matmul) to seconds of BLAS work.
"""
from __future__ import annotations

import os
import re
import sys
import gzip
import glob
import json

import numpy as np
import pandas as pd
from scipy import stats

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _limma import (fit_limma, bh, muscle_score, build_network,          # noqa: E402
                    diffusion_operator, cos_batch, ridge_baseline,
                    emp_p_two_sided, degree_matched_sets,
                    TARGETS, CONTROL_MODULES, MUSCLE_GENES)

# ---------------------------------------------------------------------------
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RAW = os.path.join(ROOT, "data", "raw")
PROC = os.path.join(ROOT, "data", "processed")
TAB = os.path.join(ROOT, "results", "tables")
FIG = os.path.join(ROOT, "results", "figures")
for d in (PROC, TAB, FIG):
    os.makedirs(d, exist_ok=True)

SEED = 20261009
rng = np.random.default_rng(SEED)

NET_HV = 3000
NET_TOPK = 12
ALPHA = 0.55
N_ITER = 80
N_PERM = 1000
STAGES = ["D1", "D7", "D14"]
INJURY = ["Atf3", "Sprr1a", "Gal", "Npy", "Vip", "Nts", "Sst", "Sox11", "Gap43"]


def log(m: str) -> None:
    print(f"[20] {m}", flush=True)


# ---------------------------------------------------------------------------
def load_counts() -> pd.DataFrame:
    """Read GSE149648 per-sample files.

    IMPORTANT — the member files are three-column, NOT two-column:

        Symbol<TAB>gene_id<TAB>A_TCI_D14_1_FPKM
        Mpz<TAB>24564<TAB>5667.92

    Column 1 is the gene symbol, column 2 is the Entrez gene id and column 3
    is the FPKM. An earlier version of this function read column 2 as the
    expression value, which silently produced a matrix of gene identifiers —
    constant across samples — and therefore zero fold changes everywhere.
    The parser now takes column 3 explicitly and asserts that it is not a
    column of identifiers.
    """
    files = sorted(glob.glob(os.path.join(RAW, "GSE149648", "*.txt.gz")))
    if not files:
        files = sorted(glob.glob(os.path.join(RAW, "GSE149648", "*.txt")))
    if not files:
        raise FileNotFoundError("GSE149648 member files absent from data/raw/GSE149648/")
    series, order, skipped, fnames = {}, [], [], {}
    for f in files:
        gsm = os.path.basename(f).split("_")[0]
        op = gzip.open if f.endswith(".gz") else open
        rows, ncol_seen = {}, set()
        try:
            with op(f, "rt") as fh:
                for line in fh:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 2:
                        continue
                    g = parts[0].strip()
                    if not g or g.lower() in ("symbol", "gene", "geneid", "gene_id"):
                        continue                      # header row
                    ncol_seen.add(len(parts))
                    # expression = LAST column; symbol = FIRST column
                    val_tok = parts[-1].strip()
                    try:
                        v = float(val_tok)
                    except ValueError:
                        continue
                    rows[g] = v
        except (EOFError, OSError) as e:
            # A truncated member (incomplete download) must not silently
            # become zeros; skip it and record it loudly.
            skipped.append((os.path.basename(f), type(e).__name__))
            log(f"  WARNING: skipped unreadable member {os.path.basename(f)} ({e})")
            continue
        if not rows:
            raise RuntimeError(f"no numeric expression parsed from {os.path.basename(f)}")
        series[gsm] = pd.Series(rows, dtype=float)
        fnames[gsm] = os.path.basename(f)
        order.append(gsm)
        log(f"    {os.path.basename(f)}: {len(rows)} symbols, "
            f"cols seen={sorted(ncol_seen)}")
    if skipped:
        log(f"  {len(skipped)} member(s) skipped -> n={len(order)}")
    counts = pd.DataFrame(series)[order]
    counts = counts[~counts.index.duplicated(keep="first")].fillna(0.0)
    # Sanity guard: a matrix of Entrez IDs would be integer-valued and nearly
    # invariant across samples; FPKM is neither.
    frac_int = float(np.isclose(counts.values % 1, 0).mean())
    cv = float(counts.replace(0, np.nan).std(axis=1).median() /
               (counts.replace(0, np.nan).mean(axis=1).median() + 1e-9))
    log(f"  matrix: {counts.shape[0]} genes x {counts.shape[1]} samples; "
        f"frac integer-valued={frac_int:.3f}; median CV={cv:.3f}")
    if frac_int > 0.95:
        raise RuntimeError("expression matrix looks like integer identifiers, "
                           "not FPKM — column parsing is wrong")
    return counts, fnames


def group_from_filename(name: str) -> str:
    """The extracted member names carry the experimental group explicitly
    (GSM4506606_A_TCI_D14_1.txt.gz, GSM4506610_A_Sham_1.txt.gz), which is more
    reliable than parsing free-text series-matrix titles.

    The day is extracted as a delimited token. A naive substring test such as
    ``"tci_d1" in name`` wrongly matches "tci_d14", which silently collapsed
    D14 into D1 in an earlier revision.
    """
    n = name.upper()
    if "SHAM" in n:
        return "Sham"
    m = re.search(r"(?:^|[_.\-])D(\d+)(?:[_.\-])", n)
    if m:
        return f"D{int(m.group(1))}"
    return "__unassigned__"


def sample_meta(counts: pd.DataFrame, fnames: dict) -> pd.DataFrame:
    meta = pd.DataFrame({"sample": counts.columns})
    meta["group"] = [group_from_filename(fnames.get(g, "")) for g in counts.columns]
    meta["file"] = [fnames.get(g, "") for g in counts.columns]
    return meta


def normalise(fpkm: pd.DataFrame):
    """GSE149648 supplies FPKM, not counts.

    FPKM is normalised for gene length and depth in principle, but the
    per-sample FPKM totals in this series are NOT equal between groups
    (sham ~0.439e6 vs day 7 ~0.461e6), and the median log fold change of a
    contrast was shifted by +0.09 to +0.12 - i.e. a 6-9% global offset. Left
    uncorrected that offset inflates the number of "differential" genes with
    small fold changes. We therefore apply a median-of-expressed-genes scaling
    (a mild, standard correction that does not assume a count model) before
    the log transform, and keep the raw FPKM total as a QC quantity.
    """
    expressed = fpkm.where(fpkm > 0)
    med = expressed.median(axis=0)
    sf = med / med.mean()
    sf = sf.replace(0, np.nan).fillna(1.0)
    scaled = fpkm.div(sf, axis=1)
    return np.log2(scaled + 1.0), fpkm.sum(axis=0)


# ---------------------------------------------------------------------------
def main() -> None:
    log("=== STAGE 1: load ===")
    counts, fnames = load_counts()
    meta = sample_meta(counts, fnames)
    if (meta["group"] == "__unassigned__").any():
        log("  ERROR: group assignment failed"); sys.exit(2)
    log(f"  groups: {meta['group'].value_counts().to_dict()}")

    log("=== STAGE 2: normalise ===")
    L, fpkm_total = normalise(counts)
    fpkm_all = counts                     # pre-filter, for contamination share
    keep = ((counts >= 1.0).sum(axis=1) >= 4) & (counts.mean(axis=1) >= 0.5)
    L, counts = L.loc[keep], counts.loc[keep]
    L.to_csv(os.path.join(PROC, "expr_log2_fpkm_GSE149648.csv"))
    meta.to_csv(os.path.join(PROC, "sample_meta_GSE149648.csv"), index=False)
    log(f"  retained {L.shape[0]} expressed genes of {counts.shape[0]}")
    log(f"  per-sample FPKM total (x1e6): "
        f"{np.round((fpkm_total/1e6).values, 3)}")

    grp = meta.set_index("sample")["group"]
    cols_by = {g: [s for s in L.columns if grp[s] == g] for g in ["Sham", "D1", "D7", "D14"]}
    log("  " + " ".join(f"{k}={len(v)}" for k, v in cols_by.items()))

    log("=== STAGE 3: confounding diagnostics ===")
    ms = muscle_score(L)
    # Absolute muscle transcript SHARE is more interpretable than a z-score,
    # because a z-score can look unremarkable when the whole group shifts.
    mgenes = [g for g in MUSCLE_GENES if g in fpkm_all.index]
    share = fpkm_all.loc[mgenes].sum(axis=0) / fpkm_all.sum(axis=0) * 100.0
    diag = pd.DataFrame({"sample": L.columns,
                         "group": [grp[s] for s in L.columns],
                         "muscle_score": ms.values,
                         "muscle_fpkm_pct": share.reindex(L.columns).values,
                         "fpkm_total_M": (fpkm_total[L.columns].values / 1e6)})
    # PCA on the expression matrix
    Xc = L.sub(L.mean(axis=1), axis=0).values.T           # samples x genes
    U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
    pc = U * S
    for k in range(4):
        diag[f"PC{k+1}"] = pc[:, k]
    diag["PC1_var_pct"] = np.nan
    var_exp = (S ** 2) / np.sum(S ** 2) * 100.0
    diag.to_csv(os.path.join(TAB, "qc_muscle_scores.csv"), index=False)
    pd.DataFrame({"PC": [f"PC{i+1}" for i in range(4)],
                  "var_pct": var_exp[:4]}).to_csv(
        os.path.join(TAB, "qc_pca_variance.csv"), index=False)
    log(f"  PCA variance explained: {np.round(var_exp[:4], 1)}")
    for g in ["Sham", "D1", "D7", "D14"]:
        v = diag.loc[diag.group == g, "muscle_score"]
        log(f"  muscle score {g}: mean={v.mean():+.3f} sd={v.std():.3f}")
    F, pF = stats.f_oneway(*[diag.loc[diag.group == g, "muscle_score"].values
                             for g in ["Sham", "D1", "D7", "D14"]])
    log(f"  muscle-score ANOVA: F={F:.3f} p={pF:.4f}")

    # Outlier screen. A MAD-based robust z is USELESS here: the muscle share
    # is extremely right-skewed (most dissections carry ~0.01% sarcomere
    # transcript, a few carry several percent), the MAD collapses towards
    # zero, and the resulting z explodes to 1e4 - flagging samples that are
    # in fact the cleanest ones. We therefore (a) work on log10 share, where
    # the distribution is roughly symmetric, and (b) apply a substantive
    # threshold of 1% of total FPKM, above which a DRG sample genuinely
    # contains skeletal muscle rather than a trace of it.
    sh = diag["muscle_fpkm_pct"].values
    lg = np.log10(np.maximum(sh, 1e-4))
    med = np.median(lg)
    mad = np.median(np.abs(lg - med))
    # 1.4826 scales the MAD to a sd under normality; if the MAD is degenerate
    # we fall back to the interquartile range so the statistic stays finite.
    scale = 1.4826 * mad if mad > 1e-9 else \
        (np.percentile(lg, 75) - np.percentile(lg, 25)) / 1.349
    rz = (lg - med) / (scale + 1e-12)
    diag["muscle_robust_z"] = rz
    diag["muscle_flag"] = np.where(sh > 1.0, "MUSCLE_CONTAMINATED", "")
    flags = diag.loc[(np.abs(rz) > 3.5) | (sh > 1.0),
                     ["sample", "group", "muscle_fpkm_pct",
                      "muscle_robust_z", "muscle_flag"]]
    diag.to_csv(os.path.join(TAB, "qc_muscle_scores.csv"), index=False)
    if len(flags):
        log("  OUTLIER FLAG (muscle share):")
        log(flags.to_string(index=False))
    else:
        log("  no sample flagged by muscle-share screen")
    MUSCLE_OUT = list(diag.loc[sh > 1.0, "sample"])
    log(f"  samples above the 1% substantive threshold: {MUSCLE_OUT}")

    C = np.corrcoef(L.values.T)
    pd.DataFrame(C, index=L.columns, columns=L.columns).to_csv(
        os.path.join(TAB, "qc_sample_correlation.csv"))
    off = C[np.triu_indices_from(C, 1)]
    log(f"  sample r range: {off.min():.3f} - {off.max():.3f}")

    log("=== STAGE 4: differential expression (limma-trend) ===")
    res_adj, res_unadj = {}, {}
    for lab in STAGES:
        case, cols = cols_by[lab], cols_by[lab] + cols_by["Sham"]
        Xc = L[cols]
        gvec = np.array([1.0] * len(cols_by[lab]) + [0.0] * len(cols_by["Sham"]))
        mc = ms[cols].values.astype(float); mc -= mc.mean()
        for tag, design in [("adj", np.column_stack([np.ones(len(cols)), gvec, mc])),
                            ("unadj", np.column_stack([np.ones(len(cols)), gvec]))]:
            r = fit_limma(Xc, design, coef_idx=1, trend=True)
            r["q"] = bh(r["P"].values)
            r["AveExpr"] = Xc.mean(axis=1).values
            r = r.sort_values("P")
            r.to_csv(os.path.join(TAB, f"DE_{lab}_vs_Sham_{tag}.csv"))
            (res_adj if tag == "adj" else res_unadj)[lab] = r
        for tag, r in [("adj", res_adj[lab]), ("unadj", res_unadj[lab])]:
            alfc = r["logFC"].abs()
            log(f"  {lab}/{tag}: q<0.05={int((r['q']<0.05).sum())} "
                f"q<0.10={int((r['q']<0.10).sum())} "
                f"P<0.01={int((r['P']<0.01).sum())} | "
                f"median|logFC|={alfc.median():.3f} "
                f"|logFC|>1={int((alfc>1).sum())} "
                f"|logFC|>2={int((alfc>2).sum())} "
                f"| median q={r['q'].median():.3f}")

    # ---- sensitivity: drop muscle-contaminated samples and refit ----
    sens_rows = []
    for lab in STAGES:
        cols = cols_by[lab] + cols_by["Sham"]
        drop = [s for s in MUSCLE_OUT if s in cols]
        if not drop:
            continue
        cols2 = [c for c in cols if c not in drop]
        g2 = np.array([1.0 if c in set(cols_by[lab]) else 0.0 for c in cols2])
        m2 = ms[cols2].values.astype(float)
        m2 -= m2.mean()
        d2 = np.column_stack([np.ones(len(cols2)), g2, m2])
        r2 = fit_limma(L[cols2], d2, coef_idx=1, trend=True)
        r2["q"] = bh(r2["P"].values)
        r1 = res_adj[lab]
        com = r1.index.intersection(r2.index)
        rho, _ = stats.spearmanr(r1.loc[com, "logFC"], r2.loc[com, "logFC"])
        s1 = set(r1.index[r1["q"] < 0.05])
        s2 = set(r2.index[r2["q"] < 0.05])
        sens_rows.append({"contrast": lab, "dropped": ";".join(drop),
                          "n_case_left": len([c for c in cols2 if c in cols_by[lab]]),
                          "n_q05_full": len(s1), "n_q05_reduced": len(s2),
                          "jaccard_q05": len(s1 & s2) / max(len(s1 | s2), 1),
                          "spearman_logFC": float(rho)})
    if sens_rows:
        sens = pd.DataFrame(sens_rows)
        sens.to_csv(os.path.join(TAB, "de_muscle_drop_sensitivity.csv"), index=False)
        log("  muscle-drop sensitivity:")
        log(sens.to_string(index=False))

    log("=== STAGE 5: positive-control markers ===")
    pc = []
    for g in INJURY + ["Cxcl13", "Il1b", "Tnf", "Ccl2", "Gfap", "Aif1"]:
        for lab in STAGES:
            if g in res_adj[lab].index:
                r = res_adj[lab].loc[g]
                pc.append({"gene": g, "contrast": lab, "logFC": round(r["logFC"], 3),
                           "P": float(f"{r['P']:.3g}"), "q": round(r["q"], 4)})
    pc = pd.DataFrame(pc)
    pc.to_csv(os.path.join(TAB, "positive_control_markers.csv"), index=False)
    log(pc.to_string(index=False))

    log("=== STAGE 6: module-level signal (targets + negative controls) ===")
    ALL_MODS = dict(TARGETS)
    ALL_MODS.update(CONTROL_MODULES)
    rows = []
    for lab in STAGES:
        r = res_adj[lab]
        bg = r["t"].values
        for mod, gs in ALL_MODS.items():
            present = [g for g in gs if g in r.index]
            if not present:
                rows.append({"contrast": lab, "module": mod, "kind": "control"
                             if mod.startswith("CTRL") else "target",
                             "n": 0, "mean_t": np.nan, "mean_logFC": np.nan,
                             "p_comp": np.nan, "genes_used": ""})
                continue
            sub = r.loc[present]
            _, pu = stats.mannwhitneyu(sub["t"].values, bg, alternative="two-sided")
            rows.append({"contrast": lab, "module": mod, "kind": "control"
                         if mod.startswith("CTRL") else "target",
                         "n": len(sub), "mean_t": round(sub["t"].mean(), 3),
                         "mean_logFC": round(sub["logFC"].mean(), 3),
                         "p_comp": pu, "genes_used": ";".join(present)})
    mod = pd.DataFrame(rows)
    ok = mod["p_comp"].notna()
    mod.loc[ok, "q_comp"] = bh(mod.loc[ok, "p_comp"].values)
    mod = mod.sort_values(["contrast", "p_comp"])
    mod.to_csv(os.path.join(TAB, "module_level_signal.csv"), index=False)
    log(mod[["contrast", "module", "n", "mean_t", "mean_logFC",
             "p_comp", "q_comp"]].to_string(index=False))
    log("\n  --- negative-control modules (do non-target sets move as much?) ---")
    ctrl = mod[mod.kind == "control"][["contrast", "module", "n", "mean_t",
                                       "mean_logFC", "q_comp"]]
    log(ctrl.to_string(index=False))

    log("=== STAGE 7: virtual perturbation ===")
    # Every gene of every module is forced onto the graph. Selecting nodes by
    # variance alone drops most analgesic target genes (ion channels and
    # opioid receptors are moderately expressed with unremarkable
    # across-sample variance), which would leave modules such as GABAergic
    # represented by a single node - or dropped entirely.
    MUST = sorted({g for gs in ALL_MODS.values() for g in gs} & set(L.index))
    genes, W = build_network(L, hv=NET_HV, topk=NET_TOPK, must_include=MUST)
    log(f"  network: {len(genes)} genes (top-{NET_HV} HV + {len(MUST)} forced "
        f"module genes), topk={NET_TOPK}, mean|W|={np.abs(W).mean():.4f}")
    n_missing = {m: len([g for g in gs if g in genes])
                 for m, gs in ALL_MODS.items()}
    log("  genes retained per module: " +
        ", ".join(f"{m}={v}/{len(ALL_MODS[m])}" for m, v in n_missing.items()))
    P = diffusion_operator(W, alpha=ALPHA, n_iter=N_ITER)
    log(f"  diffusion operator built; ||P||_F={np.linalg.norm(P):.2f}")

    ALL_MODS = dict(TARGETS)
    ALL_MODS.update(CONTROL_MODULES)
    pert = []
    for lab in STAGES:
        cols = cols_by[lab] + cols_by["Sham"]
        Lc = L[cols]                      # ridge baseline must see only this
        pain_axis = (L[cols_by["Sham"]].mean(axis=1)
                     - L[cols_by[lab]].mean(axis=1)).reindex(genes).values
        pa_c = pain_axis - pain_axis.mean()
        for modname, tlist in ALL_MODS.items():
            present = [g for g in tlist if g in genes]
            if not present:
                continue
            idx = np.array([genes.index(g) for g in present])
            # ORACLE ANCHOR: the best achievable reversal for a seed of this
            # size, obtained by perturbing the k genes that already lie
            # furthest along the pain axis with the correct sign. Without it a
            # cosine of 0.05 has no interpretable scale.
            k = len(present)
            oi = np.argsort(-np.abs(pa_c))[:k]
            oseed = np.zeros(len(genes))
            oseed[oi] = np.sign(pa_c[oi])
            rev_oracle = float(cos_batch((P @ oseed).reshape(-1, 1), pain_axis)[0])
            for ptype, sign in [("knockout", -1.0), ("overexpression", 1.0)]:
                seed = np.zeros(len(genes)); seed[idx] = sign
                d_net = P @ seed
                s_net = float(cos_batch(d_net.reshape(-1, 1), pain_axis)[0])
                s_add = float(cos_batch(seed.reshape(-1, 1), pain_axis)[0])
                s_ridge = ridge_baseline(Lc, genes, present, sign, pain_axis,
                                         loso=True)
                # Degree-matched null. Uniform random gene sets are NOT a
                # valid null here: hub genes propagate further on the graph, so
                # a hub-rich module looks "significant" against a uniform null
                # for a purely topological reason.
                msets = degree_matched_sets(W, idx, N_PERM, rng, n_bins=10)
                Ind = np.zeros((len(genes), N_PERM))
                for kk in range(N_PERM):
                    Ind[msets[kk], kk] = sign
                sp = cos_batch(P @ Ind, pain_axis)
                p_emp = emp_p_two_sided(s_net, sp)
                # uniform (unmatched) null retained for reporting only
                IndU = np.zeros((len(genes), N_PERM))
                for kk in range(N_PERM):
                    IndU[rng.choice(len(genes), len(present), replace=False), kk] = sign
                spu = cos_batch(P @ IndU, pain_axis)
                p_unmatched = emp_p_two_sided(s_net, spu)
                pert.append({"contrast": lab, "module": modname,
                             "kind": "control" if modname.startswith("CTRL") else "target",
                             "n_targets": len(present),
                             "genes_used": ";".join(present),
                             "perturbation": ptype, "reversal_network": s_net,
                             "reversal_additive": s_add, "reversal_ridge": s_ridge,
                             "reversal_oracle": rev_oracle,
                             "pct_of_oracle": 100.0 * abs(s_net) / (abs(rev_oracle) + 1e-12),
                             "cos2": s_net ** 2,
                             "delta_norm": float(np.linalg.norm(d_net - d_net.mean())),
                             "perm_mean": float(sp.mean()), "perm_sd": float(sp.std()),
                             "z_vs_perm": (s_net - sp.mean()) / (sp.std() + 1e-12),
                             "p_empirical": p_emp,
                             "p_empirical_unmatched": p_unmatched})
    # ---- leave-one-sample-out stability ----
    # A single animal drives a reversal score at n=4 per group; the range
    # across the 8 leave-one-out refits says how much of the score is carried
    # by one sample. delta = P @ seed does not depend on the axis, so only
    # the axis is recomputed - 8 fits per module, not 8 diffusions.
    loo = {}
    for lab in STAGES:
        case_c, ctrl_c = cols_by[lab], cols_by["Sham"]
        for modname, tlist in ALL_MODS.items():
            present = [g for g in tlist if g in genes]
            if not present:
                continue
            idx = np.array([genes.index(g) for g in present])
            seed = np.zeros(len(genes))
            seed[idx] = -1.0
            d = P @ seed
            vals = []
            for s in case_c + ctrl_c:
                ca = [c for c in case_c if c != s]
                co = [c for c in ctrl_c if c != s]
                if len(ca) < 2 or len(co) < 2:
                    continue
                ax = (L[co].mean(axis=1) - L[ca].mean(axis=1)).reindex(genes).values
                vals.append(float(cos_batch(d.reshape(-1, 1), ax)[0]))
            loo[(lab, modname)] = (float(np.min(vals)), float(np.max(vals)),
                                   float(np.std(vals)))
    for r in pert:
        k = (r["contrast"], r["module"])
        if k in loo:
            r["reversal_loo_min"], r["reversal_loo_max"], r["reversal_loo_sd"] = loo[k]

    pert = pd.DataFrame(pert)
    qparts = [pd.Series(bh(pert.loc[pert.contrast == lab, "p_empirical"].values),
                        index=pert.index[pert.contrast == lab]) for lab in STAGES]
    pert["q_empirical"] = pd.concat(qparts)
    pert = pert.sort_values(["contrast", "perturbation", "p_empirical"])
    pert.to_csv(os.path.join(TAB, "virtual_perturbation_scores.csv"), index=False)
    log("\n  --- top reversals by two-sided empirical p (D7) ---")
    show = pert[(pert.contrast == "D7")].sort_values("p_empirical").head(20)
    log(show[["module", "kind", "n_targets", "perturbation", "reversal_network",
              "z_vs_perm", "pct_of_oracle", "p_empirical",
              "q_empirical"]].to_string(index=False))

    # Is the network doing anything the additive model does not? If the two
    # rankings are near-identical the graph is decorative, which is exactly
    # what the 2025 Nat Methods benchmark found for far more complex models.
    rho_na, _ = stats.spearmanr(pert["reversal_network"], pert["reversal_additive"])
    log(f"\n  Spearman(network reversal, additive reversal) = {rho_na:.3f} "
        f"(n={len(pert)}) -- if ~1 the graph adds no ranking information")
    with open(os.path.join(TAB, "head_to_head_baselines.json"), "w") as fh:
        json.dump({"spearman_network_vs_additive": float(rho_na),
                   "n_conditions": int(len(pert)),
                   "mean_abs_network": float(pert["reversal_network"].abs().mean()),
                   "mean_abs_additive": float(pert["reversal_additive"].abs().mean()),
                   "mean_abs_ridge": float(pert["reversal_ridge"].abs().mean()),
                   "mean_pct_of_oracle": float(pert["pct_of_oracle"].mean())}, fh, indent=2)

    log("=== STAGE 8: figures ===")
    make_figures(L, grp, diag, res_adj, mod, pert)

    summary = {
        "dataset": "GSE149648", "platform": "GPL24782",
        "n_samples": int(L.shape[1]), "n_genes": int(L.shape[0]),
        "groups": {g: len(v) for g, v in cols_by.items()},
        "muscle_ANOVA": {"F": round(float(F), 3), "p": float(pF)},
        "n_sig_q05_adj": {k: int((v["q"] < 0.05).sum()) for k, v in res_adj.items()},
        "n_sig_q10_adj": {k: int((v["q"] < 0.10).sum()) for k, v in res_adj.items()},
        "n_sig_q05_unadj": {k: int((v["q"] < 0.05).sum()) for k, v in res_unadj.items()},
        "network": {"genes": len(genes), "topk": NET_TOPK, "alpha": ALPHA},
        "n_perm": N_PERM, "seed": SEED,
    }
    with open(os.path.join(TAB, "run_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    log(json.dumps(summary, indent=2))
    log("DONE")


def make_figures(L, grp, diag, results, mod, pert) -> None:
    order = ["Sham", "D1", "D7", "D14"]
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
    data = [diag.loc[diag.group == g, "muscle_score"].values for g in order]
    ax[0].boxplot(data, tick_labels=order)
    for i, d in enumerate(data):
        ax[0].scatter(np.full(len(d), i + 1), d, s=20, zorder=3)
    ax[0].set_ylabel("Skeletal-muscle contamination score (z)")
    ax[0].set_title("A. Muscle contamination by group")
    C = np.corrcoef(L.values.T)
    im = ax[1].imshow(C, cmap="viridis", vmin=0.7, vmax=1.0)
    ax[1].set_xticks(range(len(L.columns)))
    ax[1].set_xticklabels([grp[s] for s in L.columns], rotation=90, fontsize=7)
    ax[1].set_yticks(range(len(L.columns)))
    ax[1].set_yticklabels([grp[s] for s in L.columns], fontsize=7)
    ax[1].set_title("B. Sample-sample Pearson r")
    plt.colorbar(im, ax=ax[1], fraction=0.046)
    plt.tight_layout(); plt.savefig(os.path.join(FIG, "Fig1_QC.png"), dpi=200); plt.close()

    fig, ax = plt.subplots(1, 3, figsize=(14, 4.4))
    for i, lab in enumerate(STAGES):
        r = results[lab]
        x, y = r["logFC"].values, -np.log10(np.clip(r["P"].values, 1e-300, 1))
        sig = r["q"].values < 0.05
        ax[i].scatter(x[~sig], y[~sig], s=5, c="#bdbdbd", alpha=0.5)
        ax[i].scatter(x[sig], y[sig], s=9, c="#d62728")
        ax[i].axhline(-np.log10(0.05), ls="--", c="k", lw=0.7)
        ax[i].set_xlabel("log2 fold change"); ax[i].set_ylabel("-log10 P")
        ax[i].set_title(f"{lab} vs Sham (q<0.05: {int(sig.sum())})")
        for g2, row in r.head(6).iterrows():
            ax[i].annotate(g2, (row["logFC"], -np.log10(max(row["P"], 1e-300))), fontsize=6.5)
    plt.tight_layout(); plt.savefig(os.path.join(FIG, "Fig2_DE_volcano.png"), dpi=200); plt.close()

    piv = mod.pivot(index="module", columns="contrast", values="mean_t")[STAGES].dropna()
    fig, ax = plt.subplots(figsize=(8.5, 7))
    im = ax.imshow(piv.values, cmap="RdBu_r", vmin=-4, vmax=4, aspect="auto")
    ax.set_xticks(range(3)); ax.set_xticklabels(STAGES)
    ax.set_yticks(range(len(piv))); ax.set_yticklabels(piv.index, fontsize=8)
    plt.colorbar(im, ax=ax, fraction=0.046, label="mean moderated t")
    ax.set_title("Module-level differential signal across cancer-pain stages")
    for i in range(len(piv)):
        for j in range(3):
            ax.text(j, i, f"{piv.values[i, j]:.1f}", ha="center", va="center", fontsize=6.5)
    plt.tight_layout(); plt.savefig(os.path.join(FIG, "Fig3_module_heatmap.png"), dpi=200); plt.close()

    ko = pert[pert.perturbation == "knockout"]
    if not ko.empty:
        piv = ko.pivot(index="module", columns="contrast", values="reversal_network")[STAGES].dropna()
        fig, ax = plt.subplots(figsize=(8.5, 7))
        v = max(0.2, float(np.nanmax(np.abs(piv.values))))
        im = ax.imshow(piv.values, cmap="RdYlBu", vmin=-v, vmax=v, aspect="auto")
        ax.set_xticks(range(3)); ax.set_xticklabels(STAGES)
        ax.set_yticks(range(len(piv))); ax.set_yticklabels(piv.index, fontsize=8)
        plt.colorbar(im, ax=ax, fraction=0.046, label="reversal score (cosine)")
        ax.set_title("In-silico knockout: reversal of the cancer-pain state")
        for i in range(len(piv)):
            for j in range(3):
                ax.text(j, i, f"{piv.values[i, j]:+.2f}", ha="center", va="center", fontsize=6.5)
        plt.tight_layout(); plt.savefig(os.path.join(FIG, "Fig4_virtual_ko.png"), dpi=200); plt.close()

    sub = ko[ko.contrast == "D14"].sort_values("reversal_network")
    if not sub.empty:
        fig, ax = plt.subplots(figsize=(9.5, 4.6))
        x = np.arange(len(sub))
        ax.bar(x - 0.24, sub["reversal_network"], width=0.23, label="Network diffusion")
        ax.bar(x, sub["reversal_additive"], width=0.23, label="Additive (no propagation)")
        ax.bar(x + 0.24, sub["reversal_ridge"], width=0.23, label="Ridge-linear")
        ax.set_xticks(x); ax.set_xticklabels(sub["module"], rotation=90, fontsize=7)
        ax.axhline(0, c="k", lw=0.7); ax.set_ylabel("reversal score")
        ax.set_title("D14: network model vs baselines (head-to-head)")
        ax.legend(fontsize=8)
        plt.tight_layout(); plt.savefig(os.path.join(FIG, "Fig5_baseline_headtohead.png"), dpi=200); plt.close()


if __name__ == "__main__":
    main()
