#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
_limma.py — shared statistical / perturbation primitives.

Single source of truth so that the discovery pipeline (20) and the
validation pipeline (21) execute byte-identical code. Duplicating these
routines across scripts is exactly how version drift between a manuscript's
numbers and its code is introduced.

Exports
-------
fit_limma(X, design, coef_idx, trend)  moderated t (Smyth 2004) + intensity trend
bh(p)                                  Benjamini-Hochberg q values
muscle_score(L)                        skeletal-muscle contamination score
build_network(L, hv, topk)             k-NN co-expression graph
diffusion_operator(W, alpha, n_iter)   closed-ish form of the diffusion map
cos_batch(D, b)                        column-wise cosine similarity
ridge_baseline(...)                    ridge-linear baseline perturbation
TARGETS                                analgesic target modules
"""
from __future__ import annotations

import math
import numpy as np
import pandas as pd
from scipy import stats
from scipy.special import digamma, polygamma

# statsmodels is imported at top level ON PURPOSE. A previous revision wrapped
# it in `try: import ... except Exception:` with a silent fallback to a
# constant prior variance. That fallback is not a benign degradation: it
# silently turns limma-trend into limma-without-trend, changes every moderated
# t statistic in the study, and leaves no trace in the output. A missing
# dependency must fail loudly.
try:
    import statsmodels.api as sm
except ImportError as _exc:                                   # pragma: no cover
    raise SystemExit(
        "statsmodels is required (limma-trend lowess prior variance).\n"
        "  pip install -r requirements.txt\n"
        f"  original error: {_exc}"
    ) from _exc

# ---------------------------------------------------------------------------
# Empirical-Bayes moderated t
# ---------------------------------------------------------------------------
def _solve_prior_df(var_e: float) -> float:
    """Solve polygamma(1, df/2) = var_e for df (monotone decreasing)."""
    f = lambda x: polygamma(1, x / 2.0) - var_e
    lo, hi = 1e-4, 1e5
    for _ in range(200):
        mid = math.sqrt(lo * hi)
        if f(mid) > 0:
            lo = mid          # polygamma too large -> need larger df
        else:
            hi = mid
    return math.sqrt(lo * hi)


def squeeze_var(s2: np.ndarray, df_res: float):
    """Method-of-moments prior variance (Smyth 2004, Appendix)."""
    ok = np.isfinite(s2) & (s2 > 0)
    z = np.log(s2[ok])
    e = z - digamma(df_res / 2.0) + math.log(df_res / 2.0)
    ebar = e.mean()
    n = len(e)
    var_e = ((e - ebar) ** 2).sum() / max(n - 1, 1) - polygamma(1, df_res / 2.0)
    if var_e <= 0:
        return np.full_like(s2, math.exp(ebar)), float("inf")
    df_prior = _solve_prior_df(var_e)
    s20 = math.exp(ebar + digamma(df_prior / 2.0) - math.log(df_prior / 2.0))
    return (df_prior * s20 + df_res * s2) / (df_prior + df_res), df_prior


def fit_limma(X: pd.DataFrame, design: np.ndarray, coef_idx: int, trend: bool = True,
              frac: float = 0.4):
    """Moderated t-test for one coefficient of a linear model.

    trend=True applies limma-trend: a lowess fit of log(s2) on mean intensity
    gives a gene-specific prior variance, which matters a lot for RNA-seq
    where the mean-variance relationship is strong.

    `frac` is the lowess span, exposed as an argument so that the sensitivity
    analysis (23_sensitivity.py) can vary it without editing this file.
    """
    Y = np.asarray(X.values, dtype=float)            # genes x samples
    n, p = design.shape
    df_res = n - p
    if df_res < 1:
        raise ValueError("design leaves no residual degrees of freedom")

    beta, *_ = np.linalg.lstsq(design, Y.T, rcond=None)
    resid = Y.T - design @ beta
    s2 = (resid ** 2).sum(axis=0) / df_res
    XtX_inv = np.linalg.pinv(design.T @ design)
    unscaled = np.sqrt(XtX_inv[coef_idx, coef_idx])
    lfc = beta[coef_idx]

    if trend:
        Amean = Y.mean(axis=1)
        ok = np.isfinite(s2) & (s2 > 0) & np.isfinite(Amean)
        s20_trend = np.full_like(s2, np.nan)
        s2_scaled = np.full_like(s2, np.nan)
        base = np.log(s2[ok])
        fitv = sm.nonparametric.lowess(base, Amean[ok], frac=frac, return_sorted=False)
        s20_trend[ok] = np.exp(fitv)
        s2_scaled[ok] = np.exp(base - fitv)
        ok2 = np.isfinite(s2_scaled) & (s2_scaled > 0)
        s2p, df_prior = squeeze_var(s2_scaled[ok2], df_res)
        s2_post = np.full_like(s2, np.nan)
        s2_post[ok2] = s2p * s20_trend[ok2]
        fallback = float(np.nanmedian(s2_post))
        s2_post = np.where(np.isfinite(s2_post), s2_post, fallback)
    else:
        s2_post, df_prior = squeeze_var(s2, df_res)

    se = np.sqrt(np.maximum(s2_post, 1e-12)) * unscaled
    t = lfc / np.maximum(se, 1e-12)
    df_total = df_res + (df_prior if np.isfinite(df_prior) else 1e6)
    pv = 2.0 * stats.t.sf(np.abs(t), df_total)
    pv = np.where(np.isfinite(pv), pv, 1.0)
    return pd.DataFrame({"logFC": lfc, "t": t, "P": pv, "s2_post": s2_post},
                        index=X.index)


def bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    n = len(p)
    order = np.argsort(p)
    q = p[order] * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(q, 0, 1)
    return out


# ---------------------------------------------------------------------------
# Confounding: skeletal-muscle contamination of DRG dissections
# ---------------------------------------------------------------------------
MUSCLE_GENES = ["Acta1", "Ckm", "Myl1", "Myh1", "Myh2", "Tnni2", "Tnnc2", "Tnnt3",
                "Des", "Mybph", "Actn3", "Mybpc1", "Myoz1", "Tnnc1", "Mylpf",
                "Myh4", "Atp2a1", "Casq1", "Mb", "Actc1"]


def muscle_score(L: pd.DataFrame) -> pd.Series:
    present = [g for g in MUSCLE_GENES if g in L.index]
    if not present:
        return pd.Series(0.0, index=L.columns)
    Z = L.loc[present]
    Z = Z.sub(Z.mean(axis=1), axis=0).div(Z.std(axis=1).replace(0, np.nan), axis=0)
    return Z.mean(axis=0).fillna(0.0)


# ---------------------------------------------------------------------------
# Analgesic target modules (rat gene symbols)
# ---------------------------------------------------------------------------
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

# ---------------------------------------------------------------------------
# Negative-control modules.
#
# These are deliberately NOT analgesic targets. They exist to test whether the
# "neuronal modules are down-regulated" result is specific to analgesic
# mechanisms or is instead a compartment effect: if a set of generic neuronal
# structural genes falls by as much as the Na+/K+/GABA/opioid modules do, then
# those modules are not five independent mechanistic findings - they are five
# probes of one underlying shift in neuronal transcript share.
# ---------------------------------------------------------------------------
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


# ---------------------------------------------------------------------------
# Virtual perturbation machinery
# ---------------------------------------------------------------------------
def build_network(L: pd.DataFrame, hv: int = 3000, topk: int = 12,
                  must_include=None, genes=None):
    """Sparse symmetrised k-NN graph over the highest-variance genes.

    `genes` pins the node set explicitly (used by the cross-aetiology
    comparison, where every contrast must be propagated over the SAME node
    set or the intersection of four independent top-variance selections
    collapses to a few hundred genes).

    The diagonal is set to ZERO, not to -1. An earlier revision used -1 to
    "exclude self", but the selection is on |C| and |-1| = 1 is the row
    maximum, so the self-pair was always selected and the value written was
    the signed -1. That produced a graph with a negative self-loop at every
    node: the normalised operator then became a sign-flipping walk rather
    than a diffusion, and every perturbation result downstream was wrong.

    must_include: gene symbols that MUST appear as nodes even if they are not
    high-variance. This is not cosmetic. Analgesic target genes (Scn9a,
    Oprm1, Cacna2d1, ...) are moderately expressed and their across-sample
    variance is unremarkable, so a pure top-variance selection silently drops
    most of them - the modules then get perturbed on a graph that does not
    contain the genes being perturbed, and the modules with few surviving
    genes are dropped from the results altogether. Any perturbation study that
    selects its graph by variance alone and then perturbs a predefined gene
    set has this bug.
    """
    if genes is not None:
        genes = [g for g in genes if g in L.index]
    else:
        var = L.var(axis=1).sort_values(ascending=False)
        genes = list(var.index[:hv])
        if must_include:
            extra = [g for g in must_include if g in L.index and g not in set(genes)]
            genes = genes + list(dict.fromkeys(extra))
    M = L.loc[genes].values
    C = np.nan_to_num(np.corrcoef(M))
    np.fill_diagonal(C, 0.0)                 # smallest |C| -> never selected
    idx = np.argpartition(-np.abs(C), topk, axis=1)[:, :topk]
    W = np.zeros_like(C)
    rows = np.arange(len(genes))[:, None]
    W[rows, idx] = C[rows, idx]
    W = (W + W.T) / 2.0
    np.fill_diagonal(W, 0.0)                 # defensive: no self-loops
    return genes, W


def diffusion_operator(W: np.ndarray, alpha: float = 0.55, n_iter: int = 80,
                       exact: bool = False) -> np.ndarray:
    """P with delta_inf = P @ seed for
       delta_{t+1} = (1-a) seed + a S delta_t,  S = D^-1/2 W D^-1/2.

    Because the map is linear, P is computed once and reused for every real
    and permuted perturbation; |a*S| <= a < 1 guarantees convergence.

    exact=True evaluates the closed form P = (1-a)(I - a S)^-1 by a single
    linear solve (O(n^3) once) instead of n_iter dense matmuls (O(n_iter n^3)).
    I - a S has condition number <= (1+a)/(1-a) = 3.4 at a = 0.55, so the
    solve is numerically safe. 99_unit_tests.py asserts the two routes agree
    to 1e-10, and the default remains the iteration so that the deposited
    results are produced by the same code path as in earlier revisions.
    """
    d = np.abs(W).sum(axis=1)
    d[d == 0] = 1.0
    Dm = np.diag(1.0 / np.sqrt(d))
    S = Dm @ W @ Dm
    n = W.shape[0]
    if exact:
        return (1.0 - alpha) * np.linalg.solve(np.eye(n) - alpha * S, np.eye(n))
    P = np.eye(n) * (1.0 - alpha)
    I = np.eye(n) * (1.0 - alpha)
    for _ in range(n_iter):
        P = I + alpha * (S @ P)
    return P


def emp_p_two_sided(obs: float, null: np.ndarray) -> float:
    """Two-sided empirical p-value against a permutation null.

    WHY TWO-SIDED: an earlier revision chose the tail from the perturbation
    label (knockout -> lower tail, over-expression -> upper tail). That is
    wrong, because it hard-codes the assumption that knocking a gene out must
    produce a negative reversal. The scientifically meaningful event is
    "reversal" (positive score) or "anti-reversal" (negative), and the sign
    must be read off the estimate, not assumed in advance. With the
    label-dependent tail, the single strongest reversal in this study
    (glial/immune knockout, score +0.136, z = +4.31) was assigned p = 1.000
    and silently excluded, while anti-reversals were reported as findings.

    A two-sided p also gives the knock-out and over-expression arms identical
    significance, which is correct: they are exact mirror images
    (delta(-seed) = -delta(seed) and cosine is odd in its first argument), so
    they carry one piece of information, not two.
    """
    null = np.asarray(null, dtype=float)
    n = len(null)
    lo = (np.sum(null <= obs) + 1) / (n + 1)
    hi = (np.sum(null >= obs) + 1) / (n + 1)
    return float(min(1.0, 2.0 * min(lo, hi)))


def cos_batch(D: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Cosine similarity of each column of D with the reference vector b.

    Both are mean-centred first, which removes the global shift component and
    leaves the *pattern* agreement - the quantity of interest for a reversal.

    `b` may be a 1-D vector or a 2-D (genes x m) matrix, in which case each
    column is centred and normalised BY ITS OWN mean and norm and the result
    is (D.shape[1], m).

    BUG THIS FIXES: the original implementation always did `b - b.mean()` and
    `norm(bc)`, which for a 2-D b subtracts the GRAND mean over all entries
    and divides every column by the Frobenius norm of the whole matrix. Every
    column then came out scaled by a different wrong factor, so the
    cross-aetiology Delta (computed against a 70-column matrix of permutation
    axes) was wrong by up to 5x while the point estimates (computed against a
    1-D axis) were right - an inconsistency that is invisible unless the two
    are compared directly.
    """
    D = np.atleast_2d(np.asarray(D, dtype=float))
    if D.shape[0] == 1 and np.ndim(b) > 1 and b.shape[1] != 1:
        D = D.T
    Dc = D - D.mean(axis=0, keepdims=True)
    na = np.linalg.norm(Dc, axis=0)
    b = np.asarray(b, dtype=float)
    if b.ndim == 1:
        bc = b - b.mean()
        nb = np.linalg.norm(bc)
        if nb < 1e-12:
            return np.zeros(D.shape[1])
        return (Dc.T @ bc) / (na * nb + 1e-12)
    bc = b - b.mean(axis=0, keepdims=True)
    nb = np.linalg.norm(bc, axis=0, keepdims=True)
    out = (Dc.T @ bc) / (na[:, None] * nb + 1e-12)
    return out


# ---------------------------------------------------------------------------
# Null distributions for the perturbation score
# ---------------------------------------------------------------------------
def degree_bins(W: np.ndarray, n_bins: int = 10) -> np.ndarray:
    """Bin network nodes into `n_bins` quantile bins of weighted degree.

    WHY THIS EXISTS: a uniformly random gene set is not a valid null for a
    network-diffusion score. Hub genes sit in dense neighbourhoods, so a
    perturbation seeded on them propagates much further than one seeded on
    peripheral genes; the variance of the null is therefore driven by degree,
    not by module identity. Comparing a hub-rich target module against a
    uniform null manufactures significance. Degree matching removes that.
    """
    deg = np.abs(W).sum(axis=1)
    qs = np.unique(np.quantile(deg, np.linspace(0.0, 1.0, n_bins + 1)))
    if len(qs) < 2:
        return np.zeros(W.shape[0], dtype=int)
    return np.digitize(deg, qs[1:-1], right=False)


def degree_matched_sets(W: np.ndarray, idx: np.ndarray, n_perm: int, rng,
                        n_bins: int = 10) -> np.ndarray:
    """`n_perm` random gene index vectors, each matched to `idx` bin-by-bin.

    Returns an (n_perm, len(idx)) integer array of node indices.
    """
    b = degree_bins(W, n_bins=n_bins)
    nb = int(b.max()) + 1
    pools = [np.flatnonzero(b == v) for v in range(nb)]
    pools = [p if len(p) else np.arange(len(b)) for p in pools]
    mod_bins = b[np.asarray(idx, dtype=int)]
    out = np.empty((n_perm, len(mod_bins)), dtype=int)
    for i in range(n_perm):
        used = set()
        for j, bb in enumerate(mod_bins):
            pool = pools[bb]
            pick = int(rng.choice(pool))
            for _ in range(64):
                if pick not in used:
                    break
                pick = int(rng.choice(pool))
            used.add(pick)
            out[i, j] = pick
    return out


def all_splits(n_case: int, n_ctrl: int) -> np.ndarray:
    """Every case/control assignment of n_case + n_ctrl samples, as a boolean
    matrix of shape (n_choose_k, n_case + n_ctrl). True marks the 'case'
    member of the split.

    At n = 4 vs 4 this is C(8,4) = 70 splits, i.e. the entire label
    permutation space. That means the label-permutation test can be computed
    EXACTLY rather than by Monte Carlo — no resampling error, no seed
    dependence, and the minimum attainable p is 1/70, not 0.
    """
    from itertools import combinations
    n = n_case + n_ctrl
    combs = list(combinations(range(n), n_case))
    M = np.zeros((len(combs), n), dtype=bool)
    for r, c in enumerate(combs):
        M[r, list(c)] = True
    return M


def perm_axes(Lv: np.ndarray, splits: np.ndarray) -> np.ndarray:
    """Pain axis (mean control - mean case) for every split in `splits`.

    Lv : genes x samples matrix (already restricted to the samples of one
         contrast).  Returns genes x n_splits.
    """
    ns = Lv.shape[1]
    out = np.empty((Lv.shape[0], splits.shape[0]))
    for r in range(splits.shape[0]):
        ctrl = np.flatnonzero(~splits[r])
        case = np.flatnonzero(splits[r])
        out[:, r] = Lv[:, ctrl].mean(axis=1) - Lv[:, case].mean(axis=1)
    return out


def exact_two_sided_p(obs: float, null: np.ndarray) -> float:
    """Exact permutation p from an exhaustively enumerated null (with the
    observed value itself included, as it is one of the enumerations)."""
    null = np.asarray(null, dtype=float)
    n = len(null)
    ge = int(np.sum(np.abs(null) >= abs(obs) - 1e-12))
    return float(min(1.0, (ge + 1) / (n + 1)))


def ridge_baseline(L: pd.DataFrame, genes: list, targets: list, sign: float,
                   pain_axis: np.ndarray, l2: float = 0.1, loso: bool = True) -> float:
    """Ridge-linear baseline: regress every gene on the target genes across
    samples, then predict the response to shifting the targets by `sign`.
    No network, no propagation.

    loso=True (default) averages the regression coefficients over
    leave-one-sample-out fits. This matters: an in-sample fit on the same
    samples that define the pain axis can reach |cos| ~ 0.9 purely by
    over-fitting a handful of predictors in a small-n design, which would make
    the baseline look far stronger than the network model for a reason that
    has nothing to do with biology. LOSO keeps the baseline honest.
    """
    tidx = [genes.index(t) for t in targets if t in genes]
    if not tidx:
        return 0.0
    X = L.loc[genes].values                   # genes x samples
    # Target genes index the GENE axis, so select rows and transpose to get
    # the samples x targets design matrix. (Selecting columns here would index
    # samples, which is the bug this line previously contained.)
    T = X[tidx, :].T                          # samples x n_targets
    ns = T.shape[0]
    sd = T.std(axis=0)
    sd[sd < 1e-9] = 1.0
    Tz = (T - T.mean(axis=0)) / sd
    Yc = X - X.mean(axis=1, keepdims=True)    # genes x samples, gene-centred
    shift = np.full(len(tidx), sign)

    if loso and ns > 2:
        preds = np.zeros((X.shape[0], ns))
        for i in range(ns):
            tr = [j for j in range(ns) if j != i]
            A = Tz[tr]
            R = np.linalg.pinv(A.T @ A + l2 * np.eye(len(tidx))) @ A.T @ Yc[:, tr].T
            preds[:, i] = R.T @ shift
        pred = preds.mean(axis=1)
    else:
        R = np.linalg.pinv(Tz.T @ Tz + l2 * np.eye(len(tidx))) @ Tz.T @ Yc.T
        pred = R.T @ shift
    return float(cos_batch(pred.reshape(-1, 1), pain_axis)[0])
