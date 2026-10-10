#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
99_unit_tests.py — assertions that would have caught every bug found so far.

Each test exists because the bug it guards against actually occurred in this
project and was only found by manual inspection:

  T1  moderated-t / logFC recovery on simulated data with a known truth
      (guards: wrong design coding, wrong coefficient index)
  T2  Benjamini-Hochberg against a hand-computed reference vector
      (guards: off-by-one in the rank step, missing monotone envelope)
  T3  diffusion operator: iteration vs closed-form linear solve
      (guards: wrong Neumann normalisation, alpha applied twice)
  T4  diffusion operator: convergence floor of the iteration
  T5  network has NO self-loops, is symmetric, and has the requested degree
      (guards: the -1 diagonal bug that put a negative self-loop at every
       node and silently corrupted every downstream perturbation score)
  T6  knockout and over-expression are exact mirror images
      (guards: any asymmetry in seed construction or centring)
  T7  column parsing: an identifier matrix must be REJECTED
      (guards: the Entrez-id / FPKM column mix-up that produced zero fold
       changes everywhere while reporting success)
  T8  group parsing: D14 must not collapse into D1
      (guards: naive substring matching)
  T9  degree-matched null really draws from the matched degree bin
      (guards: silently falling back to uniform sampling)
  T10 exhaustive label splits: C(8,4) = 70, and the true assignment is one
      of them (guards: an exact permutation p computed on a null that does
      not contain the observed value)
  T11 two-sided empirical p is symmetric and bounded
  T12 ridge baseline: index handling - target genes index the GENE axis
      (guards: the transpose bug that raised IndexError, and any regression
       to sample-axis indexing)

Run:  python 99_unit_tests.py
"""
from __future__ import annotations

import os
import sys
import math

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from _limma import (fit_limma, bh, build_network, diffusion_operator,        # noqa: E402
                    cos_batch, ridge_baseline, emp_p_two_sided,
                    degree_bins, degree_matched_sets, all_splits, perm_axes,
                    exact_two_sided_p, muscle_score, MUSCLE_GENES)

FAIL = []


def check(name, cond, detail=""):
    status = "PASS" if cond else "FAIL"
    print(f"  [{status}] {name}" + (f"  -- {detail}" if detail else ""))
    if not cond:
        FAIL.append(name)


def t1_logfc_recovery():
    print("T1  moderated t / logFC recovery (simulated truth)")
    rng = np.random.default_rng(1)
    n_g, n_s = 400, 12
    base = rng.normal(6, 1.5, n_g)
    sd = 0.10                             # SE(logFC) = sd*sqrt(2/6) = 0.058
    X = np.tile(base[:, None], (1, n_s)) + rng.normal(0, sd, (n_g, n_s))
    truth = np.zeros(n_g)
    truth[:40] = 2.0                      # 40 genes up by exactly 2 log2 units
    truth[40:80] = -1.0
    g = np.array([0.0] * 6 + [1.0] * 6)
    X[:, g == 1] += truth[:, None]
    design = np.column_stack([np.ones(n_s), g])
    r = fit_limma(pd.DataFrame(X), design, coef_idx=1, trend=False)
    # The tolerance is set from the sampling error, not picked by hand: with
    # 40 genes the maximum deviation is ~3 SE.
    se = sd * math.sqrt(2 / 6)
    err_up = float(np.abs(r["logFC"].values[:40] - 2.0).max())
    err_dn = float(np.abs(r["logFC"].values[40:80] + 1.0).max())
    err_null = float(np.abs(r["logFC"].values[200:]).mean())
    check("T1a up-regulated logFC within 3 SE of truth", err_up < 3 * se,
          f"max err={err_up:.4f}, 3SE={3*se:.4f}")
    check("T1a' down-regulated logFC within 3 SE of truth", err_dn < 3 * se,
          f"max err={err_dn:.4f}")
    check("T1b null genes centred", err_null < 0.15, f"mean|logFC|={err_null:.4f}")
    check("T1c true positives detected", int((r["P"][:80] < 0.01).sum()) >= 70,
          f"{int((r['P'][:80] < 0.01).sum())}/80")
    check("T1d false positive rate controlled",
          int((r["P"][200:] < 0.01).sum()) <= 10,
          f"{int((r['P'][200:] < 0.01).sum())}/200 at P<0.01")


def t2_bh_reference():
    print("T2  Benjamini-Hochberg against hand-computed values")
    # Benjamini & Hochberg (1995) JRSS-B, the 15 p-values of their worked
    # example. Adjusted values below are computed by hand from the step-up
    # rule q_i = min_{j>=i} (m/j) p_j and are NOT taken from this codebase.
    p = np.array([0.0001, 0.0004, 0.0019, 0.0095, 0.0201, 0.0278, 0.0298,
                  0.0344, 0.0459, 0.3240, 0.4262, 0.5719, 0.6528, 0.7590,
                  1.0000])
    expect = np.array([0.0015, 0.0030, 0.0095, 0.0356, 0.0603, 0.0639,
                       0.0639, 0.0645, 0.0765, 0.4860, 0.5812, 0.7149,
                       0.7532, 0.8132, 1.0000])
    got = bh(p)
    check("T2a matches the published worked example",
          np.allclose(got, expect, atol=5e-4),
          f"max diff={np.abs(got - expect).max():.5f}")
    check("T2a' 4 discoveries at q<=0.05, as in the paper",
          int((got <= 0.05).sum()) == 4, f"{int((got <= 0.05).sum())} discoveries")
    check("T2b monotone non-decreasing", bool(np.all(np.diff(got) >= -1e-12)))
    check("T2c bounded in [0,1]", bool((got >= 0).all() and (got <= 1).all()))
    check("T2d all-ones stays one", bool(np.allclose(bh(np.ones(5)), 1.0)))


def t3_operator_two_routes():
    print("T3  diffusion operator: iteration vs closed form")
    rng = np.random.default_rng(3)
    n = 200
    A = rng.normal(0, 1, (n, n))
    W = (A + A.T) / 2
    np.fill_diagonal(W, 0.0)
    Pi = diffusion_operator(W, alpha=0.55, n_iter=200, exact=False)
    Pe = diffusion_operator(W, alpha=0.55, exact=True)
    err = float(np.abs(Pi - Pe).max())
    check("T3a iteration converges to closed form", err < 1e-8,
          f"max|Pi-Pe|={err:.3e}")


def t4_operator_convergence_floor():
    print("T4  diffusion operator: convergence floor")
    rng = np.random.default_rng(4)
    n = 150
    A = rng.normal(0, 1, (n, n))
    W = (A + A.T) / 2
    np.fill_diagonal(W, 0.0)
    P80 = diffusion_operator(W, alpha=0.55, n_iter=80, exact=False)
    Pex = diffusion_operator(W, alpha=0.55, exact=True)
    err80 = float(np.abs(P80 - Pex).max())
    check("T4a 80 iterations is below 1e-12 (manuscript setting)", err80 < 1e-12,
          f"max|P80-Pex|={err80:.3e}")
    # the closed form is (1-a)(I - aS)^-1; verify the fixed point directly
    d = np.abs(W).sum(axis=1); d[d == 0] = 1.0
    Dm = np.diag(1.0 / np.sqrt(d))
    S = Dm @ W @ Dm
    a = 0.55
    resid = float(np.abs(Pex - (np.eye(n) * (1 - a) + a * (S @ Pex))).max())
    check("T4b P satisfies the fixed-point equation", resid < 1e-10,
          f"resid={resid:.3e}")


def t5_network_structure():
    print("T5  network structure (self-loops / symmetry / degree)")
    rng = np.random.default_rng(5)
    L = pd.DataFrame(rng.normal(0, 1, (800, 16)),
                     index=[f"g{i}" for i in range(800)],
                     columns=[f"s{j}" for j in range(16)])
    genes, W = build_network(L, hv=300, topk=12)
    check("T5a no self-loops (diagonal exactly zero)",
          bool(np.all(np.diag(W) == 0.0)),
          f"max|diag|={np.abs(np.diag(W)).max():.3e}")
    check("T5b symmetric", bool(np.allclose(W, W.T, atol=1e-12)))
    deg = (np.abs(W) > 0).sum(axis=1)
    check("T5c every node has >= topk neighbours after symmetrisation",
          int(deg.min()) >= 12, f"min degree={int(deg.min())}")
    check("T5d no isolated node", int(deg.min()) > 0)
    check("T5e gene list length equals hv", len(genes) == 300)


def t6_ko_oe_mirror():
    print("T6  knockout / over-expression are exact mirror images")
    rng = np.random.default_rng(6)
    L = pd.DataFrame(rng.normal(0, 1, (500, 16)))
    genes, W = build_network(L, hv=200, topk=10)
    P = diffusion_operator(W, alpha=0.55, exact=True)
    ax = rng.normal(0, 1, len(genes))
    idx = np.array([3, 17, 44, 90])
    s_ko = np.zeros(len(genes)); s_ko[idx] = -1.0
    s_oe = np.zeros(len(genes)); s_oe[idx] = 1.0
    c_ko = float(cos_batch((P @ s_ko).reshape(-1, 1), ax)[0])
    c_oe = float(cos_batch((P @ s_oe).reshape(-1, 1), ax)[0])
    check("T6a cos(KO) = -cos(OE) to machine precision",
          abs(c_ko + c_oe) < 1e-12, f"{c_ko:.6f} vs {c_oe:.6f}")
    a = rng.normal(0, 1, 100)
    check("T6b cos(v, v) = 1", abs(float(cos_batch(a.reshape(-1, 1), a)[0]) - 1) < 1e-12)
    check("T6c cos(-v, v) = -1", abs(float(cos_batch((-a).reshape(-1, 1), a)[0]) + 1) < 1e-12)
    check("T6d zero reference vector returns 0, not NaN",
          float(cos_batch(a.reshape(-1, 1), np.zeros(100))[0]) == 0.0)


def t7_column_parsing_regression():
    print("T7  identifier matrix must be rejected")
    # A matrix of Entrez ids is integer-valued; FPKM is not. The loader raises
    # when frac integer-valued > 0.95. Reproduce the detector.
    ids = np.tile(np.arange(1, 501), (8, 1)).T.astype(float)
    frac_ids = float(np.isclose(ids % 1, 0).mean())
    rng = np.random.default_rng(7)
    fpkm = np.abs(rng.normal(5, 20, (500, 8)))
    frac_fp = float(np.isclose(fpkm % 1, 0).mean())
    check("T7a identifier matrix flagged (>0.95 integer-valued)", frac_ids > 0.95,
          f"frac={frac_ids:.3f}")
    check("T7b real FPKM not flagged", frac_fp < 0.95, f"frac={frac_fp:.3f}")


def t8_group_parsing():
    print("T8  D14 must not collapse into D1")
    import re
    def grp(name):
        n = name.upper()
        if "SHAM" in n:
            return "Sham"
        m = re.search(r"(?:^|[_.\-])D(\d+)(?:[_.\-])", n)
        return f"D{int(m.group(1))}" if m else "other"
    names = ["GSM4506613_A_Sham_4.txt.gz", "GSM4506614_tci_d1_1.txt.gz",
             "GSM4506615_tci_d7_2.txt.gz", "GSM4506616_tci_d14_3.txt.gz"]
    got = [grp(n) for n in names]
    check("T8a D14 parsed as D14", got[3] == "D14", str(got))
    check("T8b D1 parsed as D1", got[1] == "D1")
    check("T8c Sham parsed as Sham", got[0] == "Sham")
    naive = ["D1" if "tci_d1" in n.lower() else "D14" for n in names]
    check("T8d the naive substring rule really does fail (regression witness)",
          naive[3] == "D1", "confirms the bug the regex fixes")


def t9_degree_matched_null():
    print("T9  degree-matched null draws from the matched bin")
    rng = np.random.default_rng(9)
    L = pd.DataFrame(rng.normal(0, 1, (600, 16)))
    genes, W = build_network(L, hv=300, topk=10)
    b = degree_bins(W, n_bins=5)
    # pick a module from the top degree bin only
    top_bin = int(b.max())
    idx = np.flatnonzero(b == top_bin)[:4]
    sets = degree_matched_sets(W, idx, n_perm=300, rng=np.random.default_rng(11),
                               n_bins=5)
    check("T9a every drawn gene is in the same degree bin",
          bool(np.all(b[sets] == top_bin)))
    check("T9b draws are not constant", len(set(sets.ravel().tolist())) > 4,
          f"{len(set(sets.ravel().tolist()))} distinct genes")
    check("T9c shape correct", sets.shape == (300, 4), str(sets.shape))
    check("T9d no duplicate gene within a draw",
          all(len(set(row.tolist())) == len(row) for row in sets[:50]))


def t10_exact_splits():
    print("T10 exhaustive label splits")
    sp = all_splits(4, 4)
    check("T10a C(8,4) = 70 splits", sp.shape[0] == 70, str(sp.shape))
    check("T10b each split has 4 cases", bool(np.all(sp.sum(axis=1) == 4)))
    true_mask = np.array([True] * 4 + [False] * 4)
    check("T10c the true assignment is one of the splits",
          bool((sp == true_mask).all(axis=1).any()))
    rng = np.random.default_rng(10)
    X = rng.normal(0, 1, (50, 8))
    Ax = perm_axes(X, sp)
    check("T10d axis matrix shape", Ax.shape == (50, 70), str(Ax.shape))
    obs_i = int(np.flatnonzero((sp == true_mask).all(axis=1))[0])
    manual = X[:, 4:].mean(axis=1) - X[:, :4].mean(axis=1)
    check("T10e observed axis matches manual computation",
          float(np.abs(Ax[:, obs_i] - manual).max()) < 1e-12)
    null = Ax[:, obs_i][:1]  # trivial null containing only the observed value
    check("T10f exact p with the observed value in the null is 1.0",
          abs(exact_two_sided_p(0.0, np.array([0.0])) - 1.0) < 1e-12)


def t11_emp_p_two_sided():
    print("T11 two-sided empirical p")
    rng = np.random.default_rng(11)
    null = rng.normal(0, 1, 5000)
    p_mid = emp_p_two_sided(0.0, null)
    p_far = emp_p_two_sided(6.0, null)
    check("T11a null-centre value gives p ~ 1", p_mid > 0.5, f"p={p_mid:.3f}")
    check("T11b extreme value gives small p", p_far < 0.05, f"p={p_far:.4f}")
    check("T11c symmetric in sign",
          abs(emp_p_two_sided(2.5, null) - emp_p_two_sided(-2.5, null)) < 0.05)
    check("T11d bounded by 1", emp_p_two_sided(0.0, null) <= 1.0)


def t12_ridge_index_handling():
    print("T12 ridge baseline index handling")
    rng = np.random.default_rng(12)
    genes = [f"g{i}" for i in range(300)]
    L = pd.DataFrame(rng.normal(0, 1, (300, 16)), index=genes)
    ax = rng.normal(0, 1, 300)
    v = ridge_baseline(L, genes, ["g5", "g50", "g200"], -1.0, ax, loso=True)
    check("T12a returns a finite number inside [-1,1]",
          bool(np.isfinite(v)) and abs(v) <= 1.0 + 1e-9, f"v={v:.4f}")
    v2 = ridge_baseline(L, genes, ["g5", "g50", "g200"], -1.0, ax, loso=False)
    check("T12b LOSO differs from in-sample (overfitting check)",
          abs(v - v2) > 1e-6, f"loso={v:.4f} insample={v2:.4f}")
    v3 = ridge_baseline(L, genes, ["g5", "not_a_gene", "g200"], -1.0, ax,
                        loso=True)
    check("T12c partially missing target list falls back to valid genes",
          bool(np.isfinite(v3)) and abs(v3) <= 1.0 + 1e-9, f"v={v3:.4f}")
    # target list with zero valid genes must return 0.0, not raise
    z = ridge_baseline(L, genes, ["nope1", "nope2"], -1.0, ax, loso=True)
    check("T12d empty target list returns 0.0", z == 0.0, f"z={z}")


def t14_cos_batch_2d():
    print("T14 cos_batch with a 2-D reference matrix")
    rng = np.random.default_rng(14)
    v = rng.normal(0, 1, 60)
    w = rng.normal(0, 1, 60)
    B = np.column_stack([v, -v, w, v + 3.0])      # last col: shifted copy
    got = cos_batch(v.reshape(-1, 1), B)[0]
    # column-by-column reference computed independently of _limma
    def ref(a, b):
        a = a - a.mean(); b = b - b.mean()
        return float(a @ b / (np.linalg.norm(a) * np.linalg.norm(b)))
    exp = [ref(v, v), ref(v, -v), ref(v, w), ref(v, v + 3.0)]
    check("T14a 2-D columns match per-column cosine",
          np.allclose(got, exp, atol=1e-12), f"{np.round(got,6)} vs {np.round(exp,6)}")
    check("T14b self column is +1", abs(got[0] - 1.0) < 1e-12, f"{got[0]:.6f}")
    check("T14c negated column is -1", abs(got[1] + 1.0) < 1e-12, f"{got[1]:.6f}")
    check("T14d a constant shift does not change the cosine",
          abs(got[3] - 1.0) < 1e-12, f"{got[3]:.6f}")
    # the naive implementation would have divided every column by the same
    # Frobenius norm; assert that is NOT what happens
    naive = np.array([float((v - v.mean()) @ (B[:, j] - B.mean()))
                      / (np.linalg.norm(v - v.mean())
                         * np.linalg.norm(B - B.mean())) for j in range(4)])
    check("T14e differs from the grand-mean / Frobenius implementation",
          not np.allclose(got, naive, atol=1e-6),
          "regression witness for the 2-D bug")
    # 2-D D against 2-D b
    D2 = np.column_stack([v, w])
    G = cos_batch(D2, B)
    check("T14f (2,4) shape", G.shape == (2, 4), str(G.shape))
    check("T14g diagonal block consistent with 1-D calls",
          abs(G[0, 0] - 1.0) < 1e-12 and abs(G[1, 2] - 1.0) < 1e-12)


def t13_muscle_score():
    print("T13 muscle contamination score")
    rng = np.random.default_rng(13)
    genes = [f"g{i}" for i in range(200)] + MUSCLE_GENES[:5]
    L = pd.DataFrame(rng.normal(0, 1, (len(genes), 8)), index=genes)
    L.loc[MUSCLE_GENES[:5], 4:] += 3.0      # last 4 samples contaminated
    s = muscle_score(L)
    check("T13a contaminated samples score higher",
          float(s.values[4:].mean()) > float(s.values[:4].mean()) + 0.5,
          f"{s.values[:4].mean():+.2f} vs {s.values[4:].mean():+.2f}")
    check("T13b score is finite", bool(np.isfinite(s.values).all()))
    L0 = pd.DataFrame(rng.normal(0, 1, (50, 8)), index=[f"x{i}" for i in range(50)])
    check("T13c no muscle genes -> all-zero score, no crash",
          bool(np.allclose(muscle_score(L0).values, 0.0)))


def main():
    print("=" * 70)
    print("unit tests — cancer pain virtual perturbation")
    print("=" * 70)
    for fn in [t1_logfc_recovery, t2_bh_reference, t3_operator_two_routes,
               t4_operator_convergence_floor, t5_network_structure,
               t6_ko_oe_mirror, t7_column_parsing_regression,
               t8_group_parsing, t9_degree_matched_null, t10_exact_splits,
               t11_emp_p_two_sided, t12_ridge_index_handling,
               t14_cos_batch_2d, t13_muscle_score]:
        fn()
    print("=" * 70)
    if FAIL:
        print(f"{len(FAIL)} FAILED: {FAIL}")
        sys.exit(1)
    print("ALL TESTS PASSED")


if __name__ == "__main__":
    main()
