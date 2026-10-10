#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
98_check_manuscript_numbers.py
==============================
Asserts every number reported in the manuscript against the deposited
artefacts. Exits non-zero on any mismatch.

WHY THIS EXISTS: the manuscript and the analysis are two files that drift.
Three of the four fatal bugs found in this project were invisible in any
single artefact and only surfaced when two code paths were compared. A
manuscript that states "4,600 genes" while the CSV says "3,041" is the same
class of failure. Numbers here are therefore not copied by hand into prose
and then forgotten; they are re-derived from the CSV/JSON at check time and
compared with the string that appears in the manuscript file.

Run:  python 98_check_manuscript_numbers.py
"""
from __future__ import annotations

import os
import sys
import json
import math
import re

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
TAB = os.path.join(ROOT, "results", "tables")
DOC = os.path.join(ROOT, "docs", "manuscript")

FAILS = []
NCHK = 0


def chk(label, got, want, tol=0.0):
    """Compare a number, and record the outcome."""
    global NCHK
    NCHK += 1
    ok = False
    if isinstance(want, str):
        ok = (str(got) == want)
    elif want is None or (isinstance(want, float) and math.isnan(want)):
        ok = True
    else:
        ok = abs(float(got) - float(want)) <= tol
    if not ok:
        FAILS.append(f"{label}: got {got!r}, manuscript says {want!r}")
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}: {got}"
          + ("" if ok else f"  != {want}"))


MS = sys.argv[1] if len(sys.argv) > 1 else "manuscript_final.md"


def in_text(label, needle, fname=MS):
    """Assert a string appears in the manuscript (guards against silent edits
    that leave the artefacts untouched)."""
    global NCHK
    NCHK += 1
    p = os.path.join(DOC, fname)
    if not os.path.exists(p):
        FAILS.append(f"{label}: manuscript missing at {p}")
        print(f"  [FAIL] {label}: manuscript missing")
        return
    txt = open(p, encoding="utf-8").read()
    ok = needle in txt
    if not ok:
        FAILS.append(f"{label}: '{needle}' not found in {fname}")
    print(f"  [{'PASS' if ok else 'FAIL'}] text: '{needle[:60]}'")


def load(f):
    return pd.read_csv(os.path.join(TAB, f))


def jload(f):
    return json.load(open(os.path.join(TAB, f)))


def main():
    print("=" * 72)
    print("manuscript <-> artefact consistency check")
    print("=" * 72)

    # ---------------- discovery: study description ----------------
    print("\n[discovery] study description")
    s20 = jload("run_summary.json")
    chk("GSE149648 samples", s20["n_samples"], 16)
    chk("GSE149648 genes after filter", s20["n_genes"], 13893)
    for g in ["Sham", "D1", "D7", "D14"]:
        chk(f"group {g} n", s20["groups"][g], 4)
    chk("muscle ANOVA F", s20["muscle_ANOVA"]["F"], 0.637, tol=0.001)
    chk("muscle ANOVA p", s20["muscle_ANOVA"]["p"], 0.605, tol=0.001)
    chk("network genes", s20["network"]["genes"], 3042)
    chk("network topk", s20["network"]["topk"], 12)
    chk("alpha", s20["network"]["alpha"], 0.55, tol=1e-9)

    # ---------------- discovery: DE counts ----------------
    print("\n[discovery] differential expression counts (muscle-adjusted)")
    want_adj = {"D1": 0, "D7": 4600, "D14": 235}
    want_unadj = {"D1": 1, "D7": 2563, "D14": 958}
    for lab in ["D1", "D7", "D14"]:
        chk(f"{lab} adj q<0.05", s20["n_sig_q05_adj"][lab], want_adj[lab])
        chk(f"{lab} unadj q<0.05", s20["n_sig_q05_unadj"][lab], want_unadj[lab])
    chk("D7 adj q<0.10", s20["n_sig_q10_adj"]["D7"], 5867)
    chk("D7 unadj q<0.05", s20["n_sig_q05_unadj"]["D7"], 2563)

    # ---------------- discovery: muscle contamination ----------------
    print("\n[discovery] muscle contamination")
    qc = load("qc_muscle_scores.csv")
    top = qc.sort_values("muscle_fpkm_pct", ascending=False).iloc[0]
    chk("worst sample id", top["sample"], "GSM4506604")
    chk("worst sample muscle share %", top["muscle_fpkm_pct"], 6.26, tol=0.01)
    chk("only one sample above 1%", int((qc["muscle_fpkm_pct"] > 1.0).sum()), 1)
    chk("worst sample group", top["group"], "D7")

    md = load("de_muscle_drop_sensitivity.csv")
    r = md.iloc[0]
    chk("muscle-drop contrast", r["contrast"], "D7")
    chk("muscle-drop n_q05 reduced", int(r["n_q05_reduced"]), 3041)
    chk("muscle-drop jaccard", float(r["jaccard_q05"]), 0.584, tol=0.001)
    chk("muscle-drop spearman logFC", float(r["spearman_logFC"]), 0.912, tol=0.001)

    # ---------------- discovery: module-level t ----------------
    print("\n[discovery] module mean moderated t at D7")
    mod = load("module_level_signal.csv")
    m7 = mod[mod.contrast == "D7"].set_index("module")
    want_t = {"CTRL interferon": 24.923, "Glial / immune": 11.412,
              "Toll-like receptors": 6.364, "Cytokine / chemokine": 5.877,
              "CTRL housekeeping": 0.653, "GABAergic": -5.708,
              "K+ / HCN channels": -5.635, "Neurotrophin": -5.292,
              "CTRL pan-neuronal": -5.201, "CTRL erythrocyte": -4.801,
              "Glutamate / NMDA": -4.353, "Voltage-gated Na+": -4.347,
              "CTRL myelin / Schwann": -4.101, "Opioid receptors": -3.761,
              "Ca2+ channel a2d": -3.509}
    for m, v in want_t.items():
        chk(f"D7 mean t {m}", float(m7.loc[m, "mean_t"]), v, tol=0.002)
    chk("interferon q", float(m7.loc["CTRL interferon", "q_comp"]), 6.8e-5, tol=1e-5)

    # ---------------- discovery: perturbation ----------------
    print("\n[discovery] in-silico perturbation at D7")
    pert = load("virtual_perturbation_scores.csv")
    p7 = pert[pert.contrast == "D7"].set_index(["module", "perturbation"])
    want_rev = {("CTRL interferon", "knockout"): 0.228,
                ("Glial / immune", "knockout"): 0.155,
                ("CTRL pan-neuronal", "knockout"): -0.109,
                ("K+ / HCN channels", "knockout"): -0.092,
                ("Voltage-gated Na+", "knockout"): -0.090,
                ("CTRL myelin / Schwann", "knockout"): -0.090}
    for (m, pt), v in want_rev.items():
        chk(f"D7 reversal {m}/{pt}", float(p7.loc[(m, pt), "reversal_network"]),
            v, tol=0.001)
    chk("interferon pct of oracle",
        float(p7.loc[("CTRL interferon", "knockout"), "pct_of_oracle"]),
        75.4, tol=0.1)
    chk("interferon z vs null",
        float(p7.loc[("CTRL interferon", "knockout"), "z_vs_perm"]), 4.61, tol=0.01)
    # no module may be represented by fewer than 3 genes
    chk("min genes per module is >= 3", int(pert["n_targets"].min()), 3)

    h2h = jload("head_to_head_baselines.json")
    chk("Spearman(network, additive)", h2h["spearman_network_vs_additive"],
        0.991, tol=0.001)
    chk("mean |network|", h2h["mean_abs_network"], 0.057, tol=0.001)
    chk("mean |additive|", h2h["mean_abs_additive"], 0.033, tol=0.001)
    chk("mean pct of oracle", h2h["mean_pct_of_oracle"], 20.5, tol=0.1)

    # ---------------- cross-aetiology ----------------
    print("\n[cross-aetiology] design")
    s22 = jload("etiology_summary.json")
    chk("shared genes", s22["shared_genes"], 13601)
    chk("network nodes (per contrast)",
        s22["common_network_genes_per_contrast"], 3038)
    for k in ["Bone_D14", "Bone_D7", "PTX", "OXA"]:
        chk(f"{k} case n", s22["contrasts"][k]["case"], 4)
        chk(f"{k} control n", s22["contrasts"][k]["control"], 4)

    print("\n[cross-aetiology] module mean logFC Delta")
    mx = load("etiology_module_logfc_cross.csv").set_index(
        ["module", "aetiology_A", "aetiology_B"])
    want_mod = {("Glial / immune", "Bone_D7", "PTX"): (2.695, 0.0073),
                ("Glial / immune", "Bone_D14", "PTX"): (1.945, 0.0073),
                ("Cytokine / chemokine", "Bone_D14", "OXA"): (1.809, 0.0073),
                ("Cytokine / chemokine", "Bone_D7", "PTX"): (1.703, 0.0073),
                ("Toll-like receptors", "Bone_D7", "PTX"): (1.275, 0.0073),
                ("Neurotrophin", "Bone_D7", "PTX"): (-0.972, 0.0073),
                ("Ca2+ channel a2d", "Bone_D7", "PTX"): (-0.600, 0.0073),
                ("Neuropeptides", "Bone_D7", "PTX"): (-1.345, 0.0145),
                ("CTRL pan-neuronal", "Bone_D7", "PTX"): (-0.691, 0.0145),
                ("Glutamate / NMDA", "Bone_D7", "PTX"): (-0.635, 0.0162),
                ("Voltage-gated Na+", "Bone_D7", "PTX"): (-0.654, 0.0310),
                ("GABAergic", "Bone_D7", "PTX"): (-0.495, 0.0310),
                ("CTRL erythrocyte", "Bone_D7", "OXA"): (-2.483, 0.0073)}
    for (m, a, b), (d, q) in want_mod.items():
        chk(f"Delta logFC {m} {a} vs {b}",
            float(mx.loc[(m, a, b), "delta_logFC"]), d, tol=0.002)
        chk(f"q {m} {a} vs {b}",
            float(mx.loc[(m, a, b), "q_permutation"]), q, tol=0.0003)

    print("\n[cross-aetiology] canonical markers")
    kx = load("etiology_marker_cross.csv").set_index(
        ["gene", "aetiology_A", "aetiology_B"])
    want_mk = {("Cd68", "Bone_D7", "PTX"): (3.418, 0.0125),
               ("Aif1", "Bone_D7", "PTX"): (3.076, 0.0125),
               ("Cd68", "Bone_D14", "PTX"): (3.071, 0.0125),
               ("Aif1", "Bone_D14", "PTX"): (2.838, 0.0125),
               ("Gal", "Bone_D7", "PTX"): (-2.592, 0.0125),
               ("Ccl2", "Bone_D7", "OXA"): (1.391, 0.0173),
               ("Il1b", "Bone_D7", "PTX"): (1.230, 0.0182),
               ("Scn10a", "Bone_D7", "PTX"): (-0.702, 0.0260),
               ("Calca", "Bone_D7", "PTX"): (-0.770, 0.0272)}
    for (g, a, b), (d, q) in want_mk.items():
        chk(f"marker Delta {g} {a} vs {b}",
            float(kx.loc[(g, a, b), "delta_logFC"]), d, tol=0.002)
        chk(f"marker q {g} {a} vs {b}",
            float(kx.loc[(g, a, b), "q_permutation"]), q, tol=0.0003)

    print("\n[cross-aetiology] perturbation layer calibration")
    bt = load("etiology_bootstrap_differences.csv")
    pc = bt[bt.network_scope == "per_contrast"]
    ko = pc[pc.perturbation == "knockout"]
    oe = pc[pc.perturbation == "overexpression"]
    # 22 modules x 6 aetiology pairs x 2 perturbations
    chk("n tests (per contrast, raw)", len(pc), 264)
    chk("n modules", pc["module"].nunique(), 22)
    chk("n aetiology pairs", pc.groupby(["aetiology_A",
                                         "aetiology_B"]).ngroups, 6)
    # KO and OE are exact mirror images -> 132 independent tests per scope
    paired = ko.set_index(["module", "aetiology_A", "aetiology_B"])[
        "p_permutation_exact"]
    poe = oe.set_index(["module", "aetiology_A", "aetiology_B"])[
        "p_permutation_exact"]
    maxdiff = float((paired - poe).abs().max())
    chk("max|p_KO - p_OE| (must be ~0)", maxdiff, 0.0, tol=1e-9)
    chk("independent tests per scope (KO only)", len(ko), 132)
    chk("no Delta significant at q<0.05 (dedup 132)",
        int((ko.q_permutation < 0.05).sum()), 0)
    chk("min q_permutation (dedup 132)", float(ko.q_permutation.min()),
        0.165, tol=0.001)
    # randset null over the dedupulated set
    chk("randset sig (dedup per-contrast)", int((ko.q_randset < 0.05).sum()), 32)
    chk("randset sig among controls (dedup)",
        int(((ko.q_randset < 0.05) & (ko.kind == "control")).sum()), 14)
    chk("randset sig among targets (dedup)",
        int(((ko.q_randset < 0.05) & (ko.kind == "target")).sum()), 18)
    gi = pc[(pc.module == "Glial / immune") & (pc.aetiology_A == "Bone_D7")
            & (pc.aetiology_B == "PTX") & (pc.perturbation == "knockout")].iloc[0]
    chk("glial/immune Delta", float(gi["delta_observed"]), 0.274, tol=0.001)
    chk("glial/immune CI lo", float(gi["ci95_lo"]), 0.203, tol=0.001)
    chk("glial/immune CI hi", float(gi["ci95_hi"]), 0.287, tol=0.001)
    chk("glial/immune permutation p", float(gi["p_permutation_exact"]),
        0.0014, tol=0.0002)
    # the largest |z| against the random-gene-set null must be a control
    zmax = pc.loc[pc["z_vs_randset"].abs().idxmax()]
    chk("largest |z vs random sets| is a negative control",
        zmax["kind"], "control")
    chk("largest |z| module", zmax["module"], "CTRL erythrocyte")
    # bootstrap interval must contain the observed value for EVERY test
    cover = ((pc.ci95_lo <= pc.delta_observed)
             & (pc.delta_observed <= pc.ci95_hi)).mean()
    chk("bootstrap interval covers observed (all tests)", float(cover), 1.0,
        tol=1e-9)

    print("\n[independent single-cell validation] GSE292522")
    ssc = jload("sc_GSE292522_summary.json")
    chk("sc total cells", ssc["total_cells"], 55967)
    chk("sc up in B16F10 q<0.05", ssc["pseudobulk_DE"]["n_up_B16F10_q05"], 1517)
    chk("sc down in B16F10 q<0.05",
        ssc["pseudobulk_DE"]["n_down_B16F10_q05"], 703)
    chk("sc ISG up count (should be 0)", len(ssc["ISG_up_in_B16F10"]), 0)
    chk("sc myeloid fraction Con", ssc["myeloid_fraction_mean"]["Con"],
        0.00117, tol=1e-4)
    chk("sc myeloid fraction B16F10", ssc["myeloid_fraction_mean"]["B16F10"],
        0.000517, tol=1e-4)
    chk("sc myeloid MWU p (two-sided)", ssc["myeloid_MWU_p_twosided"], 0.1142857, tol=1e-3)
    chk("sc myeloid MWU p (one-sided less)", ssc["myeloid_MWU_p_less"], 0.0571429, tol=1e-3)
    chk("sc convergent spearman rho vs Walker-256 D7",
        ssc["convergent_validation_vs_bulk_D7"]["spearman_rho"], 0.098, tol=0.005)
    scde = load("sc_GSE292522_pseudobulk_DE.csv").set_index("Unnamed: 0")
    chk("sc top up gene is Sprr1a", scde.sort_values("t").iloc[-1].name, "Sprr1a")
    chk("sc Sprr1a logFC", float(scde.loc["Sprr1a", "logFC"]), 1.554, tol=0.01)

    # ---------------- marker-level DE table ----------------
    print("\n[cross-aetiology] marker logFC per contrast")
    mk = load("etiology_marker_de.csv").set_index(["gene", "contrast"])
    want_lf = {("Cd68", "Bone_D7"): 2.86, ("Aif1", "Bone_D7"): 2.18,
               ("Gfap", "Bone_D7"): 1.33, ("Ccl2", "Bone_D7"): 0.81,
               ("Il1b", "Bone_D7"): 0.80, ("Gal", "PTX"): 1.92,
               ("Atf3", "PTX"): 1.13, ("Npy", "PTX"): 1.12,
               ("Aif1", "PTX"): -0.87, ("Cd68", "PTX"): -0.58,
               ("Gal", "Bone_D7"): -0.69}
    for (g, c), v in want_lf.items():
        chk(f"logFC {g} {c}", float(mk.loc[(g, c), "logFC"]), v, tol=0.006)

    # ---------------- sensitivity grid ----------------
    print("\n[sensitivity] hyperparameter grid")
    sg = jload("sensitivity_summary.json")
    chk("grid points", sg["n_grid_points"], 75)
    chk("discovery min Spearman vs reference",
        sg["grid_discovery_min_spearman"], 0.959, tol=0.001)
    chk("discovery median Spearman vs reference",
        sg["grid_discovery_median_spearman"], 0.988, tol=0.001)
    chk("top-5 overlap minimum", sg["grid_discovery_min_top5_overlap"], 5)
    chk("Delta min Spearman vs reference",
        sg["grid_delta_min_spearman"], 0.989, tol=0.001)
    chk("Delta median Spearman vs reference",
        sg["grid_delta_median_spearman"], 0.994, tol=0.001)
    q05 = sg["de_q05_by_frac"]
    chk("DE q<0.05 range across lowess spans (min)",
        min(q05.values()), 316)
    chk("DE q<0.05 range across lowess spans (max)",
        max(q05.values()), 325)
    chk("Spearman(network, additive) cross-study space",
        sg["spearman_network_vs_additive"], 0.929, tol=0.001)

    # ---------------- negative-result quantification (T1/T2/T3/T5) ----------------
    # Expectations are DERIVED at run time: the MWU / binomial are recomputed
    # from the deposited CSV and compared against the JSON, then the JSON
    # values are required to appear in the manuscript text. Hardcoding a
    # stale number here would "protect" the error (2026-09-27 lesson).
    print("\n[negative-result quantification] discrimination + power")
    nq = jload("negative_result_quantification.json")
    t1 = nq["T1_target_vs_control"]
    t2 = nq["T2_randset_null_selectivity"]
    t5 = nq["T5_reversal_score_power"]
    # recompute T1 MWU on |z| from the deduplicated knockout universe
    from scipy.stats import mannwhitneyu as _mwu, binomtest as _bt
    bt = load("etiology_bootstrap_differences.csv")
    dk = bt[(bt.network_scope == "per_contrast") & (bt.perturbation == "knockout")]
    tgt = dk[dk.kind == "target"]
    ctl = dk[dk.kind == "control"]
    chk("T1 dedup universe size (JSON n_target+n_control)",
        t1["n_target_tests"] + t1["n_control_tests"], len(dk))
    u_abs = _mwu(ctl.z_vs_randset.abs(), tgt.z_vs_randset.abs(),
                 alternative="two-sided")
    chk("T1 MWU |z| two-sided recomputed vs JSON",
        round(float(u_abs.pvalue), 6), t1["mwu_abs_z_two_sided_p"], tol=1e-6)
    chk("T1 median |z| control recomputed vs JSON",
        round(float(ctl.z_vs_randset.abs().median()), 4),
        t1["median_abs_z_control"], tol=1e-4)
    chk("T1 median |z| target recomputed vs JSON",
        round(float(tgt.z_vs_randset.abs().median()), 4),
        t1["median_abs_z_target"], tol=1e-4)
    # T2 binomial over-representation of controls among randset-significant
    sig = dk[dk.q_randset < 0.05]
    base = len(ctl) / len(dk)
    btp = _bt(int((sig.kind == "control").sum()), len(sig), base,
              alternative="greater").pvalue
    chk("T2 randset-significant dedup total",
        t2["randset_significant_total"], len(sig))
    chk("T2 controls among randset-significant",
        t2["randset_significant_controls"], int((sig.kind == "control").sum()))
    chk("T2 binomial p recomputed vs JSON",
        round(float(btp), 5), t2["binomial_overrepresentation_p"], tol=1e-5)
    # T5 corrected power: KO/OE are exact mirrors; MDE = z_sum * bootstrap SE
    # (normal approximation, z_sum = 1.96 + 0.8416)
    se_d = (dk.ci95_hi - dk.ci95_lo) / (2 * 1.96)
    z_sum = 1.959964 + 0.841621
    mde_norm = z_sum * float(se_d.median())
    chk("T5 MDE (normal approx) recomputed vs JSON",
        round(mde_norm, 4), t5["mde_delta_units_normal_approx"], tol=1e-4)
    chk("T5 max target delta vs JSON",
        t5["max_target_delta_observed"],
        round(float(tgt.delta_observed.abs().max()), 4), tol=1e-4)
    chk("T5 median target |delta| recomputed vs JSON",
        round(float(tgt.delta_observed.abs().median()), 4),
        t5["median_target_abs_delta"], tol=1e-4)
    n_exceed = int((tgt.delta_observed.abs() > mde_norm).sum())
    chk("T5 targets exceeding MDE (normal) recomputed vs JSON",
        n_exceed, t5["n_target_tests_exceeding_mde_normal"])
    chk("T5 exact-permutation significant count",
        t5["n_exact_permutation_p_lt_0.05"],
        int((dk.p_permutation_exact < 0.05).sum()))
    # anti-selectivity is robust to dropping the interferon control
    noifn = dk[dk.module != "CTRL interferon"]
    u_noifn = _mwu(noifn[noifn.kind == "control"].z_vs_randset.abs(),
                   noifn[noifn.kind == "target"].z_vs_randset.abs(),
                   alternative="two-sided")
    chk("T1 anti-selectivity excluding interferon control (p<0.05)",
        bool(u_noifn.pvalue < 0.05), True)
    chk("T3 myeloid observed Cohen d (n-2 denominator)",
        nq["T3_single_cell_power"]["myeloid_observed_cohens_d"], 1.267, tol=0.005)
    chk("T3 gene-level MDE Cohen d",
        nq["T3_single_cell_power"]["mde_cohens_d_genes"], 2.38, tol=0.01)
    # manuscript must carry the derived numbers verbatim
    in_text("manuscript T1 MWU p", "%.4f" % t1["mwu_abs_z_two_sided_p"])
    in_text("manuscript T1 median |z| control", "%.2f" % t1["median_abs_z_control"])
    in_text("manuscript T1 median |z| target", "%.2f" % t1["median_abs_z_target"])
    in_text("manuscript T2 binomial p", "binomial p = %.4f" % t2["binomial_overrepresentation_p"])
    in_text("manuscript T5 MDE normal",
            "≈0.%.3d" % round(1000 * t5["mde_delta_units_normal_approx"]))
    in_text("manuscript T5 exact sig",
            "%d of the 132" % t5["n_exact_permutation_p_lt_0.05"])
    in_text("manuscript T5 targets above MDE (normal)",
            "%d (23 under the conservative" % n_exceed)
    in_text("manuscript target percentile",
            "%.1fth percentile"
            % nq["T4_target_rank_percentile"]["median_percentile_absz_below"])

    # ---------------- per-compartment single-cell (script 32) ----------------
    print("\n[single-cell per-compartment] GSE292522")
    ctde = load("sc_GSE292522_celltype_DE.csv")

    def ctq(cell, gene):
        r = ctde[(ctde.cell_type == cell) & (ctde.gene == gene)]
        return float(r.iloc[0]["q"]) if len(r) else float("nan")

    chk("per-compartment Ifit2 in Neuron q", ctq("Neuron", "Ifit2"),
        0.0459, tol=0.002)
    chk("per-compartment Isg15 in SatelliteGlia q",
        ctq("SatelliteGlia", "Isg15"), 0.0071, tol=0.002)
    chk("per-compartment Gal in Neuron q", ctq("Neuron", "Gal"),
        0.0001, tol=0.0005)
    chk("per-compartment Scn10a in Neuron q", ctq("Neuron", "Scn10a"),
        0.0117, tol=0.002)
    # no ISG may be up-regulated in any compartment (q < 0.05 and logFC > 0)
    ISGset = {"Ifit1", "Ifit2", "Ifit3", "Isg15", "Mx1", "Stat1", "Irf7",
              "Bst2", "Oas1a", "Oas1b"}
    isg_up = ctde[(ctde.gene.isin(ISGset)) & (ctde.logFC > 0) & (ctde.q < 0.05)]
    chk("per-compartment ISG up-regulated (should be 0)", len(isg_up), 0)

    # ---------------- Tables 2 and 3: every cell regenerated from CSVs ----------------
    print("\n[tables 2/3] cell-by-cell regeneration from deposited CSVs")
    mk = load("etiology_marker_cross.csv")
    md = load("etiology_module_logfc_cross.csv")

    def _num(s):
        s = s.replace("−", "-").replace("+", "").strip()
        try:
            return float(s)
        except ValueError:
            return None

    def parse_md_table(header, next_header):
        txt = open(MS_PATH, encoding="utf-8").read()
        sec = txt.split(header)[1].split(next_header)[0]
        rows = [l for l in sec.split("\n") if l.startswith("|")]
        return [[c.strip() for c in r.strip("|").split("|")] for r in rows]

    MS_PATH = os.path.join(DOC, MS)
    MODNAME = {"Ca²⁺ channel α2δ": "Ca2+ channel a2d",
               "Voltage-gated Na⁺": "Voltage-gated Na+",
               "K⁺ / HCN channels": "K+ / HCN channels"}
    n_cells = 0
    T2_COLS = {"Bone_D7": 1, "Bone_D14": 2, "PTX": 3, "OXA": 4}
    for row in parse_md_table("**Table 3.", "**Table 4.")[2:]:
        gene = row[0]
        for a, col in T2_COLS.items():
            cv = _num(row[col])
            if cv is None:
                continue
            r = mk[((mk.aetiology_A == a) | (mk.aetiology_B == a)) & (mk.gene == gene)]
            if not len(r):
                FAILS.append(f"Table 2 {gene} {a}: no CSV row")
                n_cells += 1
                continue
            rr = r.iloc[0]
            v = float(rr.logFC_A if rr.aetiology_A == a else rr.logFC_B)
            chk(f"Table 2 {gene} {a} logFC", round(v, 4), round(cv, 4), tol=0.005)
            n_cells += 1
    for row in parse_md_table("**Table 4.", "**Table 5.")[2:]:
        name = MODNAME.get(row[0], row[0])
        for a, col in [("Bone_D7", 2), ("Bone_D14", 3), ("PTX", 4), ("OXA", 5)]:
            cv = _num(row[col])
            if cv is None:
                continue
            r = md[(md.module == name) & ((md.aetiology_A == a) | (md.aetiology_B == a))]
            if not len(r):
                FAILS.append(f"Table 3 {name} {a}: no CSV row")
                n_cells += 1
                continue
            rr = r.iloc[0]
            v = float(rr.logFC_A if rr.aetiology_A == a else rr.logFC_B)
            chk(f"Table 3 {name} {a} logFC", round(v, 4), round(cv, 4), tol=0.005)
            n_cells += 1
    print(f"  ... {n_cells} table cells regenerated")

    # ---------------- text presence ----------------
    print(f"\n[manuscript] key statements present in {MS}")
    for s in ["4,600", "3,041", "0.991", "24.9", "+3.42", "−2.59", "6.26",
              "3/4,901", "3/71", "2/71", "0.165", "Jaccard 0.584", "0.959",
              "0.989", "264 contrast-level tests", "132 independent", "20.5",
              "0.929", "GSE292522", "Sprr1a", "model-specific", "55,967",
              "32/132", "anti-selectiv", "0.0019", "binomial p = 0.0067",
              "46 cells", "eight testable compartments", "1.27",
              "≈0.052", "counts per 10,000", "two control and four tumour"]:
        in_text("manuscript", s)

    # abstract length cap: many SCIE venues cap the abstract at 300 words
    _atxt = open(MS_PATH, encoding="utf-8").read()
    _ab = re.search(r"## Abstract\s*(.*?)\*\*Keywords", _atxt, re.S).group(1)
    _ab = re.sub(r"\*\*(Background|Methods|Results|Conclusions)\.\*\*", " ", _ab)
    _ab = re.sub(r"[*⁰¹²³⁴⁵⁻]+", " ", _ab)
    _n_ab = len([w for w in _ab.split() if any(c.isalnum() for c in w)])
    chk("abstract <= 300 words", _n_ab <= 300, True)

    # ---------------- reference DOIs + figure files ----------------
    print("\n[manuscript] reference DOIs and figure files")
    ref_block = open(MS_PATH, encoding="utf-8").read()
    ref_block = ref_block.split("## References")[1].split("## Tables")[0]
    n_refs = len([ln for ln in ref_block.split("\n") if re.match(r"^\d+\. ", ln)])
    n_doi = len([ln for ln in ref_block.split("\n") if re.match(r"^\d+\. ", ln) and "doi:" in ln])
    chk("all references carry a DOI", n_doi, n_refs)
    for fn in ("fig1.png", "fig2.png", "fig3.png"):
        chk(f"figure file present ({fn})",
            os.path.exists(os.path.join(DOC, "figures", fn)), True)

    print("=" * 72)
    print(f"{NCHK} checks, {len(FAILS)} failed")
    if FAILS:
        for f in FAILS:
            print("  FAIL:", f)
        sys.exit(1)
    print("ALL MANUSCRIPT NUMBERS VERIFIED AGAINST ARTEFACTS")


if __name__ == "__main__":
    main()
