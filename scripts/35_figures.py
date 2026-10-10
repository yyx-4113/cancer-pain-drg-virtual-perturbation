#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""35_figures.py — publication figures for the manuscript.

Fig 1  study design schematic
Fig 2  cross-aetiology divergence of canonical markers (Bone_D7 vs PTX)
Fig 3  anti-selectivity of the perturbation layer: |z| target vs control

All plotted numbers are read at run time from results/tables/*.csv.
Outputs docs/manuscript/figures/fig1..3.{png,pdf}
"""
from __future__ import annotations
import os
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAB = os.path.join(ROOT, "results", "tables")
FIG = os.path.join(ROOT, "docs", "manuscript", "figures")
os.makedirs(FIG, exist_ok=True)

plt.rcParams.update({"font.size": 9, "font.family": "DejaVu Sans",
                     "axes.linewidth": 0.8, "savefig.dpi": 600,
                     "axes.spines.top": False, "axes.spines.right": False})

UP, DN, CTL = "#c0392b", "#2471a3", "#7f8c8d"   # red up, blue down, grey control


# ------------------------------------------------------------------ Figure 1
def fig1():
    fig, ax = plt.subplots(figsize=(7.2, 3.6))
    ax.set_xlim(0, 10); ax.set_ylim(0, 6); ax.axis("off")

    def box(x, y, w, h, text, fc):
        ax.add_patch(FancyBboxPatch((x, y), w, h,
                     boxstyle="round,pad=0.06,rounding_size=0.12",
                     linewidth=0.8, edgecolor="#34495e", facecolor=fc))
        ax.text(x + w/2, y + h/2, text, ha="center", va="center", fontsize=8)

    # discovery
    box(0.15, 3.9, 2.6, 1.6,
        "GSE149648 (rat DRG)\nWalker-256 tibia\nsham, d1, d7, d14\nn = 4/group", "#fdecea")
    box(0.15, 1.9, 2.6, 1.6,
        "GSE160543 (rat DRG)\npaclitaxel / oxaliplatin\nvehicle\nn = 4/group", "#eaf2fb")
    # validation
    box(0.15, 0.15, 2.6, 1.4,
        "GSE292522 (mouse DRG)\nB16F10 vs keratinocyte\nn = 4 + 4 animals\n55,967 cells", "#eef7ee")

    # analysis
    box(3.6, 3.9, 2.6, 1.6,
        "exact permutation\nC(8,4)² = 4,900 pairs\n22 modules + markers", "#f5f5f5")
    box(3.6, 1.9, 2.6, 1.6,
        "in-silico perturbation\n17 targets vs 5 controls\ndegree-matched nulls", "#f5f5f5")
    box(3.6, 0.15, 2.6, 1.4,
        "per-animal pseudobulk\nlimma-trend\ncompartment bins", "#f5f5f5")

    # outputs
    box(7.0, 3.9, 2.85, 1.6,
        "Q1  aetiologies diverge\n(myeloid/interferon\nvs neuronal injury)", "#e8f6f3")
    box(7.0, 1.9, 2.85, 1.6,
        "Q2  no target ranking\nanti-selective vs\nnegative controls", "#fef5e7")
    box(7.0, 0.15, 2.85, 1.4,
        "interferon programme\nis model-specific", "#e8f6f3")

    for y in (4.7, 2.7, 0.85):
        ax.add_patch(FancyArrowPatch((2.75, y), (3.58, y), arrowstyle="-|>",
                     mutation_scale=9, lw=0.9, color="#34495e"))
    for y in (4.7, 2.7, 0.85):
        ax.add_patch(FancyArrowPatch((6.2, y), (6.98, y), arrowstyle="-|>",
                     mutation_scale=9, lw=0.9, color="#34495e"))
    # left labels
    for y, t in [(4.7, "discovery"), (2.7, "the layer"), (0.85, "validation")]:
        ax.text(0.15, y + 1.05, t, fontsize=7.5, style="italic", color="#566573")

    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(FIG, f"fig1.{ext}"), bbox_inches="tight")
    plt.close(fig)
    print("fig1 written")


# ------------------------------------------------------------------ Figure 2
def fig2():
    mk = pd.read_csv(os.path.join(TAB, "etiology_marker_cross.csv"))
    pair = mk[(mk.aetiology_A == "Bone_D7") & (mk.aetiology_B == "PTX")]
    genes = ["Cd68", "Aif1", "Ccl2", "Il1b", "Gal", "Calca", "Scn10a"]
    rows = []
    for g in genes:
        r = pair[pair.gene == g]
        if len(r):
            rows.append((g, float(r.iloc[0].logFC_A), float(r.iloc[0].logFC_B)))
    fig, ax = plt.subplots(figsize=(4.6, 3.4))
    ys = np.arange(len(rows))[::-1]
    for (g, a, b), y in zip(rows, ys):
        ax.plot([b, a], [y, y], color="#bdc3c7", lw=1.0, zorder=1)
        ax.scatter([a], [y], color=UP, s=34, zorder=3,
                   label="bone cancer (d7)" if y == ys[0] else None)
        ax.scatter([b], [y], color=DN, s=34, zorder=3,
                   label="paclitaxel" if y == ys[0] else None)
    ax.axvline(0, color="#34495e", lw=0.7)
    ax.set_yticks(ys); ax.set_yticklabels([r[0] for r in rows], fontsize=8.5)
    ax.set_xlabel("log$_2$ fold change vs control")
    ax.legend(frameon=False, fontsize=7.5, loc="lower right")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(FIG, f"fig2.{ext}"), bbox_inches="tight")
    plt.close(fig)
    print("fig2 written")


# ------------------------------------------------------------------ Figure 3
def fig3():
    bt = pd.read_csv(os.path.join(TAB, "etiology_bootstrap_differences.csv"))
    dk = bt[(bt.network_scope == "per_contrast") & (bt.perturbation == "knockout")].copy()
    dk["absz"] = dk.z_vs_randset.abs()
    rng = np.random.default_rng(0)
    fig, ax = plt.subplots(figsize=(4.6, 3.4))
    for i, (kind, col, lab) in enumerate([("target", DN, "analgesic target"),
                                          ("control", UP, "negative control")]):
        s = dk[dk.kind == kind].absz.values
        x = i + rng.uniform(-0.14, 0.14, len(s))
        ax.scatter(x, s, s=16, color=col, alpha=0.65, label=lab)
        ax.hlines(np.median(s), i-0.25, i+0.25, color="#2c3e50", lw=1.8)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["analgesic\ntargets", "negative\ncontrols"], fontsize=8.5)
    ax.set_ylabel("|z| vs degree-matched null")
    ax.set_ylim(-0.3, max(dk.absz.max() * 1.12, 8))
    ax.legend(frameon=False, fontsize=7.5, loc="upper left")
    nq = json.load(open(os.path.join(TAB, "negative_result_quantification.json")))
    p_mwu = nq["T1_target_vs_control"]["mwu_abs_z_two_sided_p"]
    ax.set_title("Mann–Whitney two-sided p = %.4f" % p_mwu, fontsize=8.5)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(FIG, f"fig3.{ext}"), bbox_inches="tight")
    plt.close(fig)
    print("fig3 written")


if __name__ == "__main__":
    fig1(); fig2(); fig3()
    print("figures ->", FIG)
