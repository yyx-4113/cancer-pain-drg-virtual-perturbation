"""
10_bone_end_to_end.py
=====================
骨癌痛部位端到端分析：加载 GSE249443 -> 疼痛签名 -> 共表达网络 ->
虚拟敲除/过表达 -> 逆转评分 + 置换检验。

内部验证：原论文(GSE249443)认定 Cxcl13 为骨癌痛关键基因，
本脚本会检查 Cxcl13 是否在 BCP 组显著上调，作为管线可信度的阳性对照。
"""
import os
import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RAW = os.path.join(ROOT, "data", "raw")
PROC = os.path.join(ROOT, "data", "processed")
OUT = os.path.join(ROOT, "results", "tables")
os.makedirs(PROC, exist_ok=True)
os.makedirs(OUT, exist_ok=True)

XLSX = os.path.join(RAW, "GSE249443_processed_data.xlsx")
SITE = "bone"
N_PERM = 2000
SEED = 20261009

# 治疗性方向：KO=药物抑制/阻断该靶点；OE=药物激活该靶点
TARGETS = {
    "Oprm1":    ("OE", "阿片类(mu受体激动)"),
    "Scn9a":    ("KO", "Nav1.7 钠通道阻滞剂"),
    "Scn10a":   ("KO", "Nav1.8 钠通道阻滞剂"),
    "Scn11a":   ("KO", "Nav1.9 钠通道阻滞剂"),
    "Trpv1":    ("KO", "TRPV1 拮抗剂"),
    "Trpm3":    ("KO", "TRPM3 拮抗剂"),
    "Ptgs2":    ("KO", "COX-2 抑制剂(NSAIDs)"),
    "Ptgs1":    ("KO", "COX-1 抑制剂(NSAIDs)"),
    "Cacna2d1": ("KO", "alpha2delta-1(加巴喷丁类)"),
    "Cacna2d2": ("KO", "alpha2delta-2(加巴喷丁类)"),
    "Grin2b":   ("KO", "NMDA 2B(氯胺酮)"),
    "Slc6a4":   ("KO", "SERT(SNRI 度洛西汀)"),
    "Slc6a2":   ("KO", "NET(SNRI 度洛西汀)"),
    "Kcnq2":    ("KO", "Kv7.2 channel"),
    "Htr1a":    ("OE", "5-HT1A 激动"),
}


def bh_fdr(p):
    p = np.asarray(p, float)
    ok = ~np.isnan(p)
    q = np.full_like(p, np.nan)
    pv = p[ok]
    n = len(pv)
    if n == 0:
        return q
    o = np.argsort(pv)
    adj = pv[o] * n / np.arange(1, n + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.clip(adj, 0, 1)
    q[ok] = out
    return q


def load_expression():
    df = pd.read_excel(XLSX, sheet_name="原始")
    print(f"[load] shape={df.shape}")
    fpkm_cols = [c for c in df.columns if str(c).upper().startswith("FPKM")]
    print(f"[load] FPKM columns: {fpkm_cols}")
    if "gene_name" not in df.columns:
        raise RuntimeError("gene_name column missing")

    sub = df[["gene_name"] + fpkm_cols].copy()
    sub["gene_name"] = sub["gene_name"].astype(str).str.strip()
    sub = sub[sub["gene_name"].isin(["", "nan", "NaN", "-"]) == False]
    for c in fpkm_cols:
        sub[c] = pd.to_numeric(sub[c], errors="coerce")
    # 同一基因多转录本 -> 取均值
    expr = sub.groupby("gene_name")[fpkm_cols].mean()
    expr = expr[~expr.index.duplicated()]
    print(f"[load] genes after aggregation: {expr.shape[0]}")

    # 分组：先判 sBCP(对照)，再判 BCP(病例)，避免子串误判
    ctrl = [c for c in fpkm_cols if "sBCP" in c or "sbc" in c.lower()]
    case = [c for c in fpkm_cols if c not in ctrl and "BCP" in c]
    print(f"[load] control(sBCP) n={len(ctrl)} {ctrl}")
    print(f"[load] case(BCP)     n={len(case)} {case}")
    if len(ctrl) < 2 or len(case) < 2:
        raise RuntimeError("分组样本不足")
    return expr, ctrl, case


def main():
    expr, ctrl, case = load_expression()

    # 保留在多数样本表达的基因
    keep = (expr > 0).sum(axis=1) >= max(3, expr.shape[1] - 2)
    expr = expr[keep]
    print(f"[qc] genes kept (expressed): {expr.shape[0]}")

    # log2(FPKM+1)
    L = np.log2(expr + 1.0)

    # ---------------- 1. 疼痛签名 ----------------
    a = L[ctrl].values
    b = L[case].values
    lfc = b.mean(1) - a.mean(1)
    t, p = stats.ttest_ind(b, a, axis=1, equal_var=False)
    p = np.where(np.isnan(p), 1.0, p)
    sig = pd.DataFrame({"log2FC": lfc, "t": t, "p": p, "q": bh_fdr(p)}, index=L.index)
    sig.to_csv(os.path.join(OUT, f"signature_{SITE}_GSE249443.csv"))
    print(f"\n[signature] saved. genes={len(sig)}")
    print(f"[signature] q<0.05 & |log2FC|>1 : {((sig['q']<0.05)&(sig['log2FC'].abs()>1)).sum()}")

    # ---- 内部验证：Cxcl13（原论文认定的骨癌痛关键基因）----
    print("\n[validation] 原论文关键基因检查：")
    for g in ["Cxcl13", "Atf3", "Trpv1", "Scn9a", "Ptgs2"]:
        for key in [g, g.lower(), g.upper(), g.capitalize()]:
            if key in sig.index:
                r = sig.loc[key]
                print(f"   {key:10s} log2FC={r['log2FC']:+.3f}  p={r['p']:.4g}  q={r['q']:.4g}")
                break

    print("\n[signature] Top 10 up-regulated in BCP:")
    print(sig.sort_values("log2FC", ascending=False).head(10)[["log2FC", "p", "q"]].to_string())
    print("\n[signature] Top 10 down-regulated in BCP:")
    print(sig.sort_values("log2FC").head(10)[["log2FC", "p", "q"]].to_string())

    # ---------------- 2. 共表达网络 ----------------
    TOP = 3000
    v = L.var(axis=1).sort_values(ascending=False)
    genes = v.index[:TOP]
    Z = L.loc[genes]
    Z = Z[Z.std(axis=1) > 0]
    Zn = Z.sub(Z.mean(axis=1), axis=0).div(Z.std(axis=1), axis=0).dropna()
    corr = pd.DataFrame(np.corrcoef(Zn.values), index=Zn.index, columns=Zn.index)
    print(f"\n[network] co-expression: {corr.shape[0]} genes")

    # ---------------- 3. 虚拟扰动 + 逆转评分 ----------------
    def reversal(gene, direction):
        """R = -corr(Δ, pain signature)；Δ 为虚拟扰动后的预测变化。"""
        common = corr.columns.intersection(sig.index)
        d = corr.loc[gene, common].astype(float)
        s = sig.loc[common, "log2FC"].astype(float)
        ok = d.notna() & s.notna() & np.isfinite(d) & np.isfinite(s)
        if ok.sum() < 100:
            return np.nan
        d, s = d[ok], s[ok]
        if d.std() == 0 or s.std() == 0:
            return np.nan
        delta = -d if direction == "KO" else d
        return -float(np.corrcoef(delta.values, s.values)[0, 1])

    lookup = {str(g).upper(): g for g in corr.index}
    rows = []
    for tname, (direction, desc) in TARGETS.items():
        real = lookup.get(tname.upper())
        if real is None:
            rows.append({"target": tname, "direction": direction, "desc": desc,
                         "R": np.nan, "in_network": False})
        else:
            rows.append({"target": tname, "direction": direction, "desc": desc,
                         "R": reversal(real, direction), "in_network": True})
    res = pd.DataFrame(rows)

    # ---------------- 4. 置换检验 ----------------
    found = {lookup[t.upper()] for t in TARGETS if t.upper() in lookup}
    null_genes = [g for g in corr.index if g not in found]
    rng = np.random.default_rng(SEED)
    nulls = []
    for direction in ("KO", "OE"):
        picks = rng.choice(null_genes, size=min(N_PERM // 2, len(null_genes)), replace=False)
        for g in picks:
            r = reversal(g, direction)
            if not np.isnan(r):
                nulls.append(r)
    nulls = np.array(nulls)
    mu, sd = nulls.mean(), nulls.std(ddof=1)
    print(f"\n[permutation] null n={len(nulls)} mean={mu:+.4f} sd={sd:.4f}")

    for i, r in res.iterrows():
        if not r["in_network"] or pd.isna(r["R"]):
            continue
        res.at[i, "z"] = (r["R"] - mu) / sd if sd > 0 else np.nan
        res.at[i, "p_perm"] = (np.sum(nulls >= r["R"]) + 1) / (len(nulls) + 1)
    res["q_perm"] = bh_fdr(res["p_perm"].values)
    res["site"] = SITE
    res = res.sort_values("R", ascending=False)

    out = os.path.join(OUT, f"perturbation_{SITE}.csv")
    res.to_csv(out, index=False)
    print(f"\n===== 虚拟扰动结果（骨癌痛，R 越高=越能逆转疼痛状态）=====")
    print(res[["target", "direction", "desc", "R", "z", "p_perm", "q_perm"]].to_string(index=False))
    print(f"\n[saved] {out}")

    L.to_csv(os.path.join(PROC, f"expr_log2_{SITE}.csv"))
    print(f"[saved] processed expression -> {PROC}")


if __name__ == "__main__":
    main()
