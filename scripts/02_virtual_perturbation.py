"""
02_virtual_perturbation.py
==========================
虚拟敲除 / 虚拟过表达 + 癌痛状态逆转评分 + 置换检验。

核心逻辑（可解释、可复现、带统计推断）：
  1. 从每个部位的表达矩阵学习基因共表达网络（Pearson，top variable genes）。
  2. 对给定靶点 T 做"虚拟实验"：
       - 虚拟敲除(KO)：T 表达下降 1 SD -> 下游基因预测变化 Δg = -corr(T,g)
       - 虚拟过表达(OE)：                     Δg = +corr(T,g)
     这里 corr 是数据自身学到的结构关系（线性结构传播），
     不依赖任何外部预训练模型，因此不存在"域外外推"的黑箱风险。
  3. 逆转评分 R(T) = -corr(Δ, 疼痛签名 log2FC)
     含义：扰动引起的变化若与"疼痛相对对照"的表达偏移方向相反，
           即为把细胞状态往正常方向拉回 -> 镇痛潜力越大，R 越高。
  4. 置换检验：随机抽取非靶点基因作为伪靶点，构建 R 的零分布，
     给出每个真实靶点的 z 值与经验 p 值、BH-FDR q 值。
     （这一步是防"模型打不过简单基线"批评的关键：
       真实靶点必须显著优于随机基因的零分布才算数。）

输出：results/tables/perturbation_{site}.csv  （靶点 × R/z/p/q）
"""
import os
import sys
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RAW = os.path.join(ROOT, "data", "raw")
OUT = os.path.join(ROOT, "results", "tables")
os.makedirs(OUT, exist_ok=True)

sys.path.insert(0, HERE)
from importlib.machinery import SourceFileLoader
s01 = SourceFileLoader("s01", os.path.join(HERE, "01_build_signatures.py")).load_module()

# ---------------------------------------------------------------- 靶点定义
# direction: "KO" = 该靶点被药物抑制/阻断时为治疗性方向（阿片类为激动->用 OE 模拟）
TARGETS = {
    "OPRM1":     ("OE", "opioid receptor mu (阿片类激动)"),
    "SCN9A":     ("KO", "Nav1.7 (钠通道阻滞剂)"),
    "SCN10A":    ("KO", "Nav1.8 (钠通道阻滞剂)"),
    "SCN11A":    ("KO", "Nav1.9 (钠通道阻滞剂)"),
    "TRPV1":     ("KO", "TRPV1 拮抗剂"),
    "TRPM3":     ("KO", "TRPM3 拮抗剂"),
    "PTGS2":     ("KO", "COX-2 (NSAIDs)"),
    "PTGS1":     ("KO", "COX-1 (NSAIDs)"),
    "CACNA2D1":  ("KO", "alpha2delta-1 (加巴喷丁类)"),
    "CACNA2D2":  ("KO", "alpha2delta-2 (加巴喷丁类)"),
    "GRIN2B":    ("KO", "NMDA 2B (氯胺酮)"),
    "SLC6A4":    ("KO", "SERT (度洛西汀等 SNRI)"),
    "SLC6A2":    ("KO", "NET (度洛西汀等 SNRI)"),
}

DATASETS = s01.DATASETS
N_PERM = 2000
RANDOM_SEED = 20261009
TOP_VAR_GENES = 3000


def build_coexpression(expr):
    """在 top variable genes 上计算 Pearson 相关矩阵。"""
    v = expr.var(axis=1).sort_values(ascending=False)
    genes = v.index[:TOP_VAR_GENES]
    sub = expr.loc[genes]
    sub = sub.loc[:, ~sub.columns.duplicated()]
    # 去掉全零/常量行
    sub = sub[sub.std(axis=1) > 0]
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1).replace(0, np.nan), axis=0)
    z = z.dropna()
    corr = pd.DataFrame(np.corrcoef(z.values), index=z.index, columns=z.index)
    return corr


def bh_fdr(p):
    p = np.asarray(p, float)
    ok = ~np.isnan(p)
    q = np.full_like(p, np.nan)
    pv = p[ok]
    n = len(pv)
    if n == 0:
        return q
    o = np.argsort(pv)
    adj = pv[o] * n / (np.arange(1, n + 1))
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.clip(adj, 0, 1)
    q[ok] = out
    return q


def reversal_score(corr_row, sig, direction):
    """R = -corr(Δ, pain signature)。Δ 由 direction 决定符号。"""
    common = corr_row.index.intersection(sig.index)
    if len(common) < 50:
        return np.nan
    d = corr_row.loc[common].astype(float)
    s = sig.loc[common, "log2FC"].astype(float)
    ok = d.notna() & s.notna()
    if ok.sum() < 50:
        return np.nan
    d, s = d[ok], s[ok]
    delta = -d if direction == "KO" else d
    # R = -Pearson(delta, pain_signature_log2FC)
    with np.errstate(invalid="ignore"):
        r = np.corrcoef(delta.values, s.values)[0, 1]
    return -r if not np.isnan(r) else np.nan


def main():
    rng = np.random.default_rng(RANDOM_SEED)
    all_out = []

    for cfg in DATASETS:
        gse, site = cfg["gse"], cfg["site"]
        path = os.path.join(RAW, f"{gse}_series_matrix.txt.gz")
        sig_path = os.path.join(OUT, f"signature_{site}_{gse}.csv")
        print(f"\n===== {gse} [{site}] =====")
        if not os.path.exists(path) or not os.path.exists(sig_path):
            print("  [SKIP] 缺少表达矩阵或签名，先跑 01_build_signatures.py")
            continue

        expr, meta = s01.parse_series_matrix(path)
        expr = s01.maybe_log2(expr)
        sig = pd.read_csv(sig_path, index_col=0)

        corr = build_coexpression(expr)
        print(f"  co-expression network: {corr.shape[0]} genes")

        # 基因名大小写不敏感映射（rat/mouse 常为 Scn9a 形式）
        lookup = {}
        for g in corr.index:
            lookup.setdefault(str(g).upper(), g)

        # ---- 真实靶点 ----
        rows = []
        for tname, (direction, desc) in TARGETS.items():
            real = lookup.get(tname.upper())
            if real is None:
                rows.append({"target": tname, "direction": direction, "desc": desc,
                             "R": np.nan, "z": np.nan, "p": np.nan, "found": False})
                continue
            R = reversal_score(corr.loc[real], sig, direction)
            rows.append({"target": tname, "direction": direction, "desc": desc,
                         "R": R, "found": True})
        res = pd.DataFrame(rows)

        # ---- 置换零分布 ----
        found_names = {lookup[t.upper()] for t in TARGETS if t.upper() in lookup}
        null_genes = [g for g in corr.index if g not in found_names]
        if len(null_genes) > 50:
            nulls = []
            for direction in ("KO", "OE"):
                for g in rng.choice(null_genes, size=min(N_PERM // 2, len(null_genes)), replace=False):
                    r = reversal_score(corr.loc[g], sig, direction)
                    if not np.isnan(r):
                        nulls.append(r)
            nulls = np.array(nulls)
            mu, sd = nulls.mean(), nulls.std(ddof=1)
            print(f"  permutation null: n={len(nulls)} mean={mu:+.4f} sd={sd:.4f}")
            for i, row in res.iterrows():
                if not row["found"] or pd.isna(row["R"]):
                    continue
                z = (row["R"] - mu) / sd if sd > 0 else np.nan
                p = (np.sum(nulls >= row["R"]) + 1) / (len(nulls) + 1)
                res.at[i, "z"] = z
                res.at[i, "p"] = p
            res["q"] = bh_fdr(res["p"].values)
        else:
            print("  [WARN] 零分布基因不足")

        res["site"] = site
        res["gse"] = gse
        res = res.sort_values("R", ascending=False)
        out = os.path.join(OUT, f"perturbation_{site}.csv")
        res.to_csv(out, index=False)
        print(res[["target", "direction", "R", "z", "p", "q"]].to_string(index=False))
        print(f"  saved -> {out}")
        all_out.append(res)

    if all_out:
        cat = pd.concat(all_out, ignore_index=True)
        cat.to_csv(os.path.join(OUT, "perturbation_all_sites.csv"), index=False)
        print("\n===== ALL SITES (top by R) =====")
        print(cat.dropna(subset=["R"]).sort_values("R", ascending=False)
                 .head(20)[["site", "target", "direction", "R", "z", "q"]].to_string(index=False))


if __name__ == "__main__":
    main()
