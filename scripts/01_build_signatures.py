"""
01_build_signatures.py
======================
从 GEO series matrix 构建"部位分层癌痛分子签名"。

设计要点：
  - 统一解析 GEO series_matrix.txt.gz（无需额外依赖）
  - 自动按关键词把样本分为 control vs case（疼痛/荷瘤）
  - 线性信号自动 log2 转换
  - 逐基因 t 检验 + BH-FDR，输出 log2FC / p / q
  - 每个部位单独产出签名 CSV，供后续虚拟扰动使用

关键设计决策（避免混杂）：
  为保证"不同部位"可比，比较限定在同一组织（感觉神经节/DRG）内进行：
    骨癌痛  -> GSE249443 (rat DRG, BCP vs sham)
    胰腺癌痛 -> GSE292089 (mouse ganglia, tumor vs control)
  不同组织/物种的数据集不直接互相减，只用各自内部对照做 DE。
"""
import gzip
import io
import os
import sys
import numpy as np
import pandas as pd
from scipy import stats

RAW = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data", "raw")
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results", "tables")
os.makedirs(OUT, exist_ok=True)

# 部位 -> 数据集 配置（组织、物种、对照/病例关键词）
DATASETS = [
    {
        "site": "bone",
        "gse": "GSE249443",
        "tissue": "DRG",
        "species": "rat",
        "ctrl_keys": ["sham", "control", "naive", "sbcpsham"],
        "case_keys": ["bcp", "cancer", "tumor", "tumour", "model"],
    },
    {
        "site": "visceral_pancreas",
        "gse": "GSE292089",
        "tissue": "ganglia",
        "species": "mouse",
        "ctrl_keys": ["naive", "control", "sham", "non-tumor", "nontumor", "untreated"],
        "case_keys": ["tumor", "tumour", "cancer", "pdac", "pni", "injured"],
    },
]

CTRL_KEYS_FALLBACK = ["sham", "control", "naive", "normal", "non-tumor", "nontumor"]
CASE_KEYS_FALLBACK = ["bcp", "cancer", "tumor", "tumour", "pdac", "pni", "model"]


def parse_series_matrix(path):
    """解析 GEO series matrix，返回 (expr_df, meta_df)。"""
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="ignore") as fh:
        lines = fh.readlines()

    meta_rows = {}
    table_start = None
    for i, ln in enumerate(lines):
        if ln.startswith("!Sample_title"):
            meta_rows["title"] = ln.strip().split("\t")[1:]
        elif ln.startswith("!Sample_geo_accession"):
            meta_rows["gsm"] = ln.strip().split("\t")[1:]
        elif ln.startswith("!Sample_characteristics_ch1"):
            # 可能多行，拼接
            vals = ln.strip().split("\t")[1:]
            if "characteristics" in meta_rows:
                meta_rows["characteristics"] = [
                    a + " | " + b for a, b in zip(meta_rows["characteristics"], vals)
                ]
            else:
                meta_rows["characteristics"] = vals
        elif ln.startswith("!Sample_source_name_ch1"):
            vals = ln.strip().split("\t")[1:]
            meta_rows.setdefault("source", vals)
        elif ln.startswith("!series_matrix_table_begin"):
            table_start = i + 1
            break

    if table_start is None:
        raise RuntimeError("no matrix table found")

    with io.StringIO("".join(lines[table_start:])) as buf:
        buf.seek(0)
        expr = pd.read_csv(buf, sep="\t", comment="!", index_col=0)
    # 去掉可能的结尾标记行
    expr = expr[~expr.index.astype(str).str.contains("series_matrix_table_end", na=False)]
    expr = expr.apply(pd.to_numeric, errors="coerce").dropna(how="all")

    meta = pd.DataFrame(meta_rows)
    if not meta.empty and "gsm" in meta:
        meta.index = meta["gsm"]
    return expr, meta


def sample_text(meta, gsm):
    """把样本所有元数据拼成小写字符串用于关键词匹配。"""
    if meta is None or gsm not in meta.index:
        return gsm.lower()
    row = meta.loc[gsm]
    return " ".join(str(v) for v in row.values).lower()


def assign_groups(expr, meta, ctrl_keys, case_keys):
    """按关键词把列分为 control / case，返回两个列名列表。"""
    ctrl, case = [], []
    for col in expr.columns:
        txt = sample_text(meta, col)
        if any(k in txt for k in case_keys):
            case.append(col)
        elif any(k in txt for k in ctrl_keys):
            ctrl.append(col)
    return ctrl, case


def maybe_log2(expr):
    """若数据明显是线性强度，则 log2 转换。"""
    mx = np.nanmax(expr.values)
    if mx > 50:
        expr = np.log2(np.clip(expr, 0, None) + 1)
    return expr


def bh_fdr(p):
    p = np.asarray(p, dtype=float)
    ok = ~np.isnan(p)
    q = np.full_like(p, np.nan)
    pv = p[ok]
    n = len(pv)
    if n == 0:
        return q
    order = np.argsort(pv)
    ranked = pv[order]
    adj = ranked * n / (np.arange(1, n + 1))
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(adj, 0, 1)
    q[ok] = out
    return q


def de_analysis(expr, ctrl, case):
    """逐基因 t 检验，返回 log2FC / p / q。"""
    a = expr[ctrl].values.astype(float)
    b = expr[case].values.astype(float)
    lfc = np.nanmean(b, axis=1) - np.nanmean(a, axis=1)
    with np.errstate(invalid="ignore"):
        t, p = stats.ttest_ind(b, a, axis=1, equal_var=False, nan_policy="omit")
    p = np.where(np.isnan(p), 1.0, p)
    return pd.DataFrame({"log2FC": lfc, "t": t, "p": p, "q": bh_fdr(p)}, index=expr.index)


def main():
    summary = []
    for cfg in DATASETS:
        path = os.path.join(RAW, f"{cfg['gse']}_series_matrix.txt.gz")
        print(f"\n===== {cfg['gse']} [{cfg['site']} / {cfg['tissue']} / {cfg['species']}] =====")
        if not os.path.exists(path):
            print(f"  [SKIP] 未找到 {path}（可能仍在下载）")
            summary.append({**cfg, "status": "missing"})
            continue
        try:
            expr, meta = parse_series_matrix(path)
        except Exception as e:
            print(f"  [ERROR] 解析失败: {e}")
            summary.append({**cfg, "status": f"parse_error: {e}"})
            continue

        expr = maybe_log2(expr)
        print(f"  genes={expr.shape[0]}  samples={expr.shape[1]}")
        print("  samples:")
        for c in expr.columns:
            print(f"    - {c}: {sample_text(meta, c)[:110]}")

        ctrl, case = assign_groups(
            expr, meta,
            cfg.get("ctrl_keys", CTRL_KEYS_FALLBACK),
            cfg.get("case_keys", CASE_KEYS_FALLBACK),
        )
        print(f"  control n={len(ctrl)}  case n={len(case)}")
        if len(ctrl) < 2 or len(case) < 2:
            print("  [WARN] 每组样本<2，无法做统计检验，需人工指定分组")
            summary.append({**cfg, "status": "group_fail", "n_ctrl": len(ctrl), "n_case": len(case)})
            continue

        res = de_analysis(expr, ctrl, case)
        sig = res[(res["q"] < 0.05) & (res["log2FC"].abs() > 0.5)]
        out = os.path.join(OUT, f"signature_{cfg['site']}_{cfg['gse']}.csv")
        res.sort_values("p").to_csv(out)
        print(f"  DE genes (q<0.05 & |log2FC|>0.5): {len(sig)}")
        print("  top up / down:")
        for g in sig.sort_values("log2FC").index[:3]:
            print(f"    DOWN {g}: log2FC={res.loc[g,'log2FC']:+.3f} q={res.loc[g,'q']:.3g}")
        for g in sig.sort_values("log2FC").index[-3:]:
            print(f"    UP   {g}: log2FC={res.loc[g,'log2FC']:+.3f} q={res.loc[g,'q']:.3g}")
        print(f"  saved -> {out}")
        summary.append({**cfg, "status": "ok", "n_ctrl": len(ctrl), "n_case": len(case),
                        "n_genes": expr.shape[0], "n_de": len(sig)})

    print("\n===== SUMMARY =====")
    print(pd.DataFrame(summary).to_string(index=False))


if __name__ == "__main__":
    main()
