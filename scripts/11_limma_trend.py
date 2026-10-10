"""
11_limma_trend.py
=================
limma-trend 经验贝叶斯方差收缩差异表达（Smyth 2004），Python 实现。

为什么必须用它：
  n=3 vs 3 时 Welch t 检验自由度仅 ~2-4，p 值被严重放大（本数据中
  名义 p<0.05 的基因数 382 < 随机期望 898），导致"有信号也检不出"。
  limma 通过跨基因借信息收缩方差，把 d0 个先验自由度加到每个基因上，
  是小样本表达谱的标准做法。

同时输出：
  - 样本相关性 QC（识别离群样本 / 肌肉污染线索）
  - 关键基因检查（Cxcl13、Atf3 等）
"""
import os
import numpy as np
import pandas as pd
from scipy import stats
from scipy.special import digamma, polygamma

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
PROC = os.path.join(ROOT, "data", "processed")
OUT = os.path.join(ROOT, "results", "tables")


def trigamma_inverse(x):
    """解 trigamma(y) = x 的 y（Newton 迭代，Smyth 2004）。"""
    x = np.asarray(x, float)
    y = 0.5 + 1.0 / x
    for _ in range(50):
        tri = polygamma(1, y)
        dif = tri * (1.0 - tri / x) / polygamma(2, y)
        y = y + dif
        if np.max(np.abs(dif / y)) < 1e-8:
            break
    return y


def fit_fdist(s2, df1):
    """矩估计先验尺度 s0^2 与先验自由度 d0。"""
    s2 = np.asarray(s2, float)
    ok = np.isfinite(s2) & (s2 > 0)
    z = np.log(s2[ok])
    e = z - digamma(df1 / 2.0) + np.log(df1 / 2.0)
    emean = e.mean()
    evar = e.var(ddof=1) - polygamma(1, df1 / 2.0)
    if evar <= 0:
        d0 = np.inf
        s02 = np.exp(emean)
    else:
        d0 = 2.0 * float(trigamma_inverse(evar))
        s02 = float(np.exp(emean + digamma(d0 / 2.0) - np.log(d0 / 2.0)))
    return s02, d0


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


def limma_trend_two_group(L, ctrl_cols, case_cols):
    """返回 DataFrame: log2FC / t / p / q （moderated t）。"""
    a = L[ctrl_cols].values.astype(float)
    b = L[case_cols].values.astype(float)
    n1, n2 = a.shape[1], b.shape[1]
    m1, m2 = a.mean(1), b.mean(1)
    lfc = m2 - m1

    # 合并方差
    ss = ((a - m1[:, None]) ** 2).sum(1) + ((b - m2[:, None]) ** 2).sum(1)
    df_res = n1 + n2 - 2
    s2 = ss / df_res
    # 避免 0 方差导致 log(0)
    s2 = np.where(s2 <= 0, np.nan, s2)

    s02, d0 = fit_fdist(s2[~np.isnan(s2)], df_res)
    if np.isinf(d0):
        s2_post = np.full_like(s2, s02)
        df_total = np.full_like(s2, 1e6)
    else:
        s2_post = (d0 * s02 + df_res * s2) / (d0 + df_res)
        df_total = df_res + d0

    se = np.sqrt(s2_post * (1.0 / n1 + 1.0 / n2))
    with np.errstate(invalid="ignore", divide="ignore"):
        t = lfc / se
    p = 2 * stats.t.sf(np.abs(t), df_total)
    p = np.where(np.isnan(p), 1.0, p)
    return pd.DataFrame({"log2FC": lfc, "t": t, "p": p, "q": bh_fdr(p)}, index=L.index), (s02, d0)


def main():
    site = "bone"
    L = pd.read_csv(os.path.join(PROC, f"expr_log2_{site}.csv"), index_col=0)
    ctrl = [c for c in L.columns if "sBCP" in c]
    case = [c for c in L.columns if c not in ctrl and "BCP" in c]
    print(f"[data] genes={L.shape[0]} ctrl={ctrl} case={case}")

    # ---- 样本相关性 QC ----
    C = L.corr(method="pearson")
    print("\n[QC] 样本间 Pearson 相关矩阵：")
    print(C.round(3).to_string())

    res, (s02, d0) = limma_trend_two_group(L, ctrl, case)
    print(f"\n[limma] 先验尺度 s0^2={s02:.4f}  先验自由度 d0={d0:.2f}")
    print(f"[limma] q<0.05            : {(res['q']<0.05).sum()}")
    print(f"[limma] q<0.05 &|log2FC|>1: {((res['q']<0.05)&(res['log2FC'].abs()>1)).sum()}")
    print(f"[limma] p<0.05            : {(res['p']<0.05).sum()}")
    print(f"[limma] p<0.05 &|log2FC|>1: {((res['p']<0.05)&(res['log2FC'].abs()>1)).sum()}")

    print("\n[验证] 原论文关键基因：")
    for g in ["Cxcl13", "Atf3", "Npy", "Scn9a", "Trpv1", "Ptgs2", "Oprm1", "Cacna2d1", "Trpm3"]:
        if g in res.index:
            r = res.loc[g]
            print(f"   {g:9s} log2FC={r['log2FC']:+.3f} p={r['p']:.4g} q={r['q']:.4g}")

    print("\n[limma] Top 15 up in BCP (by p):")
    print(res.sort_values("p").head(15)[["log2FC", "p", "q"]].to_string())

    out = os.path.join(OUT, f"signature_{site}_limma.csv")
    res.to_csv(out)
    print(f"\n[saved] {out}")


if __name__ == "__main__":
    main()
