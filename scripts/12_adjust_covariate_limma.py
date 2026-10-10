"""
12_adjust_covariate_limma.py
============================
识别并校正"肌肉污染"技术混杂，再做 limma 差异表达。

问题：骨癌痛数据(GSE249443)未经校正时，Top DEG 全部是骨骼肌/肌节基因
（Tnni2/Acta1/Ckm/Myl1/Myh1...，log2FC -4 ~ -5.4），说明 sBCP 组 DRG
取材混入了更多骨骼肌。这是技术混杂而非癌痛生物学。

做法：
  1. 计算每个样本的肌肉污染评分（肌节/肌肉标志基因均值）
  2. 对表达做 PCA，检验 PC1 是否承载肌肉信号并与分组纠缠
  3. 用 limma 带协变量模型：expr ~ group + covariate
     协变量取 PC1（数据驱动）或肌肉评分（先验驱动），两者都跑并比较
  4. 校正后再看 Cxcl13 / Atf3 等真实痛相关基因是否浮现
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

MUSCLE_GENES = ["Acta1", "Ckm", "Myl1", "Mylpf", "Tnni2", "Tnnt3", "Tnnc2",
                "Actn3", "Myh1", "Myh2", "Myoz1", "Cox6a2", "Sypl2", "Jsrp1",
                "Mb", "Des", "Ttn", "Mybph", "Tnni1", "Myl2", "Actc1", "Myh4"]


def trigamma_inverse(x):
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
    s2 = np.asarray(s2, float)
    ok = np.isfinite(s2) & (s2 > 0)
    z = np.log(s2[ok])
    e = z - digamma(df1 / 2.0) + np.log(df1 / 2.0)
    emean, evar = e.mean(), e.var(ddof=1) - polygamma(1, df1 / 2.0)
    if evar <= 0:
        return float(np.exp(emean)), np.inf
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


def limma_with_covariate(L, group, cov=None):
    """expr ~ group (+ covariate) 的 moderated t（针对 group 系数）。"""
    y = L.values.astype(float)                      # genes x samples
    n = y.shape[1]
    cols = [np.ones(n), group.astype(float)]
    if cov is not None:
        c = np.asarray(cov, float)
        c = (c - c.mean()) / (c.std() if c.std() > 0 else 1.0)
        cols.append(c)
    X = np.column_stack(cols)
    p_rank = X.shape[1]
    XtX_inv = np.linalg.pinv(X.T @ X)
    B = y @ X @ XtX_inv.T                            # genes x p
    resid = y - B @ X.T
    df_res = n - p_rank
    s2 = (resid ** 2).sum(1) / df_res
    s2 = np.where(s2 <= 0, np.nan, s2)
    s02, d0 = fit_fdist(s2[~np.isnan(s2)], df_res)
    if np.isinf(d0):
        s2_post, df_tot = np.full_like(s2, s02), np.full_like(s2, 1e6)
    else:
        s2_post = (d0 * s02 + df_res * s2) / (d0 + df_res)
        df_tot = np.full_like(s2, df_res + d0)
    beta_g = B[:, 1]
    se = np.sqrt(s2_post * XtX_inv[1, 1])
    with np.errstate(invalid="ignore", divide="ignore"):
        t = beta_g / se
    p = 2 * stats.t.sf(np.abs(t), df_tot)
    p = np.where(np.isnan(p), 1.0, p)
    return pd.DataFrame({"log2FC": beta_g, "t": t, "p": p, "q": bh_fdr(p)}, index=L.index), (s02, d0)


def main():
    site = "bone"
    L = pd.read_csv(os.path.join(PROC, f"expr_log2_{site}.csv"), index_col=0)
    ctrl = [c for c in L.columns if "sBCP" in c]
    case = [c for c in L.columns if c not in ctrl and "BCP" in c]
    samples = ctrl + case
    group = np.array([0] * len(ctrl) + [1] * len(case), float)
    print(f"[data] genes={L.shape[0]}  sBCP(ctrl)={ctrl}  BCP(case)={case}")

    # ---- 1. 肌肉污染评分 ----
    mg = [g for g in MUSCLE_GENES if g in L.index]
    mscore = L.loc[mg, samples].mean(0)
    print(f"\n[肌肉污染] 使用 {len(mg)} 个肌肉标志基因")
    print("  逐样本肌肉评分(log2 FPKM+1 均值):")
    for s in samples:
        tag = "CASE(BCP)" if s in case else "CTRL(sBCP)"
        print(f"    {s:14s} {tag:11s} {mscore[s]:.3f}")
    print(f"  ctrl 均值={mscore[ctrl].mean():.3f}  case 均值={mscore[case].mean():.3f}"
          f"  差={mscore[case].mean()-mscore[ctrl].mean():+.3f}")

    # ---- 2. PCA ----
    Xc = L[samples].sub(L[samples].mean(1), axis=0)
    U, S, Vt = np.linalg.svd(Xc.values, full_matrices=False)
    pcs = pd.DataFrame(Vt[:3].T, index=samples, columns=["PC1", "PC2", "PC3"])
    var = (S ** 2 / (S ** 2).sum())[:3]
    print(f"\n[PCA] 解释方差: PC1={var[0]:.1%} PC2={var[1]:.1%} PC3={var[2]:.1%}")
    print(pcs.round(3).to_string())
    # PC1 与肌肉评分、分组的相关
    for pc in ["PC1", "PC2"]:
        r_m = np.corrcoef(pcs[pc], mscore[samples])[0, 1]
        r_g = np.corrcoef(pcs[pc], group)[0, 1]
        print(f"  {pc}: corr(肌肉评分)={r_m:+.3f}  corr(分组)={r_g:+.3f}")

    # ---- 3. 三种模型比较 ----
    results = {}
    res0, (s02_0, d0_0) = limma_with_covariate(L[samples], group, cov=None)
    results["uncorrected"] = res0
    print(f"\n[model: 未校正]        s0^2={s02_0:.4f} d0={d0_0:.2f}  q<0.05: {(res0['q']<0.05).sum()}")

    resPC1, (s02_1, d0_1) = limma_with_covariate(L[samples], group, cov=pcs["PC1"].values)
    results["adj_PC1"] = resPC1
    print(f"[model: 校正PC1]       s0^2={s02_1:.4f} d0={d0_1:.2f}  q<0.05: {(resPC1['q']<0.05).sum()}")

    resMS, (s02_2, d0_2) = limma_with_covariate(L[samples], group, cov=mscore[samples].values)
    results["adj_muscle"] = resMS
    print(f"[model: 校正肌肉评分]  s0^2={s02_2:.4f} d0={d0_2:.2f}  q<0.05: {(resMS['q']<0.05).sum()}")

    # ---- 4. 关键基因（痛相关）在校正前后的变化 ----
    key = ["Cxcl13", "Atf3", "Npy", "Scn9a", "Trpv1", "Ptgs2", "Oprm1", "Cacna2d1", "Trpm3", "Gal", "Nts"]
    print("\n[痛相关关键基因] 校正前后对比：")
    print(f"{'gene':10s} " + " ".join(f"{k:>22s}" for k in results))
    print(f"{'':10s} " + " ".join(f"{'log2FC/p/q':>22s}" for k in results))
    for g in key:
        if g not in L.index:
            continue
        line = f"{g:10s} "
        for k, r in results.items():
            v = r.loc[g]
            line += f"{v['log2FC']:+.2f}/{v['p']:.3f}/{v['q']:.3f} ".rjust(23)
        print(line)

    # ---- 5. 校正后 top DEG ----
    best = results["adj_PC1"]
    print("\n[校正PC1 后] Top 15 up in BCP:")
    print(best.sort_values("p").head(15)[["log2FC", "p", "q"]].to_string())
    print("\n[校正PC1 后] Top 15 down in BCP:")
    print(best.sort_values("p").tail(15)[["log2FC", "p", "q"]].iloc[::-1].to_string())

    for k, r in results.items():
        r.to_csv(os.path.join(OUT, f"signature_{site}_{k}.csv"))
    print("\n[saved] signature_bone_{uncorrected,adj_PC1,adj_muscle}.csv")


if __name__ == "__main__":
    main()
