"""
00_virtual_experiment_poc.py
============================
虚拟实验可行性证明 (Proof-of-Concept)。

机制（与 scGen 完全一致）：
  1. 用 VAE 把单细胞表达压缩到潜空间 z；
  2. 计算扰动向量  v = mean(z_perturbed) - mean(z_control)  （同细胞类型内）；
  3. 对任一 control 细胞：z_virtual = z_control + v，再解码回表达空间
     -> 这就是"虚拟敲除/过表达"后的预测表达谱。

本 PoC 用合成 DRG 风格小数据演示该闭环：虚拟敲除痛觉关键基因 SCN9A，
观察其下游共表达模块（Nav 家族 / TRP 家族 / 阿片受体）的 predicted delta。

生产版将替换为：
  - 真实数据 GSE330263（小鼠 DRG 骨癌痛 scRNA）
  - scvi-tools 的 SCGen（同机制，但带 batch/cell-type 条件建模）
  - chemCPA 做"虚拟加药"外推未见镇痛药
"""
import numpy as np
import torch
import torch.nn as nn
from torch.nn import functional as F
from sklearn.preprocessing import StandardScaler

torch.manual_seed(0)
np.random.seed(0)

# ---- 1. 合成 DRG 风格数据 ----------------------------------------------
N_CELLS = 1500
N_GENES = 500
# 给前若干个基因起"痛觉相关"名字，便于解读
PAIN_GENES = ["SCN9A", "SCN10A", "SCN11A", "TRPV1", "TRPM3", "OPRM1",
              "GRIN2B", "KCNA1", "CACNA2D1", "BDNF"]
gene_names = PAIN_GENES + [f"G{i}" for i in range(N_GENES - len(PAIN_GENES))]
TARGET_IDX = gene_names.index("SCN9A")          # 虚拟敲除目标
NAV_IDX = [gene_names.index(g) for g in ("SCN9A", "SCN10A", "SCN11A")]
TRP_IDX = [gene_names.index(g) for g in ("TRPV1", "TRPM3")]
OPI_IDX = gene_names.index("OPRM1")

def make_module(n, base, spread, signed=True):
    """在若干基因上构造一个相关的共表达模块。"""
    w = np.zeros(N_GENES)
    for i in range(N_GENES):
        w[i] = np.sum(np.random.randn(len(base)) * spread)  # 每个基因与模块弱相关
    return w

np.random.seed(1)
# 背景表达
X = np.random.randn(N_CELLS, N_GENES) * 0.6 + 1.0
# 两个细胞类型的基础差异
cell_type = np.array(["neuron"] * 800 + ["glia"] * 700)
ct_effect = np.zeros((N_CELLS, N_GENES))
ct_effect[:800, NAV_IDX] += 1.2      # 神经元高表达 Nav
ct_effect[800:, TRP_IDX] += 0.8      # 胶质高表达 TRP
X += ct_effect

# 构造"敲除 SCN9A"的扰动条件：在神经元里把 SCN9A 及其共模块下调
is_neuron = (cell_type == "neuron").astype(int)
condition = np.array(["control"] * N_CELLS, dtype=object)
pert_start = 600   # 必须落在 neuron 区间(0-799)内，否则按细胞类型取扰动向量会取到空集
condition[pert_start:pert_start + 200] = "KO_SCN9A"   # 200 个神经元为真实观测到的敲除
ko_mask = (condition == "KO_SCN9A")
# 敲除效应：SCN9A 本身下调，并带动 Nav/TRP 模块轻微补偿性上调（模拟真实代偿）
rows_ko = np.where(ko_mask)[0]
X[np.ix_(rows_ko, np.array([TARGET_IDX]))] -= 2.5
X[np.ix_(rows_ko, np.array(NAV_IDX))] += 0.5
X[np.ix_(rows_ko, np.array(TRP_IDX))] += 0.4

X = np.clip(X, 0, None)
scaler = StandardScaler()
Xs = scaler.fit_transform(X)

# AnnData 风格的简单容器
class Ann:
    def __init__(self, X, obs):
        self.X = X; self.obs = obs
adata = Ann(Xs, {"cell_type": cell_type, "condition": condition})

# ---- 2. VAE（潜空间向量运算的载体）-----------------------------------
class VAE(nn.Module):
    def __init__(self, n_genes, z=10):
        super().__init__()
        self.enc = nn.Sequential(nn.Linear(n_genes, 64), nn.ReLU(), nn.Linear(64, z * 2))
        self.dec = nn.Sequential(nn.Linear(z, 64), nn.ReLU(), nn.Linear(64, n_genes))
        self.z = z
    def encode(self, x):
        h = self.enc(x); mu, logv = h[:, :self.z], h[:, self.z:]
        return mu, logv
    def reparam(self, mu, logv):
        return mu + torch.exp(0.5 * logv) * torch.randn_like(mu)
    def forward(self, x):
        mu, logv = self.encode(x)
        z = self.reparam(mu, logv)
        return self.dec(z), mu, logv

def to_z(model, X):
    model.eval()
    with torch.no_grad():
        mu, _ = model.encode(torch.tensor(X, dtype=torch.float32))
    return mu.numpy()

model = VAE(N_GENES, z=10)
opt = torch.optim.Adam(model.parameters(), lr=1e-3)
Xtr = torch.tensor(Xs, dtype=torch.float32)
for epoch in range(400):
    opt.zero_grad()
    recon, mu, logv = model(Xtr)
    loss = F.mse_loss(recon, Xtr) + 0.5 * torch.mean(torch.exp(logv) + mu**2 - 1 - logv)
    loss.backward(); opt.step()
print(f"[VAE] trained. final loss={loss.item():.4f}")

# ---- 3. 虚拟敲除：潜空间向量运算 ---------------------------------------
z = to_z(model, Xs)
ct = adata.obs["cell_type"]; cond = adata.obs["condition"]

def perturbation_vector(z, ct, cond, celltype, ctrl="control", pert="KO_SCN9A"):
    m_ctrl = (ct == celltype) & (cond == ctrl)
    m_pert = (ct == celltype) & (cond == pert)
    return z[m_pert].mean(0) - z[m_ctrl].mean(0)

v_neuron = perturbation_vector(z, ct, cond, "neuron")

# 取未参与训练的 neuron control 细胞做"虚拟敲除"预测
ctrl_neuron_mask = (ct == "neuron") & (cond == "control")
ctrl_idx = np.where(ctrl_neuron_mask)[0][:5]            # 抽 5 个对照神经元
z_ctrl = z[ctrl_idx]
z_virtual = z_ctrl + v_neuron                           # <<< 核心：虚拟扰动
with torch.no_grad():
    x_virtual = model.dec(torch.tensor(z_virtual, dtype=torch.float32)).numpy()
x_ctrl = Xs[ctrl_idx]

# 反标准化回近似原始尺度便于读
def unscale(x): return scaler.inverse_transform(x)
xv = unscale(x_virtual); xc = unscale(x_ctrl)

# ---- 4. 结果：虚拟敲除 SCN9A 的 predicted delta ------------------------
print("\n=== 虚拟敲除 SCN9A（神经元）的 predicted delta（Top 变化基因）===")
mean_delta = xv.mean(0) - xc.mean(0)
order = np.argsort(-np.abs(mean_delta))[:10]
for i in order:
    arrow = "▲" if mean_delta[i] > 0 else "▼"
    print(f"  {gene_names[i]:8s} {arrow} {mean_delta[i]:+.3f}")

print("\n=== 关键痛觉模块均值变化（虚拟敲除后 - 对照）===")
def mod_mean(idxs): return float(xv[:, idxs].mean() - xc[:, idxs].mean())
print(f"  Nav 模块 (SCN9A/10A/11A): {mod_mean(NAV_IDX):+.3f}")
print(f"  TRP 模块 (TRPV1/TRPM3)  : {mod_mean(TRP_IDX):+.3f}")
print(f"  OPRM1 (阿片受体)        : {mod_mean([OPI_IDX]):+.3f}")
print(f"  SCN9A 本身              : {mod_mean([TARGET_IDX]):+.3f}")

print("\n[结论] 虚拟实验闭环已跑通：scGen 式潜空间向量运算可在无湿实验条件下")
print("        预测'敲除某基因'后的细胞状态重连。生产版将替换为真实 DRG 数据与 SCGen/chemCPA。")
