# Round 3 独立评审汇总裁定 — 2026-10-10

**稿件：** `docs/manuscript/manuscript_final.md`（单作者计算再分析）
**评审模式：** 5 位独立专家，互不可见彼此意见、不可见 Round 1/2 评审历史、以首次投稿对待。
**共享 brief：** `docs/review/round3_2026-10-10/_PANEL_BRIEF.md`
**独立报告：** `A1_cancer_pain_neuroimmunology.md`、`A2_biostatistician.md`、
`A3_singlecell_methodology.md`、`A4_provenance_auditor.md`、`A5_venue_editor.md`

---

## 1. 汇总裁定

| 审稿人 | 角色 | 裁决 |
|---|---|---|
| R3-A1 | 疼痛神经免疫学 | Minor Revision |
| R3-A2 | 生物统计 / 多重检验 / 效力 | Minor Revision |
| R3-A3 | 单细胞 / pseudobulk 方法学 | Minor Revision |
| R3-A4 | 溯源 / 可重复性审计 | Minor Revision |
| R3-A5 | 期刊编辑 / 投稿就绪度 | Minor Revision |

**综合裁决：Minor Revision（5/5）。** 无 Major、无 Reject。两位审稿人独立确认阴性结果层
（T1/T2/T5）的全部头条数字与方向正确，8-case 交换 null 实现无索引/方向 bug，两个真分配
自检在运行时强制执行且通过，KO/OE 去重逻辑稳固（r = −1.000000）。同时发现两个此前
gate 未能捕捉的真实统计问题（配对 p 地板、Cohen's d 公式），以及一个事实性引文错误。

---

## 2. 编辑对每条意见的核实与处置

### 2.1 统计（A2/A4 独立命中，编辑亲验）

| # | 意见 | 编辑核实 | 判定 | 处置 |
|---|---|---|---|---|
| S1 | 配对对比的精确 p 地板写 3/71，但 reversal 端点真实地板为 **2/71** | **属实**：`etiology_bootstrap_differences.csv` 中 6 个 paired 行恰取 p = 2/71 = 0.028169。代码根源：配对 null 对 8 只 case 动物穷举 C(8,4)=70 拆分，统计量 `ra(s1) − rb(~s1)` 因控制均值固定，补 split **不构成精确镜像**，故 ge≥1 → p=(1+1)/71=2/71。marker logFC 端点反号对称，补集即镜像 → 3/71 | 真问题 | §2.7、Limitations 改为端点区分；gate 字面量补 "2/71" |
| S2 | T5 效力推导框架错配：d=5.653（t 几何 df=2）不应乘 bootstrap SE | **属实**：连贯正态近似 MDE = (1.96+0.8416)×SE = 2.802×0.0185 = **0.0518**（原 0.105 为其两倍）。且"failure cannot be attributed to power"越界——中位 target \|Δ\|=0.0313 低于 MDE，仅 41/102（保守阈值 23/102）超过 | 真问题 | §3.6 重写：给出正态近似 MDE 0.052 + t 参考 0.10，明示中位 target 低于阈值，失败 = 反选择性 + 中位 target 功效地板 |
| S3 | Cohen's d 公式错误：script 33 用 SS/7 而非 SS/(n₁+n₂−2)=SS/6 | **属实**：pooled_sd 用每组 SS/7；正确 d = 1.267（原 1.37，偏大 √(7/6)=1.080） | 真问题 | script 33 改 /6；JSON/正文/Table 5 改 1.27；gate 期望改 1.267 |
| S4 | randset q 值被平结主导（地板 2/501，24/132 在板上） | 属实（次要）：经验 p 地板 0.003992，边界依赖平结 | 部分 | 保留 500 抽样并明示为近似；结论方向不受影响（T1/T2 已在 dedup 上重算） |
| S5 | gate 约 175/287 项为硬编码字面量、33 项为同义反复 in_text | 属实 | 真问题 | T1/T2/T5/Table2/3 已改运行时派生；本次新增 DOI+图件存在性检查；残留字面量在 MEMORY 标注为已知边界 |

### 2.2 单细胞（A3/A1/A4）

| # | 意见 | 编辑核实 | 判定 | 处置 |
|---|---|---|---|---|
| C1 | "seven testable compartments" 应为 **8**（12−4=8，含 Erythrocyte，2v4 不平衡未披露） | **属实**：`sc_GSE292522_celltype_DE.csv` 含 8 腔室 | 真问题 | §3.8、Table 5 改 8，并披露红细胞 bin 2v4 不平衡 |
| C2 | "11.9% neuron purity" 无产物可溯源；"non-overlapping criterion" 不成立 | **属实**：Neuron 标记集本身含 Snap25/Syp/Tubb3；按不同判据得 3–35%，11.9% 不可复现 | 真问题 | 删除该数字与"独立非重叠判据"说法，改为"标记面板重叠、无外部独立判据" |
| C3 | "no detectable Gfap (logFC=0.000)" 用差分零证不表达属逻辑错误 | **属实**：SG bin 8 只合计仅 **7 个 Gfap UMI，5/8 动物为零**（亲验） | 真问题 | 改述为表达层面事实（7 UMI / 5 只为零） |
| C4 | 成分位移未报告：Neuron 78.6%→67.0%、Schwann 18.9%→29.3% | **属实**（亲验，均 MWU p=0.0286） | 真问题 | §3.8 补入成分位移 |
| C5 | §2.8 承诺的 leave-one-animal-out refit 在 §3.8 缺失 | 属实 | 真问题 | 见 §2.2 补做（`37_compartment_loo.py`） |
| C6 | B16F10 接种部位/时间点/对照组未说明；CD45⁺ 去除未披露 | **属实**：GSE292522 = L3–L5 DRG、皮内接种 14 天、MPEK-BL6 角质细胞对照、测序前 MagniSort 去 CD45⁺ | 真问题 | §2.8/§3.8/Discussion 补全；CD45 去除列为髓系稀缺的机制性原因 |
| C7 | 55,967 描述为"QC retained"但 script 31 无细胞 QC | 属实 | 真问题 | 改"as deposited"并可引用 depositor QC（>100 genes、<25% mt） |

### 2.3 溯源与规范（A4/A5）

| # | 意见 | 编辑核实 | 判定 | 处置 |
|---|---|---|---|---|
| P1 | ref 4 "Falk & Dickenson, BJA 2015;115(2):169–171" **不存在** | **属实**（联网核实）：真实文献为 Falk & Dickenson, *J Clin Oncol* 2014;32(16):1647–1654 | 真问题（高危） | 改为真实文献 + DOI |
| P2 | Data availability URL 返回 404 却声称"are deposited" | 属实 | 真问题 | 保留占位仓库 URL 并明示"to be archived at Zenodo upon acceptance"；作者须建仓（唯一样本外遗留） |
| P3 | 两套 logFC（de.csv 肌注校正 vs cross.csv 未校正）混用未声明 | 属实 | 真问题 | §2.3 新增披露段，说明端点各自来源及肌肉协变量不可入置换枚举的局限 |
| P4 | 摘要 457 词超限 | 属实 | 真问题 | 压至 **295 词** |
| P5 | 无图件，desk-reject 风险 | 属实 | 真问题 | 新增 Fig 1（设计示意）、Fig 2（marker 分歧）、Fig 3（反选择性），含图注与内联调用 |
| P6 | ref 10 作者列表被截断；DOI 格式不均（4/16） | 属实 | 真问题 | ref 10 补全 7 位作者；全部 16 条补 DOI（Crossref 核验） |
| P7 | 表脚注缺失；两段 >350 词 | 属实 | 部分 | 补 Table 2–5 脚注；拆分两段过长段落 |

### 2.4 生物学（A1）

| # | 意见 | 编辑核实 | 判定 | 处置 |
|---|---|---|---|---|
| B1 | §3.2 median \|logFC\| 0.073/0.158/0.119 实为全基因中位数，被误归给 q<0.05 集 | **属实**：显著集实际中位数 D7=0.340、D14=1.199 | 真问题 | 明确为"across all tested genes"，并补显著集中位数 |
| B2 | marker 端点缺 LOO | 属实 | 真问题 | 补做（`36_marker_loo.py`）：Cd68 [+3.28,+3.52]、Aif1 [+2.88,+3.21]、Gal [−2.69,−2.39]，全部符号稳定 |
| B3 | T3 gene-level MDE=0.001 被近零基因污染 | 属实 | 真问题 | script 33 限制到表达基因（\|logFC\|≥0.2），MDE 改 0.488 |
| B4 | 反选择性在剔除 IFN 模块后减弱（p 0.0019→0.019） | 属实（亲验） | 部分 | §3.6 增加分层陈述，明示反选择性不单由 IFN 控制承载 |

---

## 3. 修订后状态

- **98 号 gate：287 checks, 0 failed**（含运行时派生的 T1/T2/T5、Table 2/3 逐格再生、DOI 完整性、图件存在性）。
- 摘要 295 词；正文（至 References）约 6,250 词；14 个新产物（3 图 ×2 格式、3 个 LOO/核验 JSON）。
- 新增脚本：35（图件）、36（marker LOO）、37（compartment LOO）。

## 4. 遗留

- Data availability：作者须建立 `yyx-4113/cancer-pain-drg-virtual-perturbation` 仓库并 Zenodo 归档（唯一样本外遗留）。
- ARRIVE 性别报告项、目标期刊特定字数/关键词复核待投稿前按具体刊物处理。
