# Round 4 独立评审汇总裁定 — 2026-10-10

**稿件：** `docs/manuscript/manuscript_final.md`
**评审模式：** 5 位独立专家，互不可见、不可见任何既往评审，以首次投稿对待。
**共享 brief：** `docs/review/round4_2026-10-10/_PANEL_BRIEF.md`
**独立报告：** `A1_cancer_pain_neuroimmunology.md`、`A2_biostatistician.md`、
`A3_singlecell_methodology.md`、`A4_provenance_auditor.md`、`A5_venue_editor.md`

---

## 1. 汇总裁定

| 审稿人 | 角色 | 裁决 |
|---|---|---|
| R4-A1 | 疼痛神经免疫学 | Minor Revision |
| R4-A2 | 生物统计 / 多重检验 / 效力 | Minor Revision |
| R4-A3 | 单细胞 / pseudobulk 方法学 | Minor Revision |
| R4-A4 | 溯源 / 可重复性审计 | **Major Revision** |
| R4-A5 | 期刊编辑 / 投稿就绪度 | Minor Revision |

**综合裁决：Minor Revision（4/5）+ 1 Major（样本外）。**
A4 的唯一 Major 是 Data availability 占位仓库 URL 返回 404——这是**作者侧动作**（建仓），
不是稿件缺陷；三条审稿人独立命中同一项。除此之外无任何科学层面 Major：三位审稿人
独立重算了 T1/T2/T4/T5、精确 p 地板、Cohen's d、MDE、去重镜像与全部 Table 2/3/5
单元格，**零不可溯源、零矛盾数字**。

---

## 2. 编辑核实与处置（已全部落地）

### 2.1 真实问题（亲验后修复）

| # | 意见（来源） | 编辑核实 | 处置 |
|---|---|---|---|
| P1 | 摘要超限：A1 记 325、A5 记 325（Round-3 曾称 295） | **属实**：实测正文 310 / 含标题 314 词 | 压至 **295 正文 / 299 含标题**；gate 新增 ≤300 运行时断言 |
| P2 | §3.6 "55.7th percentile" 方向错 | **属实**：script 33 的 `pct_abs=(rank_from_top−1)/(n−1)` 是"距顶距离"非百分位；标准百分位 = **44.7**（靶点中位 \|z\| 1.0729 低于总体中位 1.3248） | 改为"just below the median (44.7th percentile)"；script 33 T4 增 `median_percentile_absz_below` 字段；gate 增运行时派生断言 |
| P3 | 红细胞腔室组数写反（A3） | **属实**：`sc_GSE292522_celltype_scores.csv` 显示红细胞 bin = 2 只 Con(S1,S3) + 4 只 B16F10(S5–S8)；正因仅 2 只对照，script 32 的 ISG 判别（需≥3/组）将其排除 | §3.8 改 "two control and four tumour" |
| P4 | "log₂(CPM+1)" 标签错（A3） | **属实**：script 31/32 均为 `×1e4`（counts-per-10k），非 CPM | §2.8 两处改 "log₂(counts per 10,000 + 1)"；修正 script 31 注释 |
| P5 | Table 4/5 正文未引用（A5） | **属实**：仅出现在自身表头 | §3.4 加 "(Table 4)"、§3.8 加 "(Table 5)" |
| P6 | ref [14] 引用乱序（A5） | **属实**：首次出现于 [15],[16] 之后 | 重编号：Baruch→[14]、Usoskin→[15]、Luger→[16]；正文首次出现顺序现为 1..16 |
| P7 | Fig 1 图注称"所有数值由 results/tables 再生"（A4/A5） | **属实**：Fig 1 是示意框图 | 改为"schematic；样本数依沉积设计" |
| P8 | §3.8 髓系占比下限 0.02% 偏低（A1） | **属实**：实测 0.0147–0.199% | 改 "0.015–0.2%" |
| P9 | Data availability 现在时"are deposited"与 404 仓库不符（A1/A4/A5） | **属实** | 改将来时"will be made available … and archived at Zenodo upon acceptance"；**建仓为作者侧唯一遗留** |
| P10 | §3.5 Atf3 用肌校正值(+1.13/+0.22) 与本节 cross 口径不一致（A1） | **属实**：cross.csv 为 +1.12/+0.31 | 改 +1.12/+0.31 |
| P11 | T1/T2 视 132 检验为独立、实则共享对比（A2） | 属实（次要） | §3.6 补非独立性提示（p 为近似，方向不受影响） |
| P12 | 无独立 Ethics 声明（A5） | 属实 | 新增 `## Ethics statement` |

### 2.2 审稿人误判（编辑不采纳）

| # | 意见 | 判定 |
|---|---|---|
| X1 | A2：artefact 键 `myeloid_MWU_p_less` 命名错误 | **误判**：script 31 用 `mannwhitneyu(B16F10, Con, alternative="less")`，键名恰为 B16F10<Con，正确 |
| X2 | A1：摘要 325 词 | 计数方式差异；实测 310/314，仍超限（已修），绝对数取编辑实测 |

### 2.3 已确认无误（保留）

- 配对地板端点依赖（logFC 3/71 vs reversal 2/71）由 A2 从原始数据独立重构确认真实（obs+comp 分别为 0.00e+00 与 +0.00745）。
- KO/OE 镜像 r = −1.000000、max\|Δp\| = 0；T1 中位 \|z\| 2.4702/1.0729、MWU p = 0.001914；T2 14/32、基率 30/132、binomial p = 0.006708；T5 MDE 0.0518、max Δ 0.2743、41/102、23/102、31/132；d = 1.2670、MDE-d 2.3808；红系最大 \|z\| 7.388；network-additive ρ = 0.991；跨模型 ρ = 0.0985。全部匹配。

---

## 3. 修订后状态

- **98 号 gate：291 checks, 0 failed**（新增摘要 ≤300、T4 百分位运行时派生、"counts per 10,000"、"two control and four tumour"）。
- 摘要 295 词（含标题 299）；正文至 References 约 6,150 词。
- 参考文献 16 条全部按首次出现顺序、全部带 DOI（Crossref 核验）。
- 图件 fig1/2/3（png+pdf）；fig3 标题 p 值改为运行时从 JSON 读取。
- 新增 `## Ethics statement`。

## 4. 遗留（均为样本外 / 期刊特定）

- **Data availability 仓库**：作者须建立 `yyx-4113/cancer-pain-drg-virtual-perturbation` 并（录用后）Zenodo 归档——唯一阻塞投稿的样本外项。
- 标题 202 字符偏长、具体期刊字数/关键词上限——投稿前按目标刊复核。
