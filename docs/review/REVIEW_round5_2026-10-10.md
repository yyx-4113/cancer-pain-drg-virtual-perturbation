# Round 5 独立评审汇总裁定 — 2026-10-10

**稿件：** `docs/manuscript/manuscript_final.md`
**评审模式：** 5 位独立专家，互不可见、不可见任何既往评审，以首次投稿对待（确认轮）。
**共享 brief：** `docs/review/round5_2026-10-10/_PANEL_BRIEF.md`
**独立报告：** `A1_cancer_pain_neuroimmunology.md`、`A2_biostatistician.md`、
`A3_singlecell_methodology.md`、`A4_provenance_auditor.md`、`A5_venue_editor.md`

---

## 1. 汇总裁定

| 审稿人 | 角色 | 裁决 |
|---|---|---|
| R5-A1 | 疼痛神经免疫学 | **Accept** |
| R5-A2 | 生物统计 / 多重检验 / 效力 | Minor Revision |
| R5-A3 | 单细胞 / pseudobulk 方法学 | Minor Revision |
| R5-A4 | 溯源 / 可重复性审计 | Minor Revision |
| R5-A5 | 期刊编辑 / 投稿就绪度 | Minor Revision |

**综合裁决：Accept（1/5）+ Minor Revision（4/5）。无 Major、无 Reject。**
A1 直接给出 Accept，无阻断性问题。其余四条 Minor 全部为一致性/格式层面，**无任何科学瑕疵**：
全部核心统计量（精确 p 地板端点依赖、T1/T2/T4/T5、Cohen's d、MDE、去重镜像、Table 3/4/5
逐格）被 5 位审稿人独立重算并复现；A4 溯源 38 个数字，**0 不可溯源、0 矛盾**；
16 条参考文献首现顺序正确、DOI 经 Crossref 全部解析、无伪造。

**结论：稿件已达可发表状态**（soundness-oriented SCIE），唯一样本外遗留见 §4。

---

## 2. 编辑核实与处置（已全部落地）

| # | 意见（来源） | 编辑核实 | 处置 |
|---|---|---|---|
| Q1 | 摘要计数超限（A5 记 310/314；编辑实测 295 白空格 / 313 连字符拆分） | **属实**（按最严连字符口径亦超 300） | 再压缩至 **272 白空格 / 290 连字符拆分** |
| Q2 | 表引用顺序为 1,4,2,3,5（A5） | **属实**：§3.4 先引表 4，晚于 §3.5 的表 2/3 | 重排表块 + 重编号：现为 1（day-7 module t）→2（calibration）→3（markers）→4（modules）→5（single-cell）；正文/图注/gate 同步 |
| Q3 | §2.3 前向引用 "Tables 3–4" 早于 Table 1（复核发现） | 属实 | 删去 §2.3 表号（改为只引 Section 3.5），表首现顺序修复为 1→5 |
| Q4 | §2.3 "adjusted value is the more conservative" 对 Cd68 反向（A1） | **属实**：adjusted 2.86 > unadjusted 2.84 | 改为"we quote the muscle-adjusted value"；并补第三张拟合表（共享基因集 13,601）与 §3.2 单研究计数（13,893）差异的说明 |
| Q5 | `etiology_marker_de.csv`(Gfap 1.328) 与 `DE_D7_vs_Sham_adj.csv`(1.294) 两个肌校正件不一致（A3） | **属实**：前者拟合于两研究共享基因集，后者拟合于单研究表达基因集 | §2.3 增加说明（同一模型、嵌套基因集，非冲突测量） |
| Q6 | T1 的 MWU p 缺非独立性提示（T2 已有）（A2） | 属实 | §3.6 T1 句补"132 效应按对比分组、非完全独立，p 为近似" |
| Q7 | "44.7th percentile" 为靶点**中位数** \|z\| 的秩百分位，非靶点占比（A2） | 属实（措辞） | 保留"just below the median (44.7th percentile)"——表述准确 |
| Q8 | Table 4 为 22 模块中按 p 取前 8 行、未声明（A2） | 属实 | Table 2（原 Table 4）脚注加"top-ranked modules; all 22 in artefacts" |
| Q9 | Data availability 404 仓库（A1/A4/A5） | 属实 | 已为将来时措辞；**建仓为作者侧唯一样本外遗留** |
| Q10 | Usoskin 引用号在 JSON 中仍为 [16]，Manuscript 为 [15]（A1） | 属实（重编号后残留） | 修 script 32 + 其 JSON/log 为 [15] |
| Q11 | §3.8 Cd68/Aif1 大效应零结果宜限定全 DRG 尺度（A3） | 属实 | 补 "at the whole-DRG pseudobulk resolution" |
| Q12 | §2.9 宜说明 gate 期望值的运行时/硬编码性质（A2） | 属实 | §2.9 补"recomputing the key perturbed quantities at run time" |

### 审稿人误判（不采纳）
- 无新误判；A5 的摘要计数（310/314）系其分词器把连字符复合词拆分所致（同轮 A5 又把 gate 的"291 checks"误读为"291 词"），以编辑实测为准。

---

## 3. 修订后状态

- **98 号 gate：291 checks, 0 failed**。
- 摘要 **272 白空格 / 290 连字符拆分**（均 ≤300）。
- 参考文献 1–16 首现顺序正确、全部带 DOI；表 1–5 首现顺序 1→5；图 1–3 齐备、图注与源文件相符。
- 声明齐全：Funding / Conflicts / CRediT / Acknowledgements / **Ethics** / Generative-AI / Data & code。
- 3 图件（png+pdf）；fig3 标题 p 值运行时读取。

## 4. 遗留（唯一样本外项 + 期刊特定项）

- **唯一阻塞项（作者侧）**：建立 `https://github.com/yyx-4113/cancer-pain-drg-virtual-perturbation` 仓库（并录用后 Zenodo 归档）。稿件措辞已改为"将提供"，故不构成虚假声明。
- 期刊特定：目标刊摘要/正文词数上限、关键词数、running head 字符上限、标题长度（202 字符偏长）——投稿前按具体刊物复核。

## 5. 评审循环终止判定

自 Round 3 起四轮结论趋势：Round 4（4 Minor + 1 样本外 Major）→ **Round 5（1 Accept + 4 Minor，0 Major）**。
本轮已无任何科学层面阻断项，全部 Minor 已处置，残留仅为作者建仓与期刊格式。
**满足"重复直到可发表"的终止条件。**
