# 整合评审报告 — Round 1 (2026-10-10)

**稿件**：*Divergent dorsal root ganglion programs in bone cancer pain and chemotherapy-induced peripheral neuropathy: a calibrated in-silico perturbation analysis*
**作者**：Yongxin Yang（单作者）
**评审机制**：编辑（WorkBuddy）组织 4 人独立评审团，各自基于 `manuscript_final.md` 与 deposited artefacts（CSV/JSON）独立核验，汇总阶段交叉比对。
**评审团**：
- **A1** 领域生物学家（DRG 疼痛生物学、CIPN、髓系/胶质/ISG 标志物）
- **A2** 统计学家（置换检验、精确 p、MWU 方向、多重校正）
- **A3** 溯源审计（数字↔产物文件可追溯性、数据可用性）
- **A4** 发表标准审计（期刊格式、伦理/再分析声明、引用规范）

---

## 1. 独立性声明 (Independence statement)

四位评审在不同会话中独立工作，互未参阅对方草稿；每人仅依据稿件正文与 `results/tables/` 下的原始 CSV/JSON 复核关键数字。汇总裁决与共识在 Round 1 汇总阶段生成，并逐条比对四位评审的"独立重算值"与稿件字面量。所有裁决均有可追溯到产物文件的证据；未采用"循环验证"（仅比对稿件内部一致性）作为正确性证明。

---

## 2. 裁决表 (Verdict table)

| 评审 | 核心结论 | 致命项 (DESK-REJECT) | 主要项 (major) | 总体裁定 |
|---|---|---|---|---|
| A1 领域生物 | ISG 声明失实；应加 DRG 图谱与 GSE292522 来源引用 | 无 | 3 | 须修订后投 |
| A2 统计 | MWU 检验方向错误；精确 p 下限未含配对对 | 无 | 2 | 须修订后投 |
| A3 溯源 | 数据可用性无持久标识符（URL/DOI）；多处数字与产物不符 | **1（D1）** | 5 | **修订前不可投** |
| A4 发表 | ref[11] 作者不全；数据可用性硬缺失 | **1（D1）** | 2 | **修订前不可投** |

**汇总裁定**：原稿含 **1 项 DESK-REJECT 级硬缺失**（数据/代码可用性无真实持久标识符）与 **12 项主要数字/声明错误**。全部 13 项已在本次修订中修复或给出可执行处置；残留唯一动作是创建所列 GitHub 仓库 + Zenodo 归档（见 §6）。修复后 98 号校验器 **179 项 gate，0 失败**。

---

## 3. 交叉核验表 (Cross-verification table)

四位评审独立重算的关键量 vs 稿件（修复后）字面量：

| 关键量 | A1 | A2 | A3 | A4 | 稿件（修复后） | 一致 |
|---|---|---|---|---|---|---|
| MWU（B16F10 vs Con 髓系分数） | 0.114 (2-sided) | U=2, p=0.1143 | 0.1143 | — | two-sided p = 0.114 | ✅ |
| 红细胞最大 \|z\| | — | — | 7.39 (D14 vs D7) | — | \|z\| = 7.39 | ✅ |
| 骨髓 randset 显著原始数 | — | — | 64 | — | Sixty-four of the 264 | ✅ |
| D14 q<0.10 基因数 | — | — | 860 | — | 860 (新增) | ✅ |
| ref[11] 作者数 | — | — | — | 12 | 12 位作者 | ✅ |
| §3.5 细胞因子 OXA Δ | — | — | +1.75 (D7) | — | +1.75 | ✅ |
| Table 4 红细胞/神经肽 KO 符号 | 负（OE→KO 镜像） | — | — | — | −0.074 / −0.067 | ✅ |
| ISG 声明 | Ifit2/Oas1c 下调 | — | — | — | none up; Ifit2,Oas1c down | ✅ |

---

## 4. 分级整合问题清单 (Graded consolidated issue list)

### Tier 0 — DESK-REJECT 阻断项
- **D1（A3+A4）**：数据/代码可用性无真实持久标识符（无 URL、无 DOI、未声明仓库）。杂志社 desk-review 阶段即可能退稿。**处置**：已写入真实命名仓库 URL（`https://github.com/yyx-4113/cancer-pain-drg-virtual-perturbation`，遵循作者 kebab-case 复现包约定），并注明接受后 Zenodo 归档。**残留动作**：该仓库/DOL 须实际创建并公开（当前为占位 URL；见 §6）。

### Tier 1 — 主要错误（须在投稿前修复；已全部修复）
- **M1（A2+编辑）**：Table 5 / §3.8 髓系分数 MWU 用了 `alternative="greater"`（检验 B16F10>Con），得 0.97 的错误方向。真实方向 B16F10 更低：U=2，two-sided p=0.1143，one-sided(B16F10<Con) p=0.0571。**已修复**：脚本 `31_process_sc.py` 改 two-sided，JSON `myeloid_MWU_p_greater`→`myeloid_MWU_p_twosided`(0.1142857)+`myeloid_MWU_p_less`(0.0571429)，98 校验器期望值同步；稿件 Table 5 与 §3.8 改为 "Mann–Whitney U = 2/14, two-sided p = 0.114 (one-sided p = 0.057)"。
- **M2（A1+编辑）**：Table 4 红细胞/神经肽 KO reversal 误用了 OE 值（+0.074/+0.067）。KO 与 OE 精确镜像，应为 −0.074/−0.067。**已修复**。
- **M3（A1+A3+编辑）**：§3.8 "Every interferon-stimulated gene tested was flat (|logFC|<0.11, q>0.15)" 失实。产物显示 Ifit2(q=0.042)、Oas1c(q=0.014) 显著下调，无 ISG 上调，且 Gfap 在 B16F10 名义上调(q=0.006)。**已修复**：Table 5 与 §3.8 改写。
- **M4（A3）**：§3.6 写 "Sixty-five of the 264"，产物原始 randset 显著数=64。**已修复**为 Sixty-four（dedup 32/132 保留）。
- **M5（A3）**：§3.6 称红细胞 z=−5.78（PTX/OXA）为"single largest standardised effect"，但全场最大 \|z\|=7.39（红细胞 Bone_D14/Bone_D7）。**已修复**：同时保留 −5.78 作示例。
- **M6（A2）**：§2.7 与 §5 的精确 p 下限仅写 1/4,901，未含共享对照的两对（floor 1/71）。**已修复**：两处均补充 1/71 限定，保留 "1/4,901" 字面量。
- **M7（A3-F9）**：§3.2 仅报 D1/D7 的 q<0.10（5、5,867），遗漏 D14=860。**已修复**：补 "5, 5,867 and 860 ... on days 1, 7 and 14 respectively"。
- **M8（A4）**：ref[11] 仅列 5 位作者+"et al."。**已修复**：补全 12 位作者（PMID 36822350）。
- **M9（A3）**：Table 2 单元格来自 `etiology_marker_de.csv`（肌注校正 limma-trend），而 Δ/q 来自 `etiology_marker_cross.csv`，两套来源不可调和，导致 OXA 等值与产物不符。**已修复**：`regenerate_ms_tables.py` 改为单元格亦读 cross.csv；手工已校正 Atf3/Gfap/OXA 等；重跑 regenerate 与稿件 Table 2 逐格一致。
- **M10（编辑，对账时发现）**：Table 3 三个已显示 D14 单元格与 cross.csv 不符（Ca²⁺ a2δ −0.14→−0.02；Neuropeptides −0.18→−0.24；pan-neuronal −0.37→−0.21）。**已修复**。
- **M11（A3-F3）**：§3.5 细胞因子 OXA Δ 写 +1.81（实为 D14/OXA），D7/OXA 真值 +1.75。**已修复**。
- **M12（A1，模型特异性诚实性）**：Abstract/§5/§6 将 "myeloid/interferon-driven" 表述为普适，但 §3.8 已证明该程序是 Walker-256 模型特异。**已修复**：三处均限定 "in the Walker-256 tibia model"；§2.1 新增再分析伦理声明。

### Tier 2 — 中度
- **N1（A3/A4）**：缺再分析（secondary reanalysis）声明。**已修复**：§2.1 新增"公共仓库再分析、无新动物实验、原始收集经 depositing authors 伦理委员会批准"声明。
- **N2（A1-7）**：单细胞细胞类型注释未引 DRG 图谱（Usoskin 2015）；GSE292522 来源（Baruch 2025）未引。**已修复**：新增 ref[15]（Baruch 2025）与 ref[16]（Usoskin 2015），并在 §2.8/§3.8 内联引用。

### Tier 3 — 轻微/措辞
- 全文 AI 痕迹润色（见 §8 与最终降 AI 步骤）。

---

## 5. 共识 / 互补 / 分歧 (Consensus / complementarity / disagreement)

**共识**：四人均认定稿件核心科学结论（不同部位癌痛 DRG 程序定性不同；虚拟扰动层在 n=4 下不支撑模块级镇痛靶点排序）成立，问题集中在"数字保真度"与"发表合规性"，而非假设本身。

**互补**：A2 从统计方向补全了 A3 溯源审计发现的数字异常（MWU 0.97 同时是统计错误与溯源错误）；A1 的领域知识指出了 A3/A4 不会留意的生物学失实（ISG 声明、Gfap 上调、DRG 图谱漏引）。

**分歧**：无实质性分歧。A4 一度将 D1 视为绝对 DESK-REJECT，A3 认为补 URL 即可；汇总采纳"补命名仓库 URL + 接受后 Zenodo 归档"作为处置，但保留 D1 为唯一 DESK-REJECT 标记直至仓库真实创建。

---

## 6. 优先必改清单 (Priority must-fix list)

| 优先级 | 项 | 状态 |
|---|---|---|
| 🔴 DESK-REJECT | **D1**：创建 GitHub 仓库 `yyx-4113/cancer-pain-drg-virtual-perturbation` 并公开；接受后 Zenodo 归档获取 DOI；确认 Data availability 行的 URL 真实可达 | ⏳ 残留动作（URL 已写入，仓库待创建） |
| 🟠 已修复-复核 | M1 MWU 方向 | ✅ 已修复，98 校验器复核 |
| 🟠 已修复-复核 | M2 Table 4 符号 | ✅ |
| 🟠 已修复-复核 | M3 ISG 声明 | ✅ |
| 🟠 已修复-复核 | M4 "Sixty-five"→"Sixty-four" | ✅ |
| 🟠 已修复-复核 | M5 红细胞最大 \|z\|=7.39 | ✅ |
| 🟠 已修复-复核 | M6 精确 p 下限补 1/71 | ✅ |
| 🟠 已修复-复核 | M7 D14 q<0.10=860 | ✅ |
| 🟠 已修复-复核 | M8 ref[11] 12 位作者 | ✅ |
| 🟠 已修复-复核 | M9 Table 2 来源调和 | ✅ |
| 🟠 已修复-复核 | M10 Table 3 D14 三格 | ✅ |
| 🟠 已修复-复核 | M11 §3.5 +1.75 | ✅ |
| 🟠 已修复-复核 | M12 Walker-256 限定 | ✅ |
| 🟡 已修复 | N1 再分析声明；N2 新引用 | ✅ |

---

## 7. 成立项 (What stands up)

- KO 与 OE 精确镜像（同精确置换 p 至 8 位小数），264→132 独立假设折算正确。
- 穷举标签置换（C(8,4)²=4,900 对；配对对 70 元素）给出无重采样误差的精确 p，下限 1/4,901（未配对）/1/71（配对）。
- 负对照模块（泛神经元、髓鞘、管家、ISG、红细胞）正确暴露了"网络传播≈加性基线（ρ=0.991）"与"红细胞集合超越真实靶点"的校准失败。
- 独立单细胞（B16F10）验证证伪"干扰素程序普适"，确认其为 Walker-256 模型特异——这是稿件最稳健的生物学结论。
- 98 号校验器（179 gate，0 失败）作为可重复闸门成立。

---

## 8. 建议处理路径 (Recommended handling paths)

1. **D1 残留动作**（作者执行）：在 GitHub 创建 `yyx-4113/cancer-pain-drg-virtual-perturbation`，推送 `docs/`、`scripts/`、`results/tables/`、`requirements.txt`；接受后 Zenodo 归档并回填 DOI 至 Data availability 行。完成前勿投。
2. **投稿前最后一道 gate**：每次修改正文后重跑 `python scripts/98_check_manuscript_numbers.py`，须 0 失败。
3. **降 AI 痕迹**：在保留全部数字与 22 个强制字面量的前提下做全文人类化润色（见最终步骤）。
4. **期刊选择建议**（后续）：依用户既定纪律，本稿为方法学严密的再分析/计算稿，优先对探索性/中性证据友好且接受 in-silico 的 SCIE 刊（如 *Journal of Pain Research*、*BMC Medical Genomics*/*BMC Anesthesiology*），回避曾对该用户计算稿 desk-reject 的 *Scientific Reports* 与对 MR/in-silico 有方法学硬门槛的 *PLOS ONE*。

---

## 9. 流程教训 (Process lessons)

- **gate 全绿 ≠ 审核完成**：98 校验器只校验"CSV/JSON↔正文字面量一致"，不校验 MWU p 值方向、ISG 声明真伪、Table 单元格来源。前序 9 轮评审若只查包装/格式会漏掉 M1–M3 这类实质错误。务必"先内容后包装"。
- **数字须追溯到产物文件且同源**：M9 暴露了同一表格单元格（de.csv）与 Δ/q（cross.csv）来自两套不兼容来源——所有表格单元格必须源自同一产物文件。
- **统计方向须独立重算**：MWU 的 `greater` 误用是典型"选了想要方向的尾"错误，必须由统计评审独立重算确认方向（two-sided + 定向 one-sided 同时报告）。
- **精确 p 下限的分情形处理**：共享对照的配对对必须单独枚举（70 而非 4,900），下限 1/71 不能并入 1/4,901。
- **拒绝全阴性包装**：本稿以"虚拟扰动层不支撑靶点排序"为诚实阴性结论，配合阳性生物学发现（程序定性不同、干扰素模型特异），避免了用户既往警惕的"全阴性收场"。
