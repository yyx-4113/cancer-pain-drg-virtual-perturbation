# GITHUB 仓库建仓 / 归档操作手册（中文）

> 适用项目：`cancer-pain-drg-virtual-perturbation`
> 稿件：*Divergent dorsal root ganglion programs in bone cancer pain and
> chemotherapy-induced peripheral neuropathy …*（单作者 Yongxin Yang）
> 本手册记录把本地分析目录变成公开复现包仓库、打版本 tag、并在接收后归档到
> Zenodo 的完整操作。**所有命令均可直接复制运行。**

---

## 0. 前置条件

- GitHub 账号 **`yyx-4113`**（科研主账号；`yongxinyang` 与科研无关，勿用）。
- 本机已安装并登录 GitHub CLI：
  ```bash
  gh auth status        # 期望：✓ Logged in to github.com account yyx-4113，scopes 含 'repo'
  ```
  若未登录，在宿主机执行 `gh auth login`（AI 无法代填令牌）。
- 仓库名（kebab-case 主题名）：**`cancer-pain-drg-virtual-perturbation`**。

---

## 1. 本地初始化（只做一次）

在项目根目录 `cancer_pain_virtual_perturbation/` 下：

```bash
git init -b main
git config user.name  "Yongxin Yang"
git config user.email "960856791@qq.com"     # 与 GitHub 账号绑定邮箱一致
```

### 1.1 提交范围（重要）

`.gitignore` 已排除 `data/raw/` 与 `data/processed/GSE292522/`（合计约 1.4 GB，GEO 为权威来源、
可再生，不随仓库分发）。**其余产物（脚本、results/tables 下的 CSV/JSON、docs、图件）全部入库。**

- **禁止** `git add -A` / `git add .` —— 本目录可能混有临时或与本仓库无关的文件。
- 只按目录显式添加：

```bash
git add scripts results docs review_A2 \
        README.md CITATION.cff LICENSE requirements.txt \
        MANIFEST.sha256 GITHUB_DEPOSIT_SOP.md author_verification_statement.md \
        .gitignore .github
git status --short          # 复核将要提交的文件，确认无 data/raw、无 __pycache__
```

---

## 2. 首次提交

```bash
git commit -m "Release v1.0.0: cross-study DRG reanalysis + negative-control calibration"
```

> 若出现 `LF will be replaced by CRLF` 警告属正常，不影响提交。

---

## 3. 在 GitHub 建仓并推送

```bash
# 建立公开仓库（不要勾选自动生成 README/LICENSE，本仓库已有）
gh repo create yyx-4113/cancer-pain-drg-virtual-perturbation \
    --public \
    --description "Code and artefacts: divergent DRG programs in bone cancer pain vs CIPN; negative-control calibration of in-silico perturbation" \
    --source . --remote origin

git push -u origin main
```

成功标志：输出 `... main -> main`。

**连通性排错**：若 `git ls-remote origin` 失败，是代理/网络问题而非认证问题；
沙箱经 `http_proxy`/`https_proxy` 环境变量出网，**不要**手动设置 `git config http.proxy`。
真实认证失败会报 `could not read Username for 'https://github.com'`。

---

## 4. 打版本 tag 并触发 Release

仓库内置 `.github/workflows/release.yml`：**推送 `v*` tag 即自动创建 GitHub Release**
（`softprops/action-gh-release` + `generate_release_notes`），无需手动建 Release。

```bash
git tag -a v1.0.0 -m "v1.0.0 — initial public release: manuscript pipeline, artefacts, gate (98) and unit tests (99)"
git push origin v1.0.0
```

核验：

```bash
git ls-remote --tags origin          # 应列出 v1.0.0
gh release view v1.0.0               # 查看自动生成的 Release
```

---

## 5. 更新稿件的数据可用性声明

建仓推送完成后，把稿件 `docs/manuscript/manuscript_final.md` 的
**Data and code availability** 段由将来时改为现在时，并给出**实名仓库 URL + 版本 tag**：

> Analysis code, artefacts and the test suite are available at
> https://github.com/yyx-4113/cancer-pain-drg-virtual-perturbation (tag v1.0.0),
> with `requirements.txt` pinning every version used, and will be archived at Zenodo
> upon acceptance.

改完务必重跑闸门：`python scripts/98_check_manuscript_numbers.py`（期望 0 failed）。

---

## 6. 接收后归档到 Zenodo（留 DOI）

1. 在 Zenodo 用 **GitHub 登录**，于 *Settings → GitHub* 打开本仓库的开关。
2. 回 GitHub 创建**新 tag**（如 `v1.0.1` 或 `v1.1.0`）并推送 —— Zenodo 自动抓取该快照、
   生成 **版本 DOI** 与 **概念 DOI**（concept DOI，永远指向最新版）。
3. 到 Zenodo 记录编辑：补齐作者（**Yongxin Yang，不加 MD/PhD 后缀**）、机构
   （The Second Affiliated Hospital of Fujian University of Traditional Chinese Medicine）、
   ORCID `0009-0004-9698-6552`、License（MIT）、关键词。
4. 把生成的 **版本 DOI** 回填到稿件 Data-and-code-availability 段与 `CITATION.cff`
   （`doi:` 字段），再重跑闸门。

> 如需用脚本直接投递（不经 GitHub 集成），可用 `zenodo-inveniormd-deposit` 流程，
> 但推荐 GitHub 集成，保证「tag ↔ DOI」一一对应、可追溯。

---

## 7. 版本号约定

- `vX.Y.Z`：`X` = 稿件级大版本（重分析/新数据）；`Y` = 分析或产物变更；`Z` = 文档/修正。
- 每轮稿件修订完成后，若产物有变，同步推一个 tag，确保「稿件数字 ↔ 仓库 tag」可互查。

---

## 8. 常见问题

| 现象 | 原因 / 处理 |
|---|---|
| `fatal: not a git repository` | 尚未 `git init`，先执行第 1 节。 |
| `! [rejected] (fetch first)` | 本地/远端分叉，先 `git fetch && git merge`（**勿 rebase**，会因 autocrlf 报未暂存改动）。 |
| 推送体积超限 / 很慢 | 检查是否误 `git add` 了 `data/raw/`；用 `git status --short` 复核，必要时 `git rm --cached` 后重提交。 |
| Release 未生成 | 确认 tag 以 `v` 开头且 workflow 已随 `main` 推送；查看仓库 *Actions* 页。 |
| 误用 `yongxinyang` 账号 | 该账号非科研账号，务必用 `yyx-4113`。 |

---

*本手册随仓库版本升级维护。最近更新：2026-10-10（v1.0.0 建仓）。*
