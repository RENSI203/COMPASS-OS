# 00 · 候选版本冻结与审查范围 ｜ candidate baseline

**阶段**：0 / 6（顺序发布审查）
**日期**：2026-10-04
**结论**：**候选已明确**（`624184c`）｜ **可进入第 1 阶段**

> 本阶段**未** merge / tag / push / 发布。所有操作限于本地分支与本地 worktree。

---

## 1. 分支、HEAD、工作树与版本号

### 1.1 冻结前的事实

| 项 | 值 |
|---|---|
| 分支 | `feature/sample-sankey` |
| 冻结前 HEAD | `55df2c6`（Add sample-level computation-path visualization） |
| 冻结前工作树 | 14 个已跟踪文件被修改 + 29 个未跟踪条目（其中 24 项属本功能，20+ 项属**并行论文工作区**） |
| 冻结前版本号 | `pyproject.toml` `1.2.0` · `CITATION.cff` `1.2.0` · `__init__.py` `1.2.0` · `README.md` 徽标 `1.1.0`（**四者不一致**） |

### 1.2 真实发布历史（权威依据）

```
$ git tag -l                    $ git ls-remote --tags origin
v1.0.0                          v1.0.0
v1.0.1                          v1.0.1

$ GitHub Releases API
v1.0.0   COMPASS-OS v1.0.0   published=2026-10-03T13:48:37Z
v1.0.1   COMPASS-OS v1.0.1   published=2026-10-03T15:57:20Z
```

| tag | commit | 该 tag 的版本号 |
|---|---|---|
| `v1.0.0` | `3f2aead` | `1.0.0` |
| `v1.0.1` | `7270935` | `1.0.1` |

**事实：`1.1.0` 从未发布。** 本地与远端均无 `v1.1.0` tag，GitHub 上亦无对应 Release。

### 1.3 为什么曾从 1.1.0 改为 1.2.0 —— 以及为什么现在改回 1.1.0

**当初改成 1.2.0 的理由**（记录在案，本身是合理的语义化版本推理）：
渲染层集成把 `SamplePathResult.figure()` 的返回类型从 Plotly `Figure` 换成 Matplotlib
`RenderedPath`，并移除了 `write_html` / `write_image` / `include_plotlyjs` 等 Plotly 专有接口。
按语义化版本，破坏公开 API 应当升 minor ⇒ 从 `1.1.0` 升到 `1.2.0`。

**为什么这个理由在本仓库不成立**：
破坏性变更只有在**存在已发布的旧版本**时才需要一个新的版本号去区分。
而 `1.1.0` **从未发布** —— Plotly 版 `sample_path` 从未出现在任何 tag、Release 或分发产物中，
没有任何用户见过它，也不存在"从 1.1.0 迁移到 1.2.0"的受众。
`55df2c6`（引入 Plotly 版 sample_path）与本次渲染层集成属于**同一条未发布的开发线**，
应当作为**一个**版本交付。

继续使用 `1.2.0` 会造成两个具体危害：

1. **虚构发布历史**：对外暗示存在过一个 `1.1.0` 版本，而它并不存在；
2. **误导迁移说明**：`docs/SAMPLE_PATH.md` 会写出"`v1.1.x` → `v1.2.0`"的迁移表，
   读者会去找并不存在的 1.1.x 产物。

**因此候选版本定为 `1.1.0`** —— 已发布历史 `v1.0.1` 之后的**下一个未发布 minor**，
且它**同时**包含 sample_path 功能与渲染层集成。

> **不是**因为测试通过才改版本号。版本号依据是**发布历史的事实**（§1.2），
> 测试结果与版本号无关。

### 1.4 冻结后的版本一致性

| 位置 | 值 | 说明 |
|---|---|---|
| `pyproject.toml` `version` | `1.1.0` | ✅ |
| `src/compass_os/__init__.py` `__version__` | `1.1.0` | ✅ |
| `CITATION.cff` `version` | `1.1.0` | ✅ |
| `README.md` 版本徽标 | `1.1.0` | ✅ |
| `CITATION.cff` `cff-version` | `1.2.0` | ⚠️ **保持不变** —— 这是 Citation File Format 的 **schema 版本**，与软件版本无关，改动它会破坏 CFF 校验 |

四处软件版本现已一致，`cff-version` 按规范保持独立。

---

## 2. 候选包含的内容

**候选提交**：`624184c`（24 个文件，+2669 / −973）

### 2.1 代码

| 文件 | 状态 | 作用 |
|---|---|---|
| `src/compass_os/sample_path.py` | 改写 | 计算 + 结果对象；展示节点改由渲染层 `select_nodes` 统一决定 |
| `src/compass_os/sample_path_render.py` | **新增** | 圆形分层渲染层（单一静态管线，惰性导入 numpy/matplotlib） |
| `src/compass_os/_render_bridge.py` | **新增** | 内部桥接：`SamplePathResult` → `SamplePathData` |
| `src/compass_os/__init__.py` | 修改 | 版本号 → `1.1.0`（公共导出不变） |

**未新增任何运行期硬依赖**：去掉了 Plotly / Kaleido / CDN 需求；matplotlib 自 v1.0 起即为核心依赖。

### 2.2 测试

| 文件 | 状态 | 内容 |
|---|---|---|
| `tests/test_sample_path.py` | 修改 | 27 项功能与回归用例 |
| `tests/test_sample_path_integration.py` | **新增** | 6 项生产集成用例（规格 §9 的 12 项要求） |
| `tests/verify_sample_path_renderer.py` | **新增** | 附件 12 项渲染层自检（原样保留，改从包内导入） |
| `tests/run_tests.py` | 修改 | 注册新模块 |

### 2.3 配置

| 文件 | 改动 |
|---|---|
| `pyproject.toml` | 软件版本 → `1.1.0`；`plot` extra 注释改为"matplotlib 已是核心依赖，本功能无新增硬依赖" |
| `CITATION.cff` | 软件 `version` → `1.1.0`；`cff-version` 不动 |

### 2.4 文档

| 文件 | 状态 |
|---|---|
| `README.md` | §2b 按集成后的 API 与视觉规格重写 |
| `docs/SAMPLE_PATH.md` | 用户指南重写（调用、分数来源、预算、风险柱、导出、依赖、迁移、解读边界） |
| `docs/SAMPLE_PATH_INTEGRATION_VALIDATION.md` | 集成验收报告 |
| `docs/SAMPLE_PATH_REDESIGN_WIRING.md` | 渲染层接线报告（历史审计轨迹） |
| `docs/SAMPLE_PATH_VISUAL_REVIEW.md` | v1.1.x 旧渲染器的人工验收报告（历史轨迹，见 §7 说明） |

### 2.5 验收证据（图片与审计记录）

`docs/figures/sample_path_integration/` 下的 **7 个 PNG**（6 张单图 + 2×3 总览）、
`REAL_M2_M3_selection_audit.tsv`、`shared_style.json`。

> `MANIFEST.in` 只把 `docs/*.md` 收进 sdist、**不收** PNG；wheel 完全不含 `docs/`。
> 因此这些图片**只存在于 git 仓库**，用于人工验收留痕，不影响任何分发产物的体积或可用性。
> 同名 `.html` / `.pdf` / `.svg`（7.0 MB，可再生）**未纳入候选**，用一条命令即可重建：
> `python examples/sample_path_example.py --cohort GSE39582 --out <dir>`。

### 2.6 模型资产

**未改动**。冻结资产的 SHA256 与已发布 `v1.0.1` **逐字节一致**（§5.3 举证）。

### 2.7 发布检查脚本

| 文件 | 说明 |
|---|---|
| `tools/build_assets.py` | 资产哈希校验（`--verify`） |
| `tools/release_audit.py` | 发布卫生审计（检查项见 §6） |
| `tools/build_review_grid.py` | **新增**：6 张图拼 2×3 总览（只拼贴，不重绘） |

### 2.8 明确**排除**在候选之外的内容

| 排除项 | 理由 |
|---|---|
| `paper_assets/` 下 4 个已修改 + 20 个未跟踪条目 | **并行论文工作区**，与本功能无关 |
| `validation/results/release_audit.tsv` | 审计**运行产物**（内容随运行时机变化），应在发布时重新生成，不作为候选输入 |
| `docs/figures/sample_path_review/`（15 MB） | v1.1.x 旧 Plotly 渲染器的验收图，已被取代 |
| `docs/figures/sample_path_redesign/`（37 MB） | 上一轮独立渲染器的验收图，已被集成版取代 |
| `docs/figures/sample_path_integration/*.{html,pdf,svg}`（7.0 MB） | 可一条命令再生，见 §2.5 |

**举证**：本轮提交未触碰任何 `paper_assets` 文件 ——

```
$ git show --name-only --format="" 624184c | grep paper_assets   → 无输出
$ git diff --name-only 55df2c6 624184c -- paper_assets | wc -l   → 0
```

工作树中并行工作区的改动**原样保留、未被丢弃、未被提交**。

---

## 3. 候选提交

```
commit  624184c161e05026cf67695ea2e7f27e343202c4
author  RENSI203 <RENSI203@users.noreply.github.com>
date    Sun Oct 4 17:13:27 2026 +0800
subject Integrate the circular layered sample-path renderer into the package
files   24 changed, 2669 insertions(+), 973 deletions(-)
parent  55df2c6
```

创建方式：**显式列出 24 个路径** `git add <paths>`（**未使用 `git add .` / `-A`**），
暂存后校验 `git diff --cached --name-only | grep -E "paper_assets|validation/"` 为空，再提交。

**未** merge、**未** 打 tag、**未** push。

---

## 4. 从候选提交建立的 clean worktree

```
$ git worktree add --detach /tmp/gh_audit/cand 624184c
HEAD is now at 624184c Integrate the circular layered sample-path renderer into the package

$ cd /tmp/gh_audit/cand && git status --porcelain | wc -l
0                       # 与候选提交完全一致，无未提交改动
```

| 项 | 值 |
|---|---|
| worktree 路径 | `/tmp/gh_audit/cand` |
| 检出提交 | `624184c1`（`git rev-parse HEAD`） |
| 未提交改动 | **0** |
| 仓库文件数 | 202 |
| 仓库体积 | 27.0 MB |

**本阶段及其后所有构建与测试均以此 worktree 为依据。**
`paper_assets/` 在该 worktree 中出现的是**历史既有已跟踪内容**（20 个文件，来自更早的提交），
**不是**本功能引入的；本轮提交与其差异为 0（§2.8 举证）。

---

## 5. 记录：commit、清单、资产哈希、环境、产物哈希

### 5.1 候选提交与清单

* 候选提交：`624184c161e05026cf67695ea2e7f27e343202c4`
* 源码清单：**`docs/audit/CANDIDATE_MANIFEST.tsv`**（202 行 + 表头）
  每行含 `path` / `size_bytes` / `sha256`，覆盖候选提交的**全部**文件。

### 5.2 环境信息

| 项 | 值 |
|---|---|
| OS | Ubuntu 24.04.4 LTS（WSL2，kernel 6.18.33.2-microsoft-standard-WSL2） |
| 架构 | x86_64 |
| Python | 3.12.14 |
| 解释器 | `mamba_root/envs/compass/bin/python` |
| numpy | 2.5.2 |
| pandas | 2.2.3 |
| scikit-learn | 1.9.0 |
| scikit-survival | 0.28.0 |
| torch | 2.10.0+cpu |
| matplotlib | 3.10.8 |
| plotly | 已安装（**但候选不再需要**；仅作对照记录） |
| kaleido | **未安装**（候选不再需要） |
| 必需环境变量 | `MPLCONFIGDIR=/tmp/mplcfg`（`~/.config` 不可写） |

### 5.3 冻结模型资产 SHA256（与已发布 v1.0.1 对比）

| 资产 | 候选 SHA256（前 16） | 与 v1.0.1 |
|---|---|---|
| `pretrainer.pt` | `4fd174b45be816f3` | ✅ 一致 |
| `locked_M0.json` | `f368eb2e399f1616` | ✅ |
| `locked_M1.json` | `d4998ecf646aa264` | ✅ |
| `locked_M2.json` | `5aac2cd089934f62` | ✅ 一致 |
| `locked_M3.json` | `dfe17fc56c573b2c` | ✅ |
| `locked_pca_M3.npz` | `e1f20fcc20be5868` | ✅ |
| `reference_quantiles.json` | `2e8fc0472020cf8c` | ✅ |
| `qc_config.json` | `098dc59b752fcc64` | ✅ 一致 |
| `model_manifest.json` | `500c7ea2d3bffb22` | ✅ |
| `third_party/compass/utils/sankey.py` | `20ed31d9eb38f438` | ✅（上游参考脚本，只读） |

`tools/build_assets.py --verify` → **校验 13 项 ⇒ 异常 0**。

### 5.4 构建产物与 SHA256

**规范化构建命令**（关键：固定 `SOURCE_DATE_EPOCH` = 候选提交时间 `1791105207`）：

```bash
cd /tmp/gh_audit/cand
rm -rf dist build src/*.egg-info
SOURCE_DATE_EPOCH=$(git log -1 --format=%ct HEAD) \
  python -m build --no-isolation
```

| 产物 | SHA256 | 可重复性 |
|---|---|---|
| `compass_os-1.1.0-py3-none-any.whl` | `a26171d62f6e6b94cfca6782e371154bb45aa5ed082f2851e10e2aa4a1ffbca4` | ✅ **逐字节可重复** |
| `compass_os-1.1.0.tar.gz` | `99d180855bc93d428d95579a88362256e1503b5aa1af213cb42781151ecffea7` | ⚠️ **逐字节不可重复**（见 §5.6） |
| sdist **内容**哈希（解包后按 名称+内容 累计） | `6042a30bb46f94b3266ff774493063fa0e2c7ecdcf90a2fc5cf608488b1d45cd` | ✅ 可重复 |

wheel METADATA（节选）：

```
Name: compass-os          Version: 1.1.0          Requires-Python: >=3.10
Requires-Dist: numpy>=1.26, pandas>=2.0, scikit-learn>=1.3,
               scikit-survival>=0.22, matplotlib>=3.7
Requires-Dist: immuno-compass==2.5.3; extra == "upstream"
Requires-Dist: pytest>=7; extra == "test"
Requires-Dist: matplotlib>=3.7; extra == "plot"
```

**无 Plotly / Kaleido 依赖** ✅

### 5.5 ⚠️ 此前 clean-wheel 验证实际使用的源码（必须明确）

**不能写"HEAD 未变、测试通过"** —— 上一轮的 clean-wheel 验证**不是**从一个提交构建的：

| 步骤 | 上一轮（渲染层集成验收）实际做法 |
|---|---|
| 1 | `git worktree add --detach /tmp/gh_audit/cleanwt HEAD`，此时 HEAD = **`55df2c6`** |
| 2 | **手工把 13 个工作树文件 `cp` 进该 worktree**（`sample_path.py`、`sample_path_render.py`、`_render_bridge.py`、`__init__.py`、`pyproject.toml`、`CITATION.cff`、`README.md`、`docs/SAMPLE_PATH.md`、`examples/sample_path_example.py`、`tests/*`） |
| 3 | 在该 worktree 内 `python -m build` |

因此上一轮验证的状态是 **`55df2c6` + 13 个未提交文件**，
**这个状态从未作为提交存在过**，无法通过任何 commit hash 复现。
它当时记录的 wheel 版本号为 `1.2.0`，与本次候选（`1.1.0`）**不同**，产物不可比。

**本轮已修正的方法学**：候选先落成提交 `624184c`，再从**该提交**建 worktree（§4），
之后所有构建与测试都可由 commit hash 完整复现。

> 这是本阶段发现并已修复的**流程缺陷**（可复现性缺陷），不是科学或功能缺陷。

### 5.6 ⚠️ 本轮新发现的真实缺陷：构建产物可重复性

**缺陷 A（wheel 非确定性）** —— 已修复并验证

* **现象**：同一源码树连续两次 `python -m build`，wheel SHA256 不同
  （`a6ca2d7f…` vs `f55f629f…`）。
* **根因**（有举证）：两个 wheel **83 个条目名称与内容全部相同（内容差异条目 = 0）**，
  仅 **6 个 `dist-info` 条目的 ZIP 时间戳**不同 —— 这些条目由 setuptools 在构建时新生成，
  ZIP 头写入构建时刻。
* **影响**：记录的 SHA256 只能标识"某一次构建"，无法标识"某一源码状态"，
  削弱发布溯源与第三方复核。
* **修复**：构建时固定 `SOURCE_DATE_EPOCH`。
* **验证**：同一 `SOURCE_DATE_EPOCH` 连续两次构建 → SHA256 **完全一致**
  （`601508db…` 两次相同；采用提交时间戳后为 `a26171d6…`）。
* **状态**：**已修复**（规范化构建命令见 §5.4，已写入本文件作为发布要求）。

**缺陷 B（sdist 目录 mtime 非确定性）** —— 已定位，**未修复**（转第 2 阶段）

* **现象**：即使固定 `SOURCE_DATE_EPOCH`，连续两次 sdist 的 SHA256 仍不同
  （`f71c3e49…` vs `bd675c35…`）。
* **根因**（有举证）：两次 sdist **解包后 `diff -rq` 无差异**，文件清单 152 项**完全一致且顺序相同**；
  差异仅出现在 **31 个条目的 mtime**，且全部是**目录**条目与构建时重新生成的文件
  （`PKG-INFO`、`setup.cfg`、`*.egg-info/*`）—— setuptools 对这些条目使用**构建时刻**，
  不遵循 `SOURCE_DATE_EPOCH`。
* **影响**：sdist 的**字节**哈希不可复现；但**内容**可复现（可用 §5.4 的 content-sha256 校验）。
  wheel 不受影响（已可重复），而 wheel 才是 `pip install` 实际使用的产物。
* **修复方案**：属 setuptools 行为，无干净的应用层修法；候选做法是
  ① 用 content-sha256 作为 sdist 的溯源依据，② 在发布说明中标注 sdist 字节哈希的**构建特定性**。
* **状态**：**未修复，转第 2 阶段（依赖项与安装审查）**评估是否需要 pin 构建后端或改用其他打包器。
  本阶段只记录，不擅自更换打包工具。

---

## 6. `release_audit.py` 检查项职责，与绝对路径 FAIL 的裁定

### 6.1 全部检查项及其职责

`tools/release_audit.py`（145 行）**扫描工作树**（`REPO.rglob("*")`），**不读 git 索引**。
任一 FAIL → 退出码 1。

| # | check | 触发条件 | 状态 | 职责 |
|---|---|---|---|---|
| 1 | `artifact` | 路径含 `__pycache__` / `.pytest_cache` | WARN | 提醒清理字节码（`.gitignore` 已覆盖） |
| 2 | `artifact` | 路径在 `validation/inputs` / `validation/output` 下 | WARN | 提醒本地输入/运行产物不得入库 |
| 3 | `forbidden-name` | 文件名匹配 `TcgaTargetGtex` / `expression_tcga_*` / `series_matrix` / `.soft` / `GSE*.gz|txt` / `downstream_mechanism` / `analysis_out_` / `cohorts?_curve|all.txt` | **FAIL** | **不得分发原始 TCGA/GEO 数据、大型表达矩阵或下游机制结果**（许可与体积双重风险） |
| 4 | `size` | 单文件 > 100 MB | **FAIL** | GitHub 硬限制 |
| 5 | `size` | 单文件 > 50 MB | WARN | 需人工确认 |
| 6 | `absolute-path` | `.py/.md/.cff/.toml/.in/.cfg/.txt` 中某行匹配 `/home/<用户名>/`、`/mnt/<盘符>/`、`<盘符>:\\`、`projects/202608`、以及本机用户名（且该行不含 `http`） | **FAIL** | **不得把本机绝对路径带进发布物** —— 保证可移植性、避免泄露本机目录结构。白名单 `ALLOW_PATH_PATTERNS` 已豁免 `tools/`、`validation/inputs/`、`ASSET_MANIFEST.tsv`、`README.md` 等**溯源用途**位置 |
| 7 | `asset-hash` | `build_assets.py --verify` 返回非 0 | **FAIL** | 冻结模型资产必须与 `ASSET_MANIFEST.tsv` 逐字节一致 |
| 8 | `required-file` | 16 个关键文件缺失 | **FAIL** | README/LICENSE/NOTICE/CITATION/pyproject/MANIFEST.in/ASSET_MANIFEST/模型/测试等齐备 |
| 9 | `internal-content` | 路径含 `docs/_` 或 `handoff` | WARN | 提醒内部草稿/交接材料 |

### 6.2 绝对路径检查的裁定：**属于正式必需检查，不删除、不豁免**

**理由**：该检查守护的正是本阶段与第 2 阶段最关心的问题 ——
"安装后的功能不依赖仓库绝对路径"。它是可移植性的一线防线，
**不能**为了让审计变绿而把命中的路径加进白名单（那等于忽略问题）。

### 6.3 该 FAIL 的证据与本候选的关系

混合工作树（`github/`）中：

```
FAIL  absolute-path  paper_assets/audit/verify_wheel_install.py:52
      <该行内容为一个硬编码的本机 Python 解释器绝对路径>
```

> **本文件对证据做了脱敏**：不直接复制该行原文，否则本审计报告自身就会命中
> `absolute-path` 检查（实测确实命中过 —— 见 §6.5）。原始字符串可用
> `python tools/release_audit.py` 复现；其内容为
> `<家目录>/<用户名>/miniconda3/bin/python3` 形式的解释器路径。

**证据链**：

| 步骤 | 命令 / 结果 |
|---|---|
| 该文件是否属于候选？ | `git cat-file -e 624184c:paper_assets/audit/verify_wheel_install.py` → **不存在**（未跟踪文件） |
| 该文件是否随任何产物分发？ | 否 —— 未跟踪 ⇒ 不进 git、不进 sdist、不进 wheel |
| 在候选 worktree 中运行审计 | `cd /tmp/gh_audit/cand && python tools/release_audit.py` → **检查项 16：FAIL 0 / WARN 0** |
| 本轮提交是否引入该行？ | 否 —— `git show 624184c -- paper_assets` 无输出 |

**裁定**：

1. **检查必需**：保留 `absolute-path`，不移除、不放宽白名单。
2. **本候选 PASS**：以**候选 worktree 的审计结果**为准 —— `FAIL 0 / WARN 0`。
3. **该缺陷真实存在但不在本候选范围内**：它属于**并行论文工作区**的未跟踪脚本
   （`paper_assets/audit/verify_wheel_install.py`）。本轮**无权修改**该工作区
   （用户明确要求"保留现有工作树中与本功能无关的修改"），故**不修**。
   已定位到具体行，建议其负责人把硬编码解释器改为 `sys.executable`
   （该脚本用途是验证 wheel 安装，本就应当用当前解释器）。
4. **正式范围说明**（写入本文件，作为后续阶段的口径）：

   > `tools/release_audit.py` **扫描工作树**而非 git 索引。因此它的结果只有在
   > **候选提交的干净 checkout** 上运行才具权威性。在混有并行工作区未跟踪文件的
   > 开发工作树上运行，报出的 `absolute-path` / `artifact` 命中可能并不属于候选，
   > 必须逐条用 `git cat-file -e <candidate>:<path>` 与
   > `git show --name-only <candidate>` 归因后才能采信。
   >
   > **判定规则**：候选是否通过发布审计，以**候选 worktree 的审计输出**为唯一依据；
   > 开发工作树上的 FAIL 只有在能归因到候选提交时才构成候选缺陷。

### 6.5 本阶段的一个自我命中及其修复（保留记录）

**现象**：本文件初稿在 §6.3 中**逐字引用**了那条命中行（含真实的本机解释器绝对路径）
作为证据，结果**本文件自身**被 `release_audit.py` 判为 `absolute-path` FAIL：

```
$ cd <clean candidate worktree> && python tools/release_audit.py
检查项 17：FAIL 1 / WARN 0
=== FAIL ===
absolute-path   docs/audit/00_CANDIDATE_BASELINE.md:363   ...   FAIL
```

**判断**：这是**真实缺陷**（发布物中不得含本机绝对路径），且由本阶段引入。
按共同规则"发现真实缺陷后直接修复"，**不**通过把 `docs/audit/` 加进
`ALLOW_PATH_PATTERNS` 白名单来绕过（那属于"为了通过而忽略问题"）。

**修复**：把证据行**脱敏**为文字描述（保留文件、行号与结论，去掉路径字面量），
并在旁边说明原始字符串可用 `python tools/release_audit.py` 复现。
检查本身、白名单、判定标准**均未改动**。

**验证**：修复后在干净候选 worktree 重跑审计 → **FAIL 0 / WARN 0**（§7 复验记录）。

**教训（供后续阶段沿用）**：在受 `absolute-path` 扫描的 `.md` 中引用此类证据时必须脱敏；
"引用证据"不等于"可以复制机读路径"。

---

## 7. 本阶段的 PASS / FAIL / NOT VERIFIED

| 项 | 判定 | 依据 |
|---|---|---|
| 分支、HEAD、版本号已核实并一致 | **PASS** | §1（四处软件版本 `1.1.0` 一致） |
| 版本号决定基于真实发布历史 | **PASS** | §1.2 远端 tag + Releases API，`1.1.0` 确未发布 |
| 候选内容清单完整 | **PASS** | §2 + `docs/audit/CANDIDATE_MANIFEST.tsv`（202 文件） |
| 无关并行改动未混入 | **PASS** | §2.8 三条举证；暂存前后均校验 |
| 存在明确的候选提交 | **PASS** | `624184c`（§3） |
| 从候选提交建立 clean worktree | **PASS** | §4（未提交改动 0） |
| 冻结资产未改变 | **PASS** | §5.3（与 v1.0.1 逐字节一致；`--verify` 13 项异常 0） |
| 候选 worktree 通过 release_audit | **PASS** | §6.3（FAIL 0 / WARN 0） |
| 全量测试 | **PASS** | 候选 worktree：**132 / 132 PASS** |
| 渲染层自检 | **PASS** | **12 / 12 OK** |
| wheel 可重复构建 | **PASS** | §5.6 缺陷 A 修复后验证一致 |
| sdist 字节可重复构建 | **FAIL** | §5.6 缺陷 B，已定位、未修复，转第 2 阶段 |
| 候选 wheel 在全新环境安装并跑通 | **NOT VERIFIED**（本阶段） | 上一轮做过，但用的是 `55df2c6`+13 未提交文件、版本 `1.2.0`，**与本候选不可比**；将在**第 2 阶段**对 `a26171d6…` 重做 |
| 数据限制审查 | **NOT VERIFIED** | 第 1 阶段 |
| 依赖项与安装审查 | **NOT VERIFIED** | 第 2 阶段 |
| 主干功能审查 | **NOT VERIFIED** | 第 3 阶段 |
| 衍生脚本审查 | **NOT VERIFIED** | 第 4 阶段 |
| 最终结果与论文声明检查 | **NOT VERIFIED** | 第 5 阶段 |

**未发现**依赖本阶段 blocker 的后续阶段 —— 本阶段无 blocker（唯一的 FAIL 是 sdist
字节可重复性，不影响功能、安装或科学结果，且已给出替代溯源手段）。

---

## 8. 本阶段遗留与移交第 1 阶段的事项

| # | 事项 | 移交 |
|---|---|---|
| 1 | sdist 字节不可重复（缺陷 B） | 第 2 阶段：评估 pin 构建后端 / 改用其他打包器 / 正式接受并改用 content-sha256 |
| 2 | 候选 wheel 需在**全新环境**重新验证（`a26171d6…`） | 第 2 阶段 |
| 3 | `paper_assets/audit/verify_wheel_install.py:52` 的硬编码解释器 | **不属本候选**；已定位并建议负责人改用 `sys.executable`，本轮无权修改该工作区 |
| 4 | `docs/SAMPLE_PATH_VISUAL_REVIEW.md` 描述的 v1.1.x 旧渲染器已退役 | 第 5 阶段：确认历史文档不会与最终结论冲突（当前已在 `SAMPLE_PATH_INTEGRATION_VALIDATION.md` §9 标注） |
| 5 | 仓库内 `docs/figures/sample_path_review/`（15 MB）与 `sample_path_redesign/`（37 MB）为已取代的未跟踪目录 | 第 5 阶段：决定是否清理（本轮不删，保留审计轨迹） |
| 6 | 数据限制（cohort 来源、许可、可分发边界） | 第 1 阶段主审 |

---

## 9. 阶段结论

* **候选是否已明确**：**是**。候选提交 `624184c`，clean worktree `/tmp/gh_audit/cand`，
  版本 `1.1.0`，完整清单见 `docs/audit/CANDIDATE_MANIFEST.tsv`。
* **是否可进入第 1 阶段**：**可以**。无 blocker；本阶段 1 项 FAIL（sdist 字节可重复性）
  已完整定位并移交第 2 阶段，不影响第 1 阶段（数据限制审查）的开展。
* 本阶段**未** merge、**未** tag、**未** push、**未** 发布；
  工作树中与本功能无关的并行改动**原样保留**。
