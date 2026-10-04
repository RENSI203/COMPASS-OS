# 04 · 衍生脚本与结果生产链审查 ｜ derived scripts and results chain

**阶段**：4 / 6（顺序发布审查）
**日期**：2026-10-04
**候选**：Phase 3 冻结 `75ee0cc` → **本阶段末候选 `c5de83e`**（代码修复 `fcbfa46`，§0.2）
**结论**：**PASS（含 2 处已修复缺陷、5 项未验证/限制）**；无 blocker

> 本阶段未改动冻结分析方案、未为提高指标而修改任何定义。
> 所有重跑都使用**明示输入**与**既有参数**；未新增机制发现。

---

## 0. 范围、候选变更与重跑范围

### 0.1 纳入范围（14 个脚本）

| 组 | 脚本 | 是否支撑本篇软件论文声明 |
|---|---|---|
| 示例 | `examples/quick_start.py`、`minimal_example.py`、`missing_genes_example.py`、`cohort_example.py`、`sample_path_example.py` | **是**（公开 API 的可运行契约） |
| 工具 | `tools/build_assets.py`、`release_audit.py`、`build_fixtures.py`、`audit_duplicate_symbols.py`、`build_review_grid.py` | **是**（发布与溯源） |
| QC / 稳健性 | `validation/missing_gene_stress_test.py`、`run_missing_gene_validation.py`、`summarize_validation.py` | **是**（QC 阈值标定 = 论文的稳健性声明） |
| 平台 mask | `validation/extract_platform_masks.py` | **是**（外部队列覆盖度 = 论文的适用性声明） |

**未纳入**：`downstream_mechanism/**`、`scripts_dev/**` 等机制研究脚本。
理由：它们不支撑本软件论文的软件/稳健性声明，且属只读并行工作区；
**本阶段未启动任何新的机制发现**。

### 0.2 候选变更与需重跑的前序检查

发现 **2 处真实缺陷**（§4.1 / §5.2），已修复并重新冻结：

| 项 | Phase 3 | 本阶段 |
|---|---|---|
| 代码冻结点 | `75ee0cc` | **`c5de83e`**（脚本修复 `fcbfa46` + 审计记录/脱敏文档提交） |

**已重跑**（clean worktree @ `fcbfa46`，未提交改动 0；**先 build 后 test**）：

| 检查 | 结果 |
|---|---|
| 发布审计（**pristine worktree，先于 build**） | **FAIL 0 / WARN 0**（16 项） |
| 全量测试（**先 build 后 test**） | **PASS 150 / FAIL 0 / SKIP 0** |
| 渲染层自检 | **12 / 12 OK** |
| 冻结资产 | **13 项异常 0**（未变） |
| wheel / sdist SHA256 | 见 §0.3 |

> ⚠️ 本阶段**亲自复现了 Phase 0 §5.7 的顺序依赖**：在 clean worktree 中先跑测试得到
> `PASS 149 / SKIP 1`（缺 `dist/` ⇒ `test_built_artifacts_contain_runtime_assets` 跳过），
> 先 `build` 后重跑才得 `150 / 0 / 0`。该现象**已被 Phase 0 记录并固化**，此处为其第二次实证。

### 0.3 产物哈希

在 pristine worktree（检出本阶段末候选 `c5de83e`）用
`SOURCE_DATE_EPOCH=$(git log -1 --format=%ct HEAD)` 规范化构建：

| 产物 | SHA256 |
|---|---|
| `compass_os-1.1.0-py3-none-any.whl` | `b905d8415c415102ac1306a219193f9f8f0799b79ffa23b67b700cd125e2ac60` |
| `compass_os-1.1.0.tar.gz` | `5dfead16c066da7d0c715b3b874178fa0e0e17d09b745ab416b46e4c52a83655` |

> Phase 3 记录的 wheel `fba9ee1…` 已被本阶段 `src/` 之外的脚本修复**取代**（`fcbfa46` 修改了
> `examples/` 与 `validation/`，二者进 sdist）；第 2 阶段须对 `b905d84…` 在全新环境复验。

---

## 1. 脚本—输入—输出—主要结论对应关系

| 脚本 | 明确输入 | 主要输出 | 支撑的结论 | 可重建 |
|---|---|---|---|---|
| `examples/quick_start.py` | 包内 `data/example_*.tsv`（`Path(__file__)` 解析） | `examples/output/quick_start/`（14 文件：pdf/png/tsv/json/txt） | 一键入口可用 | **是** |
| `examples/minimal_example.py` | 同上（单样本） | 控制台 + 输出目录 | 单样本边界行为 | **是** |
| `examples/missing_genes_example.py` | 包内 `example_expression_partial.tsv.gz` | 控制台 + QC | 缺失基因策略可运行 | **是** |
| `examples/cohort_example.py` | 包内 5 例示例 | `cohort_risk.tsv` / KM / concept 热图 | 队列流程、冻结切点 KM | **是**（修复后） |
| `examples/sample_path_example.py` | `--cohort`（默认 GSE39582）+ `--expr-dir`（默认仓库同级和声化缓存） | `docs/figures/sample_path_integration/`（6 组 + grid + 审计 TSV） | `sample_path` 真实端到端 | **是**（需外部数据） |
| `tools/build_assets.py` | `--source-root`（仅重建时需要） | `assets/models/*` + `ASSET_MANIFEST.tsv` | 冻结资产可复现 | `--verify` **是**；重建需源项目 |
| `tools/release_audit.py` | 无（扫工作树） | `validation/results/release_audit.tsv` | 发布卫生 | **是** |
| `tools/build_fixtures.py` | `--source-root` | `tests/fixtures/*` | 测试夹具可重建 | 需源项目 |
| `tools/audit_duplicate_symbols.py` | `--source-root` | `docs/duplicate_*.tsv` | 重复基因审计 | 需源项目 |
| `tools/build_review_grid.py` | 位置参数 `<review_dir>` | 2×3 grid PNG | 人工验收总览 | **是** |
| `validation/extract_platform_masks.py` | `--source-root`（**required**）、`--pheno`（**required**） | `validation/results/platform_masks.json`、`external_coverage_summary.tsv` | 39 队列真实平台覆盖 | 需源项目 |
| `validation/missing_gene_stress_test.py` | `--pheno`、`--expr-dir`（均 required） | `stress_raw.tsv`、`stress_curves.tsv`、报告与图 | 缺失基因稳健性曲线 | **是** |
| `validation/run_missing_gene_validation.py` | `--pheno`、`--expr-dir`（均 required）+ `--seed`/`--repeats` 等 | `raw.tsv`、**`raw<tag>.design.json`**、图 | QC 阈值标定 | **是**（修复后） |
| `validation/summarize_validation.py` | glob `<out>/raw*.tsv` + `_tmp_*.npy` | `stress_cohort_summary.tsv`、`stress_pancohort_summary.tsv`、Figure A–D | 论文稳健性图与 QC 阈值 | **是** |

---

## 2. 逐项检查结果

### 2.1 输入、模型版本、配置、环境与输出路径（item 1）

| 检查 | 结果 |
|---|---|
| 参数是否显式 | 4 个脚本有 `required=True` 的必要输入；`--out` 普遍存在 | **PASS** |
| **硬编码个人目录** | **3 个工具**的 `--source-root` **默认值**为**开发机项目根绝对路径**（`audit_duplicate_symbols` / `build_assets` / `build_fixtures`；原文含用户名，此处按 §7.1 脱敏） | 见下 |
| `--verify` 路径是否受影响 | **不受影响**：实测从 `/tmp` 运行 `tools/build_assets.py --verify` → `校验：13 项 ⇒ 异常 0`（`--source-root` 只在**重建**分支使用） | **PASS** |
| 模型版本来源 | 全部经 `model_manifest.json` 的 `default_model`，**无硬编码模型名** | **PASS** |
| 环境要求 | 需要 `MPLCONFIGDIR`（`~/.config` 不可写）；文档已记 | 记录 |

**裁定**：`tools/` 的硬编码默认值被 `release_audit` 的 `ALLOW_PATH_PATTERNS`（`^tools/`）
**有意豁免**，定位为"开发用、非运行时依赖"。但作为发布脚本，**重建**功能对第三方不可用
（默认值指向作者机器）。**不修**：改成必需参数会改变既有 CLI 契约，而 `--verify`（发布审计
实际调用的路径）不依赖它。记为**限制**（§6）。

### 2.2 能否从明确输入重建；是否误读旧文件/缓存（item 2）

| 检查 | 结果 |
|---|---|
| `missing_gene_stress_test.py` | 每次运行读 `--pheno`/`--expr-dir`，输出到 `--out`（默认 `validation/output`，相对 **cwd**） | **PASS** |
| `run_missing_gene_validation.py` 续跑 | **发现真实缺陷**（§4.1） | **修复后 PASS** |
| `summarize_validation.py` 读取 | glob `raw*.tsv` → concat → **显式丢弃 `error` 行** → `drop_duplicates(keep="last")` | **PASS** |
| `summarize_validation.py` 回退分支 | 无 `_tmp` 中间件时**沿用既有** `coverage_qc_validation.tsv`（有 print 提示） | **限制**（§6） |
| `extract_platform_masks.py` | 显式，无缓存 | **PASS** |

**关键正面发现**：`load_raw` 对失败行与重复行都做了**正确处理**：

```python
if "error" in df.columns:
    bad = df[df["error"].notna()]
    if len(bad): print(f"  ⚠ {len(bad)} 条条件失败（前 3 条）："); print(...)
    df = df[df["error"].isna()].copy()
df = df.drop_duplicates(subset=["mask_type","level","repeat","strategy","cohort"], keep="last")
```

即：**失败被打印出来并从聚合中剔除**，而非静默计入分母 ⇒ 直接回应 item 6（§2.5）。

### 2.3 合并键、笛卡尔积、重复计数、错位、样本量膨胀（item 3）

| 检查 | 实测 | 判定 |
|---|---|---|
| 唯一 merge 点 | `summarize_validation.py:72`，键 = `["mask_type","level","strategy","cohort"]` | **PASS** |
| 是否 1:1 | 两侧均由**同键 groupby 聚合**后产生 ⇒ 键唯一 ⇒ **不可能笛卡尔积** | **PASS** |
| 重复计数 | `stress_cohort_summary.tsv` 共 **508 行** == 508 个唯一 `(mt,level,strategy,cohort)` 组合 | **PASS** |
| 分母口径 | `n_cohorts = sub["cohort"].nunique()`（**队列数**，不是行数 ⇒ 不被 repeats 膨胀） | **PASS** |
| 分母正确性（逐组核对） | 38/38 组的 `n_cohorts` **等于**该组中指标有限的队列数（不一致 **0** 组） | **PASS** |
| 样本量膨胀 | 每队列在 summary 中出现 2–38 次（取决于该队列参与的 mask 轴），**符合设计**而非重复 | **PASS** |

**潜在风险（未实现）**：若某队列的**全部**运行都失败，其行仍带有效 `cohort` 名但指标为 NaN；
`nunique()` 会把它计入 `n_cohorts`。**当前已提交结果中全 NaN 行 = 0**，故**未实现**。
记为**限制**（§6）。

### 2.4 风险方向、事件编码、时间单位、筛选、指标定义、汇总权重（item 4）

**独立复算验证**（真实数据，严格复现脚本 RNG 规则）：

| 量 | 独立复算 | `raw.tsv` 记录 | 差 | 判定 |
|---|---|---|---|---|
| `n_masked` | 784 | 784 | **0** | **PASS** |
| `gene_coverage` | 0.949974477 | 0.949974477 | **0** | **PASS** |
| `risk_spearman` | 0.943076923 | 0.943076923 | **0** | **PASS** |
| `group_agreement` | 0.960000000 | 0.960000000 | **0** | **PASS** |
| `risk_pearson` | 0.981183236 | 0.981183264 | **2.84e-08** | **PASS**（float32 前向差异） |
| `risk_mad` | 0.062477825 | 0.062477764 | **6.09e-08** | **PASS**（同上） |

> **容差声明**：`spearman` / `coverage` / `n_masked` / `group_agreement` 要求**精确相等**
> （纯排序/计数，无浮点前向），实测为 0；`pearson` / `mad` 依赖一次**新的**模型前向，
> 容差取 **1e-6**（float32 批次相关 BLAS 的既有量级，参见 Phase 3 §1），实测 6.1e-08。
> **未为达标而放宽容差。**

| 口径 | 代码事实 | 判定 |
|---|---|---|
| 风险方向 | `risk_spearman = corr(full_risk, pert_risk)`（**扰动前后自相关**，不是与生存的关联） | **PASS** |
| 时间单位 | `os_time_days`（天），传给 `_uno` | **PASS** |
| 事件编码 | `event = (os_event >= 0.5)` → bool | **PASS** |
| 分层切点 | `cut = median_cutoff(default_model)`（**冻结** M2 切点 = 1.0783025125247092） | **PASS** |
| `group_agreement` | `mean((fr>=cut) == (pr>=cut))` —— 用**同一冻结切点**判两侧 | **PASS** |
| 汇总权重 | pan-cohort 对**队列取中位数**（`cohort 为独立单位`），**非**样本量加权；代码注释明示 | **PASS** |
| 队列纳入 | `len(ids) < 20` 或表达文件缺失 ⇒ **打印原因并跳过**，不静默丢弃 | **PASS** |

### 2.5 随机种子、mask 生成、重复次数、失败条件、分母、纳入规则（item 5）

| 项 | 事实 | 判定 |
|---|---|---|
| 种子 | `--seed` 默认 `20261003`；`default_rng([seed, round(level*10000), rep, mask_type_code])` | **PASS** |
| mask 确定性 | 由 `(seed, level, repeat, mask_type)` **完全决定**；同一条件可逐位复现（§2.4 已实证） | **PASS** |
| mask 池 | `pool = 表达矩阵列`（非硬编码词表）；`non_sig` 池由 `G_all − G_sig` **程序计算** | **PASS** |
| 重复次数 | `--repeats` 默认 20；`level == 0` 时只跑 1 次（`range(1)`）—— 合理（零扰动无随机性） | **PASS** |
| 失败条件 | `len(ids) < 20` 跳过队列；`predict` 异常 ⇒ **写入带 `error` 的行**再 `continue` | **PASS** |
| 分母 | 队列数（§2.3） | **PASS** |
| **溯源可追溯性** | `design.json` 每次运行**被覆盖** ⇒ 只描述一次运行 | **FAIL → 已修复**（§4.1） |

### 2.6 是否选择性保留好结果 / 丢弃失败 / 把未运行条件当完成（item 6）

| 检查 | 结果 | 判定 |
|---|---|---|
| 失败是否记录 | `except` 分支写入完整行（含 `error` 文本与 `n_masked`）后 `continue` —— **不是静默跳过** | **PASS** |
| 失败是否进入聚合 | `load_raw` 显式 `df[df["error"].isna()]` 剔除，并 **print 前 3 条** | **PASS** |
| 是否把未运行当完成 | 结果表只含**实际产出**的行；`n_cohorts` 由数据推出，非预设 | **PASS** |
| 是否有"按结果挑队列" | 无：队列纳入只依赖 `len(ids) >= 20` 与文件存在性 | **PASS** |
| 失败数可查 | 已提交 `stress_cohort_summary.tsv` 全 NaN 行 = **0**（即无"成组失败"） | **PASS** |

> 需注意：`summarize` 的 print 只显示**前 3 条**失败，完整失败清单需自行查 `raw.tsv`。
> 这是**披露充分性**问题，不是选择性保留。记为**限制**（§6）。

### 2.7 图表是否直接来自可追溯结果（item 7）

| 脚本 | 读文件处 | 疑似硬编码数值 | 判定 |
|---|---|---|---|
| `summarize_validation.py` | 4（`raw*.tsv`、`external_coverage_summary.tsv`、`coverage_qc_validation.tsv`） | **0** | **PASS** |
| `missing_gene_stress_test.py` | 3 | **0** | **PASS** |
| `tools/build_review_grid.py` | 2（PNG + `_meta`/selection） | **0** | **PASS** |

Figure A–D 的输入依次为 `pan`（由 `raw*.tsv` 聚合）、`external_coverage_summary.tsv`；
**未发现**手工改分数、改排序或改数据的分支。`build_review_grid.py` 只做**拼贴**
（`imread` + `imshow`），其文件头亦声明 "只拼贴，不重绘" ⇒ **PASS**。

### 2.8 能否从仓库外/明确工作目录运行；参数与报错可用性（item 8）

**实测（cwd = `/tmp`）**：

| 脚本 | 结果 |
|---|---|
| `tools/build_assets.py --verify` | ✅ `校验：13 项 ⇒ 异常 0` |
| `tools/release_audit.py` / `build_review_grid.py` | ✅ `--help` 正常 |
| `examples/quick_start.py` | ✅ 运行并写出 14 文件（**注意：无参数解析器，`--help` 也会执行**） |
| `examples/minimal_example.py` | ✅ 且**明确提示** `risk_group=None` |
| `examples/missing_genes_example.py` | ✅ 且提示示例数据缺失 ~75 % 基因、真实区间以 stress test 为准 |
| `examples/cohort_example.py` | ❌ → **已修复**（§5.2） |
| `examples/sample_path_example.py` | ✅（需外部表达缓存，`--expr-dir` 可指定） |

---

## 3. 重跑范围（小规模、关键链路）

**原则**：只重跑足以确认指标正确性的小规模真实条件；全量昂贵分析**不重跑**（无缺陷、来源明确）。

| 重跑项 | 规模 | 目的 | 结果 |
|---|---|---|---|
| `run_missing_gene_validation.py` | 2 队列 × 25 例 × 2 repeats × `global` 轴 | 验证**指标定义**与**溯源/续跑** | ✅ 48 行产出；指标独立复算一致（§2.4） |
| 同上（同参数续跑） | — | 验证正确续跑 | ✅ "续跑：已有 48 行结果（溯源参数一致）" |
| 同上（改 `--seed` / `--max-per-cohort`） | — | 验证守卫 | ✅ 拒绝并同时显示两套参数 |
| `examples/cohort_example.py` | 5 例（包内） | 验证示例可运行 | ✅ 修复后退出码 0 |
| `examples/quick_start.py` 等 3 个 | 包内示例 | 验证示例链 | ✅ 全部退出码 0 |
| 全量回归 | 150 项 | 确认修复无回归 | ✅ **150 / 0 / 0** |

**未重跑（并说明理由）**：

| 未重跑项 | 理由 |
|---|---|
| 12 队列 × 20 repeats 全量曲线轴 | 昂贵（数小时）；来源明确（`design.json` + `raw.design.json` 现已配对）；指标定义已由小规模重跑独立验证 |
| 38 队列经验轴 | 同上 |
| `extract_platform_masks.py` 全量 | 依赖只读源项目 GEO 注释；`external_coverage_summary.tsv` 已提交且与 `qc_config.json` 的 `empirical_cohorts: 38` 自洽 |
| `build_assets.py` 重建 | 需源项目 checkpoint；`--verify` 已确认包内资产与清单一致 |

---

## 4. 缺陷 1：续跑缓存键不完整 + 溯源文件被覆盖 — **已修复**

### 4.1 现象与根因

**现象**：`run_missing_gene_validation.py` 支持断点续跑。续跑键为

```python
done = {(r.mask_type, float(r.level), int(r.repeat), r.strategy, str(r.cohort)) ...}
```

**不含** `--seed` / `--max-per-cohort` / `--repeats` / `--cohorts` / `--mask-types`。
⇒ 用**不同参数**重跑（不加 `--restart`）时，已完成的格子被**静默跳过**，
新结果与旧结果**混进同一个 `raw.tsv`**，而两者定义不同。

**加重因素**：溯源写在**单个** `out/design.json`，**每次运行被覆盖**。

**已实现的后果（有举证）**：

```
validation/results/design.json  →  cohorts: 12, n_loaded: 12, max_per_cohort: 60, repeats: 20
validation/results/stress_cohort_summary.tsv → 唯一队列: 38
stress_pancohort_summary.tsv:
   empirical   n_cohorts = 38
   global      n_cohorts = 12
   non_sig     n_cohorts = 12
   signature   n_cohorts = 12
```

⇒ **`design.json` 只记录了 12 队列的曲线轴，而结果集同时含 38 队列的经验轴。**
38 队列那一次的**队列纳入规则没有溯源文件** —— 这正是 item 5 要求的可追溯性缺口。

### 4.2 影响

| 面 | 影响 |
|---|---|
| 数值本身 | **不受影响**：每个轴内部自洽，`n_cohorts` 逐组核对正确（§2.3） |
| **溯源** | `design.json` 与结果**不对应**；若再跑一次经验轴，12 队列的纳入规则将**永久丢失** |
| 复现风险 | 用不同 seed/max_per_cohort 续跑会**静默混合定义**（latent，本已提交结果未触发） |
| QC 阈值 | 阈值由 12 队列曲线轴标定（Phase 1 §9 已记 `curve_cohorts: 12`）—— 数值有效，但溯源不足 |

### 4.3 修复

1. **参数指纹** `design_signature()`，覆盖影响数值的全部参数
   （`cohorts` / `max_per_cohort` / `repeats` / `seed` / `mask_types` / `nonsig_levels`）；
2. **配对侧车** `raw<tag>.design.json` 与 `raw<tag>.tsv` **同 tag** 写入，随 raw 一起走；
3. **续跑守卫**：侧车指纹 ≠ 本次指纹 ⇒ **拒绝续跑**，同时打印两套参数并提示用 `--restart`；
4. `out/design.json` **保留**（既有读者 `summarize_validation.py:395` 仍可读），
   但标注 `_note: "LAST RUN ONLY — ..."` 指向配对侧车。

### 4.4 验证

```
A) 同参数续跑          → "续跑：已有 48 行结果（溯源参数一致）"          ✅
B) --seed 999          → 拒绝，显示已有/本次两套参数                     ✅
C) --max-per-cohort 30 → 拒绝，显示已有/本次两套参数                     ✅
D) 侧车内容            → params 含全部 6 个参数 + raw_file + tag         ✅
```

---

## 5. 缺陷 2：`cohort_example.py` 在自己的示例数据上崩溃 — **已修复**

### 5.1 现象与根因

```
File "examples/cohort_example.py", line 45, in main
    tab["risk_group(默认模型)"] = res.risk_group[res.default_model]
TypeError: 'NoneType' object is not subscriptable
```

* `risk_group` 仅在队列 **≥ `min_cohort_for_stratification`（默认 30）** 时给出
  （`api.py:275`）。
* 而 `examples/example_data/example_expression.tsv.gz` 只有 **5 例**
  ⇒ `risk_group` **必为 None** ⇒ 崩溃是**确定性的**，与数据内容无关。

### 5.2 影响

**随包发布的示例脚本无法运行**。用户复制该示例会直接得到
`TypeError: 'NoneType' object is not subscriptable`，而不是被引导到正确的边界行为。
属 v1.0.1 引入"小队列不产生分层"后的**示例未同步**回归 ——
同一目录下的 `minimal_example.py` **已正确处理**该情形，二者**不一致**。

### 5.3 修复

按 `minimal_example.py` 的既有模式显式处理：`risk_group is None` 时写入**占位说明**
（`"(未分层：队列 < 30 例)"`）并打印原因，明确该列**不得用于分层结论**。

### 5.4 验证

修复后从 `/tmp` 运行 → **退出码 0**，正常产出 `cohort_risk.tsv` / concept 热图，
并打印分层切点说明；KM 分支因 `os_event.sum() = 1 < 2` 按既有逻辑跳过并提示。

---

## 6. 未验证项与剩余限制

| # | 项 | 判定 | 说明 |
|---|---|---|---|
| 1 | 全量 12 队列 / 38 队列轴**未重跑** | **NOT VERIFIED（重跑）** | 昂贵；指标定义已由小规模重跑独立验证；现已可溯源 |
| 2 | `design.json` 的历史缺口 | **限制（已缓解）** | 新增配对侧车；但**已提交的 `design.json` 仍是 12 队列版**，38 队列那次的纳入规则**无法从仓库恢复**（其 `params` 未曾记录） |
| 3 | `tools/*` 的 `--source-root` 默认指向作者机器 | **限制** | `--verify` 不受影响（实测从 `/tmp` 通过）；**重建**功能对第三方不可用 |
| 4 | `summarize_validation.py` 无 `_tmp` 时沿用旧 `coverage_qc_validation.tsv` | **限制** | 有 print 提示，但可能造成**图与表不同批次** |
| 5 | 失败清单只 print 前 3 条 | **限制** | 完整失败需查 `raw.tsv`；不是选择性保留 |
| 6 | 全 NaN 行会抬高 `n_cohorts`（潜在） | **限制（未实现）** | 当前 0 行全 NaN；若未来成组失败则会虚高 |
| 7 | `examples/quick_start.py` 无参数解析（`--help` 也执行） | 记录 | 示例属性；会写入仓库内 `examples/output/` |
| 8 | 所有示例输出**固定在仓库内** `examples/output/` | 记录 | 只读 checkout 下无法运行；需可写工作副本 |
| 9 | `validation/inputs/`（pheno、表达缓存）**未纳入候选**（gitignored） | **限制** | 干净 clone **无法**直接重建验证结果，必须自备输入；`validation/README.md` 已列出来源路径 |

### 7.1 本阶段的自我命中与修复（与 Phase 0 §6.5 同类）

**现象**：本文件初稿在 §2.1 表格中**逐字引用**了工具里那条开发机绝对路径作为证据，
导致**本文件自身**被 `release_audit.py` 判为 `absolute-path` FAIL：

```
$ cd <pristine candidate worktree> && python tools/release_audit.py
检查项 17：FAIL 1 / WARN 0
=== FAIL ===
absolute-path   docs/audit/04_DERIVED_SCRIPTS.md:86   ...   FAIL
```

**这是同一类错误的第二次发生**（Phase 0 §6.5 已记录过一次并写明教训：
"在受 `absolute-path` 扫描的 `.md` 中引用此类证据时必须脱敏"）。
本次未吸取该教训，**记录在案**。

**修复**：把路径字面量改为文字描述（"开发机项目根绝对路径"），
保留文件、位置与结论。检查本身、白名单与判定标准**均未改动**。

**验证**：修复后在 pristine worktree 重跑 → **FAIL 0 / WARN 0**（16 项）。

> **给后续阶段的硬性提醒**：`docs/audit/**.md` 一律**不得**出现机读绝对路径；
> 引用证据时只写"开发机路径/用户名已被脱敏"，把原始字符串留给读者自行运行审计复现。

---

## 7. 与候选无关的并行改动

`paper_assets/**`、`downstream_mechanism/**`、`scripts_dev/**` 等属并行工作区：
本阶段**只读检查**，**未**修改、**未**纳入提交。`git status` 确认候选提交未触碰它们。

---

## 8. PASS / FAIL / NOT VERIFIED 汇总

| # | 事项 | 判定 |
|---|---|---|
| 1 | 输入/模型版本/配置/环境/输出路径明确；硬编码个人目录 | **PASS**（`--verify` 不受影响）；`tools/` 默认值 = **限制** |
| 2 | 可从明确输入重建；不误读旧文件/缓存 | **PASS（修复后）**（缺陷 1） |
| 3 | 合并键正确；无笛卡尔积/重复计数/错位/膨胀 | **PASS**（508 行 = 508 唯一键；38/38 组分母正确） |
| 4 | 风险方向/事件/时间/筛选/指标/权重一致 | **PASS**（独立复算：4 项精确相等，2 项 ≤6.1e-08） |
| 5 | 种子/mask/重复/失败/分母/纳入可追溯 | **PASS（修复后）**（缺陷 1） |
| 6 | 无选择性保留；失败被显式剔除且披露 | **PASS** |
| 7 | 图表直接来自可追溯结果，无手工改数 | **PASS**（3 脚本 0 处硬编码数值） |
| 8 | 可从仓库外运行；参数与报错可用 | **PASS（修复后）**（缺陷 2）；示例输出固定仓库内 = 限制 |

**缺陷 2 处，均已修复并验证。Blocker：无。**

---

## 9. 阶段结论

* 建立 **14 个脚本**的完整对应关系（脚本—输入—输出—结论—可重建性），
  并按用户要求**排除**机制研究脚本（未启动新机制发现）。
* 发现并修复 **2 处真实缺陷**：
  ① 续跑缓存键不完整 + `design.json` 被覆盖导致的**溯源缺口（已实际发生：12 vs 38 队列）**；
  ② 随包示例 `cohort_example.py` 在自身数据上**确定性崩溃**。
* **独立复算验证了指标定义**：`n_masked` / `gene_coverage` / `risk_spearman` /
  `group_agreement` **精确相等**，`risk_pearson` / `risk_mad` 差 ≤ **6.09e-08**（容差 1e-6，
  float32 前向差异）—— **未放宽容差**。
* 重跑限于**小规模关键链路**（2 队列 × 25 例 × 2 repeats）；全量轴**未重跑**并说明理由。
* 修复后：发布审计 **FAIL 0 / WARN 0**；全量测试 **150 / 0 / 0**（先 build 后 test）；
  冻结资产 **13/13 异常 0**。
* 本阶段**第二次实证**了 Phase 0 记录的顺序依赖（未先 build ⇒ SKIP 1）。
* **未** merge、**未** tag、**未** push、**未**发布。
