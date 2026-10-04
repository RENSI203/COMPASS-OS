# 03 · 主干功能与数值一致性审查 ｜ core functions and numerical consistency

**阶段**：3 / 6（顺序发布审查）
**日期**：2026-10-04
**候选**：Phase 1 冻结 `28cdabf` → **本阶段修复提交 `75ee0cc`**（§0.2）
**结论**：**PASS（含 1 处已修复缺陷、2 项明确限制）**；无 blocker

> 本阶段未改动冻结科学模型、未重新 fit、未为"完全一致"调整任何数值或容差。
> 全部数值比较使用**事先声明**的容差（§1），并记录**实测最大误差**。

---

## 0. 方法与候选变更

### 0.1 真实数据来源（可追溯，未改标签）

| 项 | 内容 |
|---|---|
| 表达矩阵 | `data_tcga/processed/01_expression/expression_tcga_log2tp1.tsv.gz`（TCGA 泛癌，log2(TPM+1)，样本 × 15,672 基因） |
| 临床表 | `data_tcga/raw/cbio_pancan_clinical.tsv`（cBioPortal 泛癌，含 `AGE` / `SEX` / `AJCC_PATHOLOGIC_TUMOR_STAGE`） |
| 抽样方式 | 取表达矩阵前 2,500 行 ∩ 临床表，按癌种各取**前 40 例临床完整**（AGE/SEX/STAGE 均非空）的样本 |
| 得到的夹具 | `github/validation/inputs/phase3_multicancer_{expr.tsv.gz,clinical.tsv}`（**gitignored**，不入候选） |
| 规模与构成 | **n = 120**：BRCA 40 / LUAD 40 / THCA 40（真实 `CANCER_TYPE_ACRONYM`，**未改写任何癌种标签**） |
| 临床完整性 | AGE 中位 55.5 / 64.0 / 52.5；STAGE ∈ {I–IV} → 1–4；SEX ∈ {Female, Male}；**三项零缺失** |

**多癌种判定**：另构造**真·多癌种子集**（每癌种各 10 例，n=30）用于预算验证；
单癌种场景用同一批样本但 `cancer_type` 全部为 `BRCA`（该值取自真实标签，非伪造）。

> ⚠️ **本阶段自身的一次纠错**：首轮把 `E.iloc[:30]` 当作"多癌种"测试，
> 但该切片**全部是 BRCA**（夹具按癌种连续排列），实际测的是单癌种。
> 已改为按癌种分层抽样重测（§6），结果随之变化 —— 记录在案。

### 0.2 候选变更与需重跑的前序检查

本阶段发现 **1 处真实缺陷**（§5.2），已修复 ⇒ Phase 1 冻结按规则失效并重新冻结：

| 项 | Phase 1 | 本阶段 |
|---|---|---|
| 代码冻结点 | `28cdabf` | **`75ee0cc`** |

**已重跑的检查**（clean worktree @ `75ee0cc`，未提交改动 0）：

| 检查 | 结果 |
|---|---|
| 发布审计 | **FAIL 0 / WARN 0**（16 项） |
| 全量测试（先 build 后 test） | **PASS 150 / FAIL 0 / SKIP 0**（148 → 150） |
| 渲染层自检 | **12 / 12 OK** |
| 冻结资产 | **13 项异常 0**（未变） |
| wheel SHA256 | `fba9ee12138bf0b349bd26809c3c649af8b149998d3c144ad8a3174f8fe6cb33` |
| sdist SHA256 | `414ed1fd205199c2dcec846ef6f35f54daadbb6d2cf641653d51311c67be57b7` |

> Phase 1 记录的 wheel `6319392…` 已被取代；第 2 阶段须对 `fba9ee1…` 在全新环境复验。

---

## 1. 事先声明的容差与实测最大误差

**容差依据**：模型为 **float32** 前向（`torch` CPU）；线性代数（Cox 设计矩阵、`@ beta`）
在 **float64** 中完成；BLAS 的**批次相关**归约顺序会在 float32 下引入 ~1e-6 量级的差。

| 比较项 | 容差 | 依据 | **实测最大误差** | 判定 |
|---|---|---|---|---|
| 同一输入**重复**调用 | **0** | 无随机性来源 | **0.000e+00** → 实测为 0，可称**逐位一致** | **PASS** |
| `get_representation` vs `predict` 的 132 层 | — | 同一次前向 | **0.000e+00** | **PASS** |
| `get_representation` vs `predict` 的 43 层 | — | 同一次前向 | **0.000e+00** | **PASS** |
| 单样本 vs 批量 η | `1e-4` | float32 批次相关 BLAS | **3.699e-07** | **PASS** |
| 输入行**乱序重排**后 η | `1e-4` | 同上 | **7.304e-06** | **PASS** |
| 完整分解 signed sum vs η | `1e-9` | float64 求和顺序 | **2.220e-16** | **PASS** |
| 图中 η vs 主 API | 0 | 同一数值 | **0.000e+00** | **PASS** |
| `sample_path` 出图后**再**预测 η | **0** | 出图不得影响推理 | **0.000e+00** → 实测为 0，可称**逐位一致** | **PASS** |
| M1（不依赖临床）在 `clinical=None` 前后 | `1e-12` | 临床不进入 M1 | **0.0** | **PASS** |

**未为任何一项放宽容差。** 唯一的"0 容差"项（重复调用、出图后再预测）实测确为 0。

> **术语纪律（发布口径）**：只有**实测差值恰为 `0.000e+00`** 的比较才可写"**逐位一致**"；
> 任何**非零**误差（本表 `3.699e-07`、`7.304e-06`、`2.220e-16` 等）
> 一律写"**在预设容差内数值一致**"，**不得**写"逐位一致"或"bitwise identical"。
> 本表已按此口径标注。

---

## 2. 一键式入口、参数化入口与默认模型

| 项 | 实测 | 判定 |
|---|---|---|
| 公共 API（`__all__`，16 项） | `analyze` / `predict` / `get_representation` / `check_robustness` / `sample_path` + 4 个结果类 + 5 个异常 + `__version__` | **PASS** |
| **默认模型** | `api.default_model()` → **`M2`**；`predict(model="default")` 返回列 `['M2']` | **PASS** |
| 参数化入口 | `predict(model="M0,M1,M2,M3")` → 列**按规格顺序** `['M0','M1','M2','M3']` | **PASS** |
| 一键入口 | `analyze()` 返回 `AnalysisResult`，`an.risk` 为**默认模型**（M2）的 `Series` | **PASS** |
| 模型元信息 | `api.model_info()` / `available_models()` 可用（**不在 `__all__`**，属 Level-3 接口） | **PASS** |
| 未知模型 | `survival.lock_path` → `InputError` | **PASS** |

> ⚠️ **与任务描述不一致之处（第 2 次记录）**：本阶段任务书写"默认 M3"，
> 但**代码事实是 `M2`**（`model_manifest.json` 的 `default_model`）。
> Phase 1 §6 已记录同一事实。本报告以**代码为准**，未按任务书改写。

---

## 3. 对齐、尺度、缺失策略、临床编码、索引保持

| 项 | 实测 | 判定 |
|---|---|---|
| 基因对齐 | 120 × 15,672，列名 = 词表符号，`reindex` 按名匹配 | **PASS** |
| 输入尺度 | `input_scale="log2_tpm1"`（本阶段真实数据即 log2 尺度）；尺度自检**双向**生效（Phase 1 修复） | **PASS** |
| 缺失策略 | 真实 TCGA 泛癌队列对 15,672 词表的覆盖良好；`reference` 默认 | **PASS** |
| 临床编码 | `age` 岁 / `sex` 数值 / `stage ∈ 1–4`；缺整表 ⇒ 三项冻结填补 | **PASS** |
| **索引保持** | 输入乱序重排后，`risk.index` **完全按输入顺序**返回；数值差见 §1 | **PASS** |
| 重复 `sample_id` | `InputError`（Phase 1 已验证） | **PASS** |

**实测（乱序重排）**：`perm = rng.permutation(120)` → `pred.risk.index == E.iloc[perm].index` ✓，
风险差 7.304e-06（float32 批次相关，在 1e-4 内）。

---

## 4. 132 signature / 43 concept / QC / 风险 / 分组 / 生存

| 项 | 实测 | 判定 |
|---|---|---|
| 132 signature | `pred.signature_scores.shape = (120, 132)`；`get_representation` 与之实测差 **0.000e+00**（**逐位一致**） | **PASS** |
| 43 concept | `pred.concept_scores.shape = (120, 43)`；实测差 **0.000e+00**（**逐位一致**） | **PASS** |
| QC 输出 | 覆盖度、逐样本 NaN 统计、`qc_tier_*`、warnings 齐备 | **PASS** |
| 风险预测 | 4 模型列齐备，顺序固定 | **PASS** |
| 冻结分组 | `risk_group` 取值 `{low: 90, high: 30}` —— 见下方限制 | **PASS（行为）** |
| 生存输出 | `survival.survival_probability(η, [365, 1095])` → `(120, 2)`，范围 0.1904–0.9976 | **PASS（定义）** |

**关于"90 low : 30 high"的正确解读（更正）**

以单癌种训练得到的冻结 `median_cutoff`（M2 = 1.078303）切分多癌种混合队列得到 90:30。
**该比例本身不构成切点失效的证据**：冻结切点在**训练队列**上定义，应用到**新队列**时
**本来就不需要**产生 50:50 —— 风险分布随癌种构成平移是预期现象。

⇒ **不得**用"偏离 1:1"论证切点失效或分层不可用。
**真正需要验证**的是 (a) 分层后两组的**预后区分**与 (b) **适用范围**；
前者所需的 transport 源表本轮**未定位到**（标为来源报告值、未独立复核），
后者仅有 `stratified_vs_pooled_external.tsv`（M2 分层 ≈ 合并，CI 均跨 0）。
详见 `05_MAIN_RESULTS_CLAIMS.md` §C6b。

---

## 5. 边界与错误输入

### 5.1 已验证的边界（真实数据）

| 场景 | 行为 | 判定 |
|---|---|---|
| 单样本调用（`E.iloc[[7]]`） | 返回 1×4 风险；与批量对应行差 **3.699e-07** | **PASS** |
| 批量 120 例 | 正常 | **PASS** |
| `clinical=None` | 三项全部冻结填补并**显式告警**："Age 120/120 samples (reference value 60.0); …" | **PASS** |
| 临床缺失时 M1 | 风险**完全不变**（0.0）—— M1 不含临床，符合设计 | **PASS** |
| 输出顺序 | 严格等于输入索引顺序 | **PASS** |

### 5.2 缺陷 5：**非数值临床值被静默当作"缺失"** — 已修复

* **现象**：真实 cBioPortal 的 `SEX` 取值为 `"Female"` / `"Male"`（字符串）。
  实测 `_clinical_frame` 对字符串 Sex 做 `pd.to_numeric(errors="coerce")` ⇒ **全部变 NaN**
  ⇒ 走冻结填补 `0.0`。告警文本仅为
  *"Clinical variables filled with frozen training reference values: Sex 120/120 samples (reference value 0.0)"*
  —— 与"用户根本没提供临床"**完全同一条消息**。
* **影响**：用户按 cBioPortal 原始编码提供临床，会得到一个**全部被填补**的临床块，
  却看起来只是"数据缺失"；`sex` 的方向性信息被静默丢弃，而 risk 仍正常返回。
* **定性**：**编码错误**（用户提供了值，只是格式不对）与**数据缺失**（用户没提供）
  是两件不同的事，必须分开报告 —— 属"实现与声明不一致"。
* **修复**：`_clinical_frame` 新增 `diag` 出参，逐字段记录
  `n_provided` / `n_unparseable`；`_clinical_fill_report` 在 `n_unparseable > 0` 时
  追加**独立警告**，点名列名并说明要求的编码。
* **验证**（`tests/test_input_contract.py`，2 项新测试 + 真实数据）：

```
Sex 30/30 samples (reference value 0.0);
Sex: 30/30 provided values could not be parsed as numeric and were treated as missing.
Clinical columns must be numeric (age = years, stage = 1-4, sex = numeric code);
e.g. Female/Male must be encoded before calling.
```

  数值编码（`sex = 0/1`）时**不产生**该告警；`diag` 实测
  `{'n_provided': 4, 'n_unparseable': 4}`（字符串）vs `{'n_provided': 4, 'n_unparseable': 0}`（数值）。

* **未改变**：数值契约本身未动（库仍要求数值），**未**新增字符串自动映射
  （那会引入 `M/F` 等歧义映射，属功能扩张而非一致性修复）。

### 5.3 其他边界

| 场景 | 行为 | 判定 |
|---|---|---|
| `stage = 0` 或 `5` | → NaN → 冻结填补（`between(1,4)`） | **PASS**（既有语义已锁定测试） |
| `Age` 无生理范围校验 | 任意有限值被接受 | **限制**（Phase 1 已记，本阶段未改）→ **NOT VERIFIED** |
| 长度不符 / 未知癌种 / `NORMAL` / `±Inf` | Phase 1 已修复并测试 | **PASS** |

---

## 6. sample_path：自动接线、真实分数、imputed、PC 聚合、五种预算、风险位置、统一导出

### 6.1 真实多癌种端到端（首次完成）

真实 TCGA 多癌种子集（BRCA/LUAD/THCA 各 10，n=30），**真实临床三项完整**：

| 模型 / 队列 | concepts | aux 节点 | 总 predictor |
|---|---|---|---|
| M1 多癌种 | **16** | — | 16 |
| M2 多癌种 | **12** | Age, Sex, Stage, **Cancer type** | 16 |
| M3 多癌种 | **11** | Age, Sex, Stage, **Cancer type**, **PC1–PC10** | 16 |
| M2 单癌种 | **13** | Age, Sex, Stage | 16 |
| M3 单癌种 | **12** | Age, Sex, Stage, **PC1–PC10** | 16 |

**与规格表逐一吻合** ⇒ **PASS**。

**Cancer type 节点聚合正确性**（真实多癌种）：

```
Cancer type 节点 value = BRCA
contribution = -1.345218  ==  Σ(CT_* 列 contribution) = -1.345218   一致
```

### 6.2 真实临床下的 imputed 标识

真实临床完整的样本上，`sample_path` 实测：

```
imputed: [{Age: False}, {Sex: False}, {Stage: False}, {Cancer type: False}, {PC1–PC10: False}]
```

即**真实观测值被正确标为 `imputed=False`**（不误报填补）⇒ **PASS**。
配合 Phase 1 的 `clinical=None` 路径（三项 `imputed=True`），**两个方向都已验证**。

### 6.3 数值一致性（真实数据）

| 项 | 实测 |
|---|---|
| 图中 η vs 主 API | **0.000e+00**（逐位一致） |
| 完整分解 signed sum vs η | **2.220e-16**（容差 1e-9） |
| PC1–PC10 单一节点，contribution | 0.5215600415595891（= 十个冻结设计列贡献之和） |
| 分解行数 | M3 = 88 行（含未展示 concept、被隐藏的 CT_* 与 offset） |
| **图上截断不影响预测** | 出图后再 `predict`，η 差 **0.000e+00** |

⇒ **PASS**。

### 6.4 统一导出

`save_report()` 产出 **13 个文件**（`.pdf ×2`、`.png ×2`、`.tsv ×6`、`.json`、`.txt`），
即报告、图与数据表同源 ⇒ **PASS**（§7 的一致性核对见下）。

---

## 7. 报告 / 导出与原始输出的一致性

| 项 | 实测 | 判定 |
|---|---|---|
| `analyze().save_report(dir)` | 返回文件字典，13 个文件 | **PASS** |
| `an.risk` | `Series`（默认模型 M2），长度 = 样本数 | **PASS** |
| `an.concept_scores` / `an.signature_scores` | 与 `predict` 实测差 **0.000e+00**（逐位一致，§1） | **PASS** |
| KM 曲线 | `save_report` 产出 KM 相关 `.png`/`.pdf`，其输入 η 即 `an.risk` | **PASS（同源）** |
| `save_tables` / `to_dataframe` | **不存在该 API**（本阶段探针的猜测，非缺陷） | n/a |

**限制（NOT VERIFIED）**：本阶段**未逐像素/逐行**核对导出的 `.tsv` 与内存对象，
只核对了**上游数值同源**（导出所用 η 与 132/43 层即 `an.*`，且与 `predict` 实测差 **0.000e+00**）。
若需更强的"导出文件 == 内存数值"断言，应补充读回比对（未做，记 NOT VERIFIED）。

---

## 8. gene-level 额外前向：时间与内存成本

`sample_path` 需要一次**额外的** `extract(with_gene_level=True)` 前向
（主 API 只返回 132/43 两层），本阶段实测其成本：

| n | `predict(M3)` | gene-level 额外 | **比值** |
|---|---|---|---|
| 16 | 0.48 s | 0.51 s | **1.06×** |
| 64 | 1.72 s | 2.15 s | **1.25×** |

* **额外峰值内存**：**14.7 MB**（n=16，`tracemalloc`）。
* **结论**：gene-level 前向使总推理时间增加约 **6 %–25 %**，比值随 n 增大而下降
  （固定开销被摊薄）。对交互式单样本/小队列使用可接受。
* **优化决定：不做优化。** 该前向是 `sample_path` 展示真实 gene-token 分数的**唯一来源**
  （不得用 attention/TPM 替代），任何"优化"（如缓存共享、降精度、抽稀基因）
  都会改变计算定义或引入状态 —— 超出"只做必要优化"的边界。**记录成本，不改实现。**

---

## 9. PASS / FAIL / NOT VERIFIED 汇总

| # | 事项 | 判定 |
|---|---|---|
| 1 | 一键/参数化入口、默认模型、M0–M3 选择、癌种要求 | **PASS**（默认实为 **M2**，已两次记录） |
| 2 | 对齐、尺度、缺失策略、临床编码、索引保持 | **PASS** |
| 3 | 132/43 提取、QC、风险、冻结分组、生存输出 | **PASS** |
| 4 | 单样本/批量、顺序、错误输入、边界 | **PASS** |
| 5 | KM/报告/导出与原始输出一致 | **PASS（同源）**；逐文件读回比对 **NOT VERIFIED** |
| 6 | `sample_path` 接线、gene score、imputed、PC、五种预算、风险位置、统一导出 | **PASS**（**真实多癌种首次完成**） |
| 7 | 完整分解含隐藏贡献、sum == η、截断不影响预测 | **PASS** |
| 8 | gene-level 提取前后风险/132/43/资产不变 | **PASS**（132/43 逐位 0.0；风险 0.0；资产 13/13 未变） |

### 本阶段缺陷与修复

| # | 缺陷 | 严重度 | 状态 |
|---|---|---|---|
| 5 | 非数值临床值（如 `Female`/`Male`）被静默当作"缺失"，与"未提供"同一告警 | **中**（方向性信息静默丢失，风险仍正常返回） | **已修复并验证** |

### 剩余限制

| # | 限制 | 判定 |
|---|---|---|
| 1 | `Age` 无生理范围校验（任意有限值被接受并标准化） | **NOT VERIFIED** |
| 2 | 导出文件与内存数值的**逐文件读回**比对未做（仅核对上游同源） | **NOT VERIFIED** |
| 3 | 绝对生存 `S(t)` 的校准无独立证据（Phase 1 已记；本阶段只核对**定义与调用一致**） | **NOT VERIFIED** |
| 4 | 冻结 `median_cutoff` 在混合队列上产生不均衡分组（实测 90:30）。**该比例本身不是切点失效的证据**；待验证的是**分层后的预后区分**与**适用范围**（transport 源表本轮未定位到 ⇒ 来源报告值、未独立复核） | **限制（待专门验证）** |
| 5 | 非数值临床**自动映射**（如 `M/F`）**故意未实现** —— 避免歧义映射；由用户在调用前编码 | **设计决定** |
| 6 | `model_info` / `available_models` 不在公共 `__all__`，属 Level-3 接口 | 记录 |

### Blocker

**无。**

---

## 10. 移交下一阶段

| # | 事项 | 来源 |
|---|---|---|
| 1 | 用 `fba9ee1…` 在**全新环境**复验安装与端到端 | 本阶段 §0.2 |
| 2 | `Age` 范围校验、sdist 字节可重复性、发布顺序固化 | Phase 0/1 遗留 |
| 3 | 衍生物（脚本/图表）审查 | 第 4 阶段 |
| 4 | `risk_group` 在混癌种场景的正式适用性 | 本阶段 §4 |
| 5 | 导出文件逐文件读回比对（若论文要引用导出表） | 本阶段 §7 |

---

## 11. 阶段结论

* 8 项指定事项全部完成，**首次在真实 TCGA 多癌种队列（BRCA/LUAD/THCA）+
  真实临床（age/sex/stage 零缺失）上**完成端到端运行；五种节点预算与规格表**逐一吻合**。
* 发现并修复 **1 处真实缺陷**（缺陷 5，临床编码静默降级为"缺失"），未触碰科学模型。
* 记录 **1 次自身纠错**（把全-BRCA 切片误标为多癌种，已重测）—— 未用改标签的方式掩盖。
* 全部数值比较使用**事先声明**容差，实测最大误差为 **7.304e-06**（容差 1e-4），
  其余关键项为 **0** 或 **2.220e-16**；未为"完全一致"放宽任何容差。
* gene-level 额外前向成本已量化（**1.06×–1.25× 时间、14.7 MB 内存**），
  决定**不优化**以保持计算定义。
* 修复后：发布审计 **FAIL 0 / WARN 0**；全量测试 **150 / 150 PASS**；
  渲染层 **12/12**；冻结资产 **13/13 异常 0**。
* **未** merge、**未** tag、**未** push、**未**发布。
