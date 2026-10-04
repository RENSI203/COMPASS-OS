# 01 · 数据限制与适用范围审查 ｜ data limitations and applicability

**阶段**：1 / 6（顺序发布审查）
**日期**：2026-10-04
**候选基线**：Phase 0 冻结 `624184c` → **本阶段修复提交 `fc63905`**（见 §0.2）
**结论**：**PASS（含 4 处已修复缺陷与 6 项明确未验证限制）**；无 blocker

> 本阶段**未**改动冻结科学模型、**未**重新 fit、**未**修改任何结果数值。
> 所有修复都只让**实现与既有声明一致**，未扩张任何科学声明。

---

## 0. 方法与范围

### 0.1 方法

对 10 个指定事项逐项建立 **支持输入 → 处理规则 → 验证证据 → 限制** 对照，
每条判断给出**代码 / 资产 / 文档**依据。判定采用 **PASS / FAIL / NOT VERIFIED**，
缺证据一律记 NOT VERIFIED，不得以"看起来合理"记 PASS。

实际执行：阅读全部核心源码（`preprocessing.py` / `qc.py` / `api.py` / `survival.py` /
`sample_path.py`）、冻结资产（`qc_config.json` / `model_manifest.json` /
`locked_M0–M3.json` / `cancer_codes.tsv` / `gene_vocabulary.txt`）、
`README.md` / `docs/*.md` / `validation/README.md` 与 `validation/results/*`；
并对输入边界**实际运行探针**（不是读代码推测）。

### 0.2 候选变更与需重跑的前序检查

本阶段发现 3 处**真实缺陷**（§3.1 / §4.1 / §4.2），已在授权范围内直接修复并验证
⇒ **Phase 0 的冻结按 `00_CANDIDATE_BASELINE.md` §3 的规则失效**，已重新冻结：

| 项 | Phase 0 | 本阶段 |
|---|---|---|
| 代码冻结点 | `624184c` | **`28cdabf`**（含 4 处修复与审计文档） |
| 受保护路径差异 | — | 见下表 |
| 冻结资产 | 13/13 异常 0 | **13/13 异常 0（未变）** |

**需重跑的 Phase 0 检查（已全部重跑）**：

| Phase 0 检查 | 为何需重跑 | 重跑结果 |
|---|---|---|
| 发布审计 | `src/` 改动 | **FAIL 0 / WARN 0**（16 项） |
| 全量测试 | 代码+测试改动 | **PASS 148 / FAIL 0 / SKIP 0**（132 → 148） |
| 渲染层自检 | 无关联但按规则重跑 | **12/12 OK** |
| 资产哈希 | 无关联但按规则重跑 | **13 项异常 0** |
| wheel 构建与 SHA256 | `src/` 改动 | 见 §0.3 |
| clean worktree | 提交变更 | 已从 `fc63905` 重建 |

### 0.3 重跑后的产物哈希

在 clean worktree（检出本阶段末候选 `28cdabf`）用
`SOURCE_DATE_EPOCH=1791109274`（= 该提交时间）规范化构建：

| 产物 | SHA256 |
|---|---|
| `compass_os-1.1.0-py3-none-any.whl` | `631939233a8366ee5847c6cc9501d643324d276b940736511cad9142b7c14031` |
| `compass_os-1.1.0.tar.gz` | `126a2556f16acf73722af35eb71c3ee7b4beb9054dbca8de92b603ad5a2bdd36` |

> Phase 0 记录的 wheel `8897c47a…` 对应旧代码树，**已被本阶段的 `src/` 修复取代**，
> 不可再用于追溯；第 2 阶段须对 `6319392…` 在全新环境重做安装验证。

> 规范化构建仍要求 `SOURCE_DATE_EPOCH=$(git log -1 --format=%ct HEAD)`（Phase 0 §5.4）。

### 0.4 本阶段结束时的权威复验（clean worktree @ `28cdabf`）

| 检查 | 结果 |
|---|---|
| 发布审计 | **FAIL 0 / WARN 0**（16 项） |
| 全量测试（先 build 后 test） | **PASS 148 / FAIL 0 / SKIP 0** |
| 渲染层自检 | **12 / 12 OK** |
| 冻结资产 | **13 项异常 0**（未变） |
| worktree 未提交改动 | **0** |

---

## 1. 表达矩阵：方向、标识、重复、类型、缺失、NaN/Inf、负值、空输入、样本对齐

| 项 | 支持输入 | 处理规则 | 依据 | 判定 |
|---|---|---|---|---|
| 方向 | `samples × genes` | 行=样本、列=基因；内部按列名取数 | `preprocessing.py:6-8` docstring、`_to_tpm` | **PASS** |
| 类型 | 必须 `pandas.DataFrame` | 非 DataFrame → `InputError` | `_to_tpm`：`isinstance(expression, pd.DataFrame)`；探针实测 ndarray → `InputError` | **PASS** |
| 基因标识 | 列名 = **基因符号**（HGNC symbol） | `reindex(columns=feature_names())`；未匹配列被丢弃 | `align_expression`；词表 `data/gene_vocabulary.txt`（15,672） | **PASS** |
| 样本标识 | 索引 = sample ID | 索引仅用于对齐与输出，不参与模型 | `api.py:_run` | **PASS** |
| 重复基因列 | **不支持** | `columns.duplicated()` → `InputError`（列出前 5 个） | `_to_tpm`；探针实测 `['A']` → `InputError` | **PASS** |
| 重复样本 ID | **不支持** | `index.duplicated()` → `InputError` | `_to_tpm`；探针实测 → `InputError` | **PASS** |
| 空输入 | **不支持** | `shape[0]==0 or shape[1]==0` → `InputError` | `_to_tpm`；探针实测 `(0,0)` 与 `(4,0)` 均 → `InputError` | **PASS** |
| 数值类型 | int / float / 数值字符串均可 | `astype(np.float64)` 统一转换 | `_to_tpm`；探针实测整数 dtype 与 object 数值串均接受 | **PASS** |
| 负值 | **不支持** | `finite.min() < -1e-9` → `InputError`，并提示"若为 log2 请设 `input_scale='log2_tpm1'`" | `_to_tpm`；探针实测 −0.1 → `InputError` | **PASS** |
| **±Inf** | **不支持**（本阶段修复） | 见 §3.1 | `_to_tpm`（新增 `np.isinf` 检查） | **PASS（修复后）** |
| NaN | **支持**，语义 = 缺失 | 视同缺失基因，按 `missing_gene_strategy` 处理；逐样本 NaN 另给 `per_sample` 统计 | `align_expression` docstring；`qc.gene_coverage`；`test_nan_is_still_treated_as_missing` | **PASS** |
| 样本对齐 | 无需手动对齐 | 癌种/临床以**同一索引**对齐；`_clinical_frame` 用 `reindex(idx)` | `survival.py:_clinical_frame`；`api.py:_run` | **PASS** |
| 部分基因缺失 | **支持**（核心设计） | 见 §3 | — | **PASS** |

**限制**：输入矩阵中**不在词表内的额外列会被静默丢弃**（`reindex` 语义），
不报错、不提示。这是有意的（允许用户直接传全转录组矩阵），但用户若列名写成
Ensembl ID 或别名，会得到 0 % 覆盖度的结果而不是"列名错误"提示。

---

## 2. 输入尺度的真实支持情况

**代码事实**（`preprocessing.py:52`）：

```python
INPUT_SCALES = ("tpm", "log2_tpm1")
```

| 尺度 | 声明方式 | 处理 | 判定 |
|---|---|---|---|
| 线性 TPM | `input_scale="tpm"`（默认） | 原样进入；冻结 `Datascaler` 内部做 `log2(x+1)` + MinMax | **PASS（已验证）** |
| `log2(TPM+1)` | `input_scale="log2_tpm1"` | `2**x − 1` 转回 TPM | **PASS（已验证）** |
| 其他任何字符串 | — | `InputError` | **PASS** |
| raw counts | **无声明入口** | ⚠️ **非负 counts 会被静默当作 TPM** | **FAIL（无法检测，见下）** |
| 微阵列原始强度 | **无声明入口** | ⚠️ 同上 | **FAIL（无法检测）** |
| cohort z-score | — | 含负值 ⇒ 通常触发负值 `InputError`；若恰为非负则静默通过 | **部分防护** |
| 任意已标准化矩阵 | — | 同上 | **部分防护** |

**实测证据**：非负整数 counts（0–5000）与 0–15 强度样矩阵**均被接受**，
按 TPM 处理（`test_counts_like_matrix_is_not_detectable` 锁定该事实）。

**这是限制而非缺陷**：库无法只凭数值判定尺度。数值上"非负矩阵"与 TPM 不可区分。
处理方式是**强制调用方显式声明** `input_scale`，并在文档中写明。

> **必须避免的表述**：**不得**声称"支持 counts / 微阵列强度"。
> 文档现状（`README.md:309-310, 323`）已正确写明
> *"Direct public inference currently accepts TPM or log2(TPM+1)-scale expression.
> Raw microarray intensity matrices are not accepted directly — they first need
> cross-platform harmonisation"*，与代码一致，**无需改动**。
> 唯一缺口是 **counts 未被点名**，本阶段已在 `docs/audit/` 与
> `docs/ADVANCED_USAGE.md` 补上（§11）。

---

## 3. 基因词表、对齐与缺失基因策略

### 3.1 词表与对齐 — **PASS**

* 词表 15,672 基因，来源 `data/gene_vocabulary.txt`（由 `tools/build_assets.py` 从
  checkpoint 导出，`tests/test_assets.py` 复核与 `feature_name` 一致）。
* 对齐 = `reindex(columns=genes)`：**按列名**匹配，与列顺序无关。

### 3.2 三种策略的实际行为 — **PASS**

| 策略 | 实际行为 | 依据 |
|---|---|---|
| `reference`（默认） | 缺失基因填**冻结训练期参考中位数**：`log2(TPM+1)=3.130937933922` ⟺ `TPM≈7.760043`。**绝不**用用户队列自身均值/中位数 | `reference_fill_log2()` = `median(reference_quantiles['values'])`；`model_manifest.json` 的 `missing_gene.reference_value_*`；实测 `reference_fill_tpm()=7.760043` |
| `zero` | **标准化空间掩码**：正常过冻结 scaler 后把缺失基因的**标准化值**置 0（与训练 `RandomMaskAugmentor` 同一归一化空间） | `representation._ZeroMaskedScaler`；`zero_equivalent_tpm()` 提供等价式供 golden test |
| `strict` | 任一必需基因缺失即 `MissingGenesError` | `align_expression`；实测缺 15,322 基因 → `MissingGenesError` |

**语义边界（文档已正确约束）**：`zero` **不得**表述为"与训练增强分布一致"——
训练期掩码只作用于正样本视图、概率 0.41 且随机；推理期是对真实缺失基因的确定性操作
⇒ 只能说"同一归一化输入空间"（`preprocessing.py:25-27`、`docs/MISSING_GENES.md`）。
**代码注释与文档一致** ⇒ PASS。

**`zero` 不得实现成 raw TPM=0**：60.2 % 的基因训练 `data_min_ > 0`，置 0 会落到训练分布外。
代码通过标准化空间掩码避免了这一点 ⇒ PASS。

### 3.3 coverage 必须在填补前计算 — **PASS（重点核对）**

| 检查 | 结果 |
|---|---|
| `gene_coverage` 调用点 | `align_expression` 内、`reindex`/`fillna` **之前**（`preprocessing.py:135` vs `139-148`） |
| 实测（缺 15,322 基因） | 填补前 `gene_coverage = 0.022333` |
| 填补后矩阵非 NaN 比例 | `1.000000` ← **若误用会虚报 100 %** |
| QC 实际使用的掩码 | `observed_mask`（填补前），实测 `0.022333` ✓ |
| `gsig_coverage`（与阈值同口径） | 取 `aligned.observed_mask` 的 916 个 signature 基因均值，**同样在填补前**（`api.py:144-148`） |

**结论：覆盖率不会被填补污染成 100 %。** 这条是本阶段最关键的核对项，**通过**。

### 3.4 未通过阈值时 — **PASS**

`qc_config.json` 不存在时 `coverage_warnings()` 返回空列表并追加说明
*"Coverage thresholds are not frozen in this installation … coverage is reported as
continuous values only, without grading"*（`api.py:155-159`）⇒ 不会静默假通过。

---

## 4. `cancer_type` 的必填、合法值、单/多癌种与未知癌种

| 项 | 行为 | 依据 | 判定 |
|---|---|---|---|
| 必填 | 无默认值；`None`/缺失 → `UnknownCancerTypeError`（"癌种必须逐个提供，不得推断"） | `cancer_codes_for`；探针实测含 `None` → `UnknownCancerTypeError` | **PASS** |
| 合法值 | `data/cancer_codes.tsv` 的 TCGA 缩写（34 行） | `cancer_codes()` | **PASS** |
| 未知癌种 | `UnknownCancerTypeError` 并指向码表 | `cancer_codes_for`；探针实测 `NOTACANCER` → `UnknownCancerTypeError` | **PASS** |
| **长度不符** | **本阶段修复**：原为 pandas 裸 `ValueError`，现为 `InputError` 并给出两个数值 | 见 §4.1 | **PASS（修复后）** |
| **`NORMAL`** | **本阶段修复**：原为 embedding 裸 `IndexError`，现为 `UnknownCancerTypeError`（"正常组织标记"） | 见 §4.2 | **PASS（修复后）** |
| drop-first 参照水平 | `ACC`（`compass_code=0`，`ct_column` 为空）⇒ 该样本 32 列全 0 | `ct_onehot_columns` docstring；码表 `is_reference_level=True` | **PASS** |
| 单/多癌种 | 不改变推理；仅影响 `sample_path` 展示预算（单癌种隐藏 Cancer type 节点） | `sample_path_render.select_nodes` | **PASS** |

### 4.1 缺陷 1：`cancer_type` 长度不符 → 裸 `ValueError` — **已修复**

* **现象**：`cancer_type` 长度 ≠ 样本数时抛 `ValueError: Length of values (2) does not
  match length of index (4)`（pandas），而非包装的 `InputError`。
* **根因**：`api.py:_run` 先 `pd.Series(list(cancer_type), index=aligned.matrix.index)`
  再检查长度 ⇒ pandas 先抛错，长度检查永远不执行。
* **影响**：错误信息不指向 `cancer_type`，用户难以定位；且与包内错误类型约定不一致。
* **修复**：先构造**无 index** 的 Series 并校验长度，再赋 index（保留 NaN 语义）。
* **验证**：`test_cancer_type_length_mismatch_is_input_error[1,2,5]`；实测
  → `InputError: cancer_type 长度 2 != 表达矩阵样本数 4；必须逐样本提供癌种`。

### 4.1b 缺陷 4：尺度自检只在 `analyze()` 生效，`predict()` 静默接受错误尺度 — **已修复**

* **现象**：把 log2(TPM+1) 的值声明为 `input_scale="tpm"` 送入 `predict()`，
  **不产生任何尺度提示**，直接返回看似正常的 risk；反向（把 TPM 说成 log2）同样静默。
  实测：真 log2 说成 tpm → `risk=[1.3525, 1.1873, 1.2137, 1.4168]`（无提示）。
* **根因**：尺度启发式只写在 `analysis.py:480`（`analyze()` 一键入口）内部，
  **不在共享预处理路径上**，因此 `predict()`（Level 2 直接公共推理）没有该检查。
* **影响**：README 的中文说明写有"若数值与声明尺度明显不符会给出提示"，
  但该声明对 `predict()` **不成立** —— 属**实现与声明不一致**。
  更实际的风险是：错误的尺度声明会静默产生一个看似合理的风险值。
* **修复**：把同一启发式（`max < 30` ⇒ 疑似 log2，与 `analysis.py` 既有口径一致）
  移入 `align_expression()` 共享路径，**双向**提示，且**只提示不阻断**
  （避免把可疑但合法的输入拒之门外）。
* **验证**：`test_scale_mismatch_note_fires_on_both_directions`、
  `test_scale_mismatch_note_is_a_warning_not_an_error`；实测
  真 log2 说成 tpm → 出现提示；两种**正确**声明 → **无**提示。

### 4.2 缺陷 2：`cancer_type='NORMAL'` → 裸 `IndexError` — **已修复**

* **现象**：`compass_code = -1` 直接送入 embedding ⇒ `IndexError: index out of range in self`。
* **根因**：`data/cancer_codes.tsv` 第 2 行 `NORMAL  -1  (空)  False` 被当作普通合法值；
  `cancer_codes_for` 只校验"是否在码表中"，不校验码值符号。
* **影响**：`NORMAL` 是**正常组织**标记而非肿瘤类型；用户按码表填写却得到与输入无关的报错。
* **附加影响（更隐蔽）**：`ct_onehot_columns` 会把 `NORMAL` 映射为**32 列全 0**，
  即**静默等同于 ACC 参照水平** —— 若某些路径不经过 embedding，会得到看似正常的错误结果。
* **修复**：`cancer_codes_for` 与 `ct_onehot_columns` 均拒绝负码值并给出明确说明。
* **验证**：`test_normal_is_rejected_with_clear_message` / `test_normal_rejected_in_predict`
  / `test_all_table_tumour_types_still_resolve`（确认 33 个非负码值肿瘤类型全部仍可用）。

---

## 5. Age / Sex / Stage：编码、填补、异常值与索引对齐

| 项 | 行为 | 依据 | 判定 |
|---|---|---|---|
| 列名 | `age` / `sex` / `stage`（小写） | `survival.CLINICAL_FIELDS` | **PASS** |
| 索引对齐 | `clin.index = str(...)` → `reindex(idx)`，与表达矩阵同一 `idx` | `_clinical_frame` | **PASS** |
| 非数值处理 | `pd.to_numeric(errors="coerce")` ⇒ 转不动的变 NaN ⇒ 走冻结填补 | `_clinical_frame:120` | **PASS** |
| Age 异常值 | **无范围校验**；任何有限值都接受 | 缺检查 | **限制（见下）** |
| Sex 编码 | 数值直通（冻结模型按 0/1 使用）；**不参与标准化** | `locked_M*.json` 的 `scaled_features` 不含 Sex | **PASS** |
| Stage 异常值 | **仅接受 1–4**；0 或缺缺失 → NaN → 冻结填补 | `_clinical_frame:122`：`between(1, 4)` | **PASS** |
| 缺失填补 | 冻结训练参考值 `Age=60.0`、`Sex=0.0`、`Stage=2.0`（M0–M3 相同） | `locked_M*.json` 的 `clinical_fill`；实测一致 | **PASS** |
| 缺整张临床表 | `clinical=None` ⇒ 三列全 NaN ⇒ 全部走冻结填补（不报错） | `_clinical_frame:110-113` | **PASS**（但见下） |

**限制（重要）**：冻结填补值**不是真实患者观测值**。
`sample_path` 已把逐字段 `imputed` 标注为 `(imputed)` 并进入图注与 `nodes.tsv`；
`docs/SAMPLE_PATH.md` §3 亦写明"实际被冻结填补的临床节点标注 `(imputed)`"。
**但 `predict()` 的返回值本身不区分"用户提供"与"冻结填补"**——
只有 `sample_path` 这条路径会标注。若论文声明涉及临床变量贡献，必须显式说明填补比例。

**限制**：**Age 无生理范围校验**（例如 999 或 −5 会被接受并标准化）。
Stage 有 1–4 约束，Age/Sex 没有等价约束。未验证这是否影响结论 ⇒ 记 **NOT VERIFIED**。

---

## 6. 默认模型与 M0–M3 的输入要求；单样本风险与队列分层边界

| 项 | 事实 | 依据 | 判定 |
|---|---|---|---|
| **默认模型** | **`M2`**（不是 M3） | `model_manifest.json` 的 `default_model: "M2"`；`survival.default_model()`；`api.predict(model="default")` | **PASS** |
| M0 | 临床基线：`CancerType(one-hot)+Age+Sex+Stage`，**无 COMPASS 分支** | `model_manifest.json` | **PASS** |
| M1 | 仅 43 concepts（无临床） | 同上 | **PASS** |
| M2 | 主模型：`CancerType+Age+Sex+Stage+43 concepts` | 同上 | **PASS** |
| M3 | M2 + `PC1–PC10`（锁定 PCA，作用于 `log2(TPM+1)`） | 同上 + `locked_pca_M3.npz` | **PASS** |
| 模型选择 | `model="default"` 或显式 `M0..M3`；未知 → `InputError` | `survival.lock_path` | **PASS** |
| 单样本 | `predict()` 接受 n=1；`sample_path` 支持单样本 cohort | 探针实测 `align_expression` 单样本 cov 正常；`tests/test_sample_path.py::test_single_sample_cohort` | **PASS** |
| 队列分层边界 | `risk_group` 用**冻结** `median_cutoff`（训练集风险中位数）；`risk_group_relative` 用本队列中位数（**display only**） | `survival.median_cutoff`；`README.md:281-283` | **PASS** |

**冻结中位切点**（实测，全 4 模型）：

| 模型 | `median_cutoff` |
|---|---|
| M0 | 1.085061 |
| M1 | 0.039928 |
| M2 | 1.078303 |
| M3 | 0.696889 |

**限制**：`median_cutoff` 来自**训练集（TCGA）**风险分布。在风险分布整体平移的
外部队列上，以它做 high/low 分组可能把绝大多数样本分到同一侧——代码已提供
`risk_group_relative` 作为**显示用**替代，并规定其不得用于正式分层结论
（`README.md:283`）。**该边界在文档中已写明** ⇒ PASS。

---

## 7. `sample_path`：表达单位、gene score 来源、132/43 层、阈值与颜色；展示筛选不得改变预测

| 项 | 事实 | 依据 | 判定 |
|---|---|---|---|
| 表达单位 | `expression_units` 随 `input_scale` 给出：`"TPM"` 或 `"TPM (converted from log2(TPM+1) input)"` | `sample_path.py` payload | **PASS** |
| gene score 来源 | 官方 `PreTrainer.extract(..., with_gene_level=True)` 的 `dfg` = `genesetprojector.geneset_scorer(encoder_output)[:, 2:]`（共享 `nn.Linear(32→1)`，形状 `(n, 15672)`） | `sample_path._extract_gene_level`；上游 `tune.py:289-291` | **PASS** |
| 132 / 43 层 | 取主 API `predict()` **同一次前向**的 `signature_scores` / `concept_scores`（与 η 同源） | `sample_path.py`；`docs/SAMPLE_PATH_INTEGRATION_VALIDATION.md` §3 | **PASS** |
| 禁止替代 | 未使用 attention weight / TPM / 风险分配量 / 旧 HTML contribution 作为 gene score；注意力权重**仅**用于筛选连接 | `_topology_edges` / `_signature_concept_edges` docstring | **PASS** |
| 展示阈值 | `gene_threshold=0.25`（真实 gene-score 尺度）、`max_genes≤32`、`max_signatures≤26`、`predictor_budget≤16`，各层 ≤50 | `sample_path_render.Style` | **PASS** |
| 颜色尺度 | 共享固定范围（非逐图 min-max）：`expression (0,1704.3)`、`gene_score (−2,2)`、`signature/concept (−1,1)`、`contribution (−1.5,1.5)`；超出仅饱和 | `examples/sample_path_example.py` 的 `SHARED_STYLE` | **PASS** |
| 阈值性质 | 明确为**展示参数**，非经验证科学阈值 | `docs/SAMPLE_PATH.md` §4；图注 | **PASS** |
| **展示筛选不改预测** | 筛选只在渲染层 `select_nodes` 内进行，**不进入推理**；`risk`/`rank`/`percentile` 取主 API | `sample_path_render.select_nodes`（纯选择）；`tests/test_sample_path_integration.py::test_score_extraction_does_not_change_existing_outputs` | **PASS** |
| 完整分解 | 未截断分解含未展示 concept、被隐藏 Cancer type 与 offset，`sum == η` | renderer `_check` 强制；实测 Δ ≤ 2.2e-16 | **PASS** |

---

## 8. cohort-relative rank/percentile、median cutoff、冻结 cutoff 与 S(t)

| 量 | 定义 | 依据 | 判定 |
|---|---|---|---|
| `risk` | Cox 线性预测子 η = Xβ（**相对风险尺度，不是概率**） | `survival.risk_from_design`；`README.md:276` | **PASS** |
| `rank` | `1 + #{η > η_i}`（队列内升序） | `sample_path.py` | **PASS** |
| `percentile` | `100 × (P(η<η_i) + 0.5·P(η=η_i))`，高值 = 高风险 | `sample_path.py`；写入 `rank_definition` | **PASS** |
| cohort median cutoff | `np.median(队列 η)`；`cutoff_mode` 明确标 **"visualization boundary only — NOT the frozen validated prognostic cutoff"**；`cutoff_is_frozen_validated` 恒 `False` | `sample_path.py` | **PASS** |
| 冻结 cutoff | `median_cutoff(model)` = **训练集**风险中位数（§6 表） | `survival.median_cutoff` | **PASS** |
| `S(t)` | `exp(−H0(t)·exp(risk))`，`H0` 取 lock 的 `base_time`/`base_cumhazard`（2864 点，1–11260 天） | `survival.survival_probability` | **PASS（定义）** |
| **绝对生存校准** | 模块 docstring **明确**：*"外部验证中 strict 口径的绝对风险校准明显弱于判别度 … 属于 baseline-hazard-based derived output，**不得**宣传为稳定的个体绝对风险预测"* | `survival.py:17-19`；`README.md` | **PASS（约束）** |

**关键核对：绝对生存概率的校准支持不得由风险排序性能代替。**
代码与文档**已经**正确区分：`S(t)` 被标为 derived output 并附校准限制；
`README.md` 明确 risk 是 log-relative-hazard 而非概率。**未发现把 C-index /
风险分层性能当作绝对生存校准的表述** ⇒ **PASS**。

**限制（NOT VERIFIED）**：仓库内**未见**独立的绝对风险校准证据
（如校准曲线、ICI、Brier 分解、时间相关校准）。目前关于"校准弱"的陈述来自
`survival.py` 与文档的定性说明，**不是本候选内的可复核校准实验**。
⇒ 任何论文中的绝对生存/绝对风险声明都应记 **NOT VERIFIED**，除非补充校准证据。
本阶段**不新增**校准实验（超出范围）。

---

## 9. QC 配置中的阈值与 null 项（以**当前资产**为准）

**权威来源**：`src/compass_os/assets/models/qc_config.json`（version 1.0.0，
`validation_date 2026-10-03`，12 数据集 / 720 样本）。

| 键 | 值 | 判定 |
|---|---|---|
| `recommended_global_coverage` | **0.90** | 已标定（见下） |
| `warning_global_coverage` | **0.70** | 已标定 |
| `recommended_Gsig_coverage` | **0.90** | 已标定 |
| `warning_Gsig_coverage` | **0.70** | 已标定 |
| **`signature_coverage_warning_threshold`** | **`null`** | **未标定、不予发布** |
| `threshold_kind` | `empirically_calibrated_descriptive_tiers` | 明确为**描述性分层** |
| 组合规则 | `overall = worse(global, Gsig)` | `qc.qc_tier` |

**阈值来源证据**（`validation/results/qc_threshold_audit.tsv`，
12 队列 × 覆盖率阶梯）：

| global coverage | risk_spearman | group_agreement | recommended(0.95/0.90) | warning(0.80/0.70) |
|---|---|---|---|---|
| 0.950 | 0.9863 | 0.9917 | ✅ | ✅ |
| **0.900** | **0.9593** | **0.9917** | ✅ ← 最低满足点 | ✅ |
| 0.800 | 0.9230 | 0.9833 | ❌ | ✅ |
| **0.700** | **0.8449** | **0.9917** | ❌ | ✅ ← 最低满足点 |

⇒ `0.90`/`0.70` 是**实测曲线上的最低满足点**，与 `qc_config.json` 一致
（`cutoffs_used.recommended = {risk_spearman: 0.95, group_agreement: 0.9}`；
`warning = {0.8, 0.7}`）⇒ 标定方法与资产自洽。

**必须遵守的表述纪律**：

1. `qc_config.json` 的 `note` 明确：阈值属**工程 QC 阈值，不是生物学边界**；
   曲线若无明显断点则它们是**描述性分层**而非"安全线"。
   ⇒ **不得**写成"经过验证的安全阈值"或"clinically validated cutoff"。
2. `signature_coverage_warning_threshold` 为 `null`，`notes_missing` 给出原因：
   可复核的覆盖度–误差表要求 (signature × coverage) 逐条件累积，而当前保留的结果表中
   `coverage_mean` 对 132 个 signature **恒为 1.0**（无覆盖度变异）⇒ 无法估计阈值。
   ⇒ **不得**写成已冻结/已验证；per-signature coverage 只以连续值报告。
3. `qc_config.json` 的 `version: "1.0.0"` 是**该 QC 标定资产的版本**，
   与软件版本（`1.1.0`）无关，**不得**随软件版本一起 bump。

**限制（NOT VERIFIED）**：`qc_threshold_audit.tsv` 的 `recommended_condition`
判定基于 12 个队列；`empirical_cohorts: 38` 用于经验覆盖分布。
阈值的**跨平台稳健性**（尤其 GPL570 之外的平台）未单独验证 ⇒ NOT VERIFIED。

---

## 10. 外部队列适用性、缺失模式与平台限制

### 10.1 已有证据

`validation/results/external_coverage_summary.tsv`：**39 个 GEO 外部队列**
（38 `ok` / 1 `unavailable` = METABRIC），平台分布：

| 平台 | 队列数 | 类型 |
|---|---|---|
| GPL570 | 19 | Affymetrix U133 Plus 2.0（**微阵列**） |
| GPL96 | 4 | Affymetrix U133A（微阵列） |
| 其余 12 个平台 | 各 1–2 | 均为微阵列/外显子阵列 |

实测覆盖度范围（38 个 ok 队列）：

| 指标 | 范围 | 中位 | < 0.90 | < 0.70 |
|---|---|---|---|---|
| `global_gene_coverage` | 0.6355 – 0.9803 | 0.9668 | **13 / 38** | 1 / 38 |
| `Gsig_coverage` | 0.8035 – 0.9782 | 0.9618 | **8 / 38** | 0 / 38 |

**涵义**：约 **1/3 的外部队列低于 recommended 阈值**（但仍高于 warning），
说明"缺失基因"是外部队列**常态而非例外**，且 QC 分级确实有区分度。

### 10.2 关键限制：TCGA 无法提供独立验证

`validation/README.md` 明确：

> **所有 9,430 例 TCGA 样本都已参与最终 Cox 系数拟合**
> （`_final_split.json` 的 `train = all_ids(9430)`，覆盖 random 划分的 train/val/test 全部）。
> ⇒ TCGA 内部**没有**"未参与拟合"的样本。因此主分析使用 **39 个外部队列**；
> TCGA val/test 只能作为 *in-domain perturbation stability reference*，
> **不得**表述为 independent model validation。

⇒ 任何以 TCGA 为验证集的性能声明都必须记 **NOT VERIFIED / 不成立**。
本阶段**确认文档已正确约束**这一点 ⇒ PASS（约束层面）。

### 10.3 平台与尺度限制

* 外部队列本就是**微阵列平台**，其"TPM"来自**跨平台和声化后的缓存**
  （`scripts_dev/analysis_out_harmonize/expr_log2/`，log2 尺度），
  **不是实测 TPM**。
  ⇒ 文档中 `input_scale="tpm"` 的语义必须是"模型输入尺度"，
  不能读作"实测 TPM"。`sample_path` 的 `expression_units` 已写
  `"TPM (converted from log2(TPM+1) input)"` ⇒ 与事实一致。
* **微阵列基因缺失 ≠ 测序表达缺失**：前者是"平台不含该探针"（结构性缺失），
  后者是"表达量低于检出"。39 队列的 `global_gene_coverage` 主要反映**平台探针覆盖**。
  ⇒ 缺失基因策略的验证证据**主要来自平台缺失**，**不能**直接外推到
  RNA-seq 的检出限缺失 ⇒ **NOT VERIFIED**。

### 10.4 三类适用性判定

| 类别 | 内容 |
|---|---|
| **已经验证** | ① TPM 与 log2(TPM+1) 两种声明的尺度；② `reference`/`zero`/`strict` 三种缺失策略的行为与一致性；③ 覆盖率在**填补前**计算；④ 12 队列 × 覆盖率阶梯的 QC 阈值标定；⑤ 39 队列的平台覆盖度分布；⑥ 冻结资产哈希；⑦ `sample_path` 的真实分数来源与 η 一致性 |
| **工程支持但未验证** | ① 非负 counts/强度矩阵（会被当作 TPM，无检测）；② RNA-seq 检出限型缺失；③ Age 异常值；④ 混合来源队列（跨平台拼接）；⑤ 单样本 cohort 的 percentile 语义；⑥ GPL570 以外平台的阈值稳健性 |
| **尚未验证** | ① 绝对生存 `S(t)` 的校准（无校准曲线/Brier/ICI 证据）；② TCGA 作为独立验证集（**本质上不成立**）；③ `signature_coverage_warning_threshold`（资产中为 null）；④ 非肿瘤/NORMAL 组织样本；⑤ 非 TCGA 癌种命名体系的映射 |

---

## 11. 本阶段对文档的更新

| 文件 | 改动 |
|---|---|
| `docs/audit/01_DATA_LIMITATIONS.md` | 本文件（新建） |
| `src/compass_os/preprocessing.py` | 缺陷 1（±Inf）与缺陷 4（尺度自检移入共享路径） |
| `src/compass_os/api.py` | 缺陷 2（`cancer_type` 长度先校验） |
| `tests/test_input_contract.py` | 新建，16 项输入契约回归 |
| `docs/ADVANCED_USAGE.md` | 补"输入尺度"一节：明确只支持 TPM / log2(TPM+1)；**点名 counts 与微阵列强度不被支持且无法自动检测**；说明 `NORMAL` 不是可用癌种 |
| `README.md` | §"Input requirements" 补 counts 与 `NORMAL` 两点（不改变既有科学声明） |

**未改动**：任何冻结资产、模型系数、QC 阈值、`docs/SAMPLE_PATH*` 的既有结论。

---

## 12. PASS / FAIL / NOT VERIFIED 汇总

| # | 事项 | 判定 |
|---|---|---|
| 1 | 方向/标识/重复/类型/缺失/NaN/Inf/负值/空输入/对齐 | **PASS**（Inf 缺陷已修复） |
| 2 | TPM 与 log2(TPM+1) 真实支持；其余尺度未宣称支持 | **PASS**（尺度自检缺陷已修复；counts/微阵列不可检测 = 已记录限制） |
| 3 | 词表与对齐；三种策略行为；**coverage 在填补前计算** | **PASS** |
| 4 | `cancer_type` 必填/合法值/未知/单多癌种 | **PASS**（长度与 NORMAL 两处缺陷已修复） |
| 5 | 临床编码/填补/异常值/索引对齐 | **PASS**（Age 无范围校验 = 限制，记 NOT VERIFIED） |
| 6 | 默认 M3 与 M0–M3 输入要求；单样本与队列分层边界 | **PASS**（**默认实为 M2**，已核实） |
| 7 | `sample_path` 单位/gene score/132/43/阈值/颜色；筛选不改预测 | **PASS** |
| 8 | rank/percentile/median cutoff/冻结 cutoff/S(t) 的区别 | **PASS**（定义与约束）；**绝对校准 NOT VERIFIED** |
| 9 | QC global/Gsig/signature 阈值与 null 项 | **PASS**（`signature_coverage_warning_threshold` 确为 null，未标定为已验证） |
| 10 | 外部队列适用性/缺失模式/平台限制 | **PASS**（分类完成）；RNA-seq 缺失与绝对校准 **NOT VERIFIED** |

### 本阶段缺陷与修复

| # | 缺陷 | 严重度 | 状态 |
|---|---|---|---|
| 1 | `±Inf` 静默通过输入校验，下游抛与输入无关的 sklearn 报错 | 中（可导致静默错误结果） | **已修复并验证** |
| 2 | `cancer_type` 长度不符 → 裸 `ValueError`，长度检查永不执行 | 低（错误信息质量） | **已修复并验证** |
| 3 | `cancer_type='NORMAL'` → 裸 `IndexError`；one-hot 路径**静默等同 ACC 参照水平** | 中（可得到看似正常的错误结果） | **已修复并验证** |
| 4 | 尺度自检只在 `analyze()` 生效，`predict()` 静默接受错误尺度 —— 与 README 声明不一致 | 中（可静默产生看似合理的风险值） | **已修复并验证** |

### Blocker

**无。** 不存在阻塞第 2 阶段的未解决问题。

---

## 13. 移交第 2 阶段（依赖项与安装审查）

| # | 事项 | 来源 |
|---|---|---|
| 1 | sdist 字节不可重复（Phase 0 缺陷 B） | Phase 0 §5.6 |
| 2 | 发布顺序"先 build 后 test"需固化 | Phase 0 §5.7 |
| 3 | 用 `fc63905` 重建的 wheel 需在**全新环境**复验 | 本阶段 §0.2 |
| 4 | 确认 `plot` extra 与实际依赖一致（matplotlib 已是核心依赖） | Phase 0 / 本阶段 §2 |
| 5 | **Age 无生理范围校验**是否需要在依赖/接口层加固 | 本阶段 §5 |
| 6 | `signature_coverage_warning_threshold = null` 的对外表述 | 本阶段 §9 |
| 7 | 绝对生存校准证据是否需要补充（超出本次审查的既定范围，需用户决定） | 本阶段 §8 |

---

## 14. 阶段结论

* 10 项指定事项**全部完成审查**，逐项给出代码/资产/文档依据。
* 发现并修复 **4 处真实缺陷**（均属"实现与声明不一致"，未触碰科学模型）。
* 明确记录 **6 项未验证限制**：非负 counts 不可检测、Age 无范围校验、
  RNA-seq 检出限缺失未验证、绝对生存校准无证据、TCGA 不能作独立验证、
  单 signature 阈值未标定。
* 修复后：发布审计 **FAIL 0 / WARN 0**；全量测试 **148 / 148 PASS**；
  渲染层 **12/12 OK**；冻结资产 **13 项异常 0**。
* **未** merge、**未** tag、**未** push、**未**发布。
