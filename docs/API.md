# API 参考

三个公共入口 + 三个结果容器。所有输出都附带 QC；缺失基因**绝不静默填补**。

## `predict(expression, cancer_type, clinical=None, model="default", *, missing_gene_strategy="reference", input_scale="tpm", batch_size=16, device="cpu", min_cohort_for_stratification=30) -> PredictionResult`

| 参数 | 说明 |
| --- | --- |
| `expression` | `DataFrame`，`samples × genes`，列名 = 基因符号，索引 = 样本 ID |
| `cancer_type` | 长度 = 样本数；TCGA 缩写（见 `src/compass_os/data/cancer_codes.tsv`）；**必需**，缺失/未知报错 |
| `clinical` | 可选 `DataFrame`，索引 = 样本 ID，列 `age` / `sex` / `stage`（stage ∈ 1–4）；缺失用**冻结参考值**填充并在 warnings 中报告 |
| `model` | `"default"`（按 manifest，当前 `M2`）或 `M0`/`M1`/`M2`/`M3`，或逗号分隔多个 |
| `missing_gene_strategy` | `reference`（默认）/ `zero` / `strict`，见 `docs/MISSING_GENES.md` |
| `input_scale` | `tpm`（默认）或 `log2_tpm1` |

返回 `PredictionResult`：

| 字段 | 形状 | 说明 |
| --- | --- | --- |
| `sample_ids` / `cancer_type` | list | 与输入同序 |
| `risk` / `linear_predictor` | n × 模型 | 相对风险（线性预测子 βx） |
| `risk_rank` | n × 模型 | 队列内平均秩（n ≥ 2） |
| `risk_group` | n × 模型 | `high`/`low`，切点 = **冻结 `median_cutoff`**（非队列中位数）；单样本不伪造分层 |
| `signature_scores` | n × 132 | 模型基因集层输出 |
| `concept_scores` | n × 43 | 模型 cell-pathway 层输出 |
| `qc` | `CoverageQC` | 见下 |
| `warnings` | list | 临床填充、覆盖提示等 |

方法：`default_risk`（Series）、`survival_probability(times)`（派生输出，见 README 限制）、
`to_dict()`。

## `get_representation(expression, cancer_type, missing_gene_strategy="reference", *, input_scale="tpm", batch_size=16, device="cpu", signature_detail=False) -> RepresentationResult`

只算表示、不做生存预测。返回 `signature_scores`（n × 132）、`concept_scores`（n × 43）、`qc`。
`signature_detail=True` 时在 `qc` 上附 `signature_missing_detail`（每个 signature 缺失的基因名）。

## `check_robustness(expression, cancer_type, clinical=None, strategies=("reference","zero"), *, model="default", input_scale="tpm", batch_size=16, device="cpu") -> RobustnessResult`

对同一输入跑多种缺失策略并报告差异。返回：`risk`（逐策略）、`risk_difference`、
`concept_correlation` / `signature_correlation`（逐样本 Pearson）、`concept_input_coverage`、
`gene_coverage`、`signature_gene_coverage_reference`、`robustness_flag`。

⚠ `robustness_flag` 恒为 `"continuous_only"`：**不存在单一的二元 PASS/FAIL 判定**。
覆盖度 QC 分级已按冻结阈值校准（recommended 0.90 / warning 0.70，global 与 signature 双轴，
`overall = worse(global, Gsig)`，见 `models/qc_config.json`）；而
reference-versus-zero 稳健性比较只提供连续指标，请与已校准的覆盖度分级一起解读。

## `CoverageQC`

### 覆盖度警告（P2 起）

若 ``models/qc_config.json`` 存在（由 missing-gene stress test 生成并冻结），
``predict``/``get_representation`` 会**自动**按其中的阈值追加 `warnings`：

* `global gene coverage` / `signature-gene coverage` 低于 `recommended_*` → 提示分级
  （`warning` / `below_warning`）；**overall = worse(global, Gsig)**；
* 逐 signature 覆盖度以**连续值**报告（`signature_gene_coverage`）；v1 **不发布**逐
  signature 数值阈值。

若配置不存在，包**只报告覆盖度数值、不做分级**，并明确提示阈值尚未冻结——
绝不使用先验拍定的 80 %/90 %。

### 字段

`n_required_genes`(15672) / `n_observed_genes` / `n_missing_genes` / `gene_coverage` /
`missing_gene_strategy` / `missing_gene_names` / `n_samples` / `warnings`，
以及（存在逐样本 NaN 时）`per_sample`（`n_observed_genes_sample`、`n_missing_genes_sample`、
`gene_coverage_sample`）、`signature_gene_coverage`（n × 132）、
`concept_input_coverage`（n × 43）。`to_dict()` 直接可 JSON 化。

## 辅助

`available_models()`、`default_model()`、`model_info(model="default")`。
异常：`CompassOSError`（基类）、`MissingGenesError`、`UnknownCancerTypeError`、
`AssetNotFoundError`、`InputError`。

## 资产解析

`COMPASS_OS_ROOT` 环境变量可指向包含 `models/` 的仓库根；未设置时按包位置自动推断。
正式包**不依赖任何原项目绝对路径**（`tests/test_assets.py::test_no_absolute_paths_in_package` 强制）。
