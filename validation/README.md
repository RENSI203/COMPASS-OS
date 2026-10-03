# 输入准备与运行（missing-gene stress test）

## 为什么主分析用外部队列

**所有 9,430 例 TCGA 样本都已参与最终 Cox 系数拟合**（`_final_split.json` 的
`train = all_ids(9430)`，覆盖 random 划分的 train/val/test 全部）。
⇒ TCGA 内部**没有**"未参与拟合"的样本。因此主分析使用 **39 个外部队列**；
TCGA val/test 只能作为 *in-domain perturbation stability reference*，
**不得**表述为 independent model validation。

## 输入

1. `--pheno`：TSV，列 `sample_id  cohort  cancer_type  os_time_days  os_event`
2. `--expr-dir`：目录，每队列一个 `<cohort>.tsv.gz`（样本 × 15,672，**log2(TPM+1)**），
   或单个 `.tsv.gz` 文件

原项目中对应的只读来源（供准备输入时参考，**不写进代码默认值**）：

| 用途 | 原项目路径 |
| --- | --- |
| 和声化表达（log2 尺度，39 队列） | `scripts_dev/analysis_out_harmonize/expr_log2/<GSE>.tsv.gz` |
| 队列癌种 / 平台 | `scripts_dev/analysis_out_geo_validate/batch_metrics.tsv`（列 `TCGA癌种`、`GEO数据集`） |
| 逐样本 OS | `data_geo/GEO_cohorts/<GSE>/`（series matrix pData，见 `scripts_dev/32_validate_geo_batch.py` 的解析口径） |

## 运行

```bash
python validation/missing_gene_stress_test.py \
    --pheno validation/inputs/pheno.tsv \
    --expr-dir validation/inputs/harmonized \
    --out validation/output --repeats 20
```

mask 层级 `0/5/10/20/30/40/50 %`；两类 mask 基因集合（全部 15,672 / 132 基因集涉及的 916）；
策略 `reference` 与 `zero`；每个重复内所有样本共用同一 mask 集合以保证可比。

## 产物

`stress_raw.tsv`（逐次重复）、`stress_curves.tsv`（coverage × 指标中位）、
`stress_report.md`、`fig_stress_curves.png`、`fig_stress_uno_c.png`。

## 阈值规则

**不预设** 80 %/90 % 等阈值。阈值由本实验数据确定（例如"risk 逐样本 Spearman ≥ 0.95 且
分层一致率 ≥ 0.90 的最低 coverage"作为 `recommended_coverage`），
写入 `qc_config.json`（带版本号）后再冻结。

⚠ smoke 运行产物见 `validation/output_smoke/`（**仅 5 例示例数据的接口验证，不是实验结论**）。
