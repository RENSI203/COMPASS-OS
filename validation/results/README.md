# validation/results —— 随仓库发布的验证产物

## 被 Git 跟踪（publication-ready，精简）

| 文件 | 内容 |
| --- | --- |
| `external_coverage_summary.tsv` | 39 个外部队列的 pre-imputation 覆盖度（38 个可还原） |
| `stress_pancohort_summary.tsv` | 三臂 masking 的 pan-cohort 汇总（cohort 为单位的中位/IQR/bootstrap CI） |
| `stress_cohort_summary.tsv` | 逐队列 × 条件的汇总 |
| `qc_threshold_audit.tsv` | **逐 coverage level** 的判据满足情况（阈值来源可逐行核对） |
| `coverage_qc_validation.tsv` | 逐 feature 的覆盖度 vs 扰动误差 |
| `qc_config.json` | 冻结的 QC 阈值（与 `models/qc_config.json` 相同） |
| `design.json` | 本轮实验设计（队列/样本/重复/切点） |
| `release_audit.tsv` | 发布卫生审计结果 |
| `figures/` | Figure A–E（PDF + PNG + **source TSV**） |

## 仅保留在本地（**不进入 Git**，可由 `validation/` 下脚本重建）

| 文件 | 体积 | 重建方式 |
| --- | --- | --- |
| `raw_shard*.tsv` | ≈2.9 MB | `run_missing_gene_validation.py`（逐条件逐队列 append，可续跑） |
| `platform_masks.json` | ≈5.4 MB | `extract_platform_masks.py --source-root <原项目根>` |
| `_tmp_*.npy` | — | 运行中间件 |

`validation/inputs/`（pheno 与和声化表达缓存）同样不入 Git；和声化缓存来自原项目的
`analysis_out_harmonize/expr_log2/`，只在本地以只读方式引用。
