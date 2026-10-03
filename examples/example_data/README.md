# 示例数据

来自 TCGA 的**公开**样本（5 例，barcode 见 `example_labels.tsv`），由
`tools/build_fixtures.py` 从原项目只读表达矩阵导出。

| 文件 | 内容 |
| --- | --- |
| `example_expression.tsv.gz` | 5 例 × 15,672 基因，**TPM**（完整输入） |
| `example_expression_partial.tsv.gz` | 同 5 例 × 3,918 基因（演示缺失基因策略与覆盖度 QC） |
| `example_clinical.tsv` | `age` / `sex` / `stage`（含缺失，用于演示冻结参考值填充） |
| `example_labels.tsv` | `cancer_type` / `os_time_days` / `os_event`（仅示例与 KM 演示） |

示例数据仅用于演示接口，**不得**用于任何性能结论。
