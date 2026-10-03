# Golden 测试夹具

期望值由**原项目生产实现**生成，不是本包自算（见 `golden_provenance.json`）：

| 内容 | 来源 |
| --- | --- |
| `golden_expression.tsv.gz` | `data_tcga/processed/01_expression/expression_tcga_tpm.tsv.gz`（原始 TPM，12 例） |
| `golden_clinical.tsv` | `02_labels` 的 age/sex/stage（刻意保留缺失以覆盖冻结参考值填充） |
| `golden_labels.tsv` | `02_labels` 的 cancer_type / os_time_days / os_event |
| 43 concepts | `data_tcga/processed/04_embedding/embeddings_43.tsv.gz`（**真实癌种码**） |
| 132 gene sets（真实码） | 生产同一调用重跑（`predict` + `enable_geneset_capture`） |
| 132 gene sets（占位码 5） | `scripts_dev/analysis_out_concepts/tcga_genesets.tsv.gz`（独立复核 132 层） |
| M0–M3 风险 | 原项目 `scripts/12_validate_external.py::_build_features` + lock(β/scaler)，**逐患者**用各自癌种码 |

判据：`max|Δ| ≤ 1e-5`（132 / 43 / M0–M3）。
重建：`python tools/build_fixtures.py --source-root <原项目根>`
（需原项目只读访问；重建后 `golden_provenance.json` 的哈希会更新）。

⚠ 夹具含 TCGA 公开样本的表达数值，**仅用于复现测试**。
