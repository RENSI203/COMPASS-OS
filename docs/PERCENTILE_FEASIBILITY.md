# 单样本 `risk_percentile` 可行性审计（P2 §15；**结论：可行，但需人工确认后才加入 API**）

## 现状

* 单样本只给 `risk_rank`（队列内秩）；没有 reference-based percentile。
* 队列分层的切点是**单一标量** `median_cutoff`（M0 = 1.0851，M2 = 1.0783），
  **不区分癌种**。

## 癌种间风险是否可比？——不可比（有量化证据）

| 证据 | 数值 |
| --- | --- |
| M2 的 `CT_*` 系数（32 列） | 中位 **+0.068**，范围 **−2.015 ~ +2.262**，**极差 4.278** |
| M2 的概念系数 L1 合计 | 1.838（癌种还通过 encoder 的 CANCER token 改变概念值本身） |
| 结论 | 不同癌种的线性预测子被整体平移/缩放 ⇒ **pan-cancer 混合排序没有意义** |

⇒ 若做 percentile，必须 **`model × cancer_type` 分层**，即
`risk_percentile_within_cancer`，而不是一个全局分布。

## 资产与成本

| 项 | 估计 |
| --- | --- |
| 资产规模 | 32 个 `CT_*` × 4 模型 × 101 分位 = 12,928 个 float ≈ **101 KB**（可忽略） |
| 计算成本 | TCGA 9,430 例一次前向（43 概念）+ 4 模型线性组合 ≈ **4.4 min** |
| 依赖 | 需要把 TCGA 逐患者 risk 冻成资产（当前**未落盘**） |

## 科学限制（必须在文档中写明）

1. 参考分布来自 **TCGA 训练域**，且该模型正是在这 9,430 例上拟合的
   ⇒ percentile 是 **in-sample 描述性参考**，不是 calibrated survival probability，
   也不代表外部队列的分位含义相同。
2. **小样本癌种不稳定**：TCGA 33 个癌种 n 范围 **36–1,071**，其中 **7 个 n < 100**
   （ACC 90、MESO 85、UVM 80、KICH 65、UCS 56、DLBC 47、…，最小 36）
   ⇒ 这些癌种的分位数不可靠，需设 **最小 n 规则**（例如 n ≥ 100 用自身分布，
   否则回退 pan-cancer 或做收缩）。

## 建议（待确认）

* 实现为 **独立函数** `risk_percentile(risk, model, cancer_type)`，并入
  `predict` 的返回字段 `risk_percentile_within_cancer`；
* 资产 `models/reference_risk_quantiles.tsv`（≈100 KB）+ `asset_manifest` 登记；
* 命名与文档统一写 **descriptive reference percentile**；
* 同时输出 `reference_n`（该癌种的参考样本量）与 `reference_scope`
  （`cancer_specific` / `pan_cancer_fallback`），让使用者知道分位的来源与可靠性；
* **本阶段不实施**：按 P2 约定需人工确认后再加入正式 API。
