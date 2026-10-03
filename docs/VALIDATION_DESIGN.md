# Missing-gene 验证设计（P2）

## 为什么要三类 missingness

真实数据的缺失不是均匀随机的：平台探针覆盖、基因长度/GC、以及"某些 gene set 整体不在
芯片上"都会造成**结构化缺失**。因此验证分三类：

| 类型 | mask 对象 | 语义 |
| --- | --- | --- |
| `global` | 全部 15,672 冻结 feature | 模拟"平台缺失"的随机部分 |
| `signature` | 916 个 signature 相关基因 | 控制性缺失；小 gene set 少量关键基因缺失即可见效应 |
| `empirical` | **每个队列真实平台的 pre-imputation 观测集合** | 直接 replay 真实平台缺失模式 |

每个重复内所有样本共用同一 mask 基因集合（模拟平台而非个体效应）。

## 三个比较基准

`FULL`（完整输入）为 perturbation reference；扰动有两种填补口径：

* `REFERENCE`：冻结 TCGA 参考中位数（生产默认）
* `ZERO`：训练空间掩码（标准化空间置 0）

## ⚠ 一个必须写进论文的方法学事实

**在 `empirical` 实验中，`REFERENCE` 结果是恒等映射**（risk Spearman = 1.0000、
Δ=0.0000、concept/signature r = 1.0000）。原因不是 bug，而是口径自洽：
外部队列的和声化缓存**本来就是用参考中位数填补缺失基因**的（`47` 的
`reindex(vocab).fillna(med)`），所以"mask 掉这些基因 → 再用参考中位数填回"
逐位还原了生产输入。

因此：

* `empirical × reference` = **生产口径的自洽性检查**（必须为恒等，实测为恒等 ✓）；
* `empirical × zero` = 真正有信息量的那一臂（同一真实缺失模式、换一种填补口径）；
* `global/signature × reference` 则**不是**恒等——被 mask 的基因在缓存里本来有真实观测值。

## 成本口径（如实报告）

COMPASS 前向实测 ≈ 28–33.5 ms/样本（CPU；加大 torch 线程无效，改用多进程分片）。
本轮采用**计算受限但完整**的设计：

| 项 | 设置 |
| --- | --- |
| 覆盖度曲线（global/signature） | **12 个队列**（按覆盖度等距选取，跨 6 个平台），每队列 ≤ 60 例 |
| global levels | 5/10/20/30/40/50 %，每 level **20 次重复** |
| signature levels | 1/2/5/10/20/30 %（占 916 基因），每 level 20 次重复 |
| empirical replay | **全部 38 个可还原队列**，每队列 ≤ 50 例，每策略 1 次 |
| 分片 | 4 进程 × 6 torch 线程 |
| 全网格 9,359 例 × 完整 level 网格 | 估算 > 40 h（故未采用） |

`--max-per-cohort` / `--cohorts` / `--repeats` 均可调大以做更充分的运行。

---

## 结论分级（Application Note 用词纪律）

写论文/说明文档时，结论只能落在下面四档之一，且必须与证据档位一致。

### Strongly supported（多来源一致、效应明确）

* **三类 masking 的分化结论**：non-signature-only 与 signature-only 在同等缺失比例下
  相差极大（0.50 覆盖度：0.9998 vs 0.849），说明稳健性主要挂在 signature 基因上。
* **`reference` 填补在低覆盖下比训练空间 `zero` 掩码更好地保留表征读数**：
  50% 全局覆盖时 43-concept 稳定性 0.804 vs 0.621、132-signature 0.889 vs 0.718，
  而 risk Spearman 两者相同（0.741 vs 0.741）；低覆盖真实平台掩码下同样成立
  （0.944 vs 1.000 的概念相关性，Wilcoxon p = 7.3e-12，n = 38 队列）。
* **词表与输入契约无歧义**：15,672 个 feature 符号全唯一；重复列名被拒绝而非猜测。
* **冻结资产可逐位复现**：官方 `extract()` 与生产路径 132/43 两层 `max|Δ| = 0.0`；
  对原项目存档 132 = 1.1e-16、43 = 2.4e-07、M0–M3 ≤1.1e-06。

### Supported within the tested range（只能按实测区间陈述）

* **中等缺失下风险排序高度一致**：全局覆盖 0.95/0.90/0.80 时 risk Spearman
  0.986/0.959/0.923、分层一致率 ≥0.98；真实平台缺失模式下（38 队列）Spearman = 1.000。
  必须与"覆盖率下降时排序稳定性平滑下降（50% 时 0.741）"同时报告。
* **Across the tested masking range, changes in Uno's C remained within approximately
  ±0.02, with no consistent deterioration across the evaluated cohorts.**
  （不得写成 "discrimination was unaffected"。）
* **Preservation of signature-related genes accounted for a substantial component of
  prognostic robustness**：916 个 signature 基因全保留时，删除其余 14,756 个基因的
  50% 仍使 risk Spearman = 0.9998、分层一致率 1.00、43/132 表征相关性 1.000、
  ΔUno C = +0.0008；而删除**同样比例**的 signature 基因（覆盖度 0.50）时为 0.849。
  不得写成 `non-signature genes have no contribution`（mean|Δrisk| 0.002→0.022）。
* **覆盖度 QC 具有预测性**：逐特征输入覆盖度与扰动误差的相关在 43/43 概念与
  132/132 signature 上均为负（中位 Spearman −0.68 / −0.60）。

### Engineering calibration（工程阈值，不是生物学边界）

* `recommended_global_coverage` / `warning_global_coverage` 与对应的 Gsig 阈值：
  由**实测离散覆盖度点**按预先写明的判据（risk Spearman ≥0.95/≥0.80 且
  agreement ≥0.90/≥0.70）选出，**不做插值**；曲线平滑无断点，
  因此它们是 *empirically calibrated descriptive tiers*。
* `signature_coverage_warning_threshold`：逐 signature 的工程告警线，
  依据是该 signature 的平均扰动误差超过低误差底噪的 2 倍。

### Not supported（不得声称）

* ❌ `COMPASS supports arbitrarily incomplete transcriptomes` / 对缺失基因不变。
* ❌ `zero` 模式复现训练增强分布（只能说*same normalized input space*）。
* ❌ `non-signature genes have no contribution`（对照只能说明"保留 signature 基因
  解释了稳健性的相当一部分"，不能说明其余基因无贡献）。
* ❌ 43 概念的 biological activation / suppression / causal mechanism。
* ❌ 绝对生存概率的个体校准；单样本 reference percentile（仅完成可行性审计）。
* ❌ 直接接受原始微阵列强度矩阵。

### 恒等臂的正确表述

`empirical × reference` 是**production-pipeline self-consistency control**
（和声化缓存本身即参考中位数填补，因此该臂按构造返回原输入），
**不得**作为独立稳健性证据。
