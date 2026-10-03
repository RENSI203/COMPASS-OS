# 缺失基因处理（设计依据与措辞边界）

> 本文件记录 `missing_gene_strategy` 的**源码证据**、数学定义与措辞边界。
> 结论来自对上游 COMPASS 训练代码与冻结 checkpoint 的逐行审计。

## 1. 三种策略的数学定义

记模型词表 `G`（|G| = 15,672，`feature_name` 顺序），冻结 MinMaxScaler 的逐基因
`dmin_g, dmax_g`，冻结参考中位数 `m* = 3.130937933922`（log2 尺度）。输入记为 `t_g`（TPM），
标准化值

```
s_g(t) = ( log2(t_g + 1) − dmin_g ) / ( dmax_g − dmin_g )      # checkpoint 内 Datascaler
```

### `reference`（默认）

```
â_g = m*  (log2 尺度)   ⟺   t_g = 2^{m*} − 1 ≈ 7.76
```

* 唯一、逐样本独立；
* **不含任何用户队列统计量**（不重算 mean/median）⇒ 同一样本的预测不依赖同批其他样本；
* 与项目生产外部验证口径（`scripts/12_validate_external.py`、`scripts_dev/41`）一致：
  `med = median(reference_quantiles["values"])`，实测 `= 3.130937933921814`（与 `41` 的
  `fill_median_log2` 逐位相同）。

### `zero`（训练空间掩码 / training-space masking）

```
s_g ← 0            # 在标准化空间置零
```

等价输入侧写法（仅用于验证，不作为实现）：`t_g = 2^{dmin_g} − 1`。

**源码证据（上游 COMPASS）**

| 事实 | 位置 |
| --- | --- |
| 先标准化、再进 Dataset | `compass/main.py:268,271`（`scaler.fit` → `scaler.transform`） |
| Dataset 直接转 tensor（不再变换） | `compass/dataloader/data.py:56` |
| 掩码实现：基因列置 0 | `compass/augmentor/aug.py:70-74`（`x_new[1:][mask] = 0`） |
| 掩码概率（发布 checkpoint 实测） | `mask_p_prob = 0.41`；`mask_a_prob = mask_n_prob = 0` |
| 逐基因独立 | `torch.rand(len(x[1:]))`（`aug.py:70`），对每个基因独立伯努利 |
| 逐样本独立 | `data.py:80-87` 每样本单独调用；`aug.py:78` 断言单样本向量 |
| 抖动（同空间） | `aug.py:124-125`：`torch.normal(0, jitter_p_std=0.41)` 加到同一 scaled 向量 |
| 掩码/抖动二选一 | `aug.py:182-184`（`np.random.choice([augmentor1, augmentor2])`） |
| 只作用于 positive view | `data.py:86-87,97`（anchor/negative 用 `mask_a_prob = mask_n_prob = 0`） |
| 推理不用 augmentor | `Predictor`/`Extractor` 不调用 augmentor，且 `tune.py:210` `model.eval()` |

**为什么不能实现成 `raw TPM = 0`**（实测，冻结 scaler，15,672 基因）

| 量 | 值 |
| --- | --- |
| `data_min_ > 0` 的基因（即 scaled-0 ≠ TPM-0） | **9,435 / 15,672（60.2 %）** |
| scaled-0 ⟺ `TPM = 2^{dmin_g} − 1`：等于 0 的基因 | 6,237 |
| 上式 > 0 者：中位 / 最大 TPM | 0.506 / 411.66 |
| 若填 `TPM = 0`，得到的 scaled 值分位 [0/25/50/75/100]% | **[−1.396, −0.108, −0.014, 0, 0]** |
| 其中 scaled < 0（训练分布 [0,1] 之外）的基因 | **9,435（60.2 %）** |

⇒ 直接把缺失基因写成 raw `TPM = 0`，对 60.2 % 的基因是把输入推到**训练分布之外**，
而不是复现训练期 masking。`tests/test_reproducibility.py::test_zero_differs_from_raw_tpm_zero`
把这一事实固化成测试。

### `strict`

```
若 G \ observed ≠ ∅  →  raise MissingGenesError(missing_gene_names)
```

用于 benchmark、golden reproduction 与严格推理。

## 2. 措辞边界（必须遵守）

可以说：

> the masking operation is applied in the same normalized input space used during
> COMPASS contrastive pretraining
> （训练空间掩码 / masking-based missing-gene strategy）

**不得**说：

* `identical to the training augmentation distribution`
  —— 训练期掩码只用于正样本视图、概率 0.41、随机选择；推理期是确定性地对真实缺失基因施加；
* `COMPASS supports arbitrarily incomplete transcriptomes`
* `COMPASS is invariant to missing genes`
  —— 只能在 stress test 实测覆盖区间内描述稳健性。

## 3. 覆盖度 QC 的定义

```
signature_gene_coverage_i = |observed ∩ genes(set_i)| / |genes(set_i)|        # i = 1..132
concept_input_coverage_c  = Σ_i w_{c,i} · signature_gene_coverage_i           # c = 1..43
```

* 132 个基因集共引用 **916 个唯一基因**，且 **916/916 全部在 15,672 词表内**（已核）；
* 基因集规模：min 1 / 中位 7 / max 51；
* `w` 为 **冻结 checkpoint** 的 `CellPathwayAttentionAggregator` softmax 注意力
  （非负、组内和为 1；43 组成员数 1–8，合计 132），由
  `compass_os.representation.concept_weights()` 读出，**不是本项目自造权重**；
* 变量名只能用 `input coverage`，**禁止** confidence / biological confidence /
  activation confidence —— 它只表示输入信息覆盖程度。

## 4. 阈值（已由 stress test 冻结 → `models/qc_config.json` v1.0.0，status=complete）

| 阈 | 值 | 依据 |
| --- | --- | --- |
| `recommended_global_coverage` | **0.900** | 满足「risk Spearman ≥0.95 且 agreement ≥0.90」的**最低** global 覆盖度 |
| `warning_global_coverage` | **0.700** | 满足「≥0.80 且 ≥0.70」的最低覆盖度 |
| `recommended_Gsig_coverage` | **0.90** | 同一判据，作用在 916 个 signature 基因轴上 |
| `warning_Gsig_coverage` | **0.70** | 同上 |
| 逐 signature 阈值 | **未发布（null）** | 原 0.912 无法从保留的可复核产物重现（见下），逐 signature 覆盖度改为**连续值**输出 |

**关于逐 signature 阈值**：P3.5 审计发现早先报告的 `signature_warning_rule.threshold = 0.912`
依赖一张**逐条件累积**的 (signature × coverage) 覆盖度–误差表；该中间件未随结果保留，
现存表中 `coverage_mean` 对 132 个 signature **恒为 1.0**（无覆盖度变异），无法据以估计阈值。
因此该字段**置 null 并撤下**，per-signature coverage 继续以连续值报告
（`signature_gene_coverage`，n×132），**不使用任何不可复现的数值阈值**。

**overall 组合规则**：`overall = worse(global tier, Gsig tier)`（取更差的一档）。

**三臂对照（global / signature-only / non-signature-only，各 12 队列 × 20 重复）**：
global 与 signature-only 在可比覆盖度下接近（最大差 0.012）；
**non-signature-only 臂几乎无影响**：

| 该臂自身轴上的覆盖度 | 0.95 | 0.90 | 0.80 | 0.70 | 0.60 | 0.50 |
| --- | --- | --- | --- | --- | --- | --- |
| non-signature-only：risk Spearman | 1.0000 | 1.0000 | 0.9999 | 0.9999 | 0.9998 | 0.9998 |
| non-signature-only：mean\|Δrisk\| | 0.0022 | 0.0045 | 0.0087 | 0.0133 | 0.0169 | 0.0218 |
| non-signature-only：ΔUno C | 0.0000 | 0.0000 | 0.0000 | 0.0001 | 0.0007 | 0.0008 |

即 **916 个 signature 基因全部保留**时，删掉其余 14,756 个基因的一半，
风险排序、分层、43/132 表征几乎完全不变。

✅ 允许的表述：**Preservation of signature-related genes accounted for a substantial
component of prognostic robustness.**
❌ 不允许：`non-signature genes have no contribution`——mean|Δrisk| 随缺失从 0.002 升到
0.022、ΔUno C 升到 +0.0008，说明它们并非毫无作用。

⚠ 曲线平滑、**无明显断点**，因此这些阈值是 **empirically calibrated descriptive
engineering tiers**（工程 QC 分层），**不是生物学边界、也不是"安全线"**；
`qc_config.json` 的 `threshold_kind` 已如此标注。阈值一律取**实测覆盖度水平**（不插值），
逐 level 判据满足情况见 `validation/results/qc_threshold_audit.tsv`。
