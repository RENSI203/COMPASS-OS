# 05 · 最终主要结果与论文声明核查 ｜ main results and claim–evidence audit

**阶段**：5 / 6（顺序发布审查）
**日期**：2026-10-04
**候选**：Phase 4 冻结 `c5de83e` → **本阶段修复提交见 §0.2**
**结论**：**READY WITH LIMITATIONS**（详见 `FINAL_RELEASE_READINESS.md`）

> 本阶段**只验证既有结果与声明**：未新增结果、未做新机制发现、
> **未重新挑选"表现最好"的展示条件**。所有数值均为**既有结果文件的复核**。

---

## 0. 方法与声明来源

### 0.1 声明来源（分层）

| 层 | 文件 | 是否本候选 |
|---|---|---|
| 软件声明 | `README.md`、`CITATION.cff`、`docs/*.md` | **是**（已提交） |
| 论文声明 | `paper_assets/*.md`、`MANUSCRIPT_NUMBERS.tsv` | **否**（并行工作区，**只读检查**，未修改） |
| 数值来源 | `validation/results/*.tsv`、`src/compass_os/assets/models/*` | **是** |

### 0.2 本阶段修复

| # | 问题 | 修复 |
|---|---|---|
| 1 | `CITATION.cff` 的 `version: 1.1.0` 与 `date-released: 2026-10-03`（**v1.0.1 的日期**）不匹配 | 已改为 `2026-10-04`（候选冻结日）并注明**发布时须改为实际发布日期** |

---

## 1. Claim–Evidence 总表

> 每条给出：**声明 → 精确数值/定义/分母 → 结果文件与生成脚本 → 独立复核 → 允许表述与限制**。

### C1 · 冻结模型与表示层可复现

| 项 | 内容 |
|---|---|
| **声明** | 冻结包可复现原模型开发流程的产物 |
| **精确数值** | 132 signatures `max\|Δ\| = 1.110e-16`；43 concepts `max\|Δ\| = 2.372e-07`；M0–M3 风险 `max\|Δ\| = 1.060e-06`（M0 2.220e-16 / M1 1.060e-06 / M2 1.011e-06 / M3 6.181e-07） |
| **指标定义** | 逐元素 `max\|Δ\|`，对 golden fixture 的生产产物；判据 **1e-5**（`tests/_util.py: TOL_GOLDEN`） |
| **分母** | 12 例 fixture × (132 / 43 / 4 模型) |
| **来源** | `tests/test_reproducibility.py`；fixture 来源见 `tests/fixtures/golden_provenance.json` |
| **复核** | **本阶段独立测量**：三项**与 README 声明逐位吻合**（1.110e-16 / 2.372e-07 / 1.060e-06）；判据 1e-5 → **PASS** |
| **允许表述** | ✅ "在 golden fixture 上复现生产产物，误差远低于 1e-5 判据" |
| **限制** | 仅 12 例 fixture、**strict** 策略；不等于"任意输入都可复现" |

> ⚠️ 另有一组**不同口径**的数字在 `MANUSCRIPT_NUMBERS.tsv`：
> `frozen_risk_path_max_abs_delta_M0..M3` = 1.11e-16 … 5.00e-16，那是
> **risk = Xβ 的独立 NumPy 再推导**，不是全链路复现。**两组数字不可混用**。

### C2 · 内部/外部预后性能与相对临床基线的增益

| 项 | 内容 |
|---|---|
| **可用证据** | 外部队列（24 个 M0-可评估队列）Uno's C 中位：**M0 0.5948 / M1 0.5215 / M2 0.5776 / M3 0.5831** |
| **配对增益（DL 随机效应）** | `M2 − M0 = +0.008075`，**95% CI [−0.0039, +0.0200]（跨 0）**；`M1 − M0 = −0.049457`，CI [−0.0829, −0.0160]（**M1 显著劣于临床基线**）；`M3 − M0 = +0.011638`；`M3 − M2 = +0.001681` |
| **未配对中位** | **M2 (0.5776) < M0 (0.5948)** —— 与配对估计方向相反 |
| **来源** | `paper_assets/MANUSCRIPT_NUMBERS.tsv`（并行工作区） |
| **复核** | 读取既有数值；**未**独立重算（需源项目外部结果表，不在候选内）→ **NOT VERIFIED（独立复核）** |
| **允许表述** | ⚠️ **只能**写："外部队列中，配对 DL 估计显示 M2 相对临床基线的增益为 +0.008（95% CI 跨 0，未达显著）；未配对中位数二者相当（M0 略高）。" |
| **❌ 禁止表述** | **"COMPASS 显著优于临床变量"**、**"相对临床基线有明确增量价值"**、任何把 +0.008 说成"提升"而不附 CI 的写法 |
| **限制** | CI 跨 0；未配对与配对方向不一致；24 队列（非全部 38）；**zero-refit** 口径 |

> **这是本篇最重要、也最容易被误写的声明。** 证据支持的是
> **"与临床基线相当，未证明有增量"**，不是"优于临床基线"。

### C3 · 输入覆盖度与缺失稳健性

| 项 | 内容 |
|---|---|
| **声明** | 真实公共队列的 COMPASS 基因覆盖并非总是完整 |
| **精确数值** | 38 个可还原队列（**填补前**）：global 覆盖中位 **96.68 %**（范围 **63.55–98.03 %**）；916 signature 基因覆盖中位 **96.18 %**（范围 **80.35–97.82 %**） |
| **指标定义** | 逐队列「词表中出现的基因数 / 15,672」与「916 signature 基因中观测到的比例」，**填补前**计算（Phase 1 §3.3 已实证不被填补污染） |
| **分母** | 39 队列中 `status=='ok'` 的 **38** 个（METABRIC `unavailable`） |
| **来源** | `validation/results/external_coverage_summary.tsv`；`validation/extract_platform_masks.py` |
| **复核** | **本阶段独立重算**：中位 0.9668 / 0.9618，范围 0.6355–0.9803 / 0.8035–0.9782 —— **与 README 声明吻合** ✅ |
| **允许表述** | ✅ "38 个外部队列的中位覆盖 96.7 %，最低 63.6 %；覆盖不足是常态而非例外" |
| **❌ 禁止表述** | "任意不完整转录组均适用"；"模型对缺失基因不变" |
| **限制** | README §9 已写明"稳健性只在已测试区间内成立，**不得**称模型对缺失基因不变" ✅ |

**缺失稳健性曲线**（global 随机掩码，reference 策略，12 队列）：

| 掩码比例 | 5 % | 10 % | 20 % | 30 % | 40 % | 50 % |
|---|---|---|---|---|---|---|
| risk Spearman | 0.9863 | 0.9593 | 0.9230 | 0.8449 | 0.7928 | 0.7413 |

对照臂：`non_sig`（只掩非 signature 基因）0.05→**1.0000**、0.5→**0.9998**；
`signature`（只掩 916 基因）0.01→0.9987、0.1→0.9643。
⇒ **稳健性主要来自 signature 基因的保留**（README §417 已正确表述为
"accounts for a substantial component of the robustness"）。

### C4 · reference 与 zero 的风险层 vs 表示层差异

| 层 | reference vs zero | 判定 |
|---|---|---|
| **风险层** | Δ(Spearman) ≤ **0.0003**（global 0.05–0.3；signature 0.01–0.1；non_sig 全部 ≈0） | 策略选择对**风险排序**几乎无影响 |
| **表示层（concept MAD）** | 0.05：0.00999 → 0.01572（**+57 %**）；0.30：0.0524 → 0.0796（**+52 %**） | **zero 明显更差** |
| **表示层（signature MAD）** | 0.05：0.01119 → 0.01768（**+58 %**）；0.30：0.0609 → 0.1070（**+76 %**） | **zero 明显更差** |

| 项 | 内容 |
|---|---|
| **来源** | `validation/results/stress_pancohort_summary.tsv` |
| **复核** | 本阶段直接读取既有结果文件并计算差值 |
| **允许表述** | ✅ "两种缺失策略在**风险排序**上几乎等价（Δ≤0.0003）；但在**表示层**，`zero` 的扰动幅度比 `reference` 高约 50–80 %" |
| **❌ 禁止表述** | **"两种策略等价"**（对表示层不成立）；**"zero 与 reference 可互换"** |
| **限制** | 未在真实缺失模式下比较两种策略的**下游生物学结论**差异 |

### C5 · empirical platform mask 的解释

| 项 | 内容 |
|---|---|
| **数据** | `mask_type=empirical`，`level=0.00`，`n_cohorts=38`，`risk_spearman_median=1.000000` |
| **含义** | level=0 ⇒ **未施加任何掩码**；`risk_spearman = 1.0` 是**同一输入两次前向的恒等性** |
| **复核** | 本阶段读取确认；`reference` 的 `risk_mad = 0.0` 亦印证"完全恒等" |
| **允许表述** | ✅ "作为**恒等性/流程自检**，level=0 条件下风险排序与全输入完全一致" |
| **❌ 严格禁止** | **把 `risk_spearman = 1.0` 表述为"独立泛化证据"、"外部验证通过"或"模型在真实平台上稳健"** |
| **附注（细微但重要）** | 同一 level=0 下 `zero` 的 `risk_mad = 0.0658` 而 `reference` 为 0.0：因为队列**天然缺失**的基因在两种策略下取值不同（0.0 vs 7.76 TPM）。这是**策略效应**，不是掩码效应，**不得**当作"掩码导致的风险变化" |

### C6 · QC 阈值依据、验证范围与未完成轴

| 项 | 内容 |
|---|---|
| **阈值** | `recommended_global_coverage=0.90`、`warning=0.70`；`recommended_Gsig=0.90`、`warning=0.70` |
| **依据** | 实测曲线的最低满足点：coverage 0.90 → risk_spearman 0.9593 / group_agreement 0.9917（满足 recommended 0.95/0.90）；0.70 → 0.8449 / 0.9917（满足 warning 0.80/0.70） |
| **分母** | **12 队列**（`curve_cohorts: 12`）× 20 repeats；经验轴 38 队列 |
| **来源** | `src/compass_os/assets/models/qc_config.json`（version **1.0.0**，独立于软件版本）；`validation/results/qc_threshold_audit.tsv` |
| **复核** | Phase 1 §9 已核对曲线与阈值自洽；本阶段确认资产未变 |
| **允许表述** | ✅ "经验标定的**描述性工程分级**（`threshold_kind: empirically_calibrated_descriptive_tiers`）" |
| **❌ 禁止表述** | **"经过验证的安全阈值"**、**"临床验证的 coverage 下限"**、**"低于 0.9 不可用"** |
| **未完成轴** | **`signature_coverage_warning_threshold = null`**（资产明文"未校准且**不予发布**"）⇒ 逐 signature 阈值**不存在**，只能报连续值 |
| **限制** | 阈值在 **12 队列**上标定；跨平台稳健性未单独验证（Phase 1 §9 记为 NOT VERIFIED） |

### C7 · 43 concept / 132 signature / gene score 的实际输出与语义

| 输出 | 形状 | 语义 | 来源 |
|---|---|---|---|
| 132 signature | `(n, 132)` | COMPASS 基因集注意力聚合层的**表征坐标** | `predict().signature_scores`；`docs/CONCEPT_SEMANTICS.md` |
| 43 concept | `(n, 43)` | 上层概念注意力聚合层的**表征坐标** | `predict().concept_scores` |
| gene score | `(n, 15672)` | 官方 `extract(with_gene_level=True)` 的 `dfg` = 共享 `nn.Linear(32→1)` 对基因编码的投影 | `sample_path._extract_gene_level`；上游 `tune.py:289-291` |

| 项 | 内容 |
|---|---|
| **复核** | Phase 3 §1：`get_representation` 与 `predict` 的 132/43 两层 **逐位一致（0.000e+00）** |
| **允许表述** | ✅ "43/132 是模型的**表征坐标**，可用于关联与优先级排序"；"gene score 是**模型内部分数**" |
| **❌ 严格禁止** | **"concept/signature 代表通路激活或抑制"**、**"gene score 是基因对风险的因果贡献"**、**"路径图证明因果机制"** |
| **依据** | `docs/CONCEPT_SEMANTICS.md:22` 明确"只支持 association / prioritisation；**不输出** activation / suppression / causal mechanism"；README §53/60/461 同样约束 ✅ |

### C8 · sample_path 的真实数据可用性与统一导出

| 项 | 内容 |
|---|---|
| **真实运行** | Phase 3 首次完成**真实 TCGA 多癌种 + 真实临床**（BRCA/LUAD/THCA，n=120；分层子集 n=30）端到端 |
| **五种预算** | M1 16 / M2 单 13+3 / M2 多 12+4 / M3 单 12+4 / M3 多 11+5，**全部 = 16** 且与规格表逐一吻合 |
| **风险位置** | 图中 η 与主 API **逐位一致**；完整分解 `sum vs η` = **2.220e-16** |
| **统一导出** | PNG/PDF/SVG 由**同一 Figure** 导出；HTML 内嵌**同一份 PNG 字节**（sha256 相等） |
| **允许表述** | ✅ "可在真实多癌种队列与真实临床数据上运行；导出格式同源" |
| **❌ 禁止表述** | "展示阈值 = 科学阈值"；"该图支持因果解读"；"截断后的图等价于完整模型" |

### C9 · 软件新增功能 vs 原 COMPASS 预训练模型贡献

| 归属 | 内容 |
|---|---|
| **上游 COMPASS 提供** | 预训练权重 `pretrainer.pt`、132 gene-set 与 43 concept 的注意力拓扑、gene-token 编码器；上游论文 = *Nat Med* 2026（README §645、CITATION `references`） |
| **本包新增（软件贡献）** | 冻结推理封装、Cox 头（M0–M3）与锁定 PCA、三种缺失策略、覆盖度 QC 与分级、`sample_path` 归因与渲染层、CLI/API、打包与发布审计 |
| **判定** | CITATION 的 `abstract` 明确"frozen-asset inference package … **no model training, re-fitting or tuning happens at prediction time**"；`references.notes` 明确"COMPASS-OS incorporates vendored components from COMPASS under its upstream MIT license" ✅ |
| **允许表述** | ✅ "本包是**推理/工程**贡献，模型与表征来自上游 COMPASS" |
| **❌ 禁止表述** | 把上游概念表征发现写成 COMPASS-OS 的成果；暗示本包训练/改进过模型 |

---

## 2. 六类禁止声明的逐项扫描

| # | 禁止内容 | 扫描结果 | 判定 |
|---|---|---|---|
| 1 | 任意不完整转录组均适用 | README/docs **无**此类表述；README §9 明写"**不得**称模型对缺失基因不变" | **PASS** |
| 2 | 风险排序 = 绝对生存校准 | README §278 明写 "Risk scores order samples and are **not calibrated**"；§567 "discrimination transfers better than absolute-risk calibration"；`survival.py:17-19` 同 | **PASS** |
| 3 | gene/concept 分数或路径图证明因果机制 | 仅以**否定形式**出现（`README.md:53/60/461`、`CONCEPT_SEMANTICS.md:22`、`SAMPLE_PATH.md:9`） | **PASS** |
| 4 | 免疫相关表示证明免疫治疗响应预测 | 仅在**上游论文标题**中出现（README §645/651/685、CITATION `references`），已明确归属上游；本包**无**免疫治疗响应声明 | **PASS** |
| 5 | 展示阈值 = 科学/临床验证阈值 | `README.md:208` "display parameters only — not clinically validated thresholds"；`:254`；`SAMPLE_PATH.md:128` | **PASS** |
| 6 | 未验证的平台/癌种/输入尺度已充分支持 | README §310/572/582：微阵列需自行和声化、v1 **无** `input_scale="microarray"`；`ADVANCED_USAGE.md:68-69` 同；Phase 1 §2 记录 counts 不可检测 | **PASS** |

⇒ **六类禁止声明在候选的软件声明中全部未出现。** ✅

---

## 3. 高低风险免疫程序差异等探索性声明的 outcome-conditioned selection 检查

| 检查 | 结果 |
|---|---|
| 候选内是否存在此类声明 | **未发现**（`grep -i "immune program|免疫程序|immune axis|high-risk.*immune"` → 0 命中） |
| 上游机制工作区（`downstream_mechanism/**`，**开发机路径已脱敏**）是否有此类分析 | **有**，但属**只读并行工作区**，**不在本候选内**，本阶段**未纳入**、未重跑 |
| 风险 | 若此类声明进入 Application Note，须自行核查其队列选择是否受结局条件化影响 |

**正式表述**：

> 本候选**不包含**任何"高风险 vs 低风险免疫程序差异"类声明。
> 该类分析位于并行机制工作区，**其数据来源、统计口径与队列选择未经本阶段核查**；
> **不得**因本软件升级而自动提升其证据等级（Phase 5 规则明令）。
> 若后续要写入论文，须**单独**做 outcome-conditioned selection 审查。

### 3.1 本阶段的第三次自我命中（已加机械化防护）

**现象**：本文件初稿在 §3 引用了机制工作区的**开发机绝对路径**作为证据，
导致本文件自身被 `release_audit.py` 判为 `absolute-path` FAIL。

**这是同一类错误的第三次发生**（Phase 0 §6.5、Phase 4 §7.1 已各记录一次）。
前两次都只写了"教训"，**没有机械化防护**，因此第三次仍然发生。

**处置（本次不再只写教训）**：

1. 脱敏该行（保留目录名与结论，去掉机读路径前缀）；
2. **新增回归测试** `tests/test_audit_docs_paths.py`：
   扫描 `docs/audit/*.md` 与 `*.tsv`，出现**家目录前缀**、**挂载盘前缀**、
   **Windows 盘符前缀**、**项目目录名**或**开发机用户名**即**失败**
   （检测口径与 `tools/release_audit.py` 的 `absolute-path` 一致；
   为避免本文件自身命中，此处**只描述类别、不写出字面量**）。
   把"人要不要记得脱敏"变成**CI 可执行的断言**。

---

## 4. 正文/图/表/补充材料/README/CITATION 的一致性

| 检查 | 结果 |
|---|---|
| 软件版本 | `pyproject.toml` / `__init__.py` / `CITATION.cff` / `README` 徽标 = **1.1.0** ✅ |
| 科学资产版本 | `qc_config.json` = **1.0.0**，与软件版本**独立** ✅（未误 bump） |
| `CITATION.cff` 日期 | 原为 v1.0.1 的 `2026-10-03` 配 `version: 1.1.0` ⇒ **不一致，已修复**（§0.2） |
| README 覆盖度数字 vs 数据 | 96.7 %/63.5–98.0 % 与 96.2 %/80.3–97.8 % ⇒ **独立重算完全吻合** ✅ |
| README 复现数字 vs 实测 | 1.1e-16 / 2.4e-07 / ≤1.06e-06 ⇒ **实测 1.110e-16 / 2.372e-07 / 1.060e-06，吻合** ✅ |
| 风险方向 | 全文一致："higher = higher model-estimated hazard"；`risk` 为 η，非概率 ✅ |
| 术语 | "concept/signature" 一律为表征坐标；"cutoff" 区分 frozen / cohort-relative / visualization ✅ |
| **`MANUSCRIPT_NUMBERS.tsv` 版本字段** | ⚠️ 仍为 `software_version 1.0.1` / `v1.0.1` / commit `7270935` / 旧 wheel-sdist 哈希 ⇒ **与 1.1.0 候选不一致** |

> ⚠️ **版本不一致项（第 3 处）**：`paper_assets/MANUSCRIPT_NUMBERS.tsv` 是**并行工作区**文件
> （未跟踪/已修改，**不在本候选内**）。若 Application Note 随 1.1.0 发布，
> 其 `software_version` / `software_release_tag` / `software_release_commit` /
> `wheel_sha256` / `sdist_sha256` **必须更新**。
> 本阶段**未修改该文件**（不属本候选，且并行工作区只读）。

---

## 5. PASS / FAIL / NOT VERIFIED 汇总

| # | 事项 | 判定 |
|---|---|---|
| 1 | 冻结模型与表示层可复现性 | **PASS**（三项独立实测吻合；判据 1e-5） |
| 2 | 内部/外部预后性能，特别是相对临床基线增益 | **PASS（证据存在）**；**独立复核 NOT VERIFIED**（需源项目结果表）；**表述须附 CI 跨 0** |
| 3 | 输入覆盖度与缺失稳健性 | **PASS**（覆盖度独立重算吻合） |
| 4 | reference/zero 的风险层与表示层差异 | **PASS**（风险层 Δ≤0.0003；表示层 zero 差 50–80 %） |
| 5 | empirical mask 解释 | **PASS**（已明确 level=0 为恒等性，**非**泛化证据） |
| 6 | QC 阈值依据/范围/未完成轴 | **PASS**（阈值 = 描述性工程分级；`signature` 阈值为 null） |
| 7 | 43/132/gene score 输出与语义 | **PASS**（语义边界在代码与文档双向约束） |
| 8 | sample_path 真实可用性与统一导出 | **PASS**（真实多癌种 + 真实临床；导出同源） |
| 9 | 软件 vs 上游 COMPASS 贡献区分 | **PASS**（CITATION 与 README 均明确归属） |
| — | 六类禁止声明 | **PASS**（全部未出现） |
| — | 版本/数值/方向/术语一致性 | **PASS（1 处已修复）**；`MANUSCRIPT_NUMBERS.tsv` 版本陈旧 = **并行工作区待办** |
| — | 探索性免疫声明 | **不在候选内**；**未核查**其 outcome-conditioned selection |

**本阶段缺陷 1 处（CITATION 日期），已修复。无 blocker。**

---

## 6. 剩余限制（须写入最终披露）

1. **M2 相对临床基线的增益未达显著**（+0.0081，95 % CI [−0.0039, +0.0200] 跨 0），
   且未配对中位数 **M0 略高于 M2**。
2. **M1（纯 concept）显著劣于临床基线**（−0.0495，CI [−0.0829, −0.0160]）。
3. **外部队列性能的独立复核 NOT VERIFIED**（源项目结果表不在候选内）。
4. **绝对生存 `S(t)` 无校准证据**（Phase 1/3 已记）——判别度不能替代校准。
5. **QC 阈值仅在 12 队列上标定**；跨平台稳健性未验证。
6. **`signature_coverage_warning_threshold` 未标定（null）**。
7. **RNA-seq 检出限型缺失未验证**（39 队列证据主要来自微阵列平台缺失）。
8. **非负 counts / 微阵列强度不可检测**（会被当作 TPM）。
9. **`Age` 无生理范围校验**；**非数值临床值**须调用方自行编码（现已告警）。
10. **干净 clone 无法重建验证结果**（`validation/inputs/` gitignored）。
11. **`MANUSCRIPT_NUMBERS.tsv` 版本字段陈旧**（并行工作区，需在发布前更新）。
12. **历史 `design.json` 只描述 12 队列轴**；38 队列轴的纳入规则无法从仓库恢复（Phase 4）。
