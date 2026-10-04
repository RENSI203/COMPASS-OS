# COMPASS-OS

**Pan-cancer survival prediction and mechanism-oriented transcriptomic representation
framework.**

**泛癌生存预测与机制导向转录组表征框架。**

* **Version / 版本:** 1.1.0 · **Python package:** `compass_os` · **License:** MIT
* **Repository:** `https://github.com/RENSI203/COMPASS-OS`

> **Language / 语言**：Each section is written in English first (canonical release text for
> the software publication), followed by a condensed Chinese version for biomedical
> researchers. Code, API names, parameters, file names and CLI commands are kept in English
> throughout.
> 每一节先英文（软件发布用的正式表述），其后为面向生物医学研究者的简明中文。
> 代码、API 名称、参数、文件名与命令行一律保持英文。

Built on bulk transcriptomic expression + cancer type (+ optional clinical variables), it
returns:

| | |
|---|---|
| **Prognostic risk** | a relative risk score per sample |
| **Cohort risk stratification** | high/low grouping within a submitted cohort |
| **43 concept representations** | COMPASS latent coordinates |
| **132 gene-signature representations** | gene-set-level profiles |
| **Input / robustness QC** | gene, signature and concept-level coverage |

It is a frozen-asset inference package: no training, no re-fitting, no tuning at
prediction time. The representations are hypothesis-generating features for downstream
prioritisation — this is not automatic mechanism discovery.

**中文**：输入为 bulk 转录组表达谱 + 癌种（可选临床变量），输出**预后风险**、**队列内风险分层**、
**43 个 concept 表征**、**132 个 gene-signature 表征**与**输入质控**。模型参数全部冻结，
预测阶段不训练、不重新拟合、不调参。表征只用于**下游假设优先级排序**，不是自动机制发现。

---

## 1. What is COMPASS-OS? ｜ 这是什么？

| | |
|---|---|
| **Input** | bulk expression (`samples x genes`, TPM or log2(TPM+1)), a cancer type per sample, optionally clinical variables and survival follow-up |
| **Output** | relative risk score + risk ranking/stratification · 43 concept profiles · 132 gene-signature profiles · coverage QC · a plain-text report and TSV/PDF files |
| **Model** | a frozen ridge Cox head on the COMPASS representation (default: cancer type + age/sex/stage + 43 concepts) |
| **Trained on** | TCGA pan-cancer overall survival |

Three kinds of output, deliberately kept separate:

1. **Prognostic outputs** — risk score, ranking, stratification.
2. **Mechanism-oriented representations** — 43 concepts, 132 gene signatures.
   These are **hypothesis-generating features** for downstream prioritisation; they are
   *not* evidence of pathway activation, suppression, or causal mechanism.
3. **Input / robustness QC** — gene coverage, signature coverage, concept input coverage.

**中文**：三类输出**刻意分开**，不要混为一谈——

1. **预后输出**：风险分数、队列内排序、高低危分层；
2. **机制导向表征**：43 concepts + 132 signatures，属**假设生成特征**，用于下游优先级排序，
   **不是**通路激活、抑制或因果机制的证据；
3. **输入/稳健性质控**：基因覆盖、signature 覆盖、concept 输入覆盖。

---

## 2. Quick start ｜ 快速开始

```bash
git clone https://github.com/RENSI203/COMPASS-OS.git
cd COMPASS-OS
pip install .
```

Model assets and the vendored COMPASS package ship **inside the Python package**
(`compass_os/assets/`), so no extra environment variable is needed. `COMPASS_OS_ROOT` is a
developer override only (see [Installation](#10-installation--安装)).

```python
import compass_os

result = compass_os.analyze(
    "expression.tsv",              # samples x genes, TPM
    cancer_type="LUAD",            # TCGA acronym
    clinical="clinical.tsv",       # optional: age / sex / stage
    survival="survival.tsv",       # optional: time + event columns
)

result.summary()
result.save_report("compass_results/")
```

That is the whole workflow. Real output of
[`examples/quick_start.py`](examples/quick_start.py) (5 LUAD samples):

```text
COMPASS-OS analysis
────────────────────────────────────────────────────────────────
Samples                         5
Cancer type                     LUAD
Model                           Default prognostic model (M2)

INPUT QUALITY
Gene coverage                   100.0% (15672/15672 genes)
Signature-related coverage      100.0% (of the 916 signature genes)
Concept input coverage          100.0% (median across 43 concepts)
QC                              Recommended range (global 100.0%, signature genes 100.0%;
                                thresholds global ≥ 90%, signature ≥ 90%)

SURVIVAL
Risk estimates                  Available (prognostic risk score; Cox linear predictor)
Risk stratification             Not available (cohort of 5 is too small)

BIOLOGICAL REPRESENTATIONS
43 concept profiles             Generated
132 gene-signature profiles     Generated
Interpretation                  hypothesis-generating features

ROBUSTNESS
Missing-gene strategy           Reference imputation (frozen median)
Robustness check                Not required (coverage in recommended range)

CLINICAL VARIABLES IMPUTED (frozen training reference values)
  Age: 1/5 samples → 60.0
  Sex: 1/5 samples → 0.0
  Stage: 3/5 samples → 2.0
```

The example is a **single-cancer-type** call (`cancer_type="LUAD"`), and the output above is
exactly what that call prints. To analyse a multi-cancer cohort, pass one acronym per
sample (`cancer_type=["BRCA", "LUAD", ...]`); see
[`examples/cohort_example.py`](examples/cohort_example.py). A 48-sample multi-cancer cohort
additionally reports the high/low split and a log-rank p-value.

You do **not** need to know about model variants, scaling, PCA, masking spaces, cancer
tokens or design matrices. Everything that must be fixed is fixed internally.

Runnable scripts: [`examples/quick_start.py`](examples/quick_start.py) (Level 1),
[`examples/cohort_example.py`](examples/cohort_example.py),
[`examples/missing_genes_example.py`](examples/missing_genes_example.py),
[`examples/minimal_example.py`](examples/minimal_example.py) (Level 2).

**中文**：安装用 `git clone` + `pip install .`（模型资产与 vendored COMPASS **随包分发**，
无需设置任何环境变量；`COMPASS_OS_ROOT` 仅供开发者覆盖）。准备两个文件即可——
表达矩阵（行=样本，列=基因 symbol，TPM）与癌种（TCGA 缩写）。
调用 `compass_os.analyze(...)` 后，`result.summary()` 打印上面的可读报告，
`result.save_report("compass_results/")` 输出 TSV 表格与 PDF/PNG 图件。
**使用者不需要了解** M0–M3 的区别、标准化、PCA、掩码空间、癌种 token 或设计矩阵——
这些都在包内部固定。上面展示的是 `examples/quick_start.py`（5 例 LUAD）的**真实输出**；
多癌种队列请按样本逐一传入癌种缩写，见 `examples/cohort_example.py`。

---

## 2b. Sample-level computation path ｜ 样本级计算路径归因

<div align="center"><em>sample-specific computational attribution / representation-flow
visualization — not a causal mechanism diagram, not a pathway-activation map</em></div>

```python
from compass_os import sample_path

res = sample_path(expr, cancer_type, "TCGA-19-1787", model="M3", clinical=clinical)
print(res.summary())
res.save_html("sample_path.html")     # HTML (same figure, expandable exact node table)
res.save_static("sample_path.pdf")    # .png / .pdf / .svg — same Figure
res.save("out/M3_low")                # all formats + nodes.tsv + selection.json + decomposition
```

One call does everything: input checks → frozen inference → real representation-layer scores →
**full** Cox contribution decomposition → display selection → exportable result object.
You never assemble node tables or JSON by hand.

**Figure structure** (circular layered path):

```
Gene expression → COMPASS gene score → Granular signature score
    → High-level concepts + clinical / PC → Cohort-relative risk bar
```

* **Real scores, not weights.** The gene-score column is the official
  `extract(with_gene_level=True)` gene-token scalar
  (`geneset_scorer(encoder_output)[:, 2:]`, a shared `nn.Linear(32→1)`); the 132 signature and
  43 concept columns come from the **same forward pass as η**. Attention/projector weights are
  used only to decide *which* connections to draw — never as node values.
* **Clinical / Cancer type / PC1–PC10** are aggregated from the **frozen Cox design matrix**:
  multi-column terms take the signed sum of their column contributions (PC1–PC10 is one node =
  sum of the ten frozen contributions, using the original standardisation). Fields that were
  actually imputed are labelled `(imputed)`.
* **Complete accounting.** The full, untruncated decomposition (including undisplayed concepts
  and the hidden cancer term) satisfies `sum(contribution) == η`; the renderer enforces it.
* **No flow claim.** All links share one width and colour and only express the selected model
  connectivity. Node colour encodes that layer's own score (clinical/PC nodes encode their signed
  Cox contribution, with a separate legend).

**Display budgets** (frozen; `gene_threshold=0.25` on the real gene-score scale,
`max_genes≤32`, `max_signatures≤26`, `predictor_budget≤16`, each high-dimensional column ≤ 50):

| model / cohort | high-level nodes shown |
|---|---|
| M1 | 16 concepts |
| M2, single-cancer | 13 concepts + Age + Sex + Stage |
| M2, multi-cancer | 12 concepts + Age + Sex + Stage + Cancer type |
| M3, single-cancer | 12 concepts + Age + Sex + Stage + one `PC1–PC10` node |
| M3, multi-cancer | 11 concepts + Age + Sex + Stage + Cancer type + one `PC1–PC10` node |

Cancer type is hidden for single-cancer cohorts (display only) but stays in the frozen model and
in the full decomposition. Truncation is **real** node/edge removal — never hidden labels, and
never a silent threshold relaxation to fill a layer.

**Tuning** (display parameters only — not clinically validated thresholds):

```python
from compass_os.sample_path_render import Style
style = Style(gene_threshold=0.5, max_genes=20, max_signatures=18, predictor_budget=14,
              expression_limits=(0, 1704.3), contribution_limits=(-1.5, 1.5))
m2 = sample_path(expr, ct, sid, model="M2", clinical=clinical, style=style)
m3 = sample_path(expr, ct, sid, model="M3", clinical=clinical, style=style)  # same Style ⇒ comparable
```

Shorter forms: `sample_path(..., gene_threshold=0.5, max_genes=20, max_signatures=18,
predictor_budget=14)`. Pass the **same** `Style` to M2 and M3 so canvas, node sizes, column
positions, fonts and colour scales match.

**Risk bar** — the cohort-percentile axis (high risk on top), a red→yellow gradient, a triangle
at the sample's **main-API percentile**, `Risk η / Rank / Percentile / N`, and a dashed cutoff.
`cutoff="median"` is the **current cohort's** η median — a visualization boundary, explicitly
**not** the frozen validated prognostic cutoff. It is drawn at its true empirical percentile
(not forced to 50 %), and `rank`/`percentile` are the main API's own values, never recomputed
for the figure.

**Export / dependencies.** PNG, PDF and SVG come from **one** Matplotlib `Figure`; the HTML
embeds the **identical PNG bytes** and adds an expandable table of exact node values (the main
figure in HTML is static — there is no separate interactive renderer and no Plotly/Kaleido/CDN
requirement). Matplotlib is already a core dependency (`AnalysisResult.save_report`); the
drawing layer imports it lazily, so `import compass_os`, prediction and the `sample_path`
**computation** all work without it, and an export attempt raises an explicit install hint.

`M0 contains no COMPASS concept branch and is not supported by the sample-level representation
Sankey visualization.`

Details: [`docs/SAMPLE_PATH.md`](docs/SAMPLE_PATH.md) ·
[integration validation](docs/SAMPLE_PATH_INTEGRATION_VALIDATION.md) ·
[feasibility audit](docs/SAMPLE_SANKEY_FEASIBILITY_AUDIT.md) · runnable:
[`examples/sample_path_example.py`](examples/sample_path_example.py).

**中文**：`sample_path()` 一次调用完成输入检查 → 冻结推理 → 真实表示层分数提取 →
**完整** Cox 贡献分解 → 展示节点筛选 → 可导出结果对象；**不需要**手工准备节点表或 JSON。
图分为五列：Gene expression → COMPASS gene score → Granular signature score →
High-level concepts + clinical/PC → Cohort-relative risk。要点：① gene score 列是官方
`extract(with_gene_level=True)` 的**真实 gene-token 标量分数**（共享 `nn.Linear(32→1)`），
132/43 两层与 η 同一次前向，注意力权重**只**用于决定画哪些连接；② 临床/癌种/PC 由**冻结
Cox 设计矩阵**逐列贡献聚合（多列取 signed sum，PC1–PC10 为十项之和，沿用原标准化口径），
真实填补的字段标注 `(imputed)`；③ 完整分解（含未展示 concept 与被隐藏的癌种项）满足
`sum(contribution) == η`，渲染层强制校验；④ 所有连线**等粗等色**、只表示筛选后的模型连接，
**不作流量守恒声明**；节点颜色表示该层自身分数。主要调参：`gene_threshold` / `max_genes` /
`max_signatures` / `predictor_budget` 与颜色范围，**均为展示参数，不是经验证的科学阈值**；
M2/M3 对比请传**同一份** `Style`。PNG/PDF/SVG 由**同一个** Figure 导出，HTML 内嵌**同一份
PNG 字节**并提供精确节点表。
⚠ 语义边界：这是**样本级计算路径归因**，不是因果机制图、不是通路激活图；颜色表示模型取值
高低或贡献方向，**不得**读作 pathway activated / inhibited。

---

## 3. Understanding the output ｜ 如何理解输出

| kind | object | what it is | how to use it |
|---|---|---|---|
| **Prognostic** | `result.risk` | prognostic risk score = **Cox linear predictor** η = Xβ; higher = higher model-estimated hazard | ranking samples *within* a cohort |
| | `result.risk_rank` | rank within the analysed cohort | ordering |
| | `result.risk_group` | `high`/`low` by the **frozen** pan-cancer median cutoff | group comparison inside one cohort |
| | `result.risk_group_relative` | cohort-relative median split, **display only** | used when the frozen cutoff fails to separate a cohort |
| **Representation** | `result.concept_scores` | 43 COMPASS latent coordinates | downstream prioritisation, clustering, association testing |
| | `result.signature_scores` | 132 gene-set-level scores | the more directly interpretable layer |
| **QC** | `result.qc` | coverage at gene / signature / concept level | decide how far to trust the readouts |

Two cautions that matter for interpretation:

* **The risk score is on the log-relative-hazard scale, not a probability.** `result.risk`
  is the Cox linear predictor η = Xβ (exposed also as `result.linear_predictor`); `exp(η)`
  — not `η` — is the relative hazard. Risk scores order samples and are not calibrated
  probabilities. `survival_probability(t)` exists but is a baseline-hazard-derived output
  (see [ADVANCED_USAGE](docs/ADVANCED_USAGE.md)).
* **Risk comparisons across cancer types are not meaningful.** The model contains
  cancer-type terms, so risk distributions differ by cancer type. Compare within a cohort.

**中文**：`result.risk` 是**预后风险分数**，即 **Cox 线性预测子** η = Xβ（别名
`result.linear_predictor`）；`exp(η)` 才是**相对风险比**（relative hazard），不要把 η 直接
称作 "relative risk"。数值越大表示模型估计的风险越高，只在**同一队列内**比较；
`result.risk_group` 用**冻结的泛癌中位切点**分高/低危；当该切点无法区分某个外部队列时，
`result.risk_group_relative` 给出**仅用于展示**的队列内中位切分。
两点必须注意：① **风险分数在对数风险比尺度上、不是概率**，`survival_probability(t)` 只是基线风险派生量；
② **跨癌种比较风险没有意义**（模型含癌种项，各癌种风险分布不同）。

---

## 4. Input requirements ｜ 输入要求

* `expression`: `samples x genes`, **columns are gene symbols**, index is a sample ID.
  Declare the scale with `input_scale` — `"tpm"` (default) or `"log2_tpm1"`.
  Scales are never auto-converted; if the values look inconsistent with the declared
  scale you get a note asking you to check `input_scale`.
* `cancer_type`: required, one TCGA acronym per sample (a single string applies to the
  whole cohort). Accepted values: `src/compass_os/data/cancer_codes.tsv`.
  It is used, from the same value, both as the COMPASS encoder token and as the model's
  cancer-type indicator.
* `clinical` (optional): columns `age`, `sex`, `stage` (stage 1–4). Missing values are
  filled with **frozen training reference values** and listed in the report.
* `survival` (optional): a time column and an event column; common names
  (`time`, `os_time_days`, … / `event`, `os_event`, `status`, …) are auto-detected.

**Direct public inference currently accepts TPM or log2(TPM+1)-scale expression.**
Raw microarray intensity matrices are not accepted directly — they first need
cross-platform harmonisation (see [Validation](#8-validation)).

Missing genes are handled with the frozen TCGA reference median; this never uses
statistics estimated from your own cohort, so a sample's prediction does not depend on
which other samples were submitted alongside it.

**中文**：`expression` 列名必须是**基因 symbol**；表达尺度**必须显式声明**
（`input_scale="tpm"` 默认，或 `"log2_tpm1"`），包**不会自动猜测或转换**，
若数值与声明尺度明显不符会给出提示。`cancer_type` **必填**（TCGA 缩写，
取值见 `src/compass_os/data/cancer_codes.tsv`），同一个值同时用于 COMPASS 编码器与模型癌种项。
可选 `clinical`（`age`/`sex`/`stage`，缺失值用**冻结参考值**填充并在报告中列出）与
`survival`（时间列 + 事件列，常见列名自动识别）。
**v1 仅接受 TPM 或 log2(TPM+1)**；原始微阵列强度矩阵需先自行做跨平台和声化。
缺失基因用**冻结的 TCGA 参考中位数**填补，**不使用你自己队列的统计量**，
因此单样本预测不依赖于同批提交的其他样本。

---

## 5. Advanced usage ｜ 进阶用法

For model selection (`M0`–`M3`), explicit missing-gene strategies, batch/device control,
and direct access to the representation:

```python
rep = compass_os.get_representation(expr, cancer_type, missing_gene_strategy="strict")
pred = compass_os.predict(expr, cancer_type, clinical=clin, model="M0,M1,M2,M3")
rob = compass_os.check_robustness(expr, cancer_type, strategies=("reference", "zero"))

# risk_rank and risk_group are separate guarantees (min_cohort_for_stratification=30):
#   n == 1        -> neither
#   2 <= n < 30   -> risk_rank only
#   n >= 30       -> risk_rank + frozen-cutoff high/low
pred = compass_os.predict(expr, cancer_type, min_cohort_for_stratification=10)
```

Full details, including when to use each model and what each strategy means:
[`docs/ADVANCED_USAGE.md`](docs/ADVANCED_USAGE.md).

Everything below Level 2 — the COMPASS scaler, gene vocabulary, PCA loadings, Cox
coefficients, cancer code table, reference median, signature→concept weights and baseline
hazard — is **read-only frozen assets with no public setter**. To change the model
definition, fork the repository and modify the source.

**中文**：公开接口分两层——**Level 1** 一键式 `analyze()`；**Level 2** 研究型
`predict()` / `get_representation()` / `check_robustness()`，可选模型（M0–M3）、
缺失基因策略、`input_scale`、`batch_size`/`device` 等有限参数，详见
[`docs/ADVANCED_USAGE.md`](docs/ADVANCED_USAGE.md)。
Level 2 以下的全部模型内部（COMPASS scaler、基因词表、PCA、Cox 系数、癌种码表、
参考中位数、signature→concept 权重、基线风险）都是**只读冻结资产，没有 public setter**；
要改模型定义请 fork 源码。

---

## 6. Missing genes and robustness ｜ 缺失基因与稳健性

Real public cohorts rarely contain all 15,672 COMPASS genes. Measured across 38 external
cohorts (before any imputation) the global COMPASS gene coverage was **0.635–0.980
(median 0.967)** and the 916-gene signature coverage **0.803–0.978 (median 0.962)**.

Three strategies exist (`reference` is the default used by `analyze`):

| strategy | behaviour |
|---|---|
| `reference` | fill missing genes with the frozen TCGA reference median |
| `zero` | training-space masking: the frozen scaler is applied normally, then missing genes' normalised values are set to 0 |
| `strict` | raise `MissingGenesError` if any required gene is missing |

`analyze(..., robustness="auto")` runs a second forward pass **only** when the *combined*
QC grade is not `recommended` — i.e. when **either** axis falls below its recommended tier
(`overall = worse(global tier, signature tier)`, see below). The comparison is reported in
plain language.

Coverage QC is always returned, at three levels: gene coverage, per-signature coverage
(`n x 132`) and per-concept input coverage (`n x 43`). Calibrated thresholds live in
`models/qc_config.json`, derived from a stress test — never hard-coded:

| threshold | value | basis |
|---|---|---|
| recommended global coverage | **0.90** | lowest measured coverage with risk Spearman ≥ 0.95 and group agreement ≥ 0.90 |
| warning global coverage | **0.70** | lowest measured coverage with risk Spearman ≥ 0.80 and group agreement ≥ 0.70 |
| recommended signature-gene coverage | **0.90** | same criterion, on the signature-gene axis |
| warning signature-gene coverage | **0.70** | same criterion, on the signature-gene axis |

These are **empirically calibrated descriptive engineering tiers**, not biological or
"safe" thresholds: the underlying curves are smooth with no breakpoint, and the values are
measured coverage levels (no interpolation). The overall QC grade is
**`overall = worse(global tier, signature tier)`**. Per-signature coverage
(`signature_gene_coverage`, *n* × 132) is reported as a **continuous** QC value; no
per-signature numeric threshold is shipped.

Three masking arms were compared (global / signature-only / non-signature-only,
12 cohorts × 20 repeats each):

| coverage on the arm's own axis | 0.95 | 0.90 | 0.80 | 0.70 | 0.60 | 0.50 |
|---|---|---|---|---|---|---|
| **global** — risk Spearman | 0.986 | 0.959 | 0.923 | 0.845 | 0.793 | 0.741 |
| **signature-only** — risk Spearman | 0.983 | 0.964 | 0.924 | 0.849 | — | — |
| **non-signature-only** — risk Spearman | 1.000 | 1.000 | 0.9999 | 0.9999 | 0.9998 | 0.9998 |

Preservation of the signature-related genes accounted for a substantial component of
prognostic robustness: with the 916 signature genes fully retained, removing half of the
remaining 14,756 genes left risk Spearman at 0.9998, group agreement at 1.00 and both
representation layers at a correlation of 1.000, whereas removing the same *fraction* of
the signature genes degraded risk Spearman to 0.849. The non-signature genes are not
without effect — the mean absolute risk shift grew from 0.002 to 0.022 across that range,
and ΔUno C rose to +0.0008 — so the supported statement is that preserving the signature
genes accounts for a substantial component of the robustness, not that the remaining genes
contribute nothing.

`concept_input_coverage` is *input coverage*, **not** a confidence in any biological
interpretation, and must not be renamed to "confidence".

Details and the underlying evidence: [`docs/MISSING_GENES.md`](docs/MISSING_GENES.md).

**中文**：真实公共队列很少包含全部 15,672 个 COMPASS 基因。38 个可还原外部队列
（填补之前实测）全局覆盖 **0.635–0.980（中位 0.967）**，916 个 signature 基因覆盖
**0.803–0.978（中位 0.962）**。
三种策略：`reference`（默认，用冻结 TCGA 参考中位数填补）、`zero`（训练空间掩码：
正常过冻结 scaler 后把缺失基因的标准化值置 0）、`strict`（缺任一必需基因即报错）。
`analyze(..., robustness="auto")` **只有在双轴组合分级不是 recommended 时**才额外跑一次
前向比较（即 global 轴**或** signature 轴任一低于推荐档；`overall = worse(global, Gsig)`），
并用平实语言报告结果。质控阈值放在 `models/qc_config.json`，全部由实测曲线导出、
**不是硬编码**：推荐 0.90、警告 0.70（全局轴与 signature 轴各一套，均为**实测覆盖度水平、
不插值**），属**工程 QC 分层，不是生物学"安全线"**；总分级规则为
**`overall = worse(global tier, signature tier)`**。逐 signature 覆盖
（`signature_gene_coverage`, n×132）以**连续值**报告，v1 **不发布**逐 signature 数值阈值。
三臂掩码对照（各 12 队列 × 20 次重复）见上表：**保留 916 个 signature 基因时，删掉其余
14,756 个基因的一半，风险排序几乎不变（Spearman 0.9998）**；而删除**同样比例**的
signature 基因时降到 0.849。因此可以说
"**保留 signature 相关基因解释了预后稳健性的相当一部分**"，
但**不能**说其余基因毫无作用（mean|Δrisk| 从 0.002 升到 0.022）。
`concept_input_coverage` 只是**输入覆盖度**，**不是**对任何生物学解释的置信度，不得改名为 confidence。

---

## 7. Biological interpretation ｜ 生物学解释

The 43 concepts are **learned latent coordinates** of the COMPASS model, not expression
abundance or pathway-activity scores: every gene-set score is one shared
`nn.Linear(32→1)` projection of the encoder embedding (with `ReLU` disabled, so signs are
unconstrained), and the 43 concepts are a softmax-attention convex combination of the 132
gene-set scores. A dedicated audit found **20/43 concepts reverse-encoded at the model
level**. Concept names are **semantic anchors**.

Practical guidance:

* Use the **132 gene-signature layer** as the primary route to biological interpretation;
  treat the 43 concepts as model-internal coordinates.
* Reported vocabulary: *mechanism-oriented representation*, *candidate biological
  program*, *hypothesis-generating feature*, *risk-associated representation*.
* Do **not** report "activated mechanism", "suppressed pathway" or "causal mechanism".
* Any association between concept scores and the risk score is a **model-linked
  representation association** — the default model uses those very concepts as Cox
  predictors, so this is not independent evidence. Report files label it as such.

See [`docs/CONCEPT_SEMANTICS.md`](docs/CONCEPT_SEMANTICS.md).

**中文**：43 个 concept 是模型的**学到的隐坐标**，不是表达丰度或通路活性分数
（每个 gene-set 分数是同一个 `nn.Linear(32→1)` 投影，`ReLU` 被禁用 ⇒ 符号不受约束；
43 concepts 是 132 个 gene-set 分数的 softmax 注意力凸组合）。专项审计发现
**20/43 个 concept 在模型层面是反向编码**，因此概念名只是**语义锚点**。
实践建议：① 机制解释**以 132 gene-signature 层为主线**，43 concepts 视为模型内部坐标；
② 统一使用 *mechanism-oriented representation / candidate biological program /
hypothesis-generating feature / risk-associated representation* 这类措辞；
③ **不得**写"激活的机制""被抑制的通路""因果机制"；
④ concept 分数与风险的关联属 **model-linked representation association**
（默认模型本身就把这些 concept 当作 Cox 预测项，因此**不是独立证据**），报告文件已如此标注。
详见 [`docs/CONCEPT_SEMANTICS.md`](docs/CONCEPT_SEMANTICS.md)。

---

## 8. Validation ｜ 验证

### Reproducibility ｜ 可复现性

The frozen package reproduces production outputs from the original model-development
pipeline **numerically within a pre-declared tolerance**: 132 signatures
`max|Δ| = 1.1e-16`, 43 concepts `max|Δ| = 2.4e-07`, and M0–M3 risks
`max|Δ| ≤ 1.06e-06` (acceptance criterion 1e-5). The official upstream
`PreTrainer.extract()` and the production `predict()` + capture path agree to
`max|Δ| = 0.0`, i.e. bit-for-bit for that comparison.

**中文**：冻结包**在预设容差内**数值一致地复现原模型开发流程的产物——
132 signatures `max|Δ| = 1.1e-16`、43 concepts `max|Δ| = 2.4e-07`、
M0–M3 风险 `max|Δ| ≤ 1.06e-06`（判据 1e-5）。
上述三项 `max|Δ|` **均非零**，因此属于"容差内一致"，**不是**逐位复现。
只有官方 `PreTrainer.extract()` 与生产 `predict()` 路径的比较实测为 `max|Δ| = 0.0`，
该一项可称逐位一致。

### External input availability ｜ 外部队列输入可达性

Across 38 restorable public cohorts (measured before any imputation):

| | median | range |
|---|---|---|
| global COMPASS gene coverage | **96.7 %** | 63.5 % – 98.0 % |
| signature-gene (916) coverage | **96.2 %** | 80.3 % – 97.8 % |

**中文**：38 个可还原公共队列（填补之前实测）——全局覆盖中位 **96.7 %**（63.5–98.0 %），
916 个 signature 基因覆盖中位 **96.2 %**（80.3–97.8 %）。可见真实公共数据并非总是完整。

### Missing-gene robustness ｜ 缺失基因稳健性

Across the tested masking range, changes in Uno's C remained within approximately ±0.02,
with no consistent deterioration across the evaluated cohorts. Risk *ranking* did degrade
smoothly with coverage (Spearman 0.986 / 0.959 / 0.923 / 0.845 / 0.793 / 0.741 at 95 / 90 /
80 / 70 / 60 / 50 % global coverage), so coverage and ranking stability are reported
together rather than described as "unaffected".

**中文**：在所测试的掩码区间内，**Uno's C 的变化约在 ±0.02 以内，未在各评估队列中表现出一致性下降**。
风险**排序**随覆盖下降而平滑降低（全局覆盖 95/90/80/70/60/50 % 时 Spearman 分别为
0.986 / 0.959 / 0.923 / 0.845 / 0.793 / 0.741），因此覆盖率与排序稳定性必须**同时报告**，
不得写成"不受影响"。

### Signature-related genes ｜ signature 相关基因

**Preservation of signature-related genes accounted for a substantial component of
prognostic robustness.** With all 916 signature-related genes retained, masking 50 % of the
remaining 14,756 genes preserved risk rankings (Spearman ≈ 0.9998, group agreement 1.00,
43/132 representation correlations 1.000). Masking the same *fraction* of the signature
genes themselves reduced risk Spearman to 0.849. The remaining genes are not irrelevant —
the mean absolute risk shift grew from 0.002 to 0.022 over that range — so the supported
statement is the one above, not that non-signature genes have no contribution.

**中文**：**保留 signature 相关基因解释了预后稳健性的相当一部分。**
916 个 signature 基因全部保留时，掩掉其余 14,756 个基因的 50 %，风险排序基本不变
（Spearman ≈ 0.9998、分层一致率 1.00、43/132 表征相关性 1.000）；
而掩掉**同样比例**的 signature 基因时，Spearman 降到 0.849。
其余基因并非毫无作用（mean|Δrisk| 从 0.002 升到 0.022），
因此只能下上面这句结论，**不能**说 non-signature 基因没有贡献。

### Real platform missing patterns and QC thresholds ｜ 真实平台缺失模式与 QC 阈值

Replaying the real pre-imputation platform masks of 38 cohorts (training-space masking)
left risk ranking unchanged (Spearman 1.000, group agreement 1.000) while representation
readouts shifted by ≈6 %. The corresponding `reference` arm is a **production-pipeline
self-consistency control** — it returns the input unchanged by construction, because the
harmonised cache is itself reference-imputed — and is not independent robustness evidence.
Coverage QC was predictive: the correlation between per-feature input coverage and
perturbation error was negative for 43/43 concepts and 132/132 signatures (median
Spearman −0.68 and −0.60).

Calibrated QC tiers (see section 6) and the per-level audit are in
`validation/results/qc_threshold_audit.tsv`; design and cost accounting in
[`docs/VALIDATION_DESIGN.md`](docs/VALIDATION_DESIGN.md).

**中文**：用 38 个队列**真实的填补前平台缺失模式**回放（训练空间掩码），
风险排序不变（Spearman 1.000、分层一致率 1.000），而表征读数偏移约 6 %。
其中 `reference` 那一臂是**生产流程的自洽性对照**（和声化缓存本身就是参考中位数填补，
因此按构造返回原输入），**不得**当作独立的稳健性证据。
覆盖度 QC 具备预测力：逐特征输入覆盖度与扰动误差的相关在 **43/43 concepts 与
132/132 signatures 上全部为负**（中位 Spearman −0.68 与 −0.60）。
QC 分层与逐 level 判据见 `validation/results/qc_threshold_audit.tsv`，
设计与算力口径见 [`docs/VALIDATION_DESIGN.md`](docs/VALIDATION_DESIGN.md)。

---

## 9. Limitations ｜ 已知限制

* Trained on TCGA; external performance is dataset-dependent. Discrimination transfers
  better than absolute-risk calibration.
* The frozen pan-cancer median cutoff often fails to separate an external cohort (risk
  levels shift by cancer type). `analyze` then reports a cohort-relative split for
  display only — never for cross-cohort comparison.
* Single samples are not stratified.
* Microarray input requires external harmonisation; an `input_scale="microarray"` mode does
  not exist in v1.
* Robustness to missing genes is characterised only within the tested range; do not
  describe the model as invariant to missing genes.
* Signature-coverage thresholds are not yet calibrated.

**中文**：① 模型在 TCGA 上训练，外部表现依赖数据集，**判别度的可迁移性好于绝对风险校准**；
② 冻结的泛癌中位切点常无法区分外部队列（风险水平随癌种平移），此时 `analyze` 给出
**仅用于展示**的队列内相对切分，**不得**用于跨队列比较；
③ 单样本不做分层；
④ **微阵列输入需自行和声化**，v1 没有 `input_scale="microarray"`；
⑤ 缺失稳健性只在已测试区间内成立，**不得**称模型对缺失基因不变；
⑥ 逐 signature 覆盖阈值尚未校准（以连续值报告）。

---

## 10. Installation ｜ 安装

```bash
git clone https://github.com/RENSI203/COMPASS-OS.git
cd COMPASS-OS
pip install .                       # standard install (recommended)

# developers only:
pip install -e .                    # editable install
```

All three installation modes (standard install, editable install, wheel) work identically.
**No environment variable is required**: model assets and the vendored COMPASS package are shipped as
package data inside `compass_os/assets/`, so a wheel install resolves everything from
`site-packages` without the Git repository. `COMPASS_OS_ROOT` remains available as an
optional **developer override** (point it at a directory containing `models/`, at an
`assets/` directory, or at `src/compass_os`).

Requires Python ≥ 3.10, PyTorch, pandas, numpy, scikit-learn, scikit-survival, matplotlib.
Model assets (~12 MB) and the upstream COMPASS package are shipped as **package data** in
`compass_os/assets/`; the vendored copy (`assets/third_party/compass/`, MIT, see
`assets/third_party/COMPASS_LICENSE`) is the default execution path. An installed
`immuno-compass==2.5.3` is only a fallback when the vendored copy is absent (optional extra
`compass-os[upstream]`).

Reproduce the checks:

```bash
python tests/run_tests.py      # no pytest needed; or: pytest tests/
```

`ASSET_MANIFEST.tsv` records every shipped asset (source path, size, SHA-256);
`assets/models/model_manifest.json` records model composition and the default model.

**中文**：依赖 Python ≥ 3.10、PyTorch、pandas、numpy、scikit-learn、scikit-survival、matplotlib。
模型资产与上游 COMPASS 均以 **package data** 形式放在 `compass_os/assets/` 内并随 wheel 分发
（`assets/models/` 约 12 MB；`assets/third_party/compass/` 为 vendored 上游副本，MIT，
见 `assets/third_party/COMPASS_LICENSE`），因此 **安装后不需要 Git 仓库、也不需要任何
环境变量**。vendored 副本是**默认执行路径**；仅当其缺失时才回退到已安装的
`immuno-compass==2.5.3`（可选 extra `compass-os[upstream]`）。
`COMPASS_OS_ROOT` 仅作开发者覆盖用。自检命令：`python tests/run_tests.py`（无需 pytest）。
⚠ v1.0.1 修复：此前模型通过上游 `compass.loadcompass()` 加载，其末尾
`if file.startswith(tempfile.gettempdir()): os.remove(file)` 会把**位于临时目录下**的
checkpoint 在加载后静默删除——安装在 `/tmp` 下的环境（容器 / CI）第一次预测后即失效。
现改为直接 `torch.load` 本地冻结资产，数值行为不变。
`ASSET_MANIFEST.tsv` 记录全部随包资产的来源/大小/SHA-256，
`assets/models/model_manifest.json` 记录模型组成与默认模型。

---

## 11. Citation ｜ 引用

See [`CITATION.cff`](CITATION.cff) (software authors: Jiahao Ren, Junyi Xin). If you use
this software, please cite it together with the upstream COMPASS publication.

**Upstream COMPASS model**

> Shen, W., Moon, I., Nguyen, T.H. et al. Generalizable AI predicts immunotherapy outcomes
> across cancers and treatments. *Nat Med* **32**, 3010–3022 (2026).
> https://doi.org/10.1038/s41591-026-04502-7

```bibtex
@article{shen2026compass,
  title   = {Generalizable AI predicts immunotherapy outcomes across cancers and treatments},
  author  = {Shen, W. and Moon, I. and Nguyen, T.H. and others},
  journal = {Nature Medicine},
  volume  = {32},
  pages   = {3010--3022},
  year    = {2026},
  doi     = {10.1038/s41591-026-04502-7}
}
```

**COMPASS-OS**

```bibtex
@software{compass_os_2026,
  title     = {COMPASS-OS: pan-cancer survival prediction and mechanism-oriented
               transcriptomic representation framework},
  author    = {Ren, Jiahao and Xin, Junyi},
  year      = {2026},
  version   = {1.0.1},
  license   = {MIT},
  url       = {https://github.com/RENSI203/COMPASS-OS}
}
```

**Contact ｜ 联系方式**

Questions, bug reports and collaboration enquiries:
**Jiahao Ren** — rjh2623826975@stu.njmu.edu.cn
(please also use the [issue tracker](https://github.com/RENSI203/COMPASS-OS/issues) for
reproducible bug reports)

**中文**：引用信息见 [`CITATION.cff`](CITATION.cff)（软件作者：Jiahao Ren, Junyi Xin）。
使用本软件时请**同时引用上游 COMPASS 论文**：

> Shen, W., Moon, I., Nguyen, T.H. et al. Generalizable AI predicts immunotherapy outcomes
> across cancers and treatments. *Nat Med* **32**, 3010–3022 (2026).
> https://doi.org/10.1038/s41591-026-04502-7

联系方式：**Jiahao Ren** — rjh2623826975@stu.njmu.edu.cn；
可复现的 bug 报告也欢迎走 [issue tracker](https://github.com/RENSI203/COMPASS-OS/issues)。

---

## 12. License and attribution ｜ 许可与归属

Code: [`LICENSE`](LICENSE) (MIT).

**COMPASS-OS incorporates vendored components from COMPASS under its upstream license.**
The upstream COMPASS source is redistributed unmodified in `third_party/compass/` together
with its MIT licence text at [`third_party/COMPASS_LICENSE`](third_party/COMPASS_LICENSE);
the same attribution is restated in [`NOTICE`](NOTICE);
the pretrained checkpoint (`models/pretrainer.pt`) is likewise distributed under that
upstream licence. Those components are **not** original work of the COMPASS-OS authors —
please cite the upstream COMPASS publication as well.

> Copyright (c) 2026 RENSI203 (MIT). Software citation authors: Jiahao Ren, Junyi Xin
> (see [`CITATION.cff`](CITATION.cff)).

**中文**：本项目代码采用 **MIT** 许可（见 [`LICENSE`](LICENSE)）。
**COMPASS-OS 以 vendored 形式包含来自 COMPASS 的组件，沿用其上游许可**：
上游源码原样置于 `third_party/compass/`，其 MIT 许可证文本见
[`third_party/COMPASS_LICENSE`](third_party/COMPASS_LICENSE)；预训练权重
（`models/pretrainer.pt`）同样按上游许可分发。这些组件**不是** COMPASS-OS 作者的原创工作，
请一并引用上游 COMPASS 论文。
版权主体：2026 RENSI203（MIT）；软件引用作者：Jiahao Ren, Junyi Xin。
