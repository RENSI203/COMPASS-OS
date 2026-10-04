# Sample-level computation path ｜ 样本级计算路径归因

<div align="center"><em>sample-specific computational attribution / representation-flow
visualization</em></div>

> **语义边界（务必遵守）**
> 这是**样本级计算路径归因**：展示"对这一个样本，冻结模型图里各层分别取什么值、
> 各 predictor 贡献了多少到 Cox 线性预测子 η"。
> 它**不是**因果机制图、**不是**通路激活图、**不是**生物因果网络、
> 也**不是**对最终风险的完整因果解释。
> 节点颜色表示**该层自身的真实分数**（临床 / PC 节点表示 signed Cox contribution，
> 用独立图例）；所有连线**等粗等色**，只表示筛选后的模型连接关系，
> **不编码权重或流量**，本图**不作流量守恒声明**。

---

## 1. 一键调用

```python
from compass_os import sample_path

res = sample_path(expr, cancer_type, "TCGA-19-1787", model="M3", clinical=clinical)
print(res.summary())
```

一次调用自动完成：输入检查与既有预处理 → 冻结模型推理 → 真实表示层分数提取 →
完整 Cox contribution 分解 → 展示节点与连接筛选 → 返回可导出结果对象。
**普通用户不需要准备 `SamplePathData`、节点表或 JSON**（那是内部数据接口 / 调试导出）。

`expr` / `cancer_type` / `clinical` / `input_scale` / `missing_gene_strategy` 的语义与
`compass_os.predict` 完全一致；`sample_id` 必须出现在 `expr.index` 中。

---

## 2. 图结构

```
Gene expression → COMPASS gene score → Granular signature score
    → High-level concepts + clinical / PC → Cohort-relative risk
```

| 列 | 内容 | 节点颜色表示 |
|---|---|---|
| 1 | Gene expression | 该基因的真实表达值（单位见 `expression_units`） |
| 2 | COMPASS gene score | **真实 gene-token 标量分数** |
| 3 | Granular signature score | 真实的 132 signature 表示 |
| 4 | High-level concepts + clinical / PC | concept = 真实 concept 分数；临床 / PC = signed Cox contribution（独立图例） |
| 5 | Cohort-relative risk | 连续风险柱（cohort percentile 轴） |

* **所有节点为圆形**；第 1、2 列同尺寸且最小，第 3 列更大，第 4 列最大。
* **所有连线等粗、等色、较细、浅灰**，只表示所展示的模型连接关系。
* 主图**不展示** Other genes / signatures / concepts、offset、residual 或其他记账节点；
  这些项的数值仍完整保留在 `res.decomposition`。
* 实际被冻结填补的临床节点标注 **`(imputed)`**。
* 完整名称与精确数值见 HTML 内嵌的**可展开节点表**，以及 `*.nodes.tsv`。

---

## 3. 分数来源（全部真实）

| 数据 | 来源 | 说明 |
|---|---|---|
| **gene-token 标量分数** | 官方 `PreTrainer.extract(..., with_gene_level=True)` 的 `dfg` | = `genesetprojector.geneset_scorer(encoder_output)[:, 2:]`，对每个基因位置施加**共享的** `nn.Linear(32→1)`；形状 `(n_samples, 15672)`；源码 `assets/third_party/compass/model/tune.py:289-291` |
| 132 signature | 主 API `predict()` **同一次前向**的 `signature_scores` | 与 η 严格同源 |
| 43 concept | 主 API `predict()` **同一次前向**的 `concept_scores` | 与 η 严格同源 |
| Age / Sex / Stage / Cancer type / PC1–PC10 | 冻结 Cox **设计矩阵逐列贡献**聚合 | 多编码列取 **signed sum**；PC1–PC10 为十个冻结设计列贡献之和，沿用原标准化口径 |
| η / rank / percentile / cutoff | 主 API 直接输出 | 图**不重算**排名 |
| 连接拓扑 | 冻结注意力（gene→gene-set、gene-set→concept） | **只**用于决定画哪些连接，**不**作为节点分数 |

**禁止替代**：attention / projector weight、TPM、gene→signature 的风险分配量、
旧 HTML 的 upstream contribution 都**不是** gene score。

`decomposition` 为**完整未截断**分解，含未展示 concept、被隐藏的 Cancer type 及必要的
algebraic offset，满足：

```
sum(all signed contributions) == η
```

渲染层会强制校验该恒等式（不满足即拒绝出图）。

---

## 4. 展示预算与筛选（规格冻结）

```python
from compass_os.sample_path_render import Style

style = Style(
    gene_threshold=0.25,        # 作用于**真实 gene score 尺度**；仅展示参数
    max_genes=32,               # ≤ 50
    max_signatures=26,          # ≤ 50
    predictor_budget=16,        # ≤ 16
    expression_limits=(0, 1704.3),   # 颜色范围（超出仅饱和，不改数据）
    gene_score_limits=(-2, 2),
    signature_limits=(-1, 1),
    concept_limits=(-1, 1),
    contribution_limits=(-1.5, 1.5),
)
m2 = sample_path(expr, ct, sid, model="M2", clinical=clinical, style=style)
m3 = sample_path(expr, ct, sid, model="M3", clinical=clinical, style=style)
```

也可以直接用参数覆盖：

```python
sample_path(..., gene_threshold=0.5, max_genes=20, max_signatures=18, predictor_budget=14)
```

**M2/M3 对比请传同一份 `Style`**，以保证画布、节点尺寸、列位置、字体规则与颜色尺度一致。
颜色范围是**共享的固定尺度**，不会对每张图的选中节点重新 min-max 或 z-score。

### High-level 预算表

| 模型 / 队列 | 实际展示 |
|---|---|
| M1 | 16 concepts |
| M2，单癌种 | 13 concepts + Age + Sex + Stage |
| M2，多癌种 | 12 concepts + Age + Sex + Stage + Cancer type |
| M3，单癌种 | 12 concepts + Age + Sex + Stage + 一个 `PC1–PC10` 节点 |
| M3，多癌种 | 11 concepts + Age + Sex + Stage + Cancer type + 一个 `PC1–PC10` 节点 |

* 单癌种隐藏 Cancer type **仅影响展示**；该项 contribution 仍在完整模型计算与
  `decomposition` 中。
* **真筛选**：按阈值与上限真正删除节点与连接，**不是**把多余 label 设空。
* Gene expression 与 Gene score **两列展示同一组基因**，逐名对应。
* **没有节点通过阈值时**，明确显示空层提示，**不会**偷偷降低阈值补齐节点。
* 默认阈值与颜色范围是**可配置展示参数**，**不是**经过临床验证的科学阈值。

### 筛选顺序（与渲染层一致）

1. 第四列：按 `|signed Cox contribution|` 选 top concept（名额 = `predictor_budget` − 必需临床/PC 节点）。
2. 第三列：限定为所展示 concept 的**上游** signature，按 `|signature score|` 取 top `max_signatures`。
3. 第二列：限定为所展示 signature 的**上游**基因，按 `|gene score|` 且 `>= gene_threshold`
   取 top `max_genes`；每条 signature 默认最多保留 4 条基因连接
   （`max_gene_edges_per_signature=None` 可保留全部）。
4. 被连接裁剪后**断开**的基因同时从第 1、2 列删除（不添加 Other placeholder）。

---

## 5. 风险柱

* 竖直、上**红**下**黄**、足够长；高风险在上、低风险在下（cohort percentile 轴 0–100 %）。
* 第四列所有连线收束到三角形附近。
* 三角形指向该样本**主 API 的真实 percentile** 位置。
* 标注 `Risk η = ... / Rank = ... / N / Percentile = ...%`。
* **cutoff 虚线**按其**真实经验分位**定位（`100 × (N[η<cutoff] + 0.5·N[η=cutoff]) / N`），
  **不会**一律画在 50 %。

`cutoff` 参数：

| 取值 | 含义 |
|---|---|
| `"median"`（默认） | **当前 cohort** 的 η 中位数，**仅可视化分界**，**不是**冻结的 validated cutoff |
| 数值 | 研究者给定的风险阈值（同样按真实队列分布定位） |

`res.cutoff_is_frozen_validated` 恒为 `False`：本功能**从不**把冻结预后分组 cutoff
当作展示分界。

**主 API 的 risk 尺度**：`compass_os.predict().risk` 就是 **Cox 线性预测子 η**
（不是 `exp(η)`、不是概率、不是另一种 risk scale），因此图中标注的 η 使用的正是 η。
`rank` / `percentile` 沿用主 API 既定结果：

```
rank       = 1 + #{cohort η > η_i}
percentile = 100 × (P(η < η_i) + 0.5 × P(η == η_i))      # 高值 = 高风险
```

---

## 6. 结果对象

```python
res.sample_id, res.model, res.cancer_type
res.risk              # 主 API 的 η
res.rank, res.n_samples, res.percentile
res.cutoff, res.cutoff_mode, res.cutoff_is_frozen_validated
res.decomposition     # DataFrame: feature / value / beta / contribution（完整未截断）
res.genes             # 展示的第 1、2 列基因（gene / expression / gene_score）
res.signatures        # 展示的 132 层
res.concepts          # 展示的 43 层（concept / score / contribution）
res.auxiliaries       # 展示的临床 / PC / Cancer type 节点
res.notes, res.semantics, res.summary()
```

`res.genes` / `res.signatures` / `res.concepts` / `res.auxiliaries` 由渲染层的
`select_nodes` 决定，因此**与图永远一致**。

---

## 7. 导出

```python
res.save_html("sample_path.html")     # 可独立调用，无需先保存 PNG
res.save_static("sample_path.png")    # .png / .pdf / .svg
res.save("out/M3_low")                # 一次导出全部
# -> M3_low.{png,pdf,svg,html,nodes.tsv,selection.json,decomposition.tsv}

rp = res.render()                     # RenderedPath：.figure / .selection / .cutoff
rp.close()
```

* PNG / PDF / SVG 由**同一个 Figure** 导出；**HTML 内嵌同一份 PNG 字节**并提供可展开的
  精确节点表，三种格式**不存在两套布局**。
* **不再使用**旧 Plotly Sankey 与 matplotlib fallback 两种外观；导出失败**不会**静默切换
  成另一种图形。
* HTML 主图是**静态一致版本**（没有 hover / 缩放）；精确数值请用内嵌节点表或
  `*.nodes.tsv`。

### 迁移说明（v1.0.x 之后未发布的绘图 API → v1.1.0）

| 旧写法 | 新写法 |
|---|---|
| `res.figure()` 返回 Plotly `Figure` | 返回 `RenderedPath`（`.figure` 为 Matplotlib Figure） |
| `fig.write_html(...)` / `res.save_html(p, include_plotlyjs=...)` | `res.save_html(p)`（无 plotly 参数） |
| `res.save_static(png)`（kaleido 或 fallback） | `res.save_static(png/pdf/svg)`（同一 Figure） |
| `top_genes` / `top_signatures` / `max_high_level_nodes` | `max_genes` / `max_signatures` / `predictor_budget`（旧名仍兼容） |
| `top_concepts` | 已弃用：concept 数 = `predictor_budget` − 必需临床/PC 节点 |
| `res.gene_links` / `residual_*` / `centering_*` | 已移除（属旧"流量守恒"视图）；完整数值见 `res.decomposition` |

---

## 8. 依赖

* 绘图使用 **matplotlib**，自 v1.0 起即为**核心依赖**（`AnalysisResult.save_report` 出图需要），
  本功能**没有引入新增硬依赖**：不再需要 Plotly / Kaleido / CDN。
* 渲染层**惰性导入** matplotlib：未安装时 `import compass_os`、核心预测与
  `sample_path()` 的**计算**部分照常可用；调用导出时给出明确安装提示
  （`pip install "compass-os[plot]"` 或 `pip install matplotlib`）。

---

## 9. 示例

`examples/sample_path_example.py` 在真实队列（GSE39582，COAD，n=573）上按 **M3** percentile
最接近 10 / 50 / 90 程序化选择 3 个样本，对**同一批样本**分别绘制 M2 与 M3（共用一份
`Style`），输出 6 组 PNG/PDF/SVG/HTML + 完整 decomposition + selection 审计：

```bash
MPLCONFIGDIR=/tmp/mplcfg python examples/sample_path_example.py \
    --cohort GSE39582 --out docs/figures/sample_path_integration
```

---

## 10. 如何解读（与不可解读）

**可以**：
* 该样本的表示层取值高低（第 1–3 列颜色）；
* 哪些 concept / 临床 / PC 项对该样本的 η 贡献最大、方向如何（第 4 列）；
* 该样本在**当前队列**风险分布中的相对位置（第 5 列）。

**不可以**：
* 把颜色读成 pathway **activated / inhibited**；
* 把连线读成流量、权重或因果传递；
* 把本图当作因果机制图、通路激活图或风险的完整因果解释；
* 跨队列比较 percentile（它是**队列内**相对量）；
* 把 cohort median cutoff 当作冻结的预后分组阈值；
* 把默认阈值 / 颜色范围当作经验证的科学阈值。
