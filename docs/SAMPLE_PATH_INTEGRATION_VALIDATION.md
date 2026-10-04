# Sample path 集成验收报告 ｜ integration validation

**对象**：把 `COMPASS_sample_path_redesign.zip` 的圆形分层路径图**正式集成进 COMPASS-OS 包**。
**版本**：**`1.1.0`**（候选，未发布）｜ **分支**：`feature/sample-sankey`（本地，未 merge / tag / push）。

> **版本号说明**：本功能起初在工作树中记为 `1.1.0`，集成渲染层时曾一度改为 `1.2.0`
> （理由是把 Plotly `Figure` 换成 Matplotlib `RenderedPath` 属破坏性变更）。
> 但经核对**真实发布历史**——远端 tag 仅 `v1.0.0` / `v1.0.1`，GitHub Release 亦仅 2 个——
> **`1.1.0` 从未发布**，因此不存在需要与之区分兼容性的已发布版本。
> 继续用 `1.2.0` 会虚构一个不存在的 `1.1.0` 发布并误导迁移说明；
> 故候选版本**回退为 `1.1.0`**。详见 `docs/audit/00_CANDIDATE_BASELINE.md`。

> 附件里的 M2/M3 图与 JSON 是**明确标注的合成版式示例**，本报告**不引用**它们作为验证结果。
> 下列全部结论来自**真实队列的真实推理**。

---

## 1. 修改文件与公开 API

### 新增

| 文件 | 作用 |
|---|---|
| `src/compass_os/sample_path_render.py` | 渲染层（附件的 `sankey.py` 按包内架构整理）。**只画图**：不推理、不改冻结模型、不重算风险。`numpy` / `matplotlib` 均为**函数内惰性导入** |
| `src/compass_os/_render_bridge.py` | 内部桥接：`SamplePathResult` → 渲染层 `SamplePathData`。**普通用户不接触** |
| `tests/test_sample_path_integration.py` | 生产集成测试（规格 §9 的 12 项） |
| `tests/verify_sample_path_renderer.py` | 附件 12 项 renderer checks，**原样保留**（改为从包内导入） |
| `tools/build_review_grid.py` | 6 张图拼 2×3 人工验收总览（只拼贴，不重绘） |

### 修改

| 文件 | 改动 |
|---|---|
| `src/compass_os/sample_path.py` | 展示节点筛选改由渲染层 `select_nodes` **统一决定**；`result` 改为统一导出接口；移除旧 Plotly `_build_figure` 与 matplotlib 回退 `_save_matplotlib` |
| `src/compass_os/__init__.py` | 版本 → `1.1.0`（公共导出不变） |
| `pyproject.toml` / `CITATION.cff` | 软件版本 → `1.1.0`（`CITATION.cff` 的 `cff-version` 是 CFF **schema** 版本，保持 `1.2.0` 不动）；`plot` extra 注释改为"matplotlib 已是核心依赖，本功能无新增硬依赖" |
| `examples/sample_path_example.py` | 改写为真实队列 M2/M3 端到端示例（**无硬编码绝对路径**，改用 `--expr-dir`） |
| `README.md` §2b / `docs/SAMPLE_PATH.md` | 按集成后的 API 与视觉规格重写 |
| `tests/test_sample_path.py` / `tests/run_tests.py` | 适配新语义与模块注册 |

### 公开 API

```python
from compass_os import sample_path                       # 不变

res = sample_path(expr, cancer_type, sample_id, model="M3", clinical=clinical)   # 不变
res.summary(); res.decomposition                          # 不变
res.save_html("x.html")      # 可独立调用（无需先保存 PNG）
res.save_static("x.pdf")     # .png / .pdf / .svg（同一 Figure）
res.save("out/M3")           # 全部格式 + nodes.tsv + selection.json + decomposition.tsv
res.render()                 # RenderedPath：.figure / .selection / .cutoff / .cutoff_percentile

from compass_os.sample_path_render import Style          # 新增：视觉规格
```

**用户不需要准备 `SamplePathData` / 节点表 / JSON** —— 这些只是内部接口与调试导出。

### 迁移（v1.0.x 之后未发布的绘图 API → v1.1.0）

| 旧 | 新 |
|---|---|
| `res.figure()` → Plotly `Figure` | `res.figure()` → `RenderedPath`（`.figure` 为 Matplotlib Figure） |
| `save_html(p, include_plotlyjs=...)` | `save_html(p)`（不再有 plotly 参数） |
| 静态图走 kaleido / matplotlib 回退（**两种外观**） | `save_static(png/pdf/svg)`（**同一 Figure**） |
| `top_genes` / `top_signatures` / `max_high_level_nodes` | `max_genes` / `max_signatures` / `predictor_budget`（**旧名仍兼容**） |
| `top_concepts` | 弃用：concept 数 = `predictor_budget` − 必需临床/PC 节点（给 note 提示） |
| `res.gene_links` / `residual_*` / `centering_*` | 移除（属旧"流量守恒"视图）；完整数值见 `res.decomposition` |

---

## 2. gene score 的真实提取位置、定义与张量形状

| 项 | 内容 |
|---|---|
| 代码位置 | `src/compass_os/sample_path.py::_extract_gene_level` → 调用 `PreTrainer.extract(..., with_gene_level=True)` |
| 上游定义 | `assets/third_party/compass/model/tune.py:289-291`：`gene_level_proj = genesetprojector.geneset_scorer(encoding)[:, 2:]` |
| 数学定义 | 对每个基因位置施加**共享的** `nn.Linear(32→1)`：`gene_score_g = w · e_g + b`，其中 `e_g ∈ R^32` 为基因编码，`w`/`b` 为冻结参数 |
| 张量形状 | `(n_samples, 15672)`（完整 COMPASS 词表；本队列 15,672 基因） |
| 实测分布 | `\|gene_score\|` 中位 ≈ 0.40、p90 ≈ 0.94、max ≈ 1.98（因此 `gene_threshold=0.25` 是**该尺度上**的展示阈值） |
| 官方 API | `extract()` 返回 `(dfg, dfgs, dfct)`，本功能只用其中的 `dfg` |

**未使用**任何替代量：attention/projector weight、TPM、gene→signature 风险分配量、旧 HTML 的
upstream contribution 都**没有**被当作 gene score。注意力权重**仅**用于决定画哪些连接。

---

## 3. signature / concept 与原表示层的一致性

**设计决定**：132 signature 与 43 concept 取主 API `predict()` **同一次前向**的输出
（`pred.signature_scores` / `pred.concept_scores`），因此与 η **严格同源**；
只有 gene 层需要额外一次 `with_gene_level=True` 前向（API 不返回 gene 层）。

实测（`test_score_extraction_does_not_change_existing_outputs`）：

| 比较 | 结果 |
|---|---|
| 独立提取的 132 signature vs 主 API `signature_scores` | `max\|Δ\| < 1e-5` |
| 独立提取的 43 concept vs 主 API `concept_scores` | `max\|Δ\| < 1e-5` |
| gene 层形状 | `(1, 15672)`，全部有限 |
| 同一输入重复预测的 risk | 实测差 **0.000e+00** → **逐位一致**（分数提取不改变风险） |

节点颜色 = 该层**真实分数**；concept 只按 `|signed Cox contribution|` **选**节点，
**没有**把颜色换成 contribution。

---

## 4. clinical / Cancer type / PC 聚合规则

全部来自**冻结 Cox 设计矩阵的逐列 signed contribution**（含原标准化口径）：

| 节点 | 规则 |
|---|---|
| `Age` | 该列 contribution = `β_Age · (age − mean)/scale`（Age 在 `scaled_features` 中） |
| `Sex` | 对应编码列的 signed sum（**不**参与标准化） |
| `Stage` | 同上（多列时求和） |
| `Cancer type` | **全部 32 个 `CT_*` 编码列** contribution 的 signed sum；单癌种队列**隐藏**但仍在完整分解中 |
| `PC1–PC10` | **一个**节点，值 = 十个 `PC1..PC10` 冻结设计列 contribution 的 **signed sum**（用标准化后的设计值 × 冻结 β，**不是**原始 PC × 标准化系数） |

* 渲染层**强制** `score == contribution`（临床/PC/cancer 节点），因此颜色与数值不可能不一致。
* `PC1–PC10` 节点名使用 **en dash**（与规格一致），十个分列的 contribution 完整保留在
  `nodes.tsv` 的 `details` 中。
* **imputed 标注来自真实逐样本、逐字段记录**：`_imputed_fields()` 检查该字段是否未提供 /
  整列为空 / 该样本为 NaN，**不根据最终数值猜测**。本队列无临床数据 ⇒ 三项均为
  `(imputed)`，`value` 显示真实冻结填补值 60.0 / 0.0 / 2.0。

---

## 5. η 与完整分解一致性

* `compass_os.predict().risk` **经核对即 Cox 线性预测子 η**（不是 `exp(η)`、不是概率、
  不是另一种 risk scale），因此图中标注的 η 用的就是 η。
* `rank` / `percentile` 直接沿用主 API 结果，图**不重算**排名：
  `rank = 1 + #{η > η_i}`，`percentile = 100 × (P(η<η_i) + 0.5·P(η=η_i))`（高值 = 高风险）。
* 完整分解（**含未展示 concept、被隐藏的 Cancer type、必要 offset**）满足 `sum == η`；
  渲染层 `_check()` 强制校验（不满足即拒绝出图）。

实测（本轮 6 张真实图）：

| 图 | 分解行数 | `sum(contribution)` vs η |
|---|---|---|
| low/M2 | 78 | Δ = 0.00e+00 |
| low/M3 | 88 | Δ = 0.00e+00 |
| mid/M2 | 78 | Δ = 0.00e+00 |
| mid/M3 | 88 | Δ = 0.00e+00 |
| high/M2 | 78 | Δ = 0.00e+00 |
| high/M3 | 88 | Δ = 0.00e+00 |

图中 `risk` 与主 API **逐值相等**（6/6）；`rank` 与主 API **逐值相等**（6/6）。

> **接线中修正的一处真实问题**：分解最初用适配器自己那次单样本前向的 concept 分数构建，
> 结果 `M3_low` 的 `sum` 与 η 相差 **3.8e-07**（η=0.058，相对容差最紧）。渲染层的
> `sum==η` 守卫**正确报错拒绝出图**。改为使用**主 API 的 `concept_scores`** 后残差降到 0。
> 这正是"不许静默降级"的价值所在。

---

## 6. 真实六张图与 review grid

**cohort**：GSE39582（COAD，单癌种，**n = 573**，和声化缓存，`input_scale="log2_tpm1"`）。
**选样规则**：按 **M3** percentile 最接近 **10 / 50 / 90** 程序化选择（未人工挑样本）；
M2/M3 使用**完全相同的三个样本**与**同一份 `Style`**。

| level | sample_id | M2 η | M2 rank (pct) | M3 η | M3 rank (pct) |
|---|---|---|---|---|---|
| low | `GSE39582::GSM972151` | 1.2842 | 254/573 (55.8 %) | −0.0307 | **516/573 (10.0 %)** |
| mid | `GSE39582::GSM972318` | 1.2183 | 318/573 (44.6 %) | 0.2623 | **287/573 (50.0 %)** |
| high | `GSE39582::GSM972096` | 1.6863 | 15/573 (97.5 %) | 0.5642 | **58/573 (90.0 %)** |

> ⚠️ **值得注意的真实结果**：因为本轮按规格用 **M3** percentile 选样，同一批样本在 **M2**
> 下的 percentile 差别很大（low 10.0 % → 55.8 %；high 90.0 % → 97.5 %）。
> 这不是错误，而是"同一患者在不同模型下队列相对位置不同"的真实体现；
> 换样本会破坏规格，因此**未**替换。

**共享 Style（按真实分布设定，非每图 min-max）**：
`expression_limits=(0, 1704.3)`（队列 TPM p99.5；默认 `(0,20)` 会全部饱和）、
`gene_score_limits=(−2, 2)`、`signature_limits=concept_limits=(−1, 1)`、
`contribution_limits=(−1.5, 1.5)`（M3 的 PC 聚合可达 1.50）。

### 产物路径

```
docs/figures/sample_path_integration/
├── low/{M2,M3}.{png,pdf,svg,html,nodes.tsv,selection.json,decomposition.tsv}
├── mid/{M2,M3}.{…}
├── high/{M2,M3}.{…}
├── shared_style.json
├── REAL_M2_M3_selection_audit.tsv          ← cohort/sample/risk/rank/percentile/cutoff/imputation/PC
└── sample_path_redesign_review_grid.png    ← 2×3 总览（由 6 张原 PNG 拼贴，未重绘）
```

可直接查看的主图：`low/M2.png`、`low/M3.png`、`mid/M2.png`、`mid/M3.png`、`high/M2.png`、
`high/M3.png`；总览 `sample_path_redesign_review_grid.png`；交互/表格版为同名 `.html`。

### 实际渲染节点预算

| 图 | gene 两列 | signatures | predictors |
|---|---|---|---|
| low/M2 · low/M3 | 28 · 28 | 26 | 16 |
| mid/M2 · mid/M3 | 25 · 25 | 26 | 16 |
| high/M2 · high/M3 | 28 · 29 | 26 | 16 |

两列基因**逐名对应**；predictor 恒 16（M2：13 concepts + Age/Sex/Stage；
M3：12 concepts + Age/Sex/Stage + `PC1–PC10`）。

---

## 7. 测试数量与 PASS/FAIL

| 套件 | 数量 | 结果 |
|---|---|---|
| 包内全量回归（`tests/run_tests.py`） | **132** | **PASS 132 / FAIL 0 / SKIP 0** |
| ├ 其中 `test_sample_path` | 27 | PASS |
| ├ 其中 `test_sample_path_integration`（规格 §9 的 12 项） | 6 | PASS |
| └ 其中既有回归（assets/api/analysis/qc/packaging/reproducibility…） | 99 | PASS |
| 渲染层自检（`tests/verify_sample_path_renderer.py`，附件原样保留） | **12** | **12/12 OK** |
| 冻结资产校验（`tools/build_assets.py --verify`） | 13 | **异常 0** |
| release audit（`tools/release_audit.py`） | 62 项 | FAIL 1（**非本轮**：`paper_assets/audit/verify_wheel_install.py:52` 绝对路径，属并行论文工作区），本轮文件 0 FAIL |

### 规格 §9 十二项的落点

| # | 要求 | 落点 |
|---|---|---|
| 1 | 五种 High-level 动态预算 | `test_high_level_budget_five_scenarios` + renderer check `test_five_model_budgets` |
| 2 | gene 阈值真筛选 & 两列对应 | `test_topk_truncation`、`test_selection_matches_figure`、renderer check `test_threshold_and_matching_gene_columns` |
| 3 | 单一 PC 节点与 signed sum | `test_pc_is_one_node_with_signed_sum`、renderer check `test_encoded_column_sum` |
| 4 | 临床 imputation 标注 | `test_clinical_complete_and_imputed` |
| 5 | 隐藏项仍进入完整分解 | `test_full_decomposition_sums_to_eta_with_hidden_terms` |
| 6 | 完整 contribution sum == η | 同上 + `test_full_risk_matches_main_api` + renderer check `test_signed_decomposition_must_equal_eta` |
| 7 | 分数提取前后既有输出一致 | `test_score_extraction_does_not_change_existing_outputs` |
| 8 | HTML 内嵌 == PNG 字节 | `test_html_embeds_the_same_png_bytes`、renderer check `test_html_is_same_png_and_node_values_unchanged` |
| 9 | marker / cutoff / ties / 极端 percentile | `test_risk_marker_extremes_and_ties`、`test_cutoff_median_and_numeric`、renderer check `test_risk_endpoints_and_nonmedian_cutoff`、`test_median_with_ties` |
| 10 | 重复导出不改输入、不随机 | `test_repeated_export_is_deterministic` |
| 11 | 缺绘图 extra 时核心可用 | `test_core_import_and_compute_without_matplotlib` |
| 12 | 既有公开 API 兼容 | `test_public_api_compatibility` |

---

## 8. clean-wheel 安装与依赖检查

**构建**：`git worktree add --detach /tmp/gh_audit/cleanwt HEAD`（HEAD = `55df2c6`），
只复制**本功能相关文件**（13 个），确认 worktree 内**没有** `paper_assets/` 等并行工作区文件，
再 `python -m build --no-isolation --wheel --sdist`。

| 检查 | 结果 |
|---|---|
| 产物 | `compass_os-1.1.0-py3-none-any.whl`、`compass_os-1.1.0.tar.gz`（SHA256 见 `docs/audit/00_CANDIDATE_BASELINE.md` §5） |
| wheel 条目 | 83 |
| `paper_assets` 污染 | **无** |
| 新模块是否入包 | `sample_path.py` / `sample_path_render.py` / `_render_bridge.py` **均在** |
| METADATA 依赖 | `numpy>=1.26`, `pandas>=2.0`, `scikit-learn>=1.3`, `scikit-survival>=0.22`, **`matplotlib>=3.7`**；extras：`upstream`/`test`/`plot` |
| Plotly / Kaleido | **不在任何依赖中** |

**fresh environment**：`python -m venv --system-site-packages` 新建环境，
`pip install --no-deps --no-cache-dir --force-reinstall <wheel>`，在**仓库之外**的目录运行：

| 检查 | 结果 |
|---|---|
| 包路径来自 `site-packages`、`sys.path` 不含仓库 | **True** |
| `import compass_os` | OK（`1.1.0`） |
| 核心预测 `predict(model="M2,M3")` | OK |
| 新 `sample_path(model="M3")` | OK（12 concepts + Age/Sex/Stage/PC1–PC10，88 行分解） |
| PNG / PDF / SVG 导出 | OK（649 KB / 57.6 KB / 316 KB） |
| HTML 内嵌 == PNG 字节 | **True** |
| 图中 η == 主 API、分解求和 == η | **True / True** |
| 重复导出确定性 | **True** |
| 运行时资产查找 | OK（冻结 checkpoint 从包内 `assets/` 读取，无仓库绝对路径） |
| 读取未提交的论文资产 | **无** |

### 关于 matplotlib 依赖

`matplotlib` 自 v1.0 起即为**核心依赖**（`AnalysisResult.save_report` 出图需要），
**本功能没有引入任何新增硬依赖**（反而去掉了 Plotly / Kaleido / CDN）。
渲染层仍**惰性导入** matplotlib：模拟其缺失时
`import compass_os`、核心预测、`sample_path()` 的**计算**部分全部正常，
仅导出时报出明确提示（含 `pip install "compass-os[plot]"`）。
本轮**未**全局修改 seaborn 样式、**未**关闭全部 warnings、**未**引入影响用户其它图形的全局配置
（渲染层仅在 `matplotlib.rc_context(...)` 内设置字体与 PDF/SVG 嵌入，退出即还原）。

---

## 9. 尚未解决的真实限制

1. **临床数据缺失使临床列无法体现样本间差异**：GSE39582 无临床信息，Age/Sex/Stage 三项
   **全部**冻结填补，因此这三个 contribution 对全队列**恒定**
   （实测 Age 恒 +0.032、Sex 恒 −0.000、Stage 恒 +0.927）。图中已如实标注 `(imputed)`，
   但若要用这张图讨论临床变量的作用，必须换**有临床信息**的队列。
2. **M3 的 PC 聚合贡献在三个样本上均较大**（本队列实测约 0.6–1.5），会明显压缩 concept 的
   视觉区分度；这是真实数值而非绘图放大，**未**做任何调整。
3. **gene 层与 signature/concept 层来自两次 float32 前向**：gene 层必须额外调用
   `with_gene_level=True`，与前两者相差约 `1e-7`。该差异**不影响**任何展示数值的正确性
   （risk / rank / percentile / 分解全部取主 API 口径），但若将来要求"单次前向内取全三层"，
   需要修改上游 `extract()` 的返回契约。
4. **阈值 `gene_threshold=0.25` 与颜色范围是展示参数**，**不是**经临床验证的科学阈值；
   若用于正式论文图，应基于真实的 gene-score 分布写出选择依据。
5. **HTML 主图为静态一致版本**（按规格取舍），没有 hover / 缩放；精确数值需查内嵌节点表或
   `*.nodes.tsv`。
6. **`decomposition` 的 `value` 列是设计矩阵取值**（Age 为标准化后取值，Sex/Stage 为原始
   编码值），与临床节点显示的 `value`（原始冻结填补值）口径不同；两者语义已在
   `nodes.tsv` 与文档中分别说明。
7. **本功能未处理**：多癌种队列的真实运行（只在单元测试用合成 cohort 覆盖了 12/11 concepts
   的预算分支）、以及有临床信息队列的 imputed 混合情形。
8. 仓库内**仍保留**上一轮的历史文档 `SAMPLE_SANKEY_FEASIBILITY_AUDIT.md` /
   `SAMPLE_SANKEY_VALIDATION.md` / `SAMPLE_PATH_VISUAL_REVIEW.md` /
   `SAMPLE_PATH_REDESIGN_WIRING.md`，它们描述的是 **v1.1.x 的旧渲染器**（Plotly Sankey、
   流量守恒视图）。**当前实现以本报告与 `docs/SAMPLE_PATH.md` 为准**；
   历史文档未删除，以免丢失审计轨迹。

---

## 10. 停止点

已完成：代码接入、真实数据试运行（真实 M2/M3 同样本 6 图）、导出（HTML/PNG/PDF/SVG +
完整 decomposition + selection 审计）、测试（132 + 12）、文档。
**未** merge、**未** tag、**未** push、**未**发布。
工作树中与本功能无关的 `paper_assets/` 等并行改动**未被触碰、未被提交**。
