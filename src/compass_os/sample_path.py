# -*- coding: utf-8 -*-
"""样本级**计算路径归因**可视化（sample-specific computational attribution）。

对已分析 cohort 中的**某一个样本**，取真实推理结果并绘制：

    Gene expression → COMPASS gene score → Granular signature score
        → High-level concepts + clinical / PC → Cohort-relative risk

──────────────────────── 语义边界（代码 / 图注 / 文档一律遵守）────────────────────────
这是 **sample-specific computational attribution / representation-flow visualization**，
即"样本级计算路径归因"。它**不是**因果机制图、**不是**通路激活图、**不是**生物因果网络，
也**不是**对最终风险的完整因果解释。
节点颜色表示该层**真实分数**（临床 / PC 节点表示 signed Cox contribution，用独立图例）；
所有连线**等粗等色**，只表达筛选后的模型连接关系，**不编码权重或流量**，
本图**不作流量守恒声明**。

──────────────────────────── 数据来源（全部真实，绝不臆测）────────────────────────────
* **gene-token 标量分数**：官方 ``PreTrainer.extract(..., with_gene_level=True)`` 的 ``dfg``
  = ``genesetprojector.geneset_scorer(encoder_output)[:, 2:]``——对每个基因位置施加**共享的**
  ``nn.Linear(32→1)``（源码 ``assets/third_party/compass/model/tune.py:289-291``，
  张量形状 ``(n_samples, 15672)``）。
  **不得**用 attention/projector weight、TPM、gene→signature 风险分配量或旧 HTML 的
  upstream contribution 冒充。
* **132 signature / 43 concept**：主 API ``predict()`` **同一次前向**的输出，与 η 严格同源。
* **临床 / Cancer type / PC1–PC10**：冻结 Cox 设计矩阵**逐列 signed contribution** 的聚合
  （多编码列取 signed sum；PC 为十个冻结设计列贡献之和，沿用原标准化口径）。
* **η / rank / percentile / cutoff**：主 API 直接输出；排名定义写入
  ``SamplePathData.rank_definition``，绘图层**不重算**排名。
* ``decomposition`` 为**完整未截断**分解，含未展示 concept 与被隐藏的 Cancer type：
  ``sum(contribution) == η``（渲染层强制校验）。

绘图层位于 :mod:`compass_os.sample_path_render`（matplotlib 惰性导入）。
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd

from . import survival as _surv
from .exceptions import CompassOSError, InputError
from .preprocessing import align_expression, cancer_codes_for
from .representation import load_model

__all__ = ["SamplePathResult", "sample_path", "M0_UNSUPPORTED_MESSAGE"]

M0_UNSUPPORTED_MESSAGE = (
    "M0 contains no COMPASS concept branch and is not supported by the sample-level "
    "representation Sankey visualization.")

#: 高维层默认展示上限（与渲染层 Style 默认一致；规格允许 1–50）
DEFAULT_MAX_GENES = 32
DEFAULT_MAX_SIGNATURES = 26
#: high-level predictor 层默认预算（规格冻结为 16）
DEFAULT_PREDICTOR_BUDGET = 16
#: 归因链条与主 API 的一致性**相对**容差（模型为 float32 前向）
ETA_TOL_REL = 1e-4

SEMANTICS = (
    "Sample-specific computational attribution through the frozen model graph "
    "(representation-flow visualization). Not a causal mechanism diagram, not a pathway "
    "activation map, not a biological causal network, and not a complete causal explanation "
    "of the predicted risk. Node colour encodes that layer's own real score (clinical/PC "
    "nodes encode their signed Cox contribution, with a separate legend). All links share one "
    "width and colour and only express the selected model connectivity; contribution is NOT "
    "encoded as flow and no flow-conservation claim is made. Representation colours denote "
    "high/low model representation values and must not be read as pathway activation or "
    "inhibition."
)


# --------------------------------------------------------------------------- #
@dataclass
class SamplePathResult:
    """样本级计算路径归因结果（绘图见 :meth:`render` / :meth:`save`）。"""

    sample_id: str
    model: str
    cancer_type: str
    #: 主 API 的 Cox 线性预测子 η（**不是**概率，也不是 exp(η)）
    risk: float
    rank: int
    n_samples: int
    percentile: float
    cutoff: float
    cutoff_mode: str
    cutoff_is_frozen_validated: bool
    #: **完整未截断**的逐特征 signed Cox 分解：feature / value / beta / contribution
    decomposition: pd.DataFrame
    #: 实际展示的三层——由渲染层 ``select_nodes`` 决定，故与图**永远一致**
    genes: pd.DataFrame
    signatures: pd.DataFrame
    concepts: pd.DataFrame
    #: 展示的临床 / PC / Cancer type 节点
    auxiliaries: pd.DataFrame
    eta_max_abs_diff: float = 0.0
    #: high-level predictor 节点数（concept + 临床/PC/cancer），≤ predictor_budget
    n_high_level_nodes: int = 0
    n_high_level_nodes_total: int = 0
    notes: list = field(default_factory=list)
    semantics: str = SEMANTICS
    #: 内部渲染载荷（完整真实数组 + 冻结拓扑）；普通用户无需准备
    render_payload: dict = None
    #: 实际使用的视觉规格（:class:`compass_os.sample_path_render.Style`）
    style: object = None

    # ---------------- 导出（单一渲染管线）----------------
    def render(self, style=None):
        """返回 :class:`compass_os.sample_path_render.RenderedPath`。

        PNG / PDF / SVG 由**同一个 Figure** 导出，HTML 内嵌**同一份 PNG 字节**。
        """
        from ._render_bridge import to_render_data
        from .sample_path_render import plot_sample_path
        return plot_sample_path(to_render_data(self), style or self.style)

    def figure(self, **kw):
        """返回 ``RenderedPath``（含 ``.figure`` / ``.selection`` / ``.cutoff``）。

        .. note:: 自 v1.1.0 起图形由包内统一渲染层（Matplotlib）绘制。
           旧 Plotly ``Figure``（``write_html`` / ``write_image`` / ``include_plotlyjs``）
           已移除；请改用 :meth:`save_html` / :meth:`save_static` / :meth:`save`。
        """
        return self.render(**kw)

    def save_html(self, path) -> Path:
        """**可独立调用**：导出 HTML（内嵌同一份 PNG 字节 + 可展开精确节点表）。"""
        return self.render().save_html(path)

    def save_static(self, path) -> Path:
        """导出静态图：``.png`` / ``.pdf`` / ``.svg``（与 HTML 同一 Figure）。"""
        return self.render().save_static(path)

    def save(self, prefix) -> dict:
        """一次导出 png / pdf / svg / html / nodes.tsv / selection.json / decomposition.tsv。"""
        return self.render().save(prefix)

    # ---------------- 文本摘要 ----------------
    def summary(self) -> str:
        L = [f"Sample computation path — {self.sample_id}",
             "─" * 64,
             f"{'Model':<30}{self.model}",
             f"{'Cancer type':<30}{self.cancer_type}",
             f"{'Risk (Cox linear predictor η)':<30}{self.risk:.6f}",
             f"{'Rank':<30}{self.rank} / {self.n_samples}",
             f"{'Percentile':<30}{self.percentile:.1f}%",
             f"{'Cutoff':<30}{self.cutoff:.6f} ({self.cutoff_mode})",
             "",
             f"{'Shown genes (both columns)':<30}{len(self.genes)}",
             f"{'Shown signatures':<30}{len(self.signatures)}",
             f"{'Shown concepts':<30}{len(self.concepts)}",
             f"{'Clinical / PC nodes':<30}{len(self.auxiliaries)}",
             f"{'High-level predictors':<30}{self.n_high_level_nodes}",
             f"{'Full decomposition rows':<30}{len(self.decomposition)}",
             f"{'Decomposition vs API risk':<30}{self.eta_max_abs_diff:.3e}",
             "", self.semantics]
        if self.notes:
            L += ["", "NOTES"] + ["  - " + n for n in self.notes]
        return "\n".join(L)

    def __repr__(self) -> str:  # pragma: no cover
        return (f"SamplePathResult(sample_id={self.sample_id!r}, model={self.model!r}, "
                f"risk={self.risk:.4f}, rank={self.rank}/{self.n_samples})")


# --------------------------------------------------------------------------- #
def sample_path(expression: pd.DataFrame, cancer_type: Sequence[str] | str, sample_id: str, *,
                model: str = "M2", clinical: pd.DataFrame | None = None,
                input_scale: str = "tpm", missing_gene_strategy: str = "reference",
                style=None, gene_threshold: float | None = None,
                max_genes: int | None = None, max_signatures: int | None = None,
                predictor_budget: int | None = None,
                cutoff: str | float = "median", batch_size: int = 16,
                show_cancer_type: bool | None = None, result=None,
                # ---- 兼容 v1.1.x 的旧参数名 ----
                top_genes: int | None = None, top_signatures: int | None = None,
                top_concepts: int | None = None,
                max_high_level_nodes: int | None = None) -> SamplePathResult:
    """取真实推理结果，绘制单个样本的计算路径归因图。

    参数
    ----
    expression, cancer_type, clinical, input_scale, missing_gene_strategy
        与 :func:`compass_os.predict` 完全一致（同一套对齐与冻结填补）。
    sample_id
        必须出现在 ``expression.index`` 中；否则抛 :class:`InputError`。
    model
        ``M1`` / ``M2`` / ``M3``。``M0`` 无 COMPASS concept 分支，明确不支持。
    style
        :class:`compass_os.sample_path_render.Style`。``None`` = 冻结默认规格。
        单独参数（``gene_threshold`` / ``max_genes`` / ``max_signatures`` /
        ``predictor_budget``）会覆盖其中对应字段。**对比 M2/M3 时应传同一份 style**，
        以保证画布、节点尺寸、列位置、字体与颜色尺度一致。
    gene_threshold
        作用于**真实 gene-token score** 尺度的展示阈值（默认 0.25，仅显示参数）。
    cutoff
        ``"median"`` = **当前 cohort** η 中位数（**仅可视化分界**，不是冻结 validated
        cutoff）；或研究者给定的数值阈值。
    result
        可传入既有 ``PredictionResult`` / ``AnalysisResult``，复用其 cohort 风险与 132/43
        两层（此时不再重复预测）。

    返回
    ----
    :class:`SamplePathResult`；图形经 ``.save_html()`` / ``.save_static()`` / ``.save()``
    导出，或 ``.render()`` 取 :class:`RenderedPath`。
    """
    model = str(model).upper()
    if model == "DEFAULT":
        model = _surv.default_model()
    if model == "M0":
        raise InputError(M0_UNSUPPORTED_MESSAGE)
    if model not in ("M1", "M2", "M3"):
        raise InputError(f"sample_path 仅支持 M1/M2/M3，收到 {model!r}")

    idx = [str(i) for i in expression.index]
    if sample_id not in idx:
        raise InputError(
            f"sample_id {sample_id!r} 不在表达矩阵索引中（共 {len(idx)} 个样本；"
            f"前几个：{idx[:3]}）")
    ct = ([cancer_type] * len(idx)
          if isinstance(cancer_type, str) else [str(c) for c in cancer_type])
    if len(ct) != len(idx):
        raise InputError(f"cancer_type 长度 {len(ct)} 与样本数 {len(idx)} 不一致")

    notes: list[str] = []
    # 旧参数名 → 新参数（保持向后兼容）
    if top_genes is not None and max_genes is None:
        max_genes = top_genes
    if top_signatures is not None and max_signatures is None:
        max_signatures = top_signatures
    if max_high_level_nodes is not None and predictor_budget is None:
        predictor_budget = max_high_level_nodes
    if top_concepts is not None:
        notes.append(
            "`top_concepts` is deprecated: the number of concept nodes is "
            "`predictor_budget` minus the mandatory clinical/PC nodes (per the frozen budget "
            "table). Use `predictor_budget` instead.")

    # ---- cohort 层：η / rank / percentile / cutoff（主 API 直接输出）----
    if result is not None:
        risk_series = _risk_from_result(result, model)
        concept_all = getattr(result, "concept_scores", None)
        sig_all = getattr(result, "signature_scores", None)
    else:
        import compass_os
        pred = compass_os.predict(expression, ct, clinical=clinical, model=model,
                                  missing_gene_strategy=missing_gene_strategy,
                                  input_scale=input_scale, batch_size=batch_size)
        risk_series = pred.risk[model]
        concept_all, sig_all = pred.concept_scores, pred.signature_scores
    risk_series = pd.Series(risk_series).reindex(idx).astype(float)
    if not np.isfinite(risk_series.to_numpy()).all():
        raise CompassOSError("cohort risk 含非有限值；请检查输入")

    risks = risk_series.to_numpy()
    risk_i = float(risk_series.loc[sample_id])
    rank = int((risks > risk_i).sum()) + 1
    percentile = float(100.0 * (risks < risk_i).mean() + 50.0 * (risks == risk_i).mean())

    if isinstance(cutoff, str):
        if cutoff != "median":
            raise InputError("cutoff 只能是 'median' 或数值阈值")
        cutoff_value = float(np.median(risks))
        cutoff_mode = ("cohort median of the current risk output (visualization boundary "
                       "only — NOT the frozen validated prognostic cutoff)")
    else:
        cutoff_value, cutoff_mode = float(cutoff), "user-specified risk threshold"

    # ---- 样本层：gene-token 分数（官方 API）+ 132/43（与 η 同源）----
    aligned = align_expression(expression, missing_gene_strategy, input_scale=input_scale)
    codes = cancer_codes_for(ct).to_numpy(int)
    srow = [idx.index(sample_id)]
    layers = _extract_gene_level(
        aligned.matrix.iloc[srow], codes[srow], batch_size, missing_gene_strategy,
        None if aligned.missing_mask is None else aligned.missing_mask.iloc[srow])
    gene_score = layers["gene"].iloc[0]

    def _layer_or_extract(all_df, key, label):
        if all_df is not None and sample_id in {str(i) for i in all_df.index}:
            row = all_df.loc[sample_id]
            return pd.Series(row.to_numpy(float), index=[str(c) for c in all_df.columns])
        notes.append(f"{label} taken from this sample's own extraction because the supplied "
                     "result object does not expose it.")
        return layers[key].iloc[0]

    sig_score = _layer_or_extract(sig_all, "signature", "signature layer")
    con_score = _layer_or_extract(concept_all, "concept", "concept layer")

    # ---- 完整 signed Cox 分解（与主 API 同一套函数；含隐藏项与 offset）----
    lock = _surv.load_lock(model)
    feats = list(lock["feature_names"])
    beta = np.asarray(lock["beta"], np.float64)
    scaled = list(lock["scaled_features"] or [])
    mu = np.asarray(lock["scaler_mean"] or np.zeros(len(scaled)), np.float64)
    sd = np.asarray(lock["scaler_scale"] or np.ones(len(scaled)), np.float64)
    con_row = pd.DataFrame([con_score.to_numpy(np.float64)], index=[sample_id],
                           columns=[str(c) for c in con_score.index])
    clin_row = None if clinical is None else clinical.reindex(idx).loc[[sample_id]]
    Xd = _surv.build_design_matrix(con_row, clin_row, aligned.log2.iloc[srow], model,
                                   [ct[srow[0]]])
    Z = Xd.copy()
    if scaled:
        Z[scaled] = (Xd[scaled] - mu) / sd
    contrib = Z.to_numpy(np.float64).ravel() * beta
    dec = pd.DataFrame({"feature": feats, "value": Z.to_numpy(np.float64).ravel(),
                        "beta": beta, "contribution": contrib})
    contrib_of = dict(zip(feats, contrib))
    eta_diff = abs(float(contrib.sum()) - risk_i)
    tol = ETA_TOL_REL * max(1.0, abs(risk_i))
    if eta_diff > tol:
        raise CompassOSError(
            f"内部一致性检查失败：完整分解求和 η={contrib.sum():.9f} 与主 API risk="
            f"{risk_i:.9f} 相差 {eta_diff:.3e}（容差 {tol:.3e}）；拒绝用不一致的数值绘图")

    # ---- 渲染载荷：完整真实数组 + 冻结拓扑（内部，用户无需准备）----
    imputed = _imputed_fields(clinical, idx, sample_id)
    payload = {
        "cohort_risks": [float(v) for v in risks],
        "cancer_types": sorted(set(ct)),
        "cohort_label": (f"{sorted(set(ct))[0]} cohort" if len(set(ct)) == 1
                         else "multi-cancer cohort"),
        "gene_expression": {str(g): float(aligned.matrix.iloc[srow[0]][g])
                            for g in aligned.matrix.columns},
        "gene_score": {str(g): float(v) for g, v in gene_score.items()},
        "signature_scores": {str(k): float(v) for k, v in sig_score.items()},
        "concept_scores": {str(k): float(v) for k, v in con_score.items()},
        "concept_contributions": {str(c): float(contrib_of[c]) for c in con_score.index},
        "clinical": {nm: {"contribution": float(contrib_of.get(nm, 0.0)),
                          "value": float(lock["clinical_fill"].get(nm, {}).get("fill", np.nan)),
                          "imputed": bool(imputed.get(nm, True))}
                     for nm in ("Age", "Sex", "Stage") if nm in feats},
        "pc_columns": [f for f in feats if f.startswith("PC") and f[2:].isdigit()],
        "cancer_columns": [f for f in feats if f.startswith("CT_")],
        "controller_columns": [f for f in feats if f.startswith("CT_")],   # 别名，供桥接使用
        "decomposition": [{"feature": f, "contribution": float(c)}
                          for f, c in zip(feats, contrib)],
        "gene_signature_edges": _topology_edges(load_model()),
        "signature_concept_edges": _signature_concept_edges(load_model()),
        "expression_units": ("TPM" if input_scale == "tpm"
                             else "TPM (converted from log2(TPM+1) input)"),
        "input_scale": input_scale,
        "missing_gene_strategy": missing_gene_strategy,
        "n_missing_genes": int(aligned.qc.n_missing_genes),
    }
    if show_cancer_type:
        notes.append(
            "`show_cancer_type` follows the frozen budget table: the constant cancer term is "
            "hidden for single-cancer cohorts and shown for multi-cancer cohorts. The flag was "
            "recorded but does not override that rule.")

    res = SamplePathResult(
        sample_id=sample_id, model=model, cancer_type=ct[srow[0]], risk=risk_i, rank=rank,
        n_samples=len(idx), percentile=percentile, cutoff=cutoff_value,
        cutoff_mode=cutoff_mode, cutoff_is_frozen_validated=False,
        decomposition=dec, genes=pd.DataFrame(columns=["gene", "expression", "gene_score"]),
        signatures=pd.DataFrame(columns=["signature", "score"]),
        concepts=pd.DataFrame(columns=["concept", "score", "contribution"]),
        auxiliaries=pd.DataFrame(columns=["predictor", "kind", "value", "contribution",
                                          "imputed", "detail"]),
        eta_max_abs_diff=eta_diff, notes=notes, render_payload=payload,
        style=_resolve_style(style, gene_threshold=gene_threshold, max_genes=max_genes,
                             max_signatures=max_signatures, predictor_budget=predictor_budget))
    _fill_selection(res)
    return res


# --------------------------------------------------------------------------- #
def _resolve_style(style, **over):
    """合并用户 style 与显式覆盖项；未给出的沿用渲染层默认值。"""
    from .sample_path_render import Style
    base = dict(style.__dict__) if style is not None else {}
    for k, v in over.items():
        if v is not None:
            base[k] = v
    return Style(**base)


def _fill_selection(res: SamplePathResult) -> None:
    """用渲染层的 ``select_nodes`` 决定展示节点，使 ``res.*`` 与图**永远一致**。

    ``select_nodes`` 只依赖标准库 + numpy，不会触发 matplotlib 导入。
    """
    from ._render_bridge import to_render_data
    from .sample_path_render import select_nodes

    sel = select_nodes(to_render_data(res), res.style)
    res.render_payload["selection_counts"] = dict(sel["counts"])
    res.genes = pd.DataFrame(
        [{"gene": n.name, "expression": e.score, "gene_score": n.score}
         for n, e in zip(sel["gene_score"], sel["gene_tpm"])],
        columns=["gene", "expression", "gene_score"])
    res.signatures = pd.DataFrame([{"signature": n.name, "score": n.score}
                                   for n in sel["signatures"]],
                                  columns=["signature", "score"])
    res.concepts = pd.DataFrame([{"concept": n.name, "score": n.score,
                                  "contribution": n.contribution}
                                 for n in sel["predictors"] if n.kind == "concept"],
                                columns=["concept", "score", "contribution"])
    res.auxiliaries = pd.DataFrame(
        [{"predictor": n.name, "kind": n.kind, "value": n.value,
          "contribution": n.contribution, "imputed": bool(n.imputed),
          "detail": (n.details or {}).get("column_contributions")}
         for n in sel["predictors"] if n.kind in ("clinical", "pc", "cancer")],
        columns=["predictor", "kind", "value", "contribution", "imputed", "detail"])
    res.n_high_level_nodes = int(sel["counts"]["predictors"])
    res.n_high_level_nodes_total = res.n_high_level_nodes
    if int(sel["counts"]["predictors"]) > res.style.predictor_budget:
        raise CompassOSError(
            f"high-level predictor 节点数 {sel['counts']['predictors']} 超过预算 "
            f"{res.style.predictor_budget}")


# --------------------------------------------------------------------------- #
def _extract_gene_level(aligned_tpm: pd.DataFrame, cancer_code: Sequence[int],
                        batch_size: int, strategy: str, missing_mask) -> dict:
    """官方 ``extract(with_gene_level=True)``：一次取到 gene / 132 / 43 三层。"""
    from .representation import _zero_masked

    model = load_model()
    genes = [str(g) for g in model.feature_name]
    X = aligned_tpm.reindex(columns=genes).astype(np.float64).copy()
    X.insert(0, "cancer_code", np.asarray(list(cancer_code), dtype=int))
    masked: list = []
    if strategy == "zero" and missing_mask is not None:
        mm = missing_mask.reindex(columns=genes).fillna(False).to_numpy(bool)
        masked = [g for j, g in enumerate(genes) if mm[:, j].any()]
        X.loc[:, masked] = 0.0
    with _zero_masked(model, masked):
        dfg, dfgs, dfct = model.extract(X, batch_size=batch_size, num_workers=0,
                                        with_gene_level=True)
    drop = ("CANCER", "PID")
    out = {}
    for key, df in (("gene", dfg), ("signature", dfgs), ("concept", dfct)):
        df = df.drop(columns=[c for c in drop if c in df.columns])
        df.index = [str(i) for i in df.index]
        out[key] = df
    return out


def _topology_edges(model) -> list:
    """冻结 gene→gene-set 拓扑（集合内重复基因的注意力权重**累加**）。"""
    import torch.nn.functional as F
    gsp = model.model.latentprojector.genesetprojector
    feats = [str(g) for g in model.feature_name]
    names = [str(x) for x in gsp.genesets_names]
    out = []
    for j, idxs in enumerate(gsp.genesets_indices):
        w = F.softmax(gsp.geneset_aggregator.aggregator
                      .attention_weights[f"geneset_{j}"].detach(), dim=0).numpy().ravel()
        agg: dict = {}
        for g, a in zip(idxs, w):
            agg[int(g)] = agg.get(int(g), 0.0) + float(a)
        for g, a in agg.items():
            out.append({"source": feats[g], "target": names[j], "weight": a})
    return out


def _signature_concept_edges(model) -> list:
    """冻结 gene-set→concept 拓扑（cellpathway 注意力；仅用于连接筛选）。"""
    import torch.nn.functional as F
    lp = model.model.latentprojector
    cpp, gsp = lp.cellpathwayprojector, lp.genesetprojector
    names = [str(x) for x in cpp.cellpathway_names]
    gs = [str(x) for x in gsp.genesets_names]
    out = []
    for c, idxs in enumerate(cpp.cellpathway_indices):
        w = F.softmax(cpp.cellpathway_aggregator.aggregator
                      .attention_weights[f"cellpathway_{c}"].detach(), dim=0).numpy().ravel()
        for j, b in zip(idxs, w):
            out.append({"source": gs[int(j)], "target": names[c], "weight": float(b)})
    return out


def _imputed_fields(clinical: pd.DataFrame | None, idx: list, sample_id: str) -> dict:
    """逐样本、逐字段判断临床值是否**实际被冻结填补**（不做数值猜测）。"""
    out = {}
    for nm, col in (("Age", "age"), ("Sex", "sex"), ("Stage", "stage")):
        if clinical is None or col not in getattr(clinical, "columns", []):
            out[nm] = True
            continue
        try:
            out[nm] = bool(pd.isna(clinical.reindex(idx).loc[sample_id, col]))
        except Exception:  # noqa: BLE001
            out[nm] = True
    return out


def _risk_from_result(result, model: str) -> pd.Series:
    r = getattr(result, "risk", None)
    if r is None:
        raise InputError("传入的 result 不含 risk")
    df = pd.DataFrame(r)
    if model in df.columns:
        return df[model]
    dm = getattr(result, "default_model", None)
    if dm in df.columns:
        return df[dm]
    raise InputError(f"result 中找不到模型 {model} 的 risk 列（现有：{list(df.columns)}）")
