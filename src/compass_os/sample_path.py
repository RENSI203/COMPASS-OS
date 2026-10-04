# -*- coding: utf-8 -*-
"""样本级**计算路径归因**可视化（sample-specific computational attribution）。

本模块回答的是：*对这一个样本，冻结模型图里各层分别贡献了多少到最终 Cox 线性预测子*，
以及该样本在当前 cohort 风险分布中的位置。

**语义边界（代码、图注、文档一律遵守）**：这是 **sample-specific computational
attribution / representation-flow visualization**，即"样本级计算路径归因"。
它**不是**因果机制图、**不是**通路激活图、**不是**生物因果网络，
也**不是**对最终风险的完整因果解释。图中的颜色只表示 **high / low model
representation value** 或**贡献方向（正=推高风险 / 负=降低风险）**，
**不得**读作 pathway activated / inhibited。

──────────────────────────────── 精确性 ────────────────────────────────
本模块**不复刻**任何预测逻辑，全部数值来自冻结资产与官方 API：

* 132/43 两层与 gene 层：官方 ``PreTrainer.extract(..., with_gene_level=True)``；
* cohort 风险：``compass_os.predict``（或调用方传入的既有结果）；
* 逐特征 Cox 贡献：包内 **同一套** ``survival.build_design_matrix`` + 冻结 ``beta``。

已实测的精确性（12 例 golden 夹具，reference 策略）：

* ``Σ_g a_{g,j}·gene_score_g`` vs 官方 gene-set 分数：``max|Δ| = 1.8e-07``（float32 精度）
* ``Σ_j b_{j,c}·geneset_score_j`` vs 官方 concept 分数：``max|Δ| = 3.3e-08``
* ``Σ(逐特征贡献)`` vs 主 API ``risk``：M1 ``2.2e-16`` / M2 ``4.4e-16`` / M3 ``8.9e-16``

⚠ 实现要点：``geneset_indices`` 中有 2 个 gene set 存在**重复基因索引**
（如 ``geneset_21``：15 项 / 14 唯一）。注意力权重的累加必须用 **scatter-add**
（``np.add.at``）；若用普通赋值会静默丢项并产生 ~7e-02 的误差。
"""
from __future__ import annotations

import html
from dataclasses import dataclass, field
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd

from . import survival as _surv
from ._paths import ensure_compass
from .exceptions import CompassOSError, InputError
from .preprocessing import align_expression, cancer_codes_for
from .representation import load_model

__all__ = ["SamplePathResult", "sample_path", "M0_UNSUPPORTED_MESSAGE"]

M0_UNSUPPORTED_MESSAGE = (
    "M0 contains no COMPASS concept branch and is not supported by the sample-level "
    "representation Sankey visualization.")

#: 高维层默认展示上限（gene / signature / concept 各自的 top-K）
DEFAULT_TOP = 50
#: high-level（concept + 辅助 predictor）层的**总节点上限**
DEFAULT_MAX_HIGH_LEVEL = 16
#: 归因链条与主 API 的一致性**相对**容差。冻结模型以 float32 前向，归因与 cohort 预测
#: 各自做一次前向，batch 相关的 BLAS 求和顺序不同 ⇒ η 实测差 ~1e-06…1e-05。
#: 取相对容差 1e-4（对 |η|≈3 即 ~3e-4），足以捕获任何真实的分解错误（错解会差几个量级），
#: 又不会因 float32 舍入而误报。图中**显示**的 risk/rank/percentile 一律取主 API 结果。
ETA_TOL_REL = 1e-4

SEMANTICS = (
    "Sample-specific computational attribution through the frozen model graph "
    "(representation-flow visualization). Not a causal mechanism diagram, not a pathway "
    "activation map, not a biological causal network, and not a complete causal explanation "
    "of the predicted risk. Link width is |contribution|; colour encodes contribution "
    "direction (positive = higher model-estimated risk, negative = lower). "
    "Representation colours denote high/low model representation values and must not be read "
    "as pathway activation or inhibition."
)

_POS = "#c0392b"        # 正贡献：推高风险
_NEG = "#2471a3"        # 负贡献：降低风险
_CONCEPT_C = "#7d3c98"
_AUX_C = "#566573"
_OTHER_C = "#95a5a6"
_RISK_C = "#1c2833"
_GENE_C = "#5d6d7e"


# --------------------------------------------------------------------------- #
# 官方中间量抽取
# --------------------------------------------------------------------------- #
def _extract_gene_level(aligned_tpm: pd.DataFrame, cancer_code: Sequence[int],
                        batch_size: int, strategy: str, missing_mask) -> dict:
    """调用官方 ``extract(with_gene_level=True)``，返回 gene/132/43 三层。"""
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
    drop = [c for c in ("CANCER", "PID")]
    dfg = dfg.drop(columns=[c for c in drop if c in dfg.columns])
    dfgs = dfgs.drop(columns=[c for c in drop if c in dfgs.columns])
    dfct = dfct.drop(columns=[c for c in drop if c in dfct.columns])
    dfg.index = [str(i) for i in dfg.index]
    dfgs.index = [str(i) for i in dfgs.index]
    dfct.index = [str(i) for i in dfct.index]
    return {"gene": dfg, "signature": dfgs, "concept": dfct}


def _attention_matrices(model, need_genesets: Sequence[int]):
    """构建 ``A``（gene-set → gene，scatter-add，处理重复索引）与 ``B``（concept → gene-set）。

    只构建到展示所需的 gene set，避免 (132 × 15672) 的稠密矩阵。
    """
    import torch.nn.functional as F

    lp = model.model.latentprojector
    gsp, cpp = lp.genesetprojector, lp.cellpathwayprojector
    gs_idx = [list(map(int, x)) for x in gsp.genesets_indices]
    cp_idx = [list(map(int, x)) for x in cpp.cellpathway_indices]
    n_genes = len(gsp.genesets_names) and len(model.feature_name)

    gw = gsp.geneset_aggregator.aggregator.attention_weights
    A = {j: (np.asarray(gs_idx[j]),
             F.softmax(gw[f"geneset_{j}"].detach(), dim=0).numpy().ravel())
         for j in need_genesets}
    cw = cpp.cellpathway_aggregator.aggregator.attention_weights
    B = {}
    for c, idxs in enumerate(cp_idx):
        B[c] = (np.asarray(idxs),
                F.softmax(cw[f"cellpathway_{c}"].detach(), dim=0).numpy().ravel())
    return A, B, gs_idx, cp_idx, n_genes


# --------------------------------------------------------------------------- #
# 结果对象
# --------------------------------------------------------------------------- #
@dataclass
class SamplePathResult:
    """样本级计算路径归因结果。"""

    sample_id: str
    model: str
    cancer_type: str
    #: 主 API 的 Cox 线性预测子（**不是**概率）
    risk: float
    #: 队列内排名（1 = 风险最高）
    rank: int
    n_samples: int
    #: 0–100，越大风险越高
    percentile: float
    cutoff: float
    cutoff_mode: str
    cutoff_is_frozen_validated: bool
    #: 逐特征完整 Cox 分解（feature / value / beta / contribution）
    decomposition: pd.DataFrame
    #: 展示的三层（含 contribution 列）
    genes: pd.DataFrame
    signatures: pd.DataFrame
    concepts: pd.DataFrame
    #: 直接连到 risk 的辅助 predictor（Age/Sex/Stage/PC 聚合/Cancer type）
    auxiliaries: pd.DataFrame
    #: 真实拓扑链接：gene → signature（含精确份额）
    gene_links: pd.DataFrame = None
    #: 真实拓扑链接：signature → concept
    signature_links: pd.DataFrame = None
    #: 残差聚合：其余基因 → 各展示 signature（保证守恒）
    residual_gene_links: pd.DataFrame = None
    #: 残差聚合：其余 signature → 各展示 concept（保证守恒）
    residual_signature_links: pd.DataFrame = None
    #: 冻结模型 centering 常数节点 → 各展示 concept（−β·mean/scale）
    centering_links: pd.DataFrame = None
    #: 未展示 concept 的 centering 合计
    centering_other_concepts: float = 0.0
    #: 未展示 concept 的聚合贡献（Σ β_c z_c）
    other_concepts_contribution: float = 0.0
    #: 重建 η 与主 API risk 的最大偏差（自检，用于回归测试）
    eta_max_abs_diff: float = 0.0
    #: high-level 层**predictor 节点**数（concept + 临床/PC/cancer），受
    #: ``max_high_level_nodes`` 约束。这是规格 §7 定义的 "16"。
    n_high_level_nodes: int = 0
    #: 含记账节点（Other concepts 聚合、centering 常数）的实际渲染节点数
    n_high_level_nodes_total: int = 0
    notes: list = field(default_factory=list)
    semantics: str = SEMANTICS

    # ---------------- 导出 ----------------
    def figure(self, *, title: str | None = None, width: int = 1500, height: int = 800,
               font_size: int = 12):
        """返回 Plotly ``Figure``（Sankey + 右侧连续 cohort-relative 风险柱）。"""
        return _build_figure(self, title=title, width=width, height=height,
                             font_size=font_size)

    def save_html(self, path, *, include_plotlyjs: str = "cdn", **kw) -> Path:
        p = Path(path)
        p.parent.mkdir(parents=True, exist_ok=True)
        self.figure(**kw).write_html(p, include_plotlyjs=include_plotlyjs)
        return p

    def save_static(self, path, **kw) -> Path:
        """静态图：优先 plotly+kaleido，缺失时回退 matplotlib（无需额外依赖）。"""
        p = Path(path)
        p.parent.mkdir(parents=True, exist_ok=True)
        if p.suffix.lower() in (".png", ".svg", ".pdf", ".jpg", ".jpeg", ".webp"):
            try:
                self.figure(**kw).write_image(p, scale=2)
                return p
            except Exception:  # noqa: BLE001 —— 缺 kaleido 时回退
                return _save_matplotlib(self, p)
        raise InputError(f"静态图后缀不支持：{p.suffix}")

    def summary(self) -> str:
        L = [f"Sample computation path — {self.sample_id}",
             "─" * 64,
             f"{'Model':<28}{self.model}",
             f"{'Cancer type':<28}{self.cancer_type}",
             f"{'Risk (Cox linear predictor)':<28}{self.risk:.6f}",
             f"{'Rank':<28}{self.rank} / {self.n_samples}",
             f"{'Percentile':<28}{self.percentile:.1f}%",
             f"{'Cutoff':<28}{self.cutoff:.6f} ({self.cutoff_mode})",
             "",
             f"{'Shown genes':<28}{len(self.genes)}",
             f"{'Shown signatures':<28}{len(self.signatures)}",
             f"{'Shown concepts':<28}{len(self.concepts)}",
             f"{'Auxiliary predictors':<28}{len(self.auxiliaries)}",
             f"{'High-level predictor nodes':<28}{self.n_high_level_nodes} "
             f"(limit {DEFAULT_MAX_HIGH_LEVEL})",
             f"{'  incl. bookkeeping nodes':<28}{self.n_high_level_nodes_total}",
             f"{'Other concepts (aggregated)':<28}{self.other_concepts_contribution:+.6f}",
             f"{'Rebuild vs API risk max|Δ|':<28}{self.eta_max_abs_diff:.3e}",
             "", self.semantics]
        if self.notes:
            L += ["", "NOTES"] + ["  - " + n for n in self.notes]
        return "\n".join(L)

    def __repr__(self) -> str:  # pragma: no cover
        return (f"SamplePathResult(sample_id={self.sample_id!r}, model={self.model!r}, "
                f"risk={self.risk:.4f}, rank={self.rank}/{self.n_samples})")


# --------------------------------------------------------------------------- #
# 主函数
# --------------------------------------------------------------------------- #
def sample_path(expression: pd.DataFrame, cancer_type: Sequence[str] | str, sample_id: str, *,
                model: str = "M2", clinical: pd.DataFrame | None = None,
                input_scale: str = "tpm", missing_gene_strategy: str = "reference",
                top_genes: int = DEFAULT_TOP, top_signatures: int = DEFAULT_TOP,
                top_concepts: int | None = None,
                max_high_level_nodes: int = DEFAULT_MAX_HIGH_LEVEL,
                show_cancer_type: bool | None = None,
                cutoff: str | float = "median",
                result=None, batch_size: int = 16) -> SamplePathResult:
    """把单个样本的计算路径归因成可绘制的分层结构。

    参数
    ----
    expression, cancer_type, clinical, input_scale, missing_gene_strategy
        与 :func:`compass_os.predict` 完全一致（同一套对齐/填补逻辑）。
    sample_id
        必须出现在 ``expression.index`` 中；否则抛 :class:`InputError`。
    model
        ``M1`` / ``M2`` / ``M3``。``M0`` 无 COMPASS concept 分支，明确不支持。
    top_genes / top_signatures / top_concepts
        前三层各自的展示上限（默认 50）。node 颜色与**筛选依据分离**：
        筛选按"对该样本所展示 concept 的精确贡献"，颜色才是表达/表征取值。
    max_high_level_nodes
        high-level 层（concept + 辅助 predictor）**总节点上限**，默认 16。
    show_cancer_type
        ``None`` = 自动：单癌种队列省略（并在 notes 说明），多癌种队列显示。
    cutoff
        ``"median"`` = **当前 cohort**risk 中位数（默认，仅作可视化/分组便利）；
        数值 = 研究者给定的 risk 阈值。frozen validated cutoff 不是本参数，
        见 ``cutoff_is_frozen_validated``（恒为 False）。
    result
        可传入已有的 ``PredictionResult`` / ``AnalysisResult`` 以复用 cohort 风险，
        避免重复预测。

    返回
    ----
    :class:`SamplePathResult`
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
    # ---- cohort 层：风险 / rank / percentile / cutoff（复用主 API）----
    if result is not None:
        risk_series = _risk_from_result(result, model)
        concept_all = getattr(result, "concept_scores", None)
    else:
        import compass_os
        pred = compass_os.predict(expression, ct, clinical=clinical, model=model,
                                  missing_gene_strategy=missing_gene_strategy,
                                  input_scale=input_scale, batch_size=batch_size)
        risk_series = pred.risk[model]
        concept_all = pred.concept_scores
    risk_series = pd.Series(risk_series).reindex(idx).astype(float)
    if not np.isfinite(risk_series.to_numpy()).all():
        raise CompassOSError("cohort risk 含非有限值；请检查输入")

    risks = risk_series.to_numpy()
    risk_i = float(risk_series.loc[sample_id])
    order = np.argsort(-risks, kind="stable")
    rank = int(np.where(order == idx.index(sample_id))[0][0]) + 1
    percentile = float(100.0 * (risks < risk_i).mean() + 50.0 * (risks == risk_i).mean())

    if isinstance(cutoff, str):
        if cutoff != "median":
            raise InputError("cutoff 只能是 'median' 或数值阈值")
        cutoff_value, cutoff_mode = float(np.median(risks)), "cohort median (visualization only)"
    else:
        cutoff_value, cutoff_mode = float(cutoff), "user-specified risk threshold"

    # ---- 样本层：官方 gene/132/43 ----
    aligned = align_expression(expression, missing_gene_strategy, input_scale=input_scale)
    codes = cancer_codes_for(ct).to_numpy(int)
    srow = [idx.index(sample_id)]
    layers = _extract_gene_level(aligned.matrix.iloc[srow], codes[srow], batch_size,
                                 missing_gene_strategy,
                                 None if aligned.missing_mask is None
                                 else aligned.missing_mask.iloc[srow])
    gene_score = layers["gene"].iloc[0]
    sig_score = layers["signature"].iloc[0]
    con_score = layers["concept"].iloc[0]

    # ---- 逐特征 Cox 分解（与主 API 同一套函数）----
    lock = _surv.load_lock(model)
    feats = list(lock["feature_names"])
    beta = np.asarray(lock["beta"], np.float64)
    scaled = list(lock["scaled_features"] or [])
    mu = np.asarray(lock["scaler_mean"] or np.zeros(len(scaled)), np.float64)
    sd = np.asarray(lock["scaler_scale"] or np.ones(len(scaled)), np.float64)

    con_row = pd.DataFrame([con_score.to_numpy(np.float64)], index=[sample_id],
                           columns=[str(c) for c in con_score.index])
    clin_row = None if clinical is None else clinical.reindex(idx).loc[[sample_id]]
    Xd = _surv.build_design_matrix(con_row, clin_row, aligned.log2.iloc[srow], model, [ct[srow[0]]])
    Z = Xd.copy()
    if scaled:
        Z[scaled] = (Xd[scaled] - mu) / sd
    contrib = Z.to_numpy(np.float64).ravel() * beta
    eta_rebuilt = float(contrib.sum())
    # 归因链条与主 API 各自做了一次前向（float32 下 batch 相关的 BLAS 求和顺序不同，
    # 概念层实测差 ~8e-07）。展示用的 risk / rank / percentile 一律取**主 API**
    # 的权威结果；这里只校验归因链条与它一致，容差按 float32 精度取 1e-5。
    eta_diff = abs(eta_rebuilt - risk_i)
    eta_tol = ETA_TOL_REL * max(1.0, abs(risk_i))
    if eta_diff > eta_tol:
        raise CompassOSError(
            f"内部一致性检查失败：归因链条求和 η={eta_rebuilt:.9f} 与主 API risk="
            f"{risk_i:.9f} 相差 {eta_diff:.3e}（容差 {eta_tol:.3e}）；"
            "拒绝用不一致的数值绘图")
    if eta_diff > 1e-5:
        notes.append(
            f"Attribution chain and the cohort prediction ran as two separate float32 forward "
            f"passes; their linear predictors differ by {eta_diff:.2e}. The displayed risk, "
            "rank and percentile are taken from the main API output.")

    dec = pd.DataFrame({"feature": feats, "value": Z.to_numpy(np.float64).ravel(),
                        "beta": beta, "contribution": contrib})

    def _contrib_of(name: str) -> float:
        r = dec[dec.feature == name]
        return float(r["contribution"].iloc[0]) if len(r) else 0.0

    # ---- 概念排序：|x_j β_j|（样本真实 Cox 贡献），不是 concept score ----
    con_names = [str(c) for c in con_score.index]
    con_contrib = np.array([_contrib_of(c) for c in con_names])
    con_vals = con_score.to_numpy(np.float64)

    # ---- 辅助 predictor（直接连 risk，不与上游相连）----
    aux_rows = []
    for nm in ("Age", "Sex", "Stage"):
        if nm in feats:
            aux_rows.append({"predictor": nm, "kind": "clinical",
                             "value": float(Z[nm].iloc[0]),
                             "contribution": _contrib_of(nm)})
    pc_cols = [f for f in feats if f.startswith("PC") and f[2:].isdigit()]
    if pc_cols:
        pc_detail = {c: _contrib_of(c) for c in pc_cols}
        aux_rows.append({"predictor": "PC1–PC10", "kind": "pca_aggregate",
                         "value": np.nan, "contribution": float(sum(pc_detail.values())),
                         "detail": pc_detail})

    ct_cols = [f for f in feats if f.startswith("CT_")]
    n_ct = len(set(ct))
    if show_cancer_type is None:
        show_ct = n_ct > 1
    else:
        show_ct = bool(show_cancer_type)
    if ct_cols and show_ct:
        ct_contrib = float(sum(_contrib_of(c) for c in ct_cols))
        aux_rows.append({"predictor": f"Cancer type ({ct[srow[0]]})", "kind": "cancer_type",
                         "value": np.nan, "contribution": ct_contrib})
    elif ct_cols:
        notes.append(
            "Cancer type is part of the full frozen Cox model but is constant within this "
            "single-cancer cohort and is omitted from the visualization for readability.")

    # ---- high-level 预算：concept + 辅助 ≤ max_high_level_nodes ----
    n_aux = len(aux_rows)
    budget = max(0, int(max_high_level_nodes) - n_aux)
    n_show = budget if top_concepts is None else min(int(top_concepts), budget)
    con_order = np.argsort(-np.abs(con_contrib), kind="stable")
    shown_idx = con_order[:n_show]
    hidden_idx = con_order[n_show:]
    other_contrib = float(con_contrib[hidden_idx].sum()) if len(hidden_idx) else 0.0
    if len(hidden_idx):
        notes.append(
            f"{len(hidden_idx)} concept(s) outside the top {n_show} are aggregated into a "
            "single 'Other concepts' node so the inflow to the risk node still sums to the "
            "full Cox linear predictor. Its upstream layers are intentionally not drawn.")

    # ---- 前三层：按"对所展示 concept 的精确贡献"筛选，并建立**真实成员拓扑** ----
    # 记 L_j = Σ_{c∈展示} (β_c/scale_c)·b_{j,c}：gene set j 通向所展示 concept 的系数。
    # 于是：
    #   链接 gene g → signature j 的值 = L_j · a_{g,j} · gene_score_g
    #   链接 signature j → concept c 的值 = (β_c/scale_c) · b_{j,c} · geneset_score_j
    #   链接 concept c → risk        的值 = (β_c/scale_c) · concept_c
    # 三式在**每个节点上都严格守恒**（见下方守恒自检）。
    model_obj = load_model()
    lp = model_obj.model.latentprojector
    gs_names = [str(x) for x in lp.genesetprojector.genesets_names]
    gene_names = [str(g) for g in model_obj.feature_name]
    cp_idx = [list(map(int, x)) for x in lp.cellpathwayprojector.cellpathway_indices]

    A, B, gs_idx, _cp_idx, _n = _attention_matrices(
        model_obj, sorted({j for k in shown_idx for j in cp_idx[k]}))

    c_scale = np.ones(len(con_names))
    for k, cname in enumerate(con_names):
        if cname in scaled:
            s = float(sd[scaled.index(cname)])
            c_scale[k] = s if s != 0 else 1.0
    c_beta = np.array([beta[feats.index(c)] if c in feats else 0.0 for c in con_names])

    L = {}
    for k in shown_idx:
        idxs, w = B[k]
        coef = c_beta[k] / c_scale[k]
        for j, bj in zip(idxs, w):
            L[int(j)] = L.get(int(j), 0.0) + coef * float(bj)

    sig_vals = sig_score.to_numpy(np.float64)
    sig_ids = [str(x) for x in sig_score.index]
    sig_contrib = np.array([sig_vals[j] * L.get(j, 0.0) for j in range(len(sig_ids))])
    keep_sig = np.argsort(-np.abs(sig_contrib), kind="stable")[:int(top_signatures)]
    keep_sig = [int(j) for j in keep_sig if abs(sig_contrib[j]) > 0]

    # 基因：Σ_{j∈展示} L_j·a_{g,j}·gene_score_g  →  gene_score_g × K_g
    K = np.zeros(len(gene_names))
    for j in keep_sig:
        gidx, aw = A[j]
        np.add.at(K, gidx, L.get(j, 0.0) * aw)
    gvals = gene_score.to_numpy(np.float64)
    gene_contrib = gvals * K
    keep_gene = [int(i) for i in np.argsort(-np.abs(gene_contrib), kind="stable")
                 [:int(top_genes)] if abs(gene_contrib[i]) > 0]
    keep_gene_set = set(keep_gene)

    # gene → signature 链接（仅保留两端都被展示、且份额非零的边）
    gl_rows = []
    shown_in = {int(j): 0.0 for j in keep_sig}
    for j in keep_sig:
        gidx, aw = A[j]
        lj = L.get(j, 0.0)
        if lj == 0.0:
            continue
        for g, a in zip(gidx, aw):
            g = int(g)
            v = lj * float(a) * float(gvals[g])
            if g in keep_gene_set:
                shown_in[int(j)] += v
                if v != 0:
                    gl_rows.append({"gene": gene_names[g], "signature": gs_names[j],
                                    "value": v})
    gene_links = pd.DataFrame(gl_rows, columns=["gene", "signature", "value"])

    # signature → concept 链接（真实成员关系）
    sl_rows = []
    shown_set = set(int(k) for k in shown_idx)
    for j in keep_sig:
        for c in shown_set:
            idxs = cp_idx[c]
            if j not in idxs:
                continue
            bjc = float(B[c][1][idxs.index(j)])
            v = (c_beta[c] / c_scale[c]) * bjc * float(sig_vals[j])
            if v != 0:
                sl_rows.append({"signature": gs_names[j], "concept": con_names[c], "value": v})
    signature_links = pd.DataFrame(sl_rows, columns=["signature", "concept", "value"])

    # ---- 残差聚合节点：让**每个**内部节点严格守恒 ----
    # 基因层被截到 top-K、signature 层同理；若不补残差，被截断的节点上"入 < 出"，
    # 视觉上会被误读为丢失。补两个残差源节点后：
    #   signature 入 = 展示基因 + 其余基因 = L_j·geneset_j = 出
    #   concept   入 = 展示 signature + 其余 signature = (β_c/scale_c)·concept_c = 出
    residual_gene_rows = []      # "Other genes (not shown)" → signature j
    for j in keep_sig:
        full = L.get(int(j), 0.0) * float(sig_vals[j])
        resid = full - shown_in[int(j)]
        if abs(resid) > 0:
            residual_gene_rows.append({"signature": gs_names[int(j)], "value": resid})
    # concept 的"未展示 signature"残差：入 = 展示 signature + 残差 = 未中心化的
    # (β_c/scale_c)·concept_c。**centering 常数**（−β_c·mean_c/scale_c）另由专门的
    # 常数节点提供，从而 concept 节点的出 = 中心化后的真实 Cox 贡献（所有节点守恒）。
    residual_sig_rows = []
    centering_rows = []
    for c in shown_set:
        uncentered = (c_beta[c] / c_scale[c]) * float(con_vals[c])
        shown = float(signature_links[signature_links.concept == con_names[c]]["value"].sum()) \
            if len(signature_links) else 0.0
        resid = uncentered - shown
        if abs(resid) > 0:
            residual_sig_rows.append({"concept": con_names[c], "value": resid})
        ctr = _concept_centering(model, con_names[c])
        if abs(ctr) > 0:
            centering_rows.append({"target": con_names[c], "value": ctr})
    centering_other = sum(_concept_centering(model, con_names[k]) for k in hidden_idx)
    # 覆盖率必须用**绝对值**份额：贡献有正负，符号求和可能相互抵消而 >100%。
    abs_shown = abs_full = 0.0
    for j in keep_sig:
        gidx, aw = A[j]
        lj = L.get(int(j), 0.0)
        for g, a in zip(gidx, aw):
            term = abs(lj * float(a) * float(gvals[int(g)]))
            abs_full += term
            if int(g) in keep_gene_set:
                abs_shown += term
    gene_flow_coverage = (abs_shown / abs_full) if abs_full > 0 else 1.0

    # ---- 守恒自检：**每个内部节点入 == 出**（含残差聚合节点）----
    for j in keep_sig:
        nm = gs_names[int(j)]
        inflow = float(gene_links[gene_links.signature == nm]["value"].sum()) + \
            next((r["value"] for r in residual_gene_rows if r["signature"] == nm), 0.0)
        outflow = float(signature_links[signature_links.signature == nm]["value"].sum())
        if abs(inflow - outflow) > 1e-6 * max(1.0, abs(outflow)):
            raise CompassOSError(
                f"gene-set {nm} 流量不守恒：入 {inflow:.6e} vs 出 {outflow:.6e}")
    for c in shown_set:
        nm = con_names[c]
        inflow = float(signature_links[signature_links.concept == nm]["value"].sum()) + \
            next((r["value"] for r in residual_sig_rows if r["concept"] == nm), 0.0) + \
            next((r["value"] for r in centering_rows if r["target"] == nm), 0.0)
        outflow = float(con_contrib[c])
        if abs(inflow - outflow) > 1e-6 * max(1.0, abs(outflow)):
            raise CompassOSError(
                f"concept {nm} 流量不守恒：入 {inflow:.6e} vs 出 {outflow:.6e}")

    notes.append(
        f"Gene layer shows the top {len(keep_gene)} contributors by |contribution| "
        f"({gene_flow_coverage:.1%} of the displayed gene-set inflow in absolute terms); "
        "the remainder is "
        "carried by an explicit 'Other genes (not shown)' node so that every node conserves "
        "flow. Signature and concept layers use the same convention.")
    residual_gene = pd.DataFrame(residual_gene_rows, columns=["signature", "value"])
    residual_sig = pd.DataFrame(residual_sig_rows, columns=["concept", "value"])
    centering_df = pd.DataFrame(centering_rows, columns=["target", "value"])

    genes_df = pd.DataFrame({
        "gene": [gene_names[i] for i in keep_gene],
        "expression": [float(aligned.matrix.iloc[srow[0]][gene_names[i]])
                       if gene_names[i] in aligned.matrix.columns else np.nan
                       for i in keep_gene],
        "gene_score": [float(gvals[i]) for i in keep_gene],
        "contribution": [float(gene_contrib[i]) for i in keep_gene]})
    sig_df = pd.DataFrame({
        "signature": [sig_ids[j] for j in keep_sig],
        "score": [float(sig_vals[j]) for j in keep_sig],
        "contribution": [float(sig_contrib[j]) for j in keep_sig]})
    con_df = pd.DataFrame({
        "concept": [con_names[k] for k in shown_idx],
        "score": [float(con_vals[k]) for k in shown_idx],
        "contribution": [float(con_contrib[k]) for k in shown_idx]})
    aux_cols = ["predictor", "kind", "value", "contribution", "detail"]
    aux_df = pd.DataFrame(aux_rows, columns=aux_cols)
    if "detail" not in aux_df.columns:
        aux_df["detail"] = None
    aux_df["detail"] = [d if isinstance(d, dict) else None for d in aux_df["detail"]]

    res = SamplePathResult(
        sample_id=sample_id, model=model, cancer_type=ct[srow[0]],
        risk=risk_i, rank=rank, n_samples=len(idx), percentile=percentile,
        cutoff=cutoff_value, cutoff_mode=cutoff_mode, cutoff_is_frozen_validated=False,
        decomposition=dec, genes=genes_df, signatures=sig_df, concepts=con_df,
        auxiliaries=aux_df, gene_links=gene_links, signature_links=signature_links,
        residual_gene_links=residual_gene, residual_signature_links=residual_sig,
        centering_links=centering_df, centering_other_concepts=centering_other,
        other_concepts_contribution=other_contrib,
        eta_max_abs_diff=eta_diff,
        n_high_level_nodes=int(len(con_df) + len(aux_df)),
        n_high_level_nodes_total=int(len(con_df) + len(aux_df)
                                     + (1 if len(hidden_idx) else 0)
                                     + (1 if len(centering_df) else 0)),
        notes=notes)
    _assert_budget(res, max_high_level_nodes)
    return res


# --------------------------------------------------------------------------- #
def _concept_centering(model: str, concept: str) -> float:
    """concept 的冻结 centering 常数项 ``−β_c·mean_c/scale_c``（未在 scaled_features 中则为 0）。"""
    lock = _surv.load_lock(model)
    feats = list(lock["feature_names"])
    if concept not in feats:
        return 0.0
    scaled = list(lock["scaled_features"] or [])
    if concept not in scaled:
        return 0.0
    k = scaled.index(concept)
    beta = float(np.asarray(lock["beta"], np.float64)[feats.index(concept)])
    mu = float(np.asarray(lock["scaler_mean"], np.float64)[k])
    sd = float(np.asarray(lock["scaler_scale"], np.float64)[k])
    return -beta * mu / sd if sd != 0 else 0.0


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


def _assert_budget(res: SamplePathResult, limit: int) -> None:
    """规格 §7：high-level 层 **predictor 节点**（concept + 临床/PC/cancer）≤ limit。

    记账节点（"Other concepts" 聚合、冻结模型 centering 常数）不占用该预算，
    但会在 ``n_high_level_nodes_total`` 中单独报告，并在图注中说明其性质。
    """
    if res.n_high_level_nodes > limit:
        raise CompassOSError(
            f"high-level 层 predictor 节点数 {res.n_high_level_nodes} 超过上限 {limit}")


# --------------------------------------------------------------------------- #
# Plotly 渲染
# --------------------------------------------------------------------------- #
def _build_figure(res: SamplePathResult, *, title=None, width=1500, height=800, font_size=12):
    """Sankey（真实成员拓扑）+ 右侧连续 cohort-relative 风险柱。

    ⚠ Plotly 陷阱：Sankey trace 没有 x/y 轴，``add_hline`` / ``add_vline`` 这类
    跨轴辅助函数会抛 ``PlotlyKeyError: Invalid property ... 'xaxis'``。
    因此 cutoff 线一律用 ``fig.add_shape(..., xref="x2", yref="y2")`` 直接写入 layout。
    """
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    labels: list[str] = []
    colors: list[str] = []
    node_x: list[float] = []
    node_y: list[float] = []
    src: list[int] = []
    dst: list[int] = []
    val: list[float] = []
    lcol: list[str] = []
    htxt: list[str] = []

    def add(label, color, x):
        labels.append(label); colors.append(color); node_x.append(x); node_y.append(0.0)
        return len(labels) - 1

    def link(a, b, v, tag):
        if not np.isfinite(v) or v == 0:
            return
        src.append(a); dst.append(b); val.append(abs(float(v)))
        lcol.append(_POS if v > 0 else _NEG)
        htxt.append(f"{tag}<br>contribution = {v:+.6f}")

    def ypos(k, n):
        return 1.0 - (k + 0.5) / max(n, 1)

    gl = res.gene_links if res.gene_links is not None else pd.DataFrame(
        columns=["gene", "signature", "value"])
    sl = res.signature_links if res.signature_links is not None else pd.DataFrame(
        columns=["signature", "concept", "value"])

    gene_ids = {r.gene: add(r.gene, _GENE_C, 0.01) for r in res.genes.itertuples()}
    sig_ids = {r.signature: add(r.signature, "#7fb3d5", 0.34) for r in res.signatures.itertuples()}
    con_ids = {r.concept: add(r.concept, _CONCEPT_C, 0.64) for r in res.concepts.itertuples()}
    aux_ids = []
    for r in res.auxiliaries.itertuples():
        extra = ""
        d = getattr(r, "detail", None)
        if isinstance(d, dict) and d:
            extra = "<br>" + "<br>".join(f"{kk} = {vv:+.6f}" for kk, vv in d.items())
        aux_ids.append((add(r.predictor, _AUX_C, 0.86), r, extra))
    rgl = res.residual_gene_links if res.residual_gene_links is not None else \
        pd.DataFrame(columns=["signature", "value"])
    rsl = res.residual_signature_links if res.residual_signature_links is not None else \
        pd.DataFrame(columns=["concept", "value"])
    other_genes_id = add("Other genes (not shown)", _OTHER_C, 0.01) if len(rgl) else None
    other_sigs_id = add("Other signatures (not shown)", _OTHER_C, 0.34) if len(rsl) else None
    ctr = res.centering_links if res.centering_links is not None else \
        pd.DataFrame(columns=["target", "value"])
    ctr_id = add("Frozen model centering (constant)", "#aab7b8", 0.64) if len(ctr) else None
    other_id = None
    if res.other_concepts_contribution != 0:
        other_id = add("Other concepts (aggregated; upstream not drawn)", _OTHER_C, 0.86)
    risk_id = add("Cox risk (linear predictor)", _RISK_C, 1.0)

    for name, i in gene_ids.items():
        r = res.genes[res.genes.gene == name].iloc[0]
        node_y[i] = ypos(list(gene_ids).index(name), len(gene_ids))
    for name, i in sig_ids.items():
        node_y[i] = ypos(list(sig_ids).index(name), len(sig_ids))
    for name, i in con_ids.items():
        node_y[i] = ypos(list(con_ids).index(name), len(con_ids))
    n_aux = len(aux_ids) + (1 if other_id is not None else 0)
    for k, (i, _r, _e) in enumerate(aux_ids):
        node_y[i] = ypos(k, n_aux)
    if other_id is not None:
        node_y[other_id] = ypos(len(aux_ids), n_aux)
    node_y[risk_id] = 0.5

    for r in gl.itertuples():
        if r.gene in gene_ids and r.signature in sig_ids:
            link(gene_ids[r.gene], sig_ids[r.signature], r.value,
                 f"{r.gene} → {r.signature}")
    for r in sl.itertuples():
        if r.signature in sig_ids and r.concept in con_ids:
            link(sig_ids[r.signature], con_ids[r.concept], r.value,
                 f"{r.signature} → {r.concept}")
    if other_genes_id is not None:
        for r in rgl.itertuples():
            if r.signature in sig_ids:
                link(other_genes_id, sig_ids[r.signature], r.value,
                     f"other genes → {r.signature} (residual)")
    if other_sigs_id is not None:
        for r in rsl.itertuples():
            if r.concept in con_ids:
                link(other_sigs_id, con_ids[r.concept], r.value,
                     f"other signatures → {r.concept} (residual)")
    if ctr_id is not None:
        for r in ctr.itertuples():
            if r.target in con_ids:
                link(ctr_id, con_ids[r.target], r.value,
                     f"model centering → {r.target}")
    for r in res.concepts.itertuples():
        link(con_ids[r.concept], risk_id, r.contribution, f"concept {r.concept} → risk")
    for i, r, extra in aux_ids:
        link(i, risk_id, r.contribution, f"{r.predictor} → risk{extra}")
    if other_id is not None:
        link(other_id, risk_id, res.other_concepts_contribution, "Other concepts → risk")

    fig = make_subplots(rows=1, cols=2, column_widths=[0.87, 0.13], horizontal_spacing=0.01,
                        specs=[[{"type": "sankey"}, {"type": "xy"}]])
    fig.add_trace(go.Sankey(
        arrangement="snap",
        node=dict(pad=7, thickness=13, line=dict(color="rgba(0,0,0,0.25)", width=0.5),
                  label=labels, color=colors, x=node_x, y=node_y,
                  hovertemplate="%{label}<extra></extra>"),
        link=dict(source=src, target=dst, value=val, color=lcol, customdata=htxt,
                  hovertemplate="%{customdata}<extra></extra>"),
    ), row=1, col=1)

    # ---- 右侧：连续 cohort-relative 风险柱（红=高 → 绿=低；不是绝对死亡概率）----
    fig.add_trace(go.Heatmap(z=np.linspace(0, 1, 256).reshape(-1, 1),
                             colorscale=[[0, "#c0392b"], [0.5, "#f7dc6f"], [1, "#1e8449"]],
                             showscale=False, zmin=0, zmax=1, hoverinfo="skip"), row=1, col=2)
    lo = min(float(res.decomposition["contribution"].sum()), res.risk)
    hi = max(lo + 1e-12, res.risk)
    frac = float(np.clip((res.risk - lo) / (hi - lo + 1e-12), 0.0, 1.0))
    fig.add_trace(go.Scatter(
        x=[0.5], y=[frac], mode="markers",
        marker=dict(symbol="triangle-left", size=15, color="black"),
        hovertemplate=(f"sample {res.sample_id}<br>Risk = {res.risk:.6f}"
                       f"<br>Rank = {res.rank} / {res.n_samples}"
                       f"<br>Percentile = {res.percentile:.1f}%<extra></extra>"),
        showlegend=False), row=1, col=2)
    fig.add_shape(type="line", x0=0, x1=1, y0=frac, y1=frac, xref="x2", yref="y2",
                  line=dict(color="rgba(0,0,0,0)", width=0))
    fig.update_xaxes(showticklabels=False, showgrid=False, zeroline=False, row=1, col=2)
    fig.update_yaxes(showticklabels=False, showgrid=False, zeroline=False,
                     range=[-0.02, 1.02], row=1, col=2)

    cap = ("Cutoff is the <b>cohort median of the current risk output</b> — a visualization "
           "convenience, not the frozen validated prognostic cutoff."
           if str(res.cutoff_mode).startswith("cohort median") else
           "Cutoff is a <b>user-specified risk threshold</b>, not the frozen validated "
           "prognostic cutoff.")
    ttl = title or f"Sample computation path — {res.sample_id} ({res.model})"
    fig.update_layout(
        title=dict(text=ttl, font=dict(size=16), x=0.01, xanchor="left"),
        font=dict(size=font_size), width=width, height=height,
        margin=dict(l=10, r=150, t=110, b=110), paper_bgcolor="white",
        annotations=[
            dict(text=(f"Risk = {res.risk:.4f} &nbsp;|&nbsp; Rank = {res.rank} / "
                       f"{res.n_samples} &nbsp;|&nbsp; Percentile = {res.percentile:.1f}% "
                       f"&nbsp;|&nbsp; Cutoff = {res.cutoff:.4f} ({res.cutoff_mode})"),
                 xref="paper", yref="paper", x=0.01, y=1.045, xanchor="left",
                 showarrow=False, font=dict(size=11, color="#333")),
            dict(text=cap, xref="paper", yref="paper", x=0.01, y=1.012, xanchor="left",
                 showarrow=False, font=dict(size=10, color="#777")),
            dict(text=("<b>Link width = |contribution|</b> &nbsp; "
                       "<span style='color:#c0392b'>■ positive → higher model-estimated "
                       "risk</span> &nbsp; <span style='color:#2471a3'>■ negative → lower"
                       "</span>"),
                 xref="paper", yref="paper", x=0.01, y=-0.035, xanchor="left",
                 showarrow=False, font=dict(size=10)),
            dict(text=SEMANTICS, xref="paper", yref="paper", x=0.01, y=-0.075,
                 xanchor="left", yanchor="top", showarrow=False,
                 font=dict(size=9, color="#666")),
        ])
    return fig


# --------------------------------------------------------------------------- #
def _save_matplotlib(res: SamplePathResult, path: Path) -> Path:
    """无 kaleido 时的静态回退：用 matplotlib 画同一份归因数据的层叠条形图。"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(15, 7), gridspec_kw={"width_ratios": [1, 1, 1]})
    panels = [("Top genes (contribution)", res.genes, "gene"),
              ("Top gene-signature representations", res.signatures, "signature"),
              ("High-level concepts & predictors", res.concepts, "concept")]
    for ax, (ttl, df, col) in zip(axes, panels):
        if df.empty:
            ax.set_visible(False)
            continue
        d = df.iloc[::-1]
        cols = [_POS if v > 0 else _NEG for v in d["contribution"]]
        ax.barh(d[col].astype(str), d["contribution"], color=cols)
        ax.set_title(ttl, fontsize=10)
        ax.tick_params(labelsize=7)
        ax.axvline(0, color="black", lw=0.6)
        ax.set_xlabel("Cox contribution", fontsize=8)
    if len(res.auxiliaries):
        ax = axes[2]
        for _i, r in res.auxiliaries.iterrows():
            lab = f"[aux] {r['predictor']}"
            ax.barh(lab, r["contribution"],
                    color=_POS if r["contribution"] > 0 else _NEG)
        ax.tick_params(labelsize=7)
    fig.suptitle(f"Sample computation path — {res.sample_id} ({res.model})  |  "
                 f"Risk = {res.risk:.4f}  Rank = {res.rank}/{res.n_samples}  "
                 f"Percentile = {res.percentile:.1f}%", fontsize=11)
    fig.text(0.01, 0.01, SEMANTICS, fontsize=7, color="#555", wrap=True)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig.savefig(path, dpi=170)
    plt.close(fig)
    return path
