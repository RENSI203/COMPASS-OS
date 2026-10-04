# -*- coding: utf-8 -*-
"""内部桥接：把 :class:`~compass_os.sample_path.SamplePathResult` 的真实计算结果
转成渲染层的数据契约 :class:`~compass_os.sample_path_render.SamplePathData`。

**普通用户不需要接触本模块**；它只负责"计算 → 绘图"的解耦。

数据口径（全部来自真实推理）：

* ``gene_score`` —— 官方 ``extract(with_gene_level=True)`` 的 ``dfg``
  （``geneset_scorer(encoder_output)[:, 2:]``，共享 ``nn.Linear(32→1)``）；
* ``signatures`` / ``concepts`` —— 主 API **同一次前向**的输出（与 η 严格同源）；
* 临床 / Cancer type / PC1–PC10 —— 冻结 Cox 设计矩阵逐列 signed contribution 的聚合；
* ``decomposition`` —— 完整未截断分解，``sum(contribution) == η``（渲染层强制校验）。
"""
from __future__ import annotations

from .exceptions import InputError

__all__ = ["to_render_data"]


def to_render_data(res):
    """构造渲染层 :class:`SamplePathData`（惰性导入，避免顶层硬依赖）。"""
    from .sample_path_render import Edge, Node, SamplePathData, aggregate_cox_predictor

    p: dict = res.render_payload or {}
    if not p:
        raise InputError("SamplePathResult 缺少渲染载荷；请通过 compass_os.sample_path() 生成")

    model = res.model
    dec = {r["feature"]: float(r["contribution"]) for r in p["decomposition"]}

    # 第一列：真实表达值（单位见 expression_units）
    gene_expr = [Node(name=g, score=float(v), kind="gene")
                 for g, v in p["gene_expression"].items()]
    # 第二列：真实 COMPASS gene-token score
    gene_score = [Node(name=g, score=float(v), kind="gene")
                  for g, v in p["gene_score"].items()]
    # 第三列：真实 132 signature 表示
    signatures = [Node(name=s, score=float(v), kind="signature")
                  for s, v in p["signature_scores"].items()]

    # 第四列：43 concepts —— 颜色 = 真实 concept 分数；筛选 = |signed contribution|
    predictors: list = [Node(name=c, score=float(sc), kind="concept",
                             contribution=float(p["concept_contributions"][c]))
                        for c, sc in p["concept_scores"].items()]
    # 临床：score == signed contribution（渲染层强制），并逐字段标注真实 imputed 状态
    for nm, info in p["clinical"].items():
        predictors.append(Node(
            name=nm, score=float(info["contribution"]), kind="clinical",
            contribution=float(info["contribution"]), value=info["value"],
            imputed=bool(info["imputed"]),
            details={"frozen_fill": info["value"],
                     "imputed": ("yes — frozen training reference value"
                                 if info["imputed"] else "no")}))
    # Cancer type：多编码列 signed sum；仅多癌种队列展示（渲染层按 cohort 判定）
    ct_cols = p["cancer_columns"]
    if ct_cols and len({str(c) for c in p["cancer_types"]}) > 1:
        predictors.append(aggregate_cox_predictor(
            "Cancer type", ct_cols, {c: dec.get(c, 0.0) for c in ct_cols},
            kind="cancer", value=res.cancer_type))
    # PC1–PC10：单一节点，值为十个冻结设计列贡献的 signed sum
    pc_cols = p["pc_columns"]
    if pc_cols:
        predictors.append(aggregate_cox_predictor(
            "PC1–PC10", pc_cols, {c: dec.get(c, 0.0) for c in pc_cols}, kind="pc"))

    gene_edges = [Edge(e["source"], e["target"], float(e["weight"]))
                  for e in p["gene_signature_edges"]]
    # signature→concept：只保留真实拓扑（gene set 属于该 concept）
    keep = {(e["source"], e["target"]) for e in p["signature_concept_edges"]}
    sig_edges = [Edge(s, c, float(w)) for s, c, w in
                 ((e["source"], e["target"], e["weight"])
                  for e in p["signature_concept_edges"])]

    multi = len({str(c) for c in p["cancer_types"]}) > 1
    hidden = ("" if (multi or model == "M1") else
              "Cancer type is constant, hidden here and retained in the frozen model. ")
    return SamplePathData(
        sample_id=res.sample_id, model=model,
        cohort=f"{p.get('cohort_label', 'cohort')} · N = {res.n_samples}",
        risk=float(res.risk), rank=float(res.rank), percentile=float(res.percentile),
        cohort_risks=[float(v) for v in p["cohort_risks"]],
        gene_tpm=gene_expr, gene_score=gene_score, signatures=signatures,
        predictors=predictors,
        gene_signature_edges=gene_edges, signature_concept_edges=sig_edges,
        cutoff=None,                       # None ⇒ 队列 η 中位数（仅可视化分界）
        cancer_types=[str(c) for c in p["cancer_types"]],
        expression_label="Gene expression",
        expression_units=p.get("expression_units", "TPM"),
        gene_score_units="COMPASS gene-token score (shared nn.Linear(32->1))",
        signature_units="COMPASS gene-set / signature score",
        concept_units="COMPASS concept score",
        rank_definition=("rank = 1 + #{cohort eta > eta_i}; "
                         "percentile = 100*(P(eta<eta_i) + 0.5*P(eta=eta_i)); "
                         "higher = higher risk"),
        provenance=("REAL COMPASS-OS inference (frozen pretrainer.pt + frozen Cox locks); "
                    "gene score from PreTrainer.extract(with_gene_level=True). " + hidden +
                    "Risk is cohort-relative, not survival probability."),
        # 完整分解：含未展示 concept、被隐藏的 Cancer type 与 offset；sum == η
        decomposition=[{"feature": f, "contribution": float(c)} for f, c in dec.items()],
        metadata={
            "model": model, "n_samples": int(res.n_samples),
            "cutoff_mode": res.cutoff_mode,
            "decomposition_vs_api_risk": float(res.eta_max_abs_diff),
            "high_level_predictor_nodes": int(res.n_high_level_nodes),
            "clinical_imputed": {k: bool(v["imputed"]) for k, v in p["clinical"].items()},
            "pc_aggregation": ("exact signed sum of the frozen PC design-column contributions"
                               if pc_cols else None),
            "input_scale": p.get("input_scale"),
            "missing_gene_strategy": p.get("missing_gene_strategy"),
            "n_missing_genes": p.get("n_missing_genes"),
            "notes": list(res.notes),
        },
    )
