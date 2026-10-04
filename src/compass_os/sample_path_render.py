#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""COMPASS-OS 样本计算路径渲染层（圆形分层图，单一静态渲染管线）。

本模块**只负责绘图**：不做推理、不修改冻结模型、不重算风险。
`numpy` / `matplotlib` 均为**函数内惰性导入**，因此
``import compass_os`` 与核心预测在未安装绘图 extra 时照常可用。

规范：PNG / PDF / SVG 由**同一个 Figure** 导出；HTML 内嵌**同一份 PNG 字节**
并提供可展开的精确节点表——不存在第二套绘图器，也不会在导出失败时静默切换外观。

上游来源：COMPASS_sample_path_redesign 附件（`sankey.py`，
sha256 daa89531e19db4fe154ae1add4713c6753e1d7912f87fccd92eedb0efcfae6b5），
按包内架构整理后并入；`SamplePathData` / `Node` / `Edge` / `Style` /
`plot_sample_path` / `aggregate_cox_predictor` 的语义与校验规则保持不变。

对附件版本的两处包内化改动：
1. `RenderedPath` 增加 ``save_html()`` / ``save_static()``，使 HTML 可**独立导出**
   （内部仍先生成 PNG 字节再内嵌，符合"单一渲染管线"约束）；
2. 模块 docstring/注释中文化，并把 CLI 入口改为可选的 ``main()``。
"""
from __future__ import annotations

import argparse
import base64
import csv
import html
import io
import json
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping, Sequence


@dataclass(frozen=True)
class Node:
    """score colors the node; contribution only ranks high-level concepts.

    Scores must already be expressed in the declared layer color units. No
    cohort z-scoring or normalizing a selected subset happens in this renderer.
    Clinical/PC nodes use their signed Cox contribution as their color score.
    """
    name: str
    score: float
    kind: str = "concept"  # gene / signature / concept / clinical / pc / cancer
    contribution: float | None = None
    value: Any = None
    imputed: bool = False
    details: Mapping[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class Edge:
    source: str
    target: str
    weight: float = 1.0  # selection only; never encoded as width or color


@dataclass
class SamplePathData:
    sample_id: str
    model: str
    cohort: str
    risk: float  # explicitly Cox linear predictor eta
    rank: float  # preserved from API; describe its convention in rank_definition
    percentile: float  # 0..100; higher = higher risk
    cohort_risks: Sequence[float]  # same API output/scale as risk
    gene_tpm: Sequence[Node]
    gene_score: Sequence[Node]
    signatures: Sequence[Node]
    predictors: Sequence[Node]
    gene_signature_edges: Sequence[Edge]
    signature_concept_edges: Sequence[Edge]
    cutoff: float | None = None  # None = cohort median, visualization only
    cancer_types: Sequence[str] = field(default_factory=tuple)
    expression_label: str = "Gene TPM"
    expression_units: str = "TPM"
    gene_score_units: str = "COMPASS gene score"
    signature_units: str = "COMPASS signature score"
    concept_units: str = "COMPASS concept score"
    rank_definition: str = "Rank as returned by main API"
    provenance: str = "Real model output"
    decomposition: Sequence[Mapping[str, Any]] = field(default_factory=tuple)
    metadata: Mapping[str, Any] = field(default_factory=dict)

    @classmethod
    def from_dict(cls, obj: Mapping[str, Any]) -> "SamplePathData":
        d = dict(obj)
        for key in ("gene_tpm", "gene_score", "signatures", "predictors"):
            d[key] = [Node(**n) for n in d[key]]
        for key in ("gene_signature_edges", "signature_concept_edges"):
            d[key] = [Edge(**e) for e in d[key]]
        return cls(**d)


@dataclass(frozen=True)
class Style:
    gene_threshold: float = 0.25  # abs(real gene score), NOT link attribution
    max_genes: int = 32
    max_signatures: int = 26
    predictor_budget: int = 16
    max_gene_edges_per_signature: int | None = 4
    width: float = 17.5  # inches
    height: float = 10.4
    dpi: int = 160
    gene_radius: float = 3.0
    signature_radius: float = 4.2
    predictor_radius: float = 6.0
    line_width: float = 0.45
    line_color: str = "#9bafbf"
    line_alpha: float = 0.32
    expression_limits: tuple[float, float] = (0.0, 20.0)
    gene_score_limits: tuple[float, float] = (-1.0, 1.0)
    signature_limits: tuple[float, float] = (-1.0, 1.0)
    concept_limits: tuple[float, float] = (-1.0, 1.0)
    contribution_limits: tuple[float, float] = (-1.0, 1.0)


def _finite(value: Any, name: str) -> float:
    if isinstance(value, bool):
        raise ValueError(f"{name} cannot be boolean")
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite")
    return value


def _check(data: SamplePathData, style: Style) -> None:
    if data.model not in {"M1", "M2", "M3"}:
        raise ValueError("model must be M1, M2 or M3")
    _finite(data.risk, "risk")
    if not 0 <= _finite(data.percentile, "percentile") <= 100:
        raise ValueError("percentile must be on 0..100 scale")
    if len(data.cohort_risks) == 0:
        raise ValueError("cohort_risks is required for the median cutoff")
    for r in data.cohort_risks:
        _finite(r, "cohort_risk")
    if not 1 <= _finite(data.rank, "rank") <= len(data.cohort_risks):
        raise ValueError("rank must lie within 1..N")
    if data.cutoff is not None:
        _finite(data.cutoff, "cutoff")
    if data.model in {"M2", "M3"} and not data.cancer_types:
        raise ValueError("Provide cohort cancer_types to choose the predictor budget")
    if not 1 <= style.max_genes <= 50 or not 1 <= style.max_signatures <= 50:
        raise ValueError("gene/signature limits must be 1..50")
    if not 1 <= style.predictor_budget <= 16:
        raise ValueError("predictor_budget must be 1..16")
    if style.gene_threshold < 0 or not math.isfinite(style.gene_threshold):
        raise ValueError("gene_threshold must be finite and nonnegative")
    if style.max_gene_edges_per_signature is not None and style.max_gene_edges_per_signature < 1:
        raise ValueError("max_gene_edges_per_signature must be positive or None")
    for limits in (style.expression_limits, style.gene_score_limits,
                   style.signature_limits, style.concept_limits, style.contribution_limits):
        if len(limits) != 2 or not all(math.isfinite(x) for x in limits) or limits[0] >= limits[1]:
            raise ValueError("color limits must be finite increasing pairs")
    for nodes in (data.gene_tpm, data.gene_score, data.signatures, data.predictors):
        names = [n.name for n in nodes]
        if len(names) != len(set(names)):
            raise ValueError("Duplicate node name within a layer")
        for n in nodes:
            if not n.name or n.kind not in {"gene", "signature", "concept", "clinical", "pc", "cancer"}:
                raise ValueError(f"Invalid node: {n.name!r}, {n.kind!r}")
            _finite(n.score, f"{n.name}.score")
            if n.contribution is not None:
                _finite(n.contribution, f"{n.name}.contribution")
            if n.kind in {"clinical", "pc", "cancer"}:
                if n.contribution is None or not math.isclose(n.score, n.contribution, abs_tol=1e-12):
                    raise ValueError(f"{n.name}: clinical/PC color score must equal signed contribution")
    for key, edges, sources, targets in (
        ("gene->signature", data.gene_signature_edges, data.gene_score, data.signatures),
        ("signature->concept", data.signature_concept_edges, data.signatures, data.predictors),
    ):
        source_names, target_names = {n.name for n in sources}, {n.name for n in targets}
        for e in edges:
            if e.source not in source_names or e.target not in target_names:
                raise ValueError(f"Unknown {key} endpoint: {e.source} -> {e.target}")
            _finite(e.weight, "edge.weight")
    if data.decomposition:
        if any("contribution" not in r for r in data.decomposition):
            raise ValueError("Each decomposition row needs contribution")
        total = math.fsum(_finite(r["contribution"], "decomposition") for r in data.decomposition)
        if not math.isclose(total, data.risk, rel_tol=1e-6, abs_tol=1e-7):
            raise ValueError(f"Signed decomposition sums to {total}, not eta={data.risk}")


def _bookkeeping(n: Node) -> bool:
    low = n.name.lower()
    return low.startswith("other ") or any(word in low for word in (
        "centering", "standardization offset", "bookkeeping", "residual", "not drawn", "not shown"))


def aggregate_cox_predictor(name: str, columns: Sequence[str],
                            signed_contributions: Mapping[str, float], *, kind: str,
                            value: Any = None, imputed: bool = False) -> Node:
    """Aggregate EXACT design-column contributions for Sex/Stage/CT/PC nodes.

    signed_contributions must come from the frozen Cox design matrix and beta,
    including its actual standardization. Do not pass raw value * beta unless
    that is exactly how the frozen model computes its linear predictor.
    """
    if not columns or len(set(columns)) != len(columns):
        raise ValueError("Supply nonempty, unique frozen design-matrix columns")
    if kind not in {"clinical", "pc", "cancer"}:
        raise ValueError("Aggregated predictors must have clinical/pc/cancer kind")
    values = {col: _finite(signed_contributions[col], col) for col in columns}
    total = math.fsum(values.values())
    return Node(name, total, kind, total, value, imputed,
                {"design_columns": list(columns), "column_contributions": values,
                 "aggregation": "exact signed sum of frozen Cox design-column contributions"})


def select_nodes(data: SamplePathData, style: Style = Style()) -> dict[str, Any]:
    """Deterministic node/path selection. No residual placeholders are drawn."""
    _check(data, style)
    predictors = [n for n in data.predictors if not _bookkeeping(n)]
    multi = len(set(data.cancer_types)) > 1
    extras = [n for n in predictors if n.kind in {"clinical", "pc", "cancer"}]
    if data.model == "M1":
        extras = []
    else:
        expected = {"Age", "Sex", "Stage"} | ({"PC1–PC10"} if data.model == "M3" else set())
        if multi:
            expected.add("Cancer type")
        extras = [n for n in extras if multi or n.kind != "cancer"]
        if {n.name for n in extras} != expected or len(extras) != len(expected):
            raise ValueError(f"{data.model}: require exactly these aggregated predictors: {sorted(expected)}")
        if any(n.kind != ("pc" if n.name == "PC1–PC10" else "cancer" if n.name == "Cancer type" else "clinical") for n in extras):
            raise ValueError("Incorrect clinical/PC/cancer kind")
    if style.predictor_budget < len(extras):
        raise ValueError("predictor budget is smaller than mandatory clinical/PC nodes")
    concepts = [n for n in predictors if n.kind == "concept"]
    if any(n.contribution is None for n in concepts):
        raise ValueError("Concepts require exact signed Cox contributions for top-node selection")
    concepts = sorted(concepts, key=lambda n: (-abs(n.contribution), n.name))[:style.predictor_budget-len(extras)]
    # Stable clinical order across M2/M3, after the concept nodes.
    extras.sort(key=lambda n: ["Age", "Sex", "Stage", "Cancer type", "PC1–PC10"].index(n.name))
    high = concepts + extras
    names = {n.name for n in concepts}
    se = [e for e in data.signature_concept_edges if e.target in names and e.weight != 0]
    active = {e.source for e in se}
    signatures = sorted((n for n in data.signatures if n.name in active and not _bookkeeping(n)),
                        key=lambda n: (-abs(n.score), n.name))[:style.max_signatures]
    sn = {n.name for n in signatures}
    candidates = {e.source for e in data.gene_signature_edges if e.target in sn and e.weight != 0}
    genes = sorted((n for n in data.gene_score if n.name in candidates and not _bookkeeping(n)
                    and abs(n.score) >= style.gene_threshold),
                   key=lambda n: (-abs(n.score), n.name))[:style.max_genes]
    gn = {n.name for n in genes}
    expressions = {n.name: n for n in data.gene_tpm}
    if gn - expressions.keys():
        raise ValueError(f"Missing expression values for selected genes: {sorted(gn-expressions.keys())}")
    ge = [e for e in data.gene_signature_edges if e.source in gn and e.target in sn and e.weight != 0]
    # Deduplicate topology; edge weights never change line styles.
    ge = list({(e.source, e.target): e for e in ge}.values())
    if style.max_gene_edges_per_signature is not None:
        ge = [e for s in signatures for e in sorted((e for e in ge if e.target == s.name),
              key=lambda e: (-abs(e.weight), e.source))[:style.max_gene_edges_per_signature]]
    # Genes left disconnected by edge thinning are truly removed from BOTH columns.
    connected = {e.source for e in ge}
    genes = [n for n in genes if n.name in connected]
    se = list({(e.source, e.target): e for e in se if e.source in sn}.values())
    return {"gene_tpm": [expressions[n.name] for n in genes], "gene_score": genes,
            "signatures": signatures, "predictors": high, "gene_edges": ge, "signature_edges": se,
            "hidden_cancer_type": data.model != "M1" and not multi,
            "counts": {"gene_tpm": len(genes), "gene_score": len(genes),
                       "signatures": len(signatures), "predictors": len(high)}}


@dataclass
class RenderedPath:
    figure: Any
    data: SamplePathData
    style: Style
    selection: dict[str, Any]
    cutoff: float
    cutoff_percentile: float

    # ------------------------------------------------------------------ 原子导出
    def _save_figure(self, path: Path, fmt: str) -> Path:
        """从**同一** Figure 导出 png/pdf/svg（字体与嵌入设置一致）。"""
        from matplotlib import rc_context
        path.parent.mkdir(parents=True, exist_ok=True)
        with rc_context({"svg.fonttype": "path", "pdf.fonttype": 42,
                         "font.family": "DejaVu Sans"}):
            self.figure.savefig(path, dpi=self.style.dpi, facecolor="white", format=fmt)
        return path

    def _node_rows(self, table_path: Path) -> list[str]:
        """写 nodes.tsv 并返回对应的 HTML 行。"""
        table_path.parent.mkdir(parents=True, exist_ok=True)
        rows = []
        with table_path.open("w", encoding="utf-8", newline="") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerow(["layer", "node", "score", "signed_contribution", "value", "imputed", "details"])
            for layer in ("gene_tpm", "gene_score", "signatures", "predictors"):
                for n in self.selection[layer]:
                    vals = [layer, n.name, n.score, n.contribution, n.value, n.imputed,
                            json.dumps(n.details, ensure_ascii=False)]
                    writer.writerow(vals)
                    rows.append("<tr>" + "".join(f"<td>{html.escape(str(v))}</td>" for v in vals) + "</tr>")
        return rows

    def _html_text(self, png_bytes: bytes, rows: list[str]) -> str:
        """HTML 正文：内嵌**同一份 PNG 字节** + 可展开的精确节点表。"""
        png64 = base64.b64encode(png_bytes).decode("ascii")
        heading = html.escape(f"{self.data.model} · {self.data.sample_id}")
        return f'''<!doctype html><html><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>{heading}</title>
<style>body{{margin:0;background:#f5f7fa;font:14px system-ui;color:#273447}}main{{max-width:1900px;margin:24px auto;background:white;padding:16px;border-radius:12px}}img{{width:100%;height:auto}}details{{margin:16px}}table{{border-collapse:collapse;width:100%}}td,th{{border:1px solid #e3e8ef;padding:8px;text-align:left}}.scroll{{overflow:auto}}</style>
</head><body><main><img alt="{heading}" src="data:image/png;base64,{png64}">
<details><summary>Exact visible node values and metadata</summary><div class="scroll"><table>
<tr><th>Layer</th><th>Node</th><th>Color score</th><th>Signed Cox contribution</th><th>Value</th><th>Imputed</th><th>Details</th></tr>{''.join(rows)}</table></div></details>
</main></body></html>'''

    # ------------------------------------------------------------ 独立导出接口
    def save_html(self, path: str | Path) -> Path:
        """**可独立调用**的 HTML 导出（无需先保存 PNG）。

        内部按同一渲染管线先取 PNG 字节再内嵌，因此 HTML 与 PNG 永远一致。
        """
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        import tempfile
        with tempfile.TemporaryDirectory() as td:
            png = self._save_figure(Path(td) / "figure.png", "png")
            rows = self._node_rows(Path(td) / "nodes.tsv")
            path.write_text(self._html_text(png.read_bytes(), rows), encoding="utf-8")
        return path

    def save_static(self, path: str | Path) -> Path:
        """按后缀导出静态图：``.png`` / ``.pdf`` / ``.svg``（同一 Figure）。"""
        path = Path(path)
        fmt = path.suffix.lower().lstrip(".")
        if fmt not in {"png", "pdf", "svg"}:
            raise ValueError(f"静态图后缀必须是 png/pdf/svg，收到 {path.suffix!r}")
        return self._save_figure(path, fmt)

    def save(self, prefix: str | Path) -> dict[str, Path]:
        """PNG/PDF/SVG from exactly the same figure; HTML embeds the SAME PNG bytes."""
        prefix = Path(prefix)
        prefix.parent.mkdir(parents=True, exist_ok=True)
        files = {ext: Path(str(prefix) + "." + ext) for ext in ("png", "pdf", "svg", "html")}
        for ext in ("png", "pdf", "svg"):
            self._save_figure(files[ext], ext)
        rows = self._node_rows(Path(str(prefix) + ".nodes.tsv"))
        files["nodes.tsv"] = Path(str(prefix) + ".nodes.tsv")
        files["html"].write_text(self._html_text(files["png"].read_bytes(), rows),
                                 encoding="utf-8")
        audit = {"sample": self.data.sample_id, "model": self.data.model, "risk": self.data.risk,
                 "rank": self.data.rank, "percentile": self.data.percentile, "n": len(self.data.cohort_risks),
                 "cutoff": self.cutoff, "cutoff_percentile": self.cutoff_percentile,
                 "cutoff_kind": "cohort_median_visualization_only" if self.data.cutoff is None else "provided_cutoff",
                 "counts": self.selection["counts"], "gene_threshold": self.style.gene_threshold,
                 "color_units": [self.data.expression_units, self.data.gene_score_units,
                                 self.data.signature_units, self.data.concept_units, "signed Cox contribution"],
                 "color_limits": [self.style.expression_limits, self.style.gene_score_limits,
                                  self.style.signature_limits, self.style.concept_limits, self.style.contribution_limits],
                 "hidden_cancer_type": self.selection["hidden_cancer_type"],
                 "displayed_edge_count": len(self.selection["gene_edges"])+len(self.selection["signature_edges"]),
                 "rank_definition": self.data.rank_definition, "provenance": self.data.provenance,
                 "decomposition_checked": bool(self.data.decomposition), "metadata": self.data.metadata,
                 "html_pipeline": "embedded identical PNG bytes; no alternate renderer"}
        files["selection.json"] = Path(str(prefix) + ".selection.json")
        files["selection.json"].write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8")
        if self.data.decomposition:
            files["decomposition.tsv"] = Path(str(prefix) + ".decomposition.tsv")
            fields = list(dict.fromkeys(k for r in self.data.decomposition for k in r))
            with files["decomposition.tsv"].open("w", encoding="utf-8", newline="") as f:
                w = csv.DictWriter(f, fieldnames=fields, delimiter="\t"); w.writeheader(); w.writerows(self.data.decomposition)
        return files

    def close(self) -> None:
        import matplotlib.pyplot as plt
        plt.close(self.figure)


def plot_sample_path(data: SamplePathData | Mapping[str, Any], style: Style = Style()) -> RenderedPath:
    """Draw four representation/predictor columns and a percentile risk bar."""
    if not isinstance(data, SamplePathData):
        data = SamplePathData.from_dict(data)
    selected = select_nodes(data, style)
    import numpy as np
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:  # pragma: no cover —— 仅在缺少绘图依赖时触发
        raise ImportError(
            "sample_path 绘图需要 matplotlib。请安装：pip install 'compass-os[plot]' "
            "或 pip install matplotlib。核心预测与 compass_os.sample_path() 的**计算**"
            "部分不依赖绘图库，可照常使用。") from exc
    from matplotlib.colors import LinearSegmentedColormap, Normalize
    from matplotlib.path import Path as MplPath
    from matplotlib.patches import PathPatch, Polygon, Rectangle
    from matplotlib import rc_context

    cutoff = float(np.median(data.cohort_risks)) if data.cutoff is None else float(data.cutoff)
    risks = np.asarray(data.cohort_risks, dtype=float)
    # Empirical mid-distribution position of the ACTUAL cutoff, not always 50%.
    cp = 100.0 * (float(np.sum(risks < cutoff)) + .5 * float(np.sum(risks == cutoff))) / len(risks)
    expression_cmap = LinearSegmentedColormap.from_list("expression", ["#426aba", "#fffdf6", "#d43b39"])
    score_cmap = LinearSegmentedColormap.from_list("scores", ["#fff8bd", "#ffd55b", "#ff862d", "#d72f32"])
    contribution_cmap = LinearSegmentedColormap.from_list("contribution", ["#477cba", "#fffdf6", "#d43b39"])
    risk_cmap = LinearSegmentedColormap.from_list("risk", ["#ffe66c", "#ffad42", "#e13b36"])
    with rc_context({"font.family": "DejaVu Sans", "font.size": 10}):
        fig, ax = plt.subplots(figsize=(style.width, style.height))
    fig.subplots_adjust(left=.015, right=.99, top=.99, bottom=.01)
    ax.set(xlim=(0, 175), ylim=(0, 100)); ax.set_axis_off()
    # Data coordinates give exact circles: keep equal geometry with fixed figure ratio.
    ax.set_aspect("equal", adjustable="box")
    xs = [13, 42, 80, 119]
    bottom, top = 17.0, 80.0
    layers = ["gene_tpm", "gene_score", "signatures", "predictors"]
    positions = {}
    # Arrange signatures by strongest concept target, then genes by strongest signature.
    hi_order = {n.name: i for i, n in enumerate(selected["predictors"])}
    sig_anchor = {}
    for n in selected["signatures"]:
        es = [e for e in selected["signature_edges"] if e.source == n.name]
        strongest = min(es, key=lambda e: (-abs(e.weight), hi_order[e.target])) if es else None
        sig_anchor[n.name] = hi_order[strongest.target] if strongest else 999
    selected["signatures"].sort(key=lambda n: (sig_anchor[n.name], -abs(n.score), n.name))
    si = {n.name: i for i, n in enumerate(selected["signatures"])}
    def gene_anchor(n):
        es = [e for e in selected["gene_edges"] if e.source == n.name]
        e = min(es, key=lambda e: (-abs(e.weight), si[e.target]))
        return si[e.target], -abs(n.score), n.name
    selected["gene_score"].sort(key=gene_anchor)
    exp = {n.name: n for n in selected["gene_tpm"]}
    selected["gene_tpm"] = [exp[n.name] for n in selected["gene_score"]]
    for layer, x in zip(layers, xs):
        ns = selected[layer]
        ys = np.linspace(top, bottom, len(ns)) if len(ns) > 1 else [(top+bottom)/2] * len(ns)
        positions[layer] = {n.name: (x, float(y)) for n, y in zip(ns, ys)}

    def link(a, b):
        dx = (b[0]-a[0])*.40
        path = MplPath([a, (a[0]+dx, a[1]), (b[0]-dx, b[1]), b],
                       [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4])
        ax.add_patch(PathPatch(path, fill=False, lw=style.line_width, edgecolor=style.line_color,
                               alpha=style.line_alpha, zorder=1))
    for n in selected["gene_score"]:
        link(positions["gene_tpm"][n.name], positions["gene_score"][n.name])
    for e in selected["gene_edges"]:
        link(positions["gene_score"][e.source], positions["signatures"][e.target])
    for e in selected["signature_edges"]:
        link(positions["signatures"][e.source], positions["predictors"][e.target])
    marker_y = bottom+(top-bottom)*data.percentile/100
    cutoff_y = bottom+(top-bottom)*cp/100
    # Lines converge on triangle base; its tip points toward the risk bar.
    for n in selected["predictors"]:
        link(positions["predictors"][n.name], (148, marker_y))

    cmaps = [expression_cmap, score_cmap, score_cmap, score_cmap]
    limits = [style.expression_limits, style.gene_score_limits, style.signature_limits, style.concept_limits]
    radii = [style.gene_radius, style.gene_radius, style.signature_radius, style.predictor_radius]
    units = [data.expression_units, data.gene_score_units, data.signature_units, data.concept_units]
    for i, layer in enumerate(layers):
        if not selected[layer]:
            ax.text(xs[i], (top+bottom)/2, "No nodes pass\ndisplay threshold", ha="center", va="center", color="#8493a3", fontsize=9)
        for n in selected[layer]:
            x, y = positions[layer][n.name]
            cmap, lim = (contribution_cmap, style.contribution_limits) if n.kind in {"clinical", "pc", "cancer"} else (cmaps[i], limits[i])
            # scatter marker sizes are in points, so circles cannot stretch.
            ax.scatter([x], [y], s=(2*radii[i])**2, c=[cmap(Normalize(*lim, clip=True)(n.score))],
                       edgecolors="#35465a", linewidths=.65, zorder=4)
            label = n.name.replace("_", " ") + (" (imputed)" if n.imputed else "")
            if len(label) > 30:
                label = label[:27] + "…"
            # All displayed nodes have labels; no hiding labels to claim truncation.
            ax.text(x-1.3 if i < 2 else x+1.4, y, label, ha="right" if i < 2 else "left",
                    va="center", fontsize=7.8 if i < 3 else 9,
                    color="#263447", zorder=5,
                    bbox={"facecolor": "white", "edgecolor": "none", "alpha": .88, "pad": .25})

    risk_x, risk_w = 152, 3.2
    gradient = np.linspace(0, 1, 512)[:, None]
    ax.imshow(gradient, origin="lower", aspect="auto", cmap=risk_cmap,
              extent=(risk_x, risk_x+risk_w, bottom, top), zorder=2)
    ax.add_patch(Rectangle((risk_x, bottom), risk_w, top-bottom, fill=False, lw=.8, edgecolor="#768394", zorder=3))
    ax.plot([risk_x-.6, risk_x+risk_w+.6], [cutoff_y]*2, linestyle=(0, (3, 2)), lw=1.2, color="#293a50", zorder=5)
    ax.add_patch(Polygon([(148, marker_y-1.0), (148, marker_y+1.0), (151.4, marker_y)],
                         facecolor="#263447", edgecolor="white", linewidth=.6, zorder=6))
    ax.text(risk_x+risk_w+1.4, top, "100%\nHigher risk", fontsize=8.5, va="top", color="#8e2e34")
    ax.text(risk_x+risk_w+1.4, bottom, "0%\nLower risk", fontsize=8.5, va="bottom", color="#9f7717")
    # Callouts live in a reserved region beyond the bar. Separate near the cutoff.
    risk_label_y = min(top-8, max(bottom+8, marker_y))
    cutoff_label_y = cutoff_y
    if abs(risk_label_y-cutoff_label_y) < 8:
        cutoff_label_y = min(top-13, risk_label_y+10) if risk_label_y < (top+bottom)/2 else max(bottom+13, risk_label_y-10)
    ax.annotate(f"Risk η = {data.risk:.4f}\nRank = {data.rank:g} / {len(risks)}\nPercentile = {data.percentile:.1f}%",
                xy=(risk_x+risk_w, marker_y), xytext=(158, risk_label_y), fontsize=9, va="center",
                color="#263447", arrowprops={"arrowstyle":"-", "color":"#8291a3", "lw":.6},
                bbox={"boxstyle":"round,pad=.5", "fc":"white", "ec":"#e4e9ee"}, zorder=7)
    cutoff_label = "Median cutoff" if data.cutoff is None else "Provided cutoff"
    ax.annotate(f"{cutoff_label}\n{cutoff:.4f}", xy=(risk_x+risk_w, cutoff_y), xytext=(158, cutoff_label_y),
                fontsize=8.3, va="center", color="#627185",
                arrowprops={"arrowstyle":"-", "color":"#a3afbd", "lw":.6},
                bbox={"fc":"white", "ec":"none", "pad":.2}, zorder=7)

    title = f"{data.model}  /  {data.sample_id}"
    ax.text(3, 97, title, fontsize=17, fontweight="bold", color="#263447", va="top")
    ax.text(3, 92, f"{data.cohort}  ·  N = {len(risks)}  ·  {data.provenance}", fontsize=10, color="#77869a")
    headings = [data.expression_label, "COMPASS\ngene score", "Granular signature\nscore", "High-level concepts\n+ clinical / PC", "Cohort-relative\nrisk"]
    for i, (x, heading) in enumerate(zip(xs+[153.6], headings)):
        ax.text(x, 86, heading, ha="center", va="center", fontsize=11.5, fontweight="bold", color="#35465a")
        if i < 4:
            ax.text(x, 82.5, f"{selected['counts'][layers[i]]} nodes", ha="center", fontsize=8, color="#8996a6")
    # Compact per-layer legends; scales are explicit and shared across M2/M3.
    def legend(x, y, cmap, lim, label):
        ax.imshow(np.linspace(0, 1, 100)[None, :], cmap=cmap, extent=(x, x+12, y, y+1), aspect="auto", zorder=3)
        ax.text(x, y-1.2, f"{lim[0]:g}", fontsize=7.5, color="#66768a", va="top")
        ax.text(x+12, y-1.2, f"{lim[1]:g}", fontsize=7.5, color="#66768a", ha="right", va="top")
        ax.text(x+6, y+2, label, fontsize=7.4, color="#66768a", ha="center")
    for x, cmap, lim, unit in zip(xs, cmaps, limits, units):
        legend(x-6, 10, cmap, lim, unit)
    legend(146, 10, contribution_cmap, style.contribution_limits, "Clinical / PC contribution")
    ax.text(3, 5.5, f"Selection: |gene score| ≥ {style.gene_threshold:g}; genes ≤ {style.max_genes}, signatures ≤ {style.max_signatures}, predictors ≤ {style.predictor_budget}.", fontsize=8.1, color="#77869a")
    ax.text(3, 3.4, "Node color = displayed score; all links have equal width and color and show selected model connections. Color values outside legend limits are clipped.", fontsize=8.1, color="#77869a")
    note = "Cancer type is constant, hidden here and retained in the frozen model. " if selected["hidden_cancer_type"] else ""
    ax.text(3, 1.3, note+"Risk is cohort-relative, not survival probability. Omitted nodes remain in prediction; this selected computation path is not a causal mechanism graph.", fontsize=7.5, color="#77869a")
    return RenderedPath(fig, data, style, selected, cutoff, cp)


def get_projector_weights(model: Any, concept2plot: Sequence[str] | None = None):
    """Portable replacement of the uploaded topology helper; not gene scores.

    Returns a pandas table with the original source/target/weights/group/concept
    columns. Torch and pandas are imported only if this helper is invoked.
    Supports the architecture used in the uploaded sankey.py.
    """
    import numpy as np
    import pandas as pd
    import torch
    gene_names = getattr(model, "feature_name", None)
    if gene_names is None:
        gene_names = model.pretrainer.feature_name
    projector = model.model.latentprojector
    gp, cp = projector.genesetprojector, projector.cellpathwayprojector
    gn, sn = np.asarray(gene_names), np.asarray(gp.genesets_names)
    rows = []
    for i, (concept, indices) in enumerate(zip(cp.cellpathway_names, cp.cellpathway_indices)):
        if concept2plot is not None and concept not in concept2plot:
            continue
        weights = torch.softmax(cp.cellpathway_aggregator.aggregator.attention_weights[f"cellpathway_{i}"].detach(), dim=0).cpu().numpy().reshape(-1)
        for j, weight in zip(indices, weights):
            j = int(j)
            rows.append((str(sn[j]), str(concept), float(weight), "geneset->celltype", str(concept)))
            gw = torch.softmax(gp.geneset_aggregator.aggregator.attention_weights[f"geneset_{j}"].detach(), dim=0).cpu().numpy().reshape(-1)
            for gene, w in zip(gn[gp.genesets_indices[j]], gw):
                rows.append((str(gene), str(sn[j]), float(w), "gene->geneset", str(concept)))
    return pd.DataFrame(rows, columns=["source", "target", "weights", "group", "concept"])


def plot_sankey_diagram(model=None, *, sample_data=None, style=Style(), **kwargs):
    """Migration entry point. Explicit sample scores are REQUIRED, never inferred from attention.

    Returns (selected_topology_table, matplotlib_figure). Prefer plot_sample_path
    for unified export. Old model-only calls fail with an actionable message.
    """
    if kwargs:
        raise TypeError("Legacy layout/importance arguments are unsupported; pass Style(...) and SamplePathData")
    if sample_data is None:
        raise ValueError("Provide sample_data=SamplePathData(...). Model attention weights are not sample gene scores or Cox risks.")
    r = plot_sample_path(sample_data, style)
    import pandas as pd
    rows = [(e.source, e.target, e.weight) for e in r.selection["gene_edges"]+r.selection["signature_edges"]]
    return pd.DataFrame(rows, columns=["source", "target", "weights"]), r.figure


def _demo(model: str) -> SamplePathData:
    """Deterministic synthetic LAYOUT example; deliberately labeled everywhere."""
    import numpy as np
    rng = np.random.default_rng(704)
    genes = ["CD3D", "CD3E", "CD8A", "GZMB", "PRF1", "NKG7", "IFNG", "CXCL9", "CXCL10", "STAT1", "FOXP3", "IL2RA", "CTLA4", "LAG3", "PDCD1", "MS4A1", "CD79A", "CD74", "LYZ", "CD68", "CSF1R", "C1QA", "C1QB", "COL1A1", "COL3A1", "DCN", "VWF", "PECAM1", "EPCAM", "MKI67", "TOP2A", "RAD51"]
    concepts = ["Cytotoxic T cell", "IFNγ pathway", "Macrophage", "Immune checkpoint", "B cell general", "Cell proliferation", "Genome integrity", "Endothelial", "Treg", "NK cell", "Collagen remodeling", "Epithelial", "Monocyte", "Plasma cell", "Mast", "T cell general"]
    scores = rng.uniform(.3, 1.0, len(genes))
    gs = [Node(g, float(s), "gene") for g, s in zip(genes, scores)]
    tp = [Node(g, float(t), "gene") for g, t in zip(genes, rng.uniform(0, 20, len(genes)))]
    sig = [Node(f"{c} / signature {j+1}", float(rng.uniform(-.2, 1)), "signature") for c in concepts for j in range(2)]
    co = [Node(c, float(rng.uniform(-.2, 1)), "concept", float(rng.uniform(-.2, .3))) for c in concepts]
    extra = [Node("Age", .08, "clinical", .08, 62), Node("Sex", -.035, "clinical", -.035, "Female"),
             Node("Stage", .11, "clinical", .11, "Frozen fill", True)]
    if model == "M3":
        extra.append(Node("PC1–PC10", .18, "pc", .18, details={"aggregation": "signed sum of 10 standardized Cox terms"}))
    predictors = co + extra
    decomp = [{"feature": n.name, "contribution": n.contribution} for n in predictors]
    eta = math.fsum(n.contribution for n in predictors)
    risks = list(map(float, np.linspace(eta-.4, eta+1.2, 81)))
    se = [Edge(s.name, concepts[i//2], .5) for i, s in enumerate(sig)]
    ge = [Edge(genes[(i*2+j)%len(genes)], s.name, float(.2+j*.07)) for i, s in enumerate(sig) for j in range(3)]
    return SamplePathData("LAYOUT_DEMO_NOT_PATIENT", model, "Synthetic layout example", eta, 21, 25,
                          risks, tp, gs, sig, predictors, ge, se, cancer_types=["DEMO"],
                          rank_definition="ascending risk rank; percentile=(rank-1)/(N-1)*100",
                          provenance="SYNTHETIC DATA · LAYOUT ONLY", decomposition=decomp)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    g = parser.add_mutually_exclusive_group(required=True)
    g.add_argument("--input", type=Path, help="JSON with REAL sample scores and main-API risk metadata")
    g.add_argument("--demo", action="store_true", help="synthetic layout-only M2/M3; never real inference")
    parser.add_argument("--out", type=Path, default=Path("output/sample_path"))
    parser.add_argument("--gene-threshold", type=float, default=.25)
    parser.add_argument("--max-genes", type=int, default=32)
    parser.add_argument("--max-signatures", type=int, default=26)
    parser.add_argument("--style-json", type=Path, help="shared Style values, including actual score color limits")
    args = parser.parse_args()
    opts = json.loads(args.style_json.read_text()) if args.style_json else {}
    opts.setdefault("gene_threshold", args.gene_threshold)
    opts.setdefault("max_genes", args.max_genes)
    opts.setdefault("max_signatures", args.max_signatures)
    style = Style(**opts)
    if args.demo:
        from dataclasses import asdict
        args.out.mkdir(parents=True, exist_ok=True)
        for model in ("M2", "M3"):
            d = _demo(model)
            (args.out/f"{model}_layout_input.json").write_text(json.dumps(asdict(d), ensure_ascii=False, indent=2), encoding="utf-8")
            r = plot_sample_path(d, style); r.save(args.out/f"{model}_layout_preview"); r.close()
        print(f"Synthetic layout examples saved to {args.out.resolve()}")
    else:
        data = SamplePathData.from_dict(json.loads(args.input.read_text(encoding="utf-8")))
        r = plot_sample_path(data, style)
        for kind, path in r.save(args.out).items():
            print(kind, path.resolve())
        r.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
