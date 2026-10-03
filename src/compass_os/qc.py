# -*- coding: utf-8 -*-
"""QC：基因覆盖度、signature 覆盖度、concept 输入覆盖度。

**命名纪律**：43 个概念只输出 ``concept_input_coverage``（输入信息覆盖程度的加权），
**不得**称为 confidence / biological confidence / activation confidence——它不表示
概念的生物学有效性（见 docs/CONCEPT_SEMANTICS.md）。
"""
from __future__ import annotations

import functools
from dataclasses import dataclass, field
from typing import List, Mapping, Sequence

import numpy as np
import pandas as pd

from ._paths import asset


@functools.lru_cache(maxsize=1)
def gene_sets() -> pd.DataFrame:
    """132 基因集定义（``gene_set`` / ``broad_pathway`` / ``n_genes`` / ``genes``）。"""
    df = pd.read_csv(asset("src/compass_os/data/concept_gene_sets.tsv"), sep="\t")
    df["gene_list"] = df["genes"].str.split(":")
    return df


@functools.lru_cache(maxsize=1)
def cancer_codes() -> pd.DataFrame:
    """癌种码表（``cancer_type`` / ``compass_code`` / ``ct_column`` / ``is_reference_level``）。"""
    return pd.read_csv(asset("src/compass_os/data/cancer_codes.tsv"), sep="\t")


@dataclass
class CoverageQC:
    """基因层与表示层覆盖度汇总（``to_dict()`` 可直接进返回值 / JSON）。"""

    n_required_genes: int
    n_observed_genes: int
    n_missing_genes: int
    gene_coverage: float
    missing_gene_strategy: str
    missing_gene_names: List[str] = field(default_factory=list)
    n_samples: int = 0
    per_sample: pd.DataFrame | None = None
    signature_gene_coverage: pd.DataFrame | None = None
    concept_input_coverage: pd.DataFrame | None = None
    #: 916 个 signature 基因中被真实观测到的比例（**与 QC 阈值同一口径**）
    gsig_coverage: float | None = None
    warnings: List[str] = field(default_factory=list)
    #: 组合分级：{"global":…, "Gsig":…, "overall":…, "rule": "worse(global, Gsig)"}
    tiers: dict = field(default_factory=dict)

    def to_dict(self, max_names: int = 200) -> dict:
        return {
            "n_required_genes": self.n_required_genes,
            "n_observed_genes": self.n_observed_genes,
            "n_missing_genes": self.n_missing_genes,
            "gene_coverage": self.gene_coverage,
            "missing_gene_strategy": self.missing_gene_strategy,
            "missing_gene_names": self.missing_gene_names[:max_names],
            "n_missing_gene_names_reported": min(len(self.missing_gene_names), max_names),
            "n_samples": self.n_samples,
            "gsig_coverage": self.gsig_coverage,
            "qc_tier_global": (self.tiers or {}).get("global"),
            "qc_tier_Gsig": (self.tiers or {}).get("Gsig"),
            "qc_tier_overall": (self.tiers or {}).get("overall"),
            "qc_tier_rule": (self.tiers or {}).get("rule"),
            "warnings": self.warnings,
        }


def gene_coverage(expression: pd.DataFrame, required: Sequence[str],
                  strategy: str) -> CoverageQC:
    """基因层覆盖度。

    * ``n_observed_genes``／``gene_coverage`` 为**队列层**（按列是否存在计）；
    * 输入含 sample-specific NaN 时额外给 ``per_sample``（逐样本非 NaN 计数）。
    """
    required = [str(g) for g in required]
    present = [g for g in required if g in expression.columns]
    missing = [g for g in required if g not in expression.columns]
    n_req, n_obs = len(required), len(present)
    cov = n_obs / n_req if n_req else float("nan")
    qc = CoverageQC(n_required_genes=n_req, n_observed_genes=n_obs, n_missing_genes=len(missing),
                    gene_coverage=cov, missing_gene_strategy=strategy,
                    missing_gene_names=missing, n_samples=int(expression.shape[0]))

    sub = expression[present] if present else pd.DataFrame(index=expression.index)
    nan_per_sample = sub.isna().sum(axis=1) if present else pd.Series(0, index=expression.index)
    if present and int(nan_per_sample.sum()) > 0:
        qc.per_sample = pd.DataFrame({
            "n_observed_genes_sample": n_obs - nan_per_sample.astype(int),
            "n_missing_genes_sample": (len(missing) + nan_per_sample).astype(int),
        })
        qc.per_sample["gene_coverage_sample"] = (
            qc.per_sample["n_observed_genes_sample"] / n_req)
        qc.per_sample.index.name = "sample_id"
        qc.warnings.append(
            f"Input contains sample-specific NaN: {int((nan_per_sample > 0).sum())} "
            "sample(s) affected; per-sample statistics are reported in per_sample")
    return qc


def signature_coverage(observed_mask: pd.DataFrame, *,
                       detail: bool = False) -> tuple[pd.DataFrame, pd.DataFrame | None]:
    """每个 signature 的输入基因覆盖度 ``coverage_i = |observed ∩ set_i| / |set_i|``。

    参数
    ----
    observed_mask : DataFrame(samples × genes, bool)
        True = 该基因在该样本**真实观测到**（未被填补）。列集合应为输入中出现过的基因。

    返回
    ----
    (coverage, missing_detail)
        ``coverage``：samples × 132（若所有样本缺失模式一致则等价于 132 维向量）；
        ``missing_detail``：``detail=True`` 时给出每个 signature 缺失的基因名（join 后的字符串列）。
    """
    gs = gene_sets()
    cols = set(observed_mask.columns)
    mat = np.zeros((observed_mask.shape[0], len(gs)), dtype=float)
    detail_rows = []
    for j, gl in enumerate(gs["gene_list"]):
        in_cols = [g for g in gl if g in cols]
        if in_cols:
            mat[:, j] = observed_mask[in_cols].to_numpy(bool).mean(axis=1)
        else:  # 该 signature 的基因一个都不在输入里
            mat[:, j] = 0.0
        if detail:
            miss = [g for g in gl if (g not in cols) or (not observed_mask[g].all())]
            detail_rows.append({
                "gene_set": gs["gene_set"].iloc[j],
                "n_genes": int(gs["n_genes"].iloc[j]),
                "n_missing_genes": len(miss),
                "missing_genes": ":".join(miss),
            })
    cov = pd.DataFrame(mat, index=observed_mask.index, columns=gs["gene_set"].tolist())
    cov.index.name = "sample_id"
    return cov, (pd.DataFrame(detail_rows) if detail else None)


def concept_input_coverage(sig_cov: pd.DataFrame,
                           weights: Mapping[str, tuple[Sequence[str], np.ndarray]]
                           ) -> pd.DataFrame:
    """``concept_input_coverage_c = Σ_i w_ci · signature_gene_coverage_i``。

    ``weights`` 来自**冻结 checkpoint** 的 ``CellPathwayAttentionAggregator``
    （softmax 归一、非负、组内和为 1），由 :func:`compass_os.representation.concept_weights`
    提取；不是本项目自造的权重。
    """
    cols = list(sig_cov.columns)
    out = {}
    for concept, (members, w) in weights.items():
        idx = [cols.index(m) for m in members if m in cols]
        if not idx:
            out[concept] = np.zeros(sig_cov.shape[0])
            continue
        ww = np.asarray(w, dtype=float)[[members.index(cols[i]) for i in idx]]
        ww = ww / ww.sum() if ww.sum() > 0 else ww
        out[concept] = sig_cov.iloc[:, idx].to_numpy(float) @ ww
    df = pd.DataFrame(out, index=sig_cov.index)
    df.index.name = "sample_id"
    return df


@functools.lru_cache(maxsize=1)
def qc_config() -> dict | None:
    """冻结的 QC 阈值（``models/qc_config.json``）；不存在时返回 ``None``。

    阈值由 ``validation/run_missing_gene_validation.py`` 的实测曲线确定，
    **不是**先验拍定；未冻结前 package 不给覆盖度警告（只给覆盖度数值）。
    """
    from ._paths import repo_root
    p = repo_root() / "models" / "qc_config.json"
    if not p.is_file():
        return None
    import json
    return json.loads(p.read_text(encoding="utf-8"))


TIERS = ("recommended", "warning", "below_warning")


def gsig_coverage_of(qc: CoverageQC) -> float | None:
    """用于 QC 分级的 Gsig 覆盖度 = 916 个 signature 基因中被观测到的比例。

    ⚠ 与逐 signature 覆盖度矩阵（``signature_gene_coverage``，n×132）区分：
    后者用于逐 signature 的精细检查，前者与冻结阈值同一口径（阈值就是在
    "保留了多大比例的 signature 基因"这个轴上标定的）。
    """
    if qc.gsig_coverage is not None and np.isfinite(qc.gsig_coverage):
        return float(qc.gsig_coverage)
    if qc.signature_gene_coverage is not None and qc.signature_gene_coverage.size:
        return float(np.nanmedian(qc.signature_gene_coverage.to_numpy(float)))
    return None


def _tier_for(value: float, recommended, warning) -> str | None:
    if value is None or not np.isfinite(value):
        return None
    if recommended is not None and value >= float(recommended):
        return "recommended"
    if warning is not None and value >= float(warning):
        return "warning"
    if warning is None and recommended is None:
        return None
    return "below_warning"


def qc_tier(qc: CoverageQC, cfg: dict | None = None) -> dict:
    """组合 QC 分级。

    **组合规则：overall = worse(global tier, Gsig tier)**（取更差的一档）。
    理由：两个轴衡量的是同一件事的不同侧面（缺失基因里"有没有动到 signature 基因"），
    任一侧不可靠时整体读数都不可靠；任一侧阈值缺失时该侧不参与，并在返回中标注。
    """
    cfg = cfg if cfg is not None else qc_config()
    out = {"global": None, "Gsig": None, "overall": None, "rule": "worse(global, Gsig)"}
    if not cfg:
        return out
    out["global"] = _tier_for(qc.gene_coverage, cfg.get("recommended_global_coverage"),
                              cfg.get("warning_global_coverage"))
    gsig = gsig_coverage_of(qc)
    if gsig is not None:
        out["Gsig"] = _tier_for(gsig, cfg.get("recommended_Gsig_coverage"),
                                cfg.get("warning_Gsig_coverage"))
    graded = [t for t in (out["global"], out["Gsig"]) if t is not None]
    if graded:
        # TIERS 按"由好到差"排列 ⇒ 取 index 最大者 = 最差的一档
        out["overall"] = max(graded, key=lambda t: TIERS.index(t))
    return out


def coverage_warnings(qc: CoverageQC, cfg: dict | None = None) -> list:
    """按冻结阈值生成覆盖度警告（无配置时返回空列表）。"""
    cfg = cfg if cfg is not None else qc_config()
    out = []
    if not cfg:
        return out
    g, sg = qc.gene_coverage, gsig_coverage_of(qc)
    for label, val, rec_key, warn_key in (
            ("global gene coverage", g, "recommended_global_coverage", "warning_global_coverage"),
            ("signature (Gsig) coverage", sg, "recommended_Gsig_coverage",
             "warning_Gsig_coverage")):
        rec, warn = cfg.get(rec_key), cfg.get(warn_key)
        if val is None or not np.isfinite(val):
            continue
        if rec is not None and val < rec:
            lvl = "warning" if (warn is not None and val >= warn) else "below_warning"
            out.append(f"{label} = {val:.3f} is below the recommended {rec:.3f} "
                       f"(QC grade: {lvl}; thresholds from stress test "
                       f"v{cfg.get('version', '?')})")
    rule = cfg.get("signature_coverage_warning_threshold") or {}
    thr = rule.get("value")
    if thr is not None and qc.signature_gene_coverage is not None:
        bad = int((qc.signature_gene_coverage.to_numpy(float) < float(thr)).any(axis=0).sum())
        if bad:
            out.append(f"{bad}/132 signatures have input coverage below "
                       f"{float(thr):.3f}; interpret mechanism-oriented readouts with "
                       "caution")
    return out


def coverage_summary(qc: CoverageQC) -> str:
    """一行式人类可读摘要（供日志/示例打印）。"""
    s = (f"基因覆盖 {qc.n_observed_genes}/{qc.n_required_genes}"
         f"（{qc.gene_coverage:.2%}，缺失 {qc.n_missing_genes}，策略 {qc.missing_gene_strategy}）")
    if qc.signature_gene_coverage is not None:
        s += f"；signature 覆盖中位 {qc.signature_gene_coverage.to_numpy().mean():.2%}"
    if qc.concept_input_coverage is not None:
        s += f"；concept 输入覆盖中位 {qc.concept_input_coverage.to_numpy().mean():.2%}"
    return s
