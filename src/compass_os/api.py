# -*- coding: utf-8 -*-
"""面向使用者的三个入口：:func:`get_representation` / :func:`predict` / :func:`check_robustness`。

设计要点
--------
* 用户**不需要**自己拼 COMPASS 输入：癌种单独用 ``cancer_type`` 提供，包内部完成
  ``user cancer_type → canonical cancer type → COMPASS cancer code + CT_* one-hot``，
  两处来自**同一个** canonical 值；
* ``predict`` 的默认模型由 ``models/model_manifest.json`` 决定（当前 ``M2``），
  换默认**只需改 manifest**，不需要改代码；
* 主输出是**相对风险**（``linear_predictor`` / ``relative_risk``）；绝对生存概率是
  派生输出，单独标注校准限制；
* 缺失基因策略见 :mod:`compass_os.preprocessing`；无论用哪种策略，QC 一定返回，
  **绝不静默填补**。
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable, Mapping, Sequence

import numpy as np
import pandas as pd

from . import qc as _qc
from . import survival as _surv
from .exceptions import InputError
from .preprocessing import (STRATEGIES, align_expression, cancer_codes_for, feature_names,
                            reference_fill_log2, reference_fill_tpm)

__all__ = ["RepresentationResult", "PredictionResult", "RobustnessResult",
           "get_representation", "predict", "check_robustness", "available_models",
           "default_model", "model_info"]


# --------------------------------------------------------------------------- #
# 结果容器
# --------------------------------------------------------------------------- #
@dataclass
class RepresentationResult:
    """132 gene-signature 分数 + 43 concept 分数 + QC。"""

    sample_ids: list
    cancer_type: list
    signature_scores: pd.DataFrame
    concept_scores: pd.DataFrame
    qc: _qc.CoverageQC
    warnings: list = field(default_factory=list)

    def to_dict(self, max_names: int = 200) -> dict:
        return {"sample_ids": self.sample_ids, "cancer_type": self.cancer_type,
                "signature_scores": self.signature_scores, "concept_scores": self.concept_scores,
                "qc": self.qc.to_dict(max_names), "warnings": self.warnings}


@dataclass
class PredictionResult:
    """风险预测：M0–M3 / default + 两层表示 + QC。"""

    sample_ids: list
    cancer_type: list
    default_model: str
    risk: pd.DataFrame                 # samples × (M0..M3, default)
    linear_predictor: pd.DataFrame     # 同 risk（风险 = 线性预测子）
    risk_rank: pd.DataFrame | None     # 队列内秩（n 足够时）
    risk_group: pd.DataFrame | None    # 按冻结 median_cutoff 的 high/low（n 足够时）
    signature_scores: pd.DataFrame
    concept_scores: pd.DataFrame
    qc: _qc.CoverageQC
    warnings: list = field(default_factory=list)

    @property
    def default_risk(self) -> pd.Series:
        return self.risk[self.default_model]

    def survival_probability(self, times: Sequence[float]) -> pd.DataFrame:
        """派生输出：``S_i(t)``（baseline-hazard-based；校准限制见 survival 模块）。"""
        return _surv.survival_probability(self.default_risk.to_numpy(), times,
                                          self.default_model, self.sample_ids)

    def to_dict(self, max_names: int = 200) -> dict:
        return {"sample_ids": self.sample_ids, "cancer_type": self.cancer_type,
                "default_model": self.default_model, "risk": self.risk,
                "linear_predictor": self.linear_predictor, "risk_rank": self.risk_rank,
                "risk_group": self.risk_group, "signature_scores": self.signature_scores,
                "concept_scores": self.concept_scores, "qc": self.qc.to_dict(max_names),
                "warnings": self.warnings}


@dataclass
class RobustnessResult:
    """``check_robustness`` 结果：连续指标 + 尚未校准的标志位。"""

    sample_ids: list
    gene_coverage: float
    signature_gene_coverage_reference: pd.DataFrame
    risk: dict
    risk_difference: pd.Series
    concept_correlation: dict
    signature_correlation: dict
    concept_input_coverage: dict
    robustness_flag: str | None = "continuous_only"
    note: str = (
        "不存在单一的二元 pass/fail 判定：覆盖度 QC 分级已按冻结阈值校准（recommended 0.90 / "
        "warning 0.70，global 与 signature 双轴，overall = worse(global, Gsig)），"
        "而 reference-versus-zero 稳健性比较在 v1.x 中**只提供连续指标**。"
        "请把连续的风险/表征稳定性指标与已校准的覆盖度分级一起解读。")

    def to_dict(self) -> dict:
        return {"sample_ids": self.sample_ids, "gene_coverage": self.gene_coverage,
                "risk": self.risk, "risk_difference": self.risk_difference,
                "concept_correlation": self.concept_correlation,
                "signature_correlation": self.signature_correlation,
                "concept_input_coverage": self.concept_input_coverage,
                "robustness_flag": self.robustness_flag, "note": self.note}


# --------------------------------------------------------------------------- #
# 内部：一次完整前向
# --------------------------------------------------------------------------- #
def _run(expression: pd.DataFrame, cancer_type, *, missing_gene_strategy: str,
         input_scale: str, batch_size: int, device: str) -> dict:
    genes = list(feature_names())
    aligned = align_expression(expression, missing_gene_strategy, input_scale=input_scale,
                               genes=genes)
    ct = pd.Series(list(cancer_type), index=aligned.matrix.index, name="cancer_type")
    if len(ct) != aligned.matrix.shape[0]:
        raise InputError(f"cancer_type 长度 {len(ct)} != 样本数 {aligned.matrix.shape[0]}")
    codes = cancer_codes_for(ct.tolist())

    from . import representation as _rep
    reps = _rep.compute_representation(
        aligned.matrix, codes.to_numpy(), missing_mask=aligned.missing_mask,
        strategy=missing_gene_strategy, batch_size=batch_size, device=device)

    sig, con = reps["signature_scores"], reps["concept_scores"]
    sig_cov, _ = _qc.signature_coverage(aligned.observed_mask)
    con_cov = _qc.concept_input_coverage(sig_cov, _rep.concept_weights())
    aligned.qc.signature_gene_coverage = sig_cov
    aligned.qc.concept_input_coverage = con_cov
    # 与冻结阈值同口径：916 个 signature 基因里被真实观测到的比例
    from compass_os.qc import gene_sets as _gs
    _sig_genes = {g for gl in _gs()["gene_list"] for g in gl}
    _present = [g for g in _sig_genes if g in aligned.observed_mask.columns]
    aligned.qc.gsig_coverage = (
        float(aligned.observed_mask[_present].to_numpy(bool).mean())
        if _present else 0.0)
    if aligned.qc.gene_coverage < 1.0:
        aligned.qc.warnings.append(
            f"{aligned.qc.n_missing_genes} missing genes were handled with the "
            f"'{missing_gene_strategy}' strategy; see also signature_gene_coverage and "
            "concept_input_coverage")
    cfg = _qc.qc_config()
    aligned.qc.tiers = _qc.qc_tier(aligned.qc, cfg)
    aligned.qc.warnings.extend(_qc.coverage_warnings(aligned.qc, cfg))
    if cfg is None:
        aligned.qc.warnings.append(
            "Coverage thresholds are not frozen in this installation "
            "(models/qc_config.json missing): coverage is reported as continuous values "
            "only, without grading")
    return {"aligned": aligned, "cancer_type": ct, "concept_scores": con,
            "signature_scores": sig, "qc": aligned.qc}


# --------------------------------------------------------------------------- #
# 公共入口
# --------------------------------------------------------------------------- #
def get_representation(expression: pd.DataFrame, cancer_type: Iterable,
                       missing_gene_strategy: str = "reference", *,
                       input_scale: str = "tpm", batch_size: int = 16, device: str = "cpu",
                       signature_detail: bool = False) -> RepresentationResult:
    """计算 132 gene-signature 分数与 43 concept 分数（不做生存预测）。

    参数
    ----
    expression : DataFrame
        ``samples × genes``，列名 = 基因符号，索引 = 样本 ID。
    cancer_type : 长度 = 样本数
        TCGA 缩写（如 ``"LUAD"``）；**必需**，缺失即报错。
    missing_gene_strategy : {"reference", "zero", "strict"}
        见 :mod:`compass_os.preprocessing`。
    input_scale : {"tpm", "log2_tpm1"}
        输入表达尺度；模型内部只接受 TPM，``log2_tpm1`` 会被转回。

    返回
    ----
    :class:`RepresentationResult`
        ``signature_scores``（n × 132）、``concept_scores``（n × 43）、``qc``。
    """
    out = _run(expression, cancer_type, missing_gene_strategy=missing_gene_strategy,
               input_scale=input_scale, batch_size=batch_size, device=device)
    res = RepresentationResult(sample_ids=[str(i) for i in out["aligned"].matrix.index],
                               cancer_type=out["cancer_type"].tolist(),
                               signature_scores=out["signature_scores"],
                               concept_scores=out["concept_scores"],
                               qc=out["qc"], warnings=list(out["qc"].warnings))
    if signature_detail:
        _, detail = _qc.signature_coverage(out["aligned"].observed_mask, detail=True)
        res.warnings.append("signature 缺失明细见 qc.signature_missing_detail")
        res.qc.__dict__["signature_missing_detail"] = detail
    return res


def predict(expression: pd.DataFrame, cancer_type: Iterable,
            clinical: pd.DataFrame | None = None, model: str = "default", *,
            missing_gene_strategy: str = "reference", input_scale: str = "tpm",
            batch_size: int = 16, device: str = "cpu",
            min_cohort_for_stratification: int = 30) -> PredictionResult:
    """预测 M0–M3（以及 ``default``）的**预后风险分数**（Cox 线性预测子）+ 132/43 表征 + QC。

    参数
    ----
    clinical : DataFrame, 可选
        ``index = sample_id``，列名 ``age`` / ``sex`` / ``stage``（stage ∈ 1–4）。
        缺失的临床字段用**训练期冻结的参考值**填充（Age 60 / Sex 0 / Stage 2），
        并在 ``qc`` 中记录被填充的字段与个数——不重新估计任何统计量。
    model : str
        ``"default"``（按 manifest，当前 ``M2``）或显式 ``M0``/``M1``/``M2``/``M3``；
        也可传逗号分隔的多个（如 ``"M1,M2"``）。

    返回
    ----
    :class:`PredictionResult`
        ``risk`` / ``linear_predictor``（samples × 模型；**Cox 线性预测子 η = Xβ**，
        数值越大表示模型估计的风险越高；``exp(η)`` 才是相对风险比 relative hazard）、
        ``risk_rank``、``risk_group``、
        ``signature_scores``、``concept_scores``、``qc``。

    备注
    ----
    单样本**不产生** high/low 分层（``risk_group=None``）：分层依赖队列内中位数。
    队列 ≥ ``min_cohort_for_stratification`` 例时才给 ``risk_group``，
    切点用冻结的 ``median_cutoff``（训练集中位风险），不是用户队列中位数。
    """
    out = _run(expression, cancer_type, missing_gene_strategy=missing_gene_strategy,
               input_scale=input_scale, batch_size=batch_size, device=device)
    aligned, ct, con = out["aligned"], out["cancer_type"], out["concept_scores"]
    sample_ids = [str(i) for i in aligned.matrix.index]

    models = _parse_models(model)
    want_default = _surv.default_model()
    if "default" in models:
        models = [m for m in models if m != "default"]
        if want_default not in models:
            models.append(want_default)

    warnings = list(out["qc"].warnings)
    clin_filled = _clinical_fill_report(clinical, sample_ids, models)
    if clin_filled:
        warnings.append("Clinical variables filled with frozen training reference values: "
                        + "; ".join(clin_filled))

    risk = {}
    for mk in models:
        X = _surv.build_design_matrix(con, clinical, aligned.log2, mk, ct.tolist())
        risk[mk] = _surv.risk_from_design(X, mk)
    risk_df = pd.DataFrame(risk, index=sample_ids)
    risk_df = risk_df[[m for m in ("M0", "M1", "M2", "M3") if m in risk_df.columns]
                      + [c for c in risk_df.columns if c not in ("M0", "M1", "M2", "M3")]]

    # 显式请求 "default" 时额外给 default 列；默认行为已包含 want_default
    n = len(sample_ids)
    # v1.0.1 fix：风险**排序**与 **high/low 分层**明确分开，且与文档一致：
    #   n == 1                    → risk_rank = None, risk_group = None
    #   2 <= n <  min_cohort      → risk_rank 可用,  risk_group = None
    #   n >= min_cohort           → risk_rank 可用,  risk_group = 冻结切点 high/low
    rank = None
    if n >= 2:
        rank = risk_df.rank(method="average").astype(float)
    group = None
    if n >= min_cohort_for_stratification:
        cut = _surv.median_cutoff(want_default)
        group = pd.DataFrame(
            {mk: np.where(risk_df[mk].to_numpy() >= cut, "high", "low") for mk in risk_df.columns},
            index=sample_ids)
    elif n >= 2:
        warnings.append(
            f"Cohort has {n} samples (< min_cohort_for_stratification="
            f"{min_cohort_for_stratification}): risk_rank is provided, but risk_group is not "
            "assigned (within-cohort high/low comparison would be unreliable at this size)")
    else:
        warnings.append("Single sample: neither risk_rank nor risk_group is defined")

    return PredictionResult(sample_ids=sample_ids, cancer_type=ct.tolist(),
                            default_model=want_default, risk=risk_df,
                            linear_predictor=risk_df.copy(), risk_rank=rank, risk_group=group,
                            signature_scores=out["signature_scores"],
                            concept_scores=con, qc=out["qc"], warnings=warnings)


def check_robustness(expression: pd.DataFrame, cancer_type: Iterable,
                     clinical: pd.DataFrame | None = None,
                     strategies: Sequence[str] = ("reference", "zero"), *,
                     model: str = "default", input_scale: str = "tpm",
                     batch_size: int = 16, device: str = "cpu") -> RobustnessResult:
    """对同一输入分别用多种缺失基因策略预测，报告差异（**连续指标，无阈值判定**）。

    返回 ``risk_reference`` / ``risk_zero`` 形式的逐策略风险、``risk_difference``、
    concept/signature 的逐样本相关，以及 ``robustness_flag="continuous_only"``。
    """
    strategies = tuple(strategies)
    bad = [s for s in strategies if s not in STRATEGIES]
    if bad:
        raise InputError(f"未知策略 {bad}；可用 {STRATEGIES}")
    if len(strategies) < 2:
        raise InputError("check_robustness 至少需要两种策略才能比较")

    runs = {s: _run(expression, cancer_type, missing_gene_strategy=s, input_scale=input_scale,
                    batch_size=batch_size, device=device)
            for s in strategies}
    mk = _surv.default_model() if model == "default" else str(model)
    risks, sigs, cons = {}, {}, {}
    for s, out in runs.items():
        X = _surv.build_design_matrix(out["concept_scores"], clinical, out["aligned"].log2,
                                      mk, out["cancer_type"].tolist())
        risks[s] = _surv.risk_from_design(X, mk)
        sigs[s], cons[s] = out["signature_scores"], out["concept_scores"]

    base = strategies[0]
    risk_diff = pd.Series(risks[base] - risks[strategies[1]],
                          index=runs[base]["aligned"].matrix.index, name=f"{base}_minus_"
                          f"{strategies[1]}")
    cc, sc = {}, {}
    for s in strategies[1:]:
        cc[f"{base}_vs_{s}"] = _row_corr(cons[base], cons[s])
        sc[f"{base}_vs_{s}"] = _row_corr(sigs[base], sigs[s])
    return RobustnessResult(
        sample_ids=[str(i) for i in runs[base]["aligned"].matrix.index],
        gene_coverage=float(runs[base]["qc"].gene_coverage),
        signature_gene_coverage_reference=runs[base]["qc"].signature_gene_coverage,
        risk={s: pd.Series(v, index=runs[base]["aligned"].matrix.index) for s, v in risks.items()},
        risk_difference=risk_diff, concept_correlation=cc, signature_correlation=sc,
        concept_input_coverage={s: runs[s]["qc"].concept_input_coverage for s in strategies},
    )


def _row_corr(a: pd.DataFrame, b: pd.DataFrame) -> pd.Series:
    """逐样本相关（先按共同行/列对齐）。"""
    idx = [i for i in a.index if i in set(b.index)]
    cols = [c for c in a.columns if c in set(b.columns)]
    A, B = a.loc[idx, cols].to_numpy(float), b.loc[idx, cols].to_numpy(float)
    out = []
    for i in range(A.shape[0]):
        x, y = A[i], B[i]
        ok = np.isfinite(x) & np.isfinite(y)
        out.append(float(np.corrcoef(x[ok], y[ok])[0, 1]) if ok.sum() > 2 and
                   np.std(x[ok]) > 0 and np.std(y[ok]) > 0 else np.nan)
    return pd.Series(out, index=idx, name="pearson_r")


def _parse_models(model) -> list:
    if isinstance(model, str):
        items = [m.strip() for m in model.split(",") if m.strip()]
    else:
        items = [str(m).strip() for m in model]
    if not items:
        raise InputError("model 不能为空")
    avail = set(_surv.available_models())
    bad = [m for m in items if m not in avail and m != "default"]
    if bad:
        raise InputError(f"未知模型 {bad}；可用：{sorted(avail)} 或 'default'")
    return items


def _clinical_fill_report(clinical, sample_ids, models) -> list:
    """报告被冻结参考值填充的临床字段与个数（避免静默填补）。"""
    if not any(any(f in _surv.load_lock(m)["feature_names"] for f in ("Age", "Sex", "Stage"))
               for m in models):
        return []
    clin = _surv._clinical_frame(clinical, sample_ids)
    out = []
    for f, nm in (("age", "Age"), ("sex", "Sex"), ("stage", "Stage")):
        k = int(clin[f].isna().sum())
        if k:
            out.append(f"{nm} {k}/{len(sample_ids)} samples (reference value "
                   f"{_surv.load_lock(models[0])['clinical_fill'][nm]['fill']})")
    return out


def available_models() -> tuple:
    """可用模型名（``M0``/``M1``/``M2``/``M3``），顺序为 manifest 顺序。"""
    return _surv.available_models()


def default_model() -> str:
    """当前默认模型（由 ``models/model_manifest.json`` 决定）。"""
    return _surv.default_model()


def model_info(model: str = "default") -> dict:
    """模型元信息（标签、特征组成、描述、锁定件路径）。"""
    mm = _surv.model_manifest()
    key = _surv.default_model() if model == "default" else str(model)
    if key not in mm["models"]:
        raise InputError(f"未知模型 {model!r}；可用：{sorted(mm['models'])}")
    info = dict(mm["models"][key])
    info["name"] = key
    info["is_default"] = (key == _surv.default_model())
    info["missing_gene"] = mm.get("missing_gene", {})
    info["reference_fill_log2"] = reference_fill_log2()
    info["reference_fill_tpm"] = reference_fill_tpm()
    return info
