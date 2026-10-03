# -*- coding: utf-8 -*-
"""Level 1 —— 一键式研究者工作流 :func:`analyze` 与结果对象 :class:`AnalysisResult`。

设计原则（见 README「Public API」）
----------------------------------
**Level 1（本模块）**：普通生物医学/生信研究者只需提供表达矩阵、癌种，可选临床与随访，
即可拿到「风险预测 + 43 concept / 132 signature 表征 + 输入质控 + 人类可读报告」。

内部**固定**、不暴露为参数：

* 模型 = ``models/model_manifest.json`` 的 ``default_model``（当前 M2）
* 缺失基因策略 = ``reference``（冻结 TCGA 参考中位数）
* QC 阈值 = 冻结的 ``models/qc_config.json``

**Level 2**：:func:`compass_os.predict` / :func:`get_representation` /
:func:`compass_os.check_robustness`（见 docs/ADVANCED_USAGE.md）。

措辞纪律
--------
43 concept / 132 signature 是 **mechanism-oriented representation /
hypothesis-generating feature**，用于 **downstream prioritization**。
本模块**不输出** activation / suppression / causal mechanism，也**不**把
「concept 与风险相关」表述为独立的下游证据：default 模型的 Cox 设计矩阵本身
就包含 43 concepts，此类关联一律标注为
**model-linked representation association**。
"""
from __future__ import annotations

import json
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import numpy as np
import pandas as pd

from . import qc as _qc
from .exceptions import InputError
from .preprocessing import reference_fill_log2, reference_fill_tpm

__all__ = ["analyze", "AnalysisResult", "load_table"]

_TIME_KEYS = ("time", "os_time_days", "os_time", "time_days", "followup_days", "duration")
_EVENT_KEYS = ("event", "os_event", "status", "death", "event_observed")


# --------------------------------------------------------------------------- #
# 输入载入
# --------------------------------------------------------------------------- #
def load_table(obj, *, what: str = "table") -> pd.DataFrame:
    """接受「文件路径」或 ``DataFrame``；路径按分隔符自动识别（TSV/CSV，支持 .gz）。"""
    if isinstance(obj, pd.DataFrame):
        df = obj.copy()
    elif isinstance(obj, (str, Path)):
        p = Path(obj)
        if not p.is_file():
            raise InputError(f"{what} 文件不存在：{p}")
        sep = "," if p.suffix.lower() in (".csv",) or p.name.lower().endswith(".csv.gz") else "\t"
        df = pd.read_csv(p, sep=sep, index_col=0, low_memory=False)
    else:
        raise InputError(f"{what} 必须是 DataFrame 或文件路径，收到 {type(obj).__name__}")
    df.columns = [str(c).strip() for c in df.columns]
    df.index = [str(i).strip() for i in df.index]
    return df


def _norm_cancer_type(cancer_type, index: Sequence[str]) -> pd.Series:
    """癌种输入：单字符串（全队列同一癌种）或与样本等长的序列。"""
    idx = [str(i) for i in index]
    if isinstance(cancer_type, str):
        return pd.Series([cancer_type.strip()] * len(idx), index=idx, name="cancer_type")
    if isinstance(cancer_type, Mapping):
        s = pd.Series({str(k): str(v) for k, v in cancer_type.items()})
        miss = [i for i in idx if i not in s.index]
        if miss:
            raise InputError(f"cancer_type 映射缺少 {len(miss)} 个样本的癌种")
        return s.reindex(idx).rename("cancer_type")
    seq = list(cancer_type)
    if len(seq) != len(idx):
        raise InputError(f"cancer_type 长度 {len(seq)} != 样本数 {len(idx)}")
    return pd.Series([str(s).strip() for s in seq], index=idx, name="cancer_type")


def _norm_survival(survival, index: Sequence[str]) -> pd.DataFrame | None:
    """随访数据：DataFrame/路径，自动识别 time 与 event 列名。"""
    if survival is None:
        return None
    df = load_table(survival, what="survival") if not isinstance(survival, Mapping) else \
        pd.DataFrame(survival)
    idx = [str(i) for i in index]
    low = {str(c).strip().lower(): c for c in df.columns}
    tcol = next((low[k] for k in _TIME_KEYS if k in low), None)
    ecol = next((low[k] for k in _EVENT_KEYS if k in low), None)
    if tcol is None or ecol is None:
        raise InputError(
            "survival 需要时间列与事件列；可识别的时间列名 "
            f"{list(_TIME_KEYS)}，事件列名 {list(_EVENT_KEYS)}；"
            f"实际列为 {list(df.columns)}")
    df = df.rename(columns={tcol: "time", ecol: "event"})
    sub = df.reindex(idx)
    sub["time"] = pd.to_numeric(sub["time"], errors="coerce")
    sub["event"] = pd.to_numeric(sub["event"], errors="coerce")
    if sub["time"].notna().sum() == 0:
        raise InputError("survival 的 time 列没有可用数值")
    return sub[["time", "event"]]


# --------------------------------------------------------------------------- #
# 结果对象
# --------------------------------------------------------------------------- #
@dataclass
class AnalysisResult:
    """Level 1 结果：风险、表征、质控，以及 ``summary()`` / ``save_report()``。"""

    cancer_type: list
    sample_ids: list
    model: str
    risk: pd.Series
    risk_rank: pd.Series | None
    risk_group: pd.Series | None            # 冻结切点（pan-cancer 训练中位）分层
    concept_scores: pd.DataFrame
    signature_scores: pd.DataFrame
    qc: _qc.CoverageQC
    #: 队列相对分层（**仅用于展示**；冻结切点未能区分该队列时才给出，绝不用于跨队列比较）
    risk_group_relative: pd.Series | None = None
    clinical_imputed: dict = field(default_factory=dict)
    survival: pd.DataFrame | None = None
    survival_stats: dict | None = None
    robustness: dict | None = None
    warnings: list = field(default_factory=list)
    notes: list = field(default_factory=list)
    _prediction: object = None      # 内部：完整 PredictionResult（不鼓励直接使用）

    # ---------------- 人类可读摘要 ----------------
    def summary(self, width: int = 64) -> str:
        """打印并返回人类可读摘要。"""
        txt = self._summary_text(width)
        print(txt)
        return txt

    def _summary_text(self, width: int = 64) -> str:
        q = self.qc
        L = ["COMPASS-OS analysis", "─" * width]
        cts = sorted(set(self.cancer_type))
        ct_txt = cts[0] if len(cts) == 1 else f"{len(cts)} cancer types ({', '.join(cts[:4])}" + \
            (", …)" if len(cts) > 4 else ")")
        L += [f"{'Samples':<32}{len(self.sample_ids)}",
              f"{'Cancer type':<32}{ct_txt}",
              f"{'Model':<32}Default prognostic model ({self.model})", ""]

        L.append("INPUT QUALITY")
        L.append(f"{'Gene coverage':<32}{q.gene_coverage:.1%} "
                 f"({q.n_observed_genes}/{q.n_required_genes} genes)")
        if q.gsig_coverage is not None:
            L.append(f"{'Signature-related coverage':<32}{q.gsig_coverage:.1%} "
                     "(of the 916 signature genes)")
        if q.concept_input_coverage is not None:
            con_cov = float(np.nanmedian(q.concept_input_coverage.to_numpy(float)))
            L.append(f"{'Concept input coverage':<32}{con_cov:.1%} (median across 43 concepts)")
        L.append(f"{'QC':<32}{self._qc_verdict()}")
        L.append("")

        L.append("SURVIVAL")
        L.append(f"{'Risk estimates':<32}Available (relative risk, linear predictor)")
        if self.risk_group is not None:
            n_hi = int((self.risk_group == "high").sum())
            n_lo = int((self.risk_group == "low").sum())
            L += [f"{'Risk stratification':<32}Available (frozen cutoff)",
                  f"{'High risk':<32}{n_hi}", f"{'Low risk':<32}{n_lo}"]
        elif self.risk_group_relative is not None:
            n_hi = int((self.risk_group_relative == "high").sum())
            n_lo = int((self.risk_group_relative == "low").sum())
            L += [f"{'Risk stratification':<32}Cohort-relative median split "
                  f"(frozen cutoff did not separate this cohort)",
                  f"{'Higher-risk half':<32}{n_hi}", f"{'Lower-risk half':<32}{n_lo}"]
        elif len(self.sample_ids) == 1:
            L.append(f"{'Risk stratification':<32}Not applicable (single sample)")
        else:
            L.append(f"{'Risk stratification':<32}Not available "
                     f"(cohort of {len(self.sample_ids)} is too small)")
        if self.survival_stats:
            s = self.survival_stats
            if s.get("logrank_p") is not None:
                L.append(f"{'Log-rank p (risk groups)':<32}{s['logrank_p']:.3g}")
            if s.get("hr_per_sd") is not None:
                lo, hi = s.get("hr_ci", (np.nan, np.nan))
                ci = f" [{lo:.2f}, {hi:.2f}]" if np.isfinite(lo) and np.isfinite(hi) else ""
                L.append(f"{'Hazard ratio per SD of risk':<32}{s['hr_per_sd']:.2f}{ci}")
        L.append("")

        L.append("BIOLOGICAL REPRESENTATIONS")
        L.append(f"{'43 concept profiles':<32}Generated")
        L.append(f"{'132 gene-signature profiles':<32}Generated")
        L.append(f"{'Interpretation':<32}hypothesis-generating features")
        L.append("")

        L.append("ROBUSTNESS")
        L.append(f"{'Missing-gene strategy':<32}Reference imputation (frozen median)")
        if self.robustness:
            for w in self.robustness.get("messages", []):
                L.append(f"{'Robustness check':<32}{w}")
        else:
            L.append(f"{'Robustness check':<32}Not required (coverage in recommended range)"
                     if self._coverage_ok() else f"{'Robustness check':<32}Not run")

        if self.clinical_imputed:
            L += ["", "CLINICAL VARIABLES IMPUTED (frozen training reference values)"]
            for k, v in self.clinical_imputed.items():
                L.append(f"  {k}: {v['n_missing']}/{len(self.sample_ids)} samples "
                         f"→ {v['fill']}")
        other = [w for w in self.warnings
                 if "Clinical variables filled" not in w
                 and "missing genes were handled" not in w]
        if other:
            L += ["", "NOTES"]
            L += ["  - " + w for w in other]
        L += ["", "─" * width,
              "43 concepts / 132 signatures are mechanism-oriented biological representations",
              "(hypothesis-generating features) for downstream prioritisation. They are not",
              "evidence of pathway activation, suppression or causal mechanism. Any association",
              "between concept scores and the risk score is a model-linked representation",
              "association: the default model uses these concepts as Cox predictors."]
        return "\n".join(L)

    def _qc_verdict(self) -> str:
        """把冻结阈值转成研究者能读懂的一句话（overall = worse(global, Gsig)）。"""
        cfg = _qc.qc_config()
        if not cfg:
            return (f"{self.qc.gene_coverage:.1%} observed "
                    "(QC thresholds not frozen in this installation)")
        tier = (self.qc.tiers or {}).get("overall")
        rec_g = cfg.get("recommended_global_coverage")
        rec_s = cfg.get("recommended_Gsig_coverage")
        detail = f"global {self.qc.gene_coverage:.1%}"
        if self.qc.gsig_coverage is not None:
            detail += f", signature genes {self.qc.gsig_coverage:.1%}"
        if tier == "recommended":
            return (f"Recommended range ({detail}; thresholds global ≥ {rec_g:.0%}, "
                    f"signature ≥ {rec_s:.0%})")
        if tier == "warning":
            return (f"Usable with caution ({detail}): risk ranking is expected to remain "
                    "reliable, representation-level readouts less so")
        if tier == "below_warning":
            return (f"Below the calibrated warning level ({detail}): treat risk ranking as "
                    "approximate and avoid quantitative interpretation of representation "
                    "scores")
        return f"{detail} (no graded threshold applied)"

    def _coverage_ok(self) -> bool:
        cfg = _qc.qc_config()
        rec = (cfg or {}).get("recommended_global_coverage")
        return rec is not None and self.qc.gene_coverage >= rec

    # ---------------- 报告落盘 ----------------
    def save_report(self, output_dir, *, make_plots: bool = True) -> dict:
        """把结果写入目录（TSV + 纯文本 + 可选 PDF/PNG）。返回文件清单。"""
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        written = {}

        def _tsv(name: str, df: pd.DataFrame):
            p = out / name
            df.to_csv(p, sep="\t")
            written[name] = str(p)

        (out / "summary.txt").write_text(self._summary_text() + "\n", encoding="utf-8")
        written["summary.txt"] = str(out / "summary.txt")

        pred = pd.DataFrame({"cancer_type": self.cancer_type,
                             "risk": self.risk.to_numpy(float)})
        pred["risk_rank"] = (self.risk_rank.to_numpy(float) if self.risk_rank is not None
                             else np.arange(1, len(pred) + 1))
        if self.risk_group is not None:
            pred["risk_group"] = self.risk_group.to_numpy()
            pred["risk_group_definition"] = "frozen_median_cutoff"
        elif self.risk_group_relative is not None:
            pred["risk_group"] = self.risk_group_relative.to_numpy()
            pred["risk_group_definition"] = "cohort_relative_median_display_only"
        else:
            pred["risk_group"] = "not_assessed"
            pred["risk_group_definition"] = "not_assessed"
        pred.index.name = "sample_id"
        _tsv("prediction.tsv", pred)

        qrow = self.qc.to_dict()
        qd = pd.DataFrame({"value": [qrow["n_required_genes"], qrow["n_observed_genes"],
                                     qrow["n_missing_genes"], qrow["gene_coverage"],
                                     qrow["missing_gene_strategy"], len(self.sample_ids)]},
                          index=["n_required_genes", "n_observed_genes", "n_missing_genes",
                                 "gene_coverage", "missing_gene_strategy", "n_samples"])
        qd.index.name = "metric"
        _tsv("qc.tsv", qd)
        if self.qc.signature_gene_coverage is not None:
            _tsv("signature_gene_coverage.tsv",
                 self.qc.signature_gene_coverage.T.rename_axis("gene_set"))
        if self.qc.concept_input_coverage is not None:
            _tsv("concept_input_coverage.tsv",
                 self.qc.concept_input_coverage.T.rename_axis("concept"))
        if self.qc.per_sample is not None:
            _tsv("qc_per_sample.tsv", self.qc.per_sample)

        _tsv("concept_scores.tsv", self.concept_scores)
        _tsv("signature_scores.tsv", self.signature_scores)

        grp_any = self.risk_group if self.risk_group is not None else self.risk_group_relative
        if grp_any is not None and len(self.sample_ids) >= 10:
            comp = (self.concept_scores.assign(_g=grp_any.to_numpy())
                    .groupby("_g").median().T)
            comp.index.name = "concept"
            comp["note"] = "model-linked representation association (concepts are Cox predictors)"
            _tsv("representation_by_risk_group.tsv", comp)

        if self.survival is not None and self.survival_stats:
            _tsv("survival_statistics.tsv",
                 pd.DataFrame({"value": list(self.survival_stats.values())},
                              index=list(self.survival_stats.keys())))
            surv = self.survival.copy()
            surv["risk"] = self.risk.to_numpy(float)
            if self.risk_group is not None:
                surv["risk_group"] = self.risk_group.to_numpy()
            elif self.risk_group_relative is not None:
                surv["risk_group"] = self.risk_group_relative.to_numpy()
                surv["risk_group_definition"] = "cohort_relative_median_display_only"
            _tsv("survival_inputs.tsv", surv)

        if self.robustness:
            rs = pd.DataFrame(self.robustness.get("table", {}))
            if len(rs):
                _tsv("robustness_summary.tsv", rs)

        if make_plots:
            written.update(self._plots(out))

        (out / "report_manifest.json").write_text(
            json.dumps({"files": written, "model": self.model,
                        "n_samples": len(self.sample_ids)},
                       ensure_ascii=False, indent=1), encoding="utf-8")
        written["report_manifest.json"] = str(out / "report_manifest.json")
        return written

    def _plots(self, out: Path) -> dict:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from . import plotting

        made = {}

        def _save(fig, stem):
            for ext in ("pdf", "png"):
                fig.savefig(out / f"{stem}.{ext}", bbox_inches="tight")
                made[f"{stem}.{ext}"] = str(out / f"{stem}.{ext}")
            plt.close(fig)

        fig, ax = plt.subplots(figsize=(4.6, 3.4), dpi=160)
        plotting.risk_distribution(self.risk.to_numpy(float), ax=ax,
                                   title=f"Risk score distribution ({self.model})")
        _save(fig, "risk_distribution")

        grp_plot = self.risk_group if self.risk_group is not None else self.risk_group_relative
        if self.survival is not None and grp_plot is not None and \
                int(self.survival["event"].fillna(0).sum()) >= 2:
            fig, ax = plt.subplots(figsize=(4.6, 4.0), dpi=160)
            try:
                title = ("Kaplan-Meier by risk group (frozen cutoff)"
                         if self.risk_group is not None else
                         "Kaplan-Meier by cohort-relative median split "
                         "(display only, not the frozen cutoff)")
                plotting.km_plot(self.risk.to_numpy(float),
                                 self.survival["time"].to_numpy(float),
                                 self.survival["event"].fillna(0).to_numpy(bool),
                                 group=grp_plot.to_numpy(), ax=ax, title=title)
                _save(fig, "km_curve")
            except Exception:  # 事件过少等
                plt.close(fig)

        if len(self.sample_ids) >= 3:
            fig, ax = plt.subplots(figsize=(6.0, 2.6), dpi=160)
            plotting.representation_profile(self.concept_scores, top=43, ax=ax,
                                            title="43 concept profiles (latent coordinates)")
            _save(fig, "concept_profiles")

        if self.robustness and self.robustness.get("plot_data"):
            fig, ax = plt.subplots(figsize=(4.8, 3.4), dpi=160)
            d = self.robustness["plot_data"]
            ax.bar(list(d.keys()), list(d.values()), color="#1f77b4")
            ax.set_ylabel("stability (vs full input)", fontsize=8)
            ax.set_ylim(0, 1.02)
            ax.tick_params(labelsize=8)
            ax.set_title("Robustness check: reference vs zero handling", fontsize=8)
            ax.figure.autofmt_xdate(rotation=20)
            _save(fig, "robustness_plot")
        return made


# --------------------------------------------------------------------------- #
# Level 1 入口
# --------------------------------------------------------------------------- #
def analyze(expression, cancer_type, clinical=None, survival=None, *,
            input_scale: str = "tpm", output_dir=None,
            robustness: str | bool = "auto") -> AnalysisResult:
    """一键式分析：表达（+ 癌种，可选临床与随访）→ 风险 + 表征 + 质控 + 报告。

    参数
    ----
    expression : DataFrame 或路径
        ``samples × genes``，列名 = gene symbol。**TPM**（默认）或 ``log2(TPM+1)``。
    cancer_type : str 或序列/映射
        TCGA 缩写（如 ``"LUAD"``）。单一字符串表示全队列同一癌种。
    clinical : DataFrame 或路径，可选
        列名 ``age`` / ``sex`` / ``stage``；缺失值用**训练期冻结参考值**填充并在报告中列出。
    survival : DataFrame 或路径，可选
        含时间列与事件列（自动识别 ``time``/``os_time_days``… 与 ``event``/``os_status``…）。
        提供后自动做 KM / log-rank / 连续风险 Cox。
    input_scale : {"tpm", "log2_tpm1"}
        表达尺度，**必须显式声明**；本函数不会自动猜测。
    output_dir : 路径，可选
        提供则等价于额外调用一次 :meth:`AnalysisResult.save_report`。
    robustness : {"auto", True, False}
        ``"auto"``（默认）：仅当基因覆盖低于**冻结的** recommended 阈值时，
        才运行一次 reference vs zero 的稳健性比较；否则不额外做前向计算。

    返回
    ----
    :class:`AnalysisResult`（``summary()`` / ``save_report()``）
    """
    from . import api as _api

    expr = load_table(expression, what="expression")
    idx = list(expr.index)
    ct = _norm_cancer_type(cancer_type, idx)
    clin = load_table(clinical, what="clinical") if clinical is not None else None
    surv = _norm_survival(survival, idx)

    if input_scale == "tpm":
        mx = float(np.nanmax(expr.to_numpy(float))) if expr.size else 0.0
        if mx and mx < 30:
            # 不做自动转换，只提示（可能是 log2 或已标准化的矩阵）
            pass

    pred = _api.predict(expr, ct.tolist(), clinical=clin, model="default")
    q = pred.qc
    cfg = _qc.qc_config()
    rec = (cfg or {}).get("recommended_global_coverage")

    rob = None
    if robustness is True or (robustness == "auto" and rec is not None
                              and q.gene_coverage < rec):
        rob = _robustness_block(expr, ct, clin, input_scale, pred)
    notes = []
    if robustness == "auto" and rec is None:
        notes.append("QC thresholds are not frozen in this installation; "
                     "coverage is reported as continuous values only.")
    if q.signature_gene_coverage is not None and not (cfg or {}).get("recommended_Gsig_coverage"):
        notes.append("A calibrated signature-coverage threshold is not available yet "
                     "(signature-targeted calibration pending); signature coverage is "
                     "reported as a continuous value.")
    if input_scale == "tpm" and expr.size and float(np.nanmax(expr.to_numpy(float))) < 30:
        notes.append("The expression values look small for linear TPM; if they are "
                     "log2(TPM+1), re-run with input_scale='log2_tpm1'.")

    frozen = (pred.risk_group[pred.default_model]
              if pred.risk_group is not None and len(idx) >= 30 else None)
    relative = None
    if _degenerate(frozen, min_n=5) and len(idx) >= 20:
        relative = pd.Series(np.where(pred.default_risk.to_numpy(float) >=
                                      float(np.median(pred.default_risk.to_numpy(float))),
                                      "high", "low"), index=idx, name="risk_group_relative")
        notes.append(
            "The frozen pan-cancer median cutoff does not separate this cohort "
            f"(all/most samples fall in one group). A cohort-relative median split is "
            "provided for display only; it is NOT the frozen cutoff and must not be "
            "compared across cohorts.")
        frozen = None

    res = AnalysisResult(
        cancer_type=ct.tolist(), sample_ids=idx, model=pred.default_model,
        risk=pred.default_risk, risk_rank=(pred.risk_rank[pred.default_model]
                                            if pred.risk_rank is not None else None),
        risk_group=frozen, risk_group_relative=relative,
        concept_scores=pred.concept_scores, signature_scores=pred.signature_scores,
        qc=q, survival=surv, robustness=rob, warnings=list(pred.warnings), notes=notes,
        clinical_imputed=_clinical_imputation(clin, idx, pred.default_model),
        _prediction=pred)
    if surv is not None:
        res.survival_stats = _survival_stats(res)
    if output_dir is not None:
        res.save_report(output_dir)
    return res


def _robustness_block(expr, ct, clin, input_scale, pred) -> dict:
    """低覆盖时自动比较 reference vs zero，并转成人类可读信息。"""
    from .api import check_robustness
    r = check_robustness(expr, ct.tolist(), clinical=clin,
                         strategies=("reference", "zero"), input_scale=input_scale)
    cs = float(np.nanmedian(r.concept_correlation["reference_vs_zero"].to_numpy(float)))
    ss = float(np.nanmedian(r.signature_correlation["reference_vs_zero"].to_numpy(float)))
    base = pred.default_risk.to_numpy(float)
    pert = r.risk["zero"].to_numpy(float)
    from scipy.stats import spearmanr
    sp = float(spearmanr(base, pert)[0]) if len(base) > 2 else float("nan")
    agree = float(np.mean((base >= np.median(base)) == (pert >= np.median(pert))))
    msgs = [
        f"risk ranking Spearman = {sp:.3f} (reference vs zero)",
        f"risk-group agreement = {agree:.0%}",
        f"median correlation: 43 concepts = {cs:.3f}, 132 signatures = {ss:.3f}",
    ]
    table = {
        "comparison": ["reference vs zero"],
        "risk_spearman": [sp], "risk_group_agreement": [agree],
        "concept_median_r": [cs], "signature_median_r": [ss],
        "gene_coverage": [float(pred.qc.gene_coverage)],
    }
    return {"messages": msgs, "table": table,
            "plot_data": {"risk ranking": sp, "risk group agreement": agree,
                          "43 concepts": cs, "132 signatures": ss}}


def _clinical_imputation(clinical, idx, model: str) -> dict:
    """统计各临床字段被冻结参考值填充的样本数（供报告显示，不进 warning 文本）。"""
    from . import survival as _surv
    lock = _surv.load_lock(model)
    feats = lock["feature_names"]
    out = {}
    if not any(f in feats for f in ("Age", "Sex", "Stage")):
        return out
    clin = _surv._clinical_frame(clinical, idx)
    for raw, nm in (("age", "Age"), ("sex", "Sex"), ("stage", "Stage")):
        if nm not in feats:
            continue
        k = int(clin[raw].isna().sum())
        if k:
            out[nm] = {"n_missing": k, "fill": lock["clinical_fill"][nm]["fill"]}
    return out


def _degenerate(group: pd.Series | None, min_frac: float = 0.10, min_n: int = 5) -> bool:
    """冻结切点是否把整个队列分到同一组（外部队列常见：风险整体低于 pan-cancer 训练中位）。"""
    if group is None:
        return True
    n = len(group)
    k = int((group == "high").sum())
    return min(k, n - k) < max(min_n, int(np.ceil(min_frac * n)))


def _survival_stats(res: AnalysisResult) -> dict:
    """队列随访分析：log-rank（按冻结切点分层）+ 连续风险的 Cox HR。"""
    from scipy.stats import chi2
    s = res.survival
    t = s["time"].to_numpy(float)
    e = s["event"].fillna(0).to_numpy(bool)
    stats: dict = {"n_samples_with_survival": int(np.isfinite(t).sum()),
                   "n_events": int(e.sum())}
    grp = res.risk_group if res.risk_group is not None else res.risk_group_relative
    if grp is not None and e.sum() >= 2 and len(set(grp)) > 1:
        try:
            from .plotting import logrank
            chi2v, p = logrank(t, e, grp.to_numpy())
            stats["logrank_chi2"] = chi2v
            stats["logrank_p"] = p
            stats["logrank_stratification"] = (
                "frozen median cutoff" if res.risk_group is not None
                else "cohort-relative median (display only; frozen cutoff did not separate "
                     "this cohort)")
        except Exception as exc:  # noqa: BLE001
            stats["logrank_p"] = None
            stats["logrank_reason"] = f"{type(exc).__name__}: {exc}"
    if e.sum() >= 5:
        try:
            from sksurv.linear_model import CoxPHSurvivalAnalysis
            from sksurv.util import Surv
            r = res.risk.to_numpy(float)
            sd = float(np.nanstd(r, ddof=1)) or 1.0
            y = Surv.from_arrays(e, t)
            m = CoxPHSurvivalAnalysis(alpha=1e-9).fit((r / sd).reshape(-1, 1), y)
            hr = float(np.exp(m.coef_[0]))
            stats["hr_per_sd"] = hr
            se = float(getattr(m, "standard_error_", np.array([np.nan]))[0])
            if np.isfinite(se):
                stats["hr_ci"] = (float(np.exp(m.coef_[0] - 1.96 * se)),
                                  float(np.exp(m.coef_[0] + 1.96 * se)))
                stats["hr_p"] = float(2 * (1 - chi2.cdf((m.coef_[0] / se) ** 2, 1)))
        except Exception:
            pass
    return stats
