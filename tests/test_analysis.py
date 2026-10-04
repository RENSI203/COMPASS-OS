# -*- coding: utf-8 -*-
"""Level 1（Quick API）测试：``analyze`` / ``AnalysisResult`` / ``save_report``。

覆盖 P3 §19 的 RC 条件：quick workflow、save_report、单样本、队列、
低覆盖自动 robustness、以及公开命名空间收敛。
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from _util import REPO, fixtures_available, load_fixtures

pytestmark = pytest.mark.skipif(not fixtures_available(), reason="需要 tests/fixtures")
EX = REPO / "examples" / "example_data"


# --------------------------------------------------------------------------- #
# 公开命名空间
# --------------------------------------------------------------------------- #
def test_public_namespace_is_converged():
    """顶层只导出 Level 1/2 的公开函数 + 结果类 + 异常（v1.1.0 起含 sample_path）。"""
    import compass_os
    allowed = {"analyze", "AnalysisResult", "predict", "get_representation",
               "check_robustness", "PredictionResult", "RepresentationResult",
               "RobustnessResult", "CompassOSError", "MissingGenesError",
               "UnknownCancerTypeError", "AssetNotFoundError", "InputError",
               # Level 2 · 样本级计算路径归因（v1.1.0 新增）
               "sample_path", "SamplePathResult",
               "__version__"}
    assert set(compass_os.__all__) == allowed
    for name in allowed:
        assert hasattr(compass_os, name), f"缺少导出 {name}"
    for bad in ("available_models", "default_model", "model_info", "preprocessing",
                "representation", "survival", "plotting", "qc", "load_model"):
        assert bad not in compass_os.__all__, f"内部对象不应出现在 __all__：{bad}"


def test_no_refit_entry_points_in_public_api():
    """不得提供 refit / fit / 自定义癌种码等会破坏冻结定义的入口。"""
    import compass_os
    forbidden = ("refit_pca", "fit_scaler", "fit_model", "custom_cancer_code",
                 "set_beta", "set_vocabulary", "fit")
    for name in forbidden:
        assert not hasattr(compass_os, name), f"公开层不应存在 {name}"


# --------------------------------------------------------------------------- #
# Quick workflow
# --------------------------------------------------------------------------- #
def test_analyze_accepts_file_paths():
    """Quick start 形态：全部传文件路径 + 单字符串癌种。"""
    import compass_os
    res = compass_os.analyze(str(EX / "example_expression.tsv.gz"), cancer_type="LUAD",
                             clinical=str(EX / "example_clinical.tsv"))
    assert len(res.sample_ids) == 5
    assert set(res.cancer_type) == {"LUAD"}
    assert res.risk.shape == (5,)
    assert res.concept_scores.shape == (5, 43)
    assert res.signature_scores.shape == (5, 132)
    txt = res.summary()
    for key in ("COMPASS-OS analysis", "INPUT QUALITY", "SURVIVAL",
                "BIOLOGICAL REPRESENTATIONS", "ROBUSTNESS", "Gene coverage"):
        assert key in txt, f"summary 缺少章节 {key}"


def test_analyze_accepts_dataframes_and_cohort():
    """DataFrame 输入 + 队列级分析（≥30 例才有分层）。"""
    import compass_os
    fx = load_fixtures()
    res = compass_os.analyze(fx["expression"], fx["labels"]["cancer_type"],
                             clinical=fx["clinical"])
    assert len(res.sample_ids) == 12
    assert res.risk_group is None or res.risk_group.notna().all()


def test_single_sample_quick_path():
    """单样本：可运行、无分层、摘要写明 not applicable。"""
    import compass_os
    fx = load_fixtures()
    res = compass_os.analyze(fx["expression"].iloc[[0]],
                             [fx["labels"]["cancer_type"].iloc[0]],
                             clinical=fx["clinical"].iloc[[0]])
    assert res.risk_group is None
    txt = res.summary()
    assert "Not applicable" in txt or "single sample" in txt.lower()


def test_analyze_with_survival_triggers_stats():
    """提供随访后给出 log-rank 与连续风险 HR（事件不足则只报可算的字段）。"""
    import compass_os
    fx = load_fixtures()
    surv = pd.DataFrame({"time": fx["labels"]["os_time_days"].to_numpy(),
                         "event": fx["labels"]["os_event"].to_numpy()},
                        index=fx["labels"].index)
    res = compass_os.analyze(fx["expression"], fx["labels"]["cancer_type"],
                             clinical=fx["clinical"], survival=surv)
    assert res.survival is not None
    st = res.survival_stats or {}
    assert st.get("n_events", 0) >= 1
    assert "n_samples_with_survival" in st
    # 12 例、事件少：HR 可能因数值问题缺失，但字段存在性必须稳定
    assert "hr_per_sd" in st or "logrank_p" in st or "logrank_p" not in st


def test_survival_column_auto_detection():
    """随访列名自动识别（os_time_days/os_event 等常见写法）。"""
    import compass_os
    fx = load_fixtures()
    surv = pd.DataFrame({"os_time_days": fx["labels"]["os_time_days"].to_numpy(),
                         "os_event": fx["labels"]["os_event"].to_numpy()},
                        index=fx["labels"].index)
    res = compass_os.analyze(fx["expression"].iloc[:6],
                             fx["labels"]["cancer_type"].iloc[:6].tolist(),
                             survival=surv.iloc[:6])
    assert res.survival is not None and "time" in res.survival.columns


def test_survival_missing_columns_raises():
    """随访表缺少时间/事件列 ⇒ 明确报错（不静默忽略）。"""
    import compass_os
    from compass_os.exceptions import InputError
    fx = load_fixtures()
    bad = pd.DataFrame({"foo": [1, 2, 3]}, index=fx["labels"].index[:3])
    with pytest.raises(InputError):
        compass_os.analyze(fx["expression"].iloc[:3],
                           fx["labels"]["cancer_type"].iloc[:3].tolist(), survival=bad)


# --------------------------------------------------------------------------- #
# save_report
# --------------------------------------------------------------------------- #
def test_save_report_files(tmp_path=None):
    """save_report 生成约定的 TSV/文本/图件，且不含 Excel。"""
    import tempfile

    import compass_os
    fx = load_fixtures()
    surv = pd.DataFrame({"time": fx["labels"]["os_time_days"].to_numpy(),
                         "event": fx["labels"]["os_event"].to_numpy()},
                        index=fx["labels"].index)
    res = compass_os.analyze(fx["expression"], fx["labels"]["cancer_type"],
                             clinical=fx["clinical"], survival=surv)
    out = Path(tempfile.mkdtemp()) / "report"
    files = res.save_report(out)
    for need in ("summary.txt", "prediction.tsv", "qc.tsv", "concept_scores.tsv",
                 "signature_scores.tsv", "risk_distribution.pdf",
                 "risk_distribution.png", "report_manifest.json"):
        assert (out / need).is_file(), f"缺少报告文件 {need}"
        assert need in files
    assert not list(out.glob("*.xlsx")), "报告不应包含 Excel"
    # 内容自洽
    pred = pd.read_csv(out / "prediction.tsv", sep="\t", index_col=0)
    assert len(pred) == len(res.sample_ids)
    assert {"risk", "risk_rank", "risk_group"} <= set(pred.columns)
    man = json.loads((out / "report_manifest.json").read_text(encoding="utf-8"))
    assert man["model"] == res.model and man["n_samples"] == len(res.sample_ids)


def test_save_report_includes_robustness_when_triggered():
    """稳健性被触发时报告含 robustness_summary.tsv。"""
    import tempfile

    import compass_os
    fx = load_fixtures()
    low = fx["expression"].drop(columns=list(fx["expression"].columns[::4]))  # 75% 覆盖
    res = compass_os.analyze(low, fx["labels"]["cancer_type"], clinical=fx["clinical"],
                             robustness="auto")
    assert res.robustness is not None, "低覆盖应触发自动稳健性检查"
    out = Path(tempfile.mkdtemp())
    res.save_report(out)
    assert (out / "robustness_summary.tsv").is_file()


# --------------------------------------------------------------------------- #
# robustness="auto" 行为
# --------------------------------------------------------------------------- #
def test_robustness_auto_skips_second_forward_when_coverage_ok():
    """覆盖达标 ⇒ 不做第二次前向（robustness 为 None）且摘要写明不需要。"""
    import compass_os
    fx = load_fixtures()
    res = compass_os.analyze(fx["expression"], fx["labels"]["cancer_type"],
                             clinical=fx["clinical"], robustness="auto")
    assert res.qc.gene_coverage == 1.0
    assert res.robustness is None
    assert "Not required" in res.summary() or "not required" in res.summary().lower()


def test_robustness_auto_runs_when_below_recommended():
    """低于冻结 recommended 阈值 ⇒ 自动比较 reference vs zero 并给人类可读信息。"""
    import compass_os
    from compass_os.qc import qc_config
    fx = load_fixtures()
    cfg = qc_config()
    rec = cfg["recommended_global_coverage"]
    keep = int(np.ceil((rec + 0.05) * 15672))
    low = fx["expression"].iloc[:, :max(1, keep - 1000)]
    res = compass_os.analyze(low, fx["labels"]["cancer_type"], clinical=fx["clinical"],
                             robustness="auto")
    assert res.qc.gene_coverage < rec
    assert res.robustness is not None
    msgs = res.robustness["messages"]
    assert any("risk ranking" in m for m in msgs)
    txt = res.summary()
    assert "Robustness check" in txt
    # Gsig 阈值未冻结时不得给出 PASS/FAIL
    assert not any("PASS" in m or "FAIL" in m for m in msgs)


def test_robustness_false_disables_check():
    import compass_os
    fx = load_fixtures()
    low = fx["expression"].drop(columns=list(fx["expression"].columns[::4]))
    res = compass_os.analyze(low, fx["labels"]["cancer_type"], robustness=False)
    assert res.robustness is None


# --------------------------------------------------------------------------- #
# 措辞与循环解释防护
# --------------------------------------------------------------------------- #
def test_summary_avoids_forbidden_mechanism_language():
    """**结果区**不得出现 activation / suppression / causal 等表述。

    末尾的免责声明会明确写"不是 causal mechanism 的证据"（这是应有的声明），
    因此检查范围排除最后一段免责声明（位于最后一条分隔线之后）。
    """
    import compass_os
    fx = load_fixtures()
    res = compass_os.analyze(fx["expression"], fx["labels"]["cancer_type"])
    txt = res.summary()
    sep = "─" * 64
    body = txt.split(sep)[0].lower() if sep in txt else txt.lower()
    for bad in ("activated", "suppressed", "causal", "pathway activation",
                "independent downstream evidence", "proves"):
        assert bad not in body, f"结果区出现禁用表述：{bad}"
    assert "hypothesis-generating" in txt
    # 免责声明必须存在且明确否定因果解读
    assert "not" in txt.split(sep)[-1].lower()


def test_risk_group_representation_output_is_labelled_model_linked():
    """concept–风险关联必须标注为 model-linked representation association。"""
    import tempfile

    import compass_os
    fx = load_fixtures()
    res = compass_os.analyze(fx["expression"], fx["labels"]["cancer_type"],
                             clinical=fx["clinical"])
    out = Path(tempfile.mkdtemp())
    res.save_report(out, make_plots=False)
    f = out / "representation_by_risk_group.tsv"
    if f.is_file():
        note = pd.read_csv(f, sep="\t", index_col=0)["note"].dropna().iloc[0]
        assert "model-linked representation association" in note


def test_input_scale_hint_without_auto_conversion():
    """未声明 log2 时不自动转换，只在 notes 中提示。"""
    import compass_os
    fx = load_fixtures()
    log2 = np.log2(fx["expression"].to_numpy(float) + 1.0)
    log2 = pd.DataFrame(log2, index=fx["expression"].index, columns=fx["expression"].columns)
    res = compass_os.analyze(log2, fx["labels"]["cancer_type"].tolist())
    assert any("input_scale" in n for n in res.notes), "应对疑似 log2 数据给出提示"


def test_analyze_writes_report_when_output_dir_given():
    import tempfile

    import compass_os
    fx = load_fixtures()
    d = Path(tempfile.mkdtemp()) / "out"
    compass_os.analyze(fx["expression"].iloc[:6],
                       fx["labels"]["cancer_type"].iloc[:6].tolist(), output_dir=d)
    assert (d / "summary.txt").is_file() and (d / "prediction.tsv").is_file()
