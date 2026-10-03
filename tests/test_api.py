# -*- coding: utf-8 -*-
"""API 行为测试：输入校验、QC 字段、单样本不伪造分层、癌种码口径、robustness 接口。"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from _util import fixtures_available, load_fixtures

pytestmark = pytest.mark.skipif(not fixtures_available(), reason="需要 tests/fixtures")


def test_qc_fields_present():
    """每个预测都返回完整 QC（绝不静默填补）。"""
    import compass_os
    fx = load_fixtures()
    sub = fx["expression"].drop(columns=list(fx["expression"].columns[::10]))
    res = compass_os.predict(sub, fx["labels"]["cancer_type"], clinical=fx["clinical"],
                             model="M2", missing_gene_strategy="reference")
    d = res.qc.to_dict()
    for k in ("n_required_genes", "n_observed_genes", "n_missing_genes", "gene_coverage",
              "missing_gene_strategy", "missing_gene_names", "n_samples"):
        assert k in d, f"QC 缺字段 {k}"
    assert d["n_required_genes"] == 15672
    assert d["n_missing_genes"] == len(list(fx["expression"].columns[::10]))
    assert d["missing_gene_strategy"] == "reference"
    assert res.qc.signature_gene_coverage is not None
    assert res.qc.signature_gene_coverage.shape[1] == 132
    assert res.qc.concept_input_coverage is not None
    assert res.qc.concept_input_coverage.shape[1] == 43


def test_signature_coverage_bounds_and_naming():
    """signature/concept 覆盖度 ∈ [0,1]；concept 层只能叫 input_coverage。"""
    import compass_os
    fx = load_fixtures()
    drop = list(fx["expression"].columns[::5])
    res = compass_os.predict(fx["expression"].drop(columns=drop),
                             fx["labels"]["cancer_type"], clinical=fx["clinical"], model="M2")
    sig = res.qc.signature_gene_coverage.to_numpy()
    con = res.qc.concept_input_coverage.to_numpy()
    assert np.nanmin(sig) >= 0 and np.nanmax(sig) <= 1
    assert np.nanmin(con) >= 0 and np.nanmax(con) <= 1
    # 完全缺失的基因集 ⇒ 覆盖 0；全观测 ⇒ 1
    gs = compass_os.qc.gene_sets()
    for j, gl in enumerate(gs["gene_list"]):
        expected = 1.0 - len([g for g in gl if g in set(drop)]) / len(gl)
        assert abs(sig[0, j] - expected) < 1e-12


def test_single_sample_has_no_stratification():
    """单样本不得伪造 cohort high/low；仍给出风险与表示。"""
    import compass_os
    fx = load_fixtures()
    one = fx["expression"].iloc[[0]]
    res = compass_os.predict(one, [fx["labels"]["cancer_type"].iloc[0]],
                             clinical=fx["clinical"].iloc[[0]], model="M2")
    assert res.risk.shape == (1, 1)
    assert res.risk_group is None or res.risk_group.shape[0] == 1
    assert res.signature_scores.shape == (1, 132)
    assert res.concept_scores.shape == (1, 43)
    assert np.isfinite(res.default_risk.iloc[0])


def test_cohort_risk_rank_and_group():
    """队列 ≥ min_cohort_for_stratification：给出秩与按**冻结切点**的 high/low。"""
    import compass_os
    from compass_os.survival import median_cutoff
    fx = load_fixtures()
    res = compass_os.predict(fx["expression"], fx["labels"]["cancer_type"],
                             clinical=fx["clinical"], model="M2",
                             min_cohort_for_stratification=2)
    assert res.risk_rank.shape == res.risk.shape
    cut = median_cutoff("M2")
    expect = np.where(res.risk["M2"].to_numpy() >= cut, "high", "low")
    assert (res.risk_group["M2"].to_numpy() == expect).all()


def test_cancer_type_required_and_validated():
    """癌种必需且必须在冻结码表内；不得静默推断。"""
    import compass_os
    from compass_os.exceptions import UnknownCancerTypeError
    fx = load_fixtures()
    with pytest.raises(UnknownCancerTypeError):
        compass_os.predict(fx["expression"], ["NOT_A_CANCER"] * len(fx["expression"]))
    with pytest.raises(UnknownCancerTypeError):
        compass_os.predict(fx["expression"], [None] * len(fx["expression"]))
    with pytest.raises(Exception):
        compass_os.predict(fx["expression"], ["LUAD"])  # 长度不匹配


def test_real_cancer_code_used_not_placeholder():
    """癌种码必须真实：同一表达、不同癌种 ⇒ concept/风险不同（占位码会让它们相同）。"""
    import compass_os
    fx = load_fixtures()
    expr = fx["expression"].iloc[[0]]
    row = fx["clinical"].iloc[[0]]
    a = compass_os.get_representation(expr, ["LUAD"])
    b = compass_os.get_representation(expr, ["BRCA"])
    d = float(np.abs(a.concept_scores.to_numpy() - b.concept_scores.to_numpy()).max())
    assert d > 0, "癌种不同但概念完全相同 ⇒ 疑似使用了占位癌种码"


def test_input_scale_handling():
    """``log2_tpm1`` 输入与等价 TPM 输入结果一致；负值/重复列被拒绝。"""
    import compass_os
    from compass_os.exceptions import InputError
    fx = load_fixtures()
    expr = fx["expression"].iloc[:3]
    log2 = np.log2(expr.to_numpy() + 1.0)
    log2 = pd.DataFrame(log2, index=expr.index, columns=expr.columns)
    r_tpm = compass_os.predict(expr, fx["labels"]["cancer_type"].iloc[:3],
                               clinical=fx["clinical"].iloc[:3], model="M2")
    r_log = compass_os.predict(log2, fx["labels"]["cancer_type"].iloc[:3],
                               clinical=fx["clinical"].iloc[:3], model="M2",
                               input_scale="log2_tpm1")
    d = float(np.abs(r_tpm.risk["M2"].to_numpy() - r_log.risk["M2"].to_numpy()).max())
    assert d < 1e-6, f"TPM 与 log2 输入口径不一致：max|Δrisk|={d:.3e}"

    bad = expr.copy()
    bad.iloc[0, 0] = -1.0
    with pytest.raises(InputError):
        compass_os.predict(bad, fx["labels"]["cancer_type"].iloc[:3])
    dup = expr.copy()                      # 真正制造重复列名：把第 2 列改名成与第 1 列相同
    dup.columns = [dup.columns[0], dup.columns[0]] + list(dup.columns[2:])
    assert dup.columns.duplicated().any(), "测试自身未造出重复列"
    with pytest.raises(InputError):
        compass_os.predict(dup, fx["labels"]["cancer_type"].iloc[:3])


def test_clinical_imputation_reported():
    """临床缺失用冻结参考值填充，并在 warnings 中报告（不静默）。"""
    import compass_os
    fx = load_fixtures()
    no_clin = compass_os.predict(fx["expression"], fx["labels"]["cancer_type"],
                                 clinical=None, model="M2")
    assert any("Clinical variables filled" in w for w in no_clin.warnings)
    assert all("reference value" in w for w in no_clin.warnings
               if "Clinical variables filled" in w)


def test_check_robustness_contract():
    """check_robustness 返回连续指标；无二元 pass/fail（v1.0.1 起为 continuous_only）。"""
    import compass_os
    fx = load_fixtures()
    r = compass_os.check_robustness(fx["expression"].iloc[:6],
                                    fx["labels"]["cancer_type"].iloc[:6],
                                    fx["clinical"].iloc[:6], model="M2")
    assert r.robustness_flag == "continuous_only"
    assert "Coverage QC tiers are calibrated" in r.note or "校准" in r.note
    assert set(r.risk) == {"reference", "zero"}
    assert len(r.risk_difference) == 6
    assert "reference_vs_zero" in r.concept_correlation
    assert "reference_vs_zero" in r.signature_correlation
    assert 0.0 <= r.gene_coverage <= 1.0
    assert r.concept_input_coverage["reference"].shape == (6, 43)


def test_model_info_and_manifest_default():
    """model_info 暴露特征组成；默认模型来自 manifest。"""
    from compass_os.api import default_model, model_info
    info = model_info("default")
    assert info["name"] == default_model() == "M2"
    assert info["is_default"] is True
    assert "43 COMPASS concepts" in info["features"]
    m3 = model_info("M3")
    assert m3["is_default"] is False
    assert "PC1-PC10" in m3["features"]


def test_survival_probability_is_monotone_and_bounded():
    """S(t) ∈ (0,1] 且随 t 单调不增（派生输出；校准限制见文档）。"""
    import compass_os
    fx = load_fixtures()
    res = compass_os.predict(fx["expression"].iloc[:4], fx["labels"]["cancer_type"].iloc[:4],
                             clinical=fx["clinical"].iloc[:4], model="M2")
    S = res.survival_probability([0, 365, 730, 1825])
    v = S.to_numpy()
    assert np.all(v > 0) and np.all(v <= 1.0)
    assert np.all(np.diff(v, axis=1) <= 1e-12), "生存曲线非单调不增"
    assert np.allclose(v[:, 0], 1.0)
