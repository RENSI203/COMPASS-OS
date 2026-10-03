# -*- coding: utf-8 -*-
"""Golden 复现测试：132 / 43 / M0–M3 必须与原项目生产结果一致（max|Δ| ≤ 1e-5）。

期望值来源见 ``tests/fixtures/golden_provenance.json``：
43 concepts = 原项目 ``04`` 产物；132 = 生产同一调用重跑；风险 = 原项目
``12_validate_external.py::_build_features`` + lock(β/scaler)。
"""
from __future__ import annotations

import numpy as np
import pytest

from _util import TOL_GOLDEN, TOL_ZERO_EQ, fixtures_available, load_fixtures, max_abs_diff

pytestmark = pytest.mark.skipif(not fixtures_available(), reason="需要 tests/fixtures（见 tools/build_fixtures.py）")


def _rep(**kw):
    import compass_os
    fx = load_fixtures()
    return compass_os.get_representation(fx["expression"], fx["labels"]["cancer_type"], **kw)


# --------------------------------------------------------------------------- #
# A. 表示层
# --------------------------------------------------------------------------- #
def test_signature_scores_132_match_production():
    """132 gene-signature 分数（真实癌种码）对生产重跑结果 max|Δ| ≤ 1e-5。"""
    fx = load_fixtures()
    got = _rep(missing_gene_strategy="strict").signature_scores
    d = max_abs_diff(got, fx["signature_scores"])
    assert d <= TOL_GOLDEN, f"132 signatures max|Δ|={d:.3e} > {TOL_GOLDEN:g}"


def test_concept_scores_43_match_archive():
    """43 concept 分数对原项目 ``04`` 档案（真实癌种码）max|Δ| ≤ 1e-5。"""
    fx = load_fixtures()
    got = _rep(missing_gene_strategy="strict").concept_scores
    d = max_abs_diff(got, fx["concept_scores"])
    assert d <= TOL_GOLDEN, f"43 concepts max|Δ|={d:.3e} > {TOL_GOLDEN:g}"


def test_signature_scores_code5_match_archive():
    """132 层独立复核：占位癌种码 5 下对原项目 ``41`` 档案 max|Δ| ≤ 1e-5。

    本包默认用真实癌种码；此测试直接驱动底层 ``compute_representation`` 传入 code=5，
    以验证 132 层的**列序与取值**与原项目实现一致（与癌种码口径无关的部分）。
    """
    import compass_os.representation as rep
    from compass_os.preprocessing import align_expression
    fx = load_fixtures()
    aligned = align_expression(fx["expression"], "strict")
    codes = np.full(aligned.matrix.shape[0], 5, dtype=int)
    out = rep.compute_representation(aligned.matrix, codes, strategy="strict")
    d = max_abs_diff(out["signature_scores"], fx["signature_scores_code5"])
    assert d <= TOL_GOLDEN, f"132 signatures（code5）max|Δ|={d:.3e} > {TOL_GOLDEN:g}"


def test_concept_aggregation_identity():
    """结构一致性：132 经**冻结 softmax 权重**聚合应还原 43（float32 精度）。

    这是 132 层的独立正确性证据（列序/取值错了就还原不出 43）。
    """
    import pandas as pd
    from compass_os.representation import concept_weights
    r = _rep(missing_gene_strategy="strict")
    w = concept_weights()
    sig = r.signature_scores
    recon = {}
    for concept, (members, ww) in w.items():
        cols = [m for m in members if m in sig.columns]
        if not cols:
            continue
        wsel = np.asarray(ww)[[members.index(c) for c in cols]]
        wsel = wsel / wsel.sum()
        recon[concept] = sig[cols].to_numpy(float) @ wsel
    rec = pd.DataFrame(recon, index=sig.index)
    d = max_abs_diff(rec, r.concept_scores)
    assert d <= 1e-6, f"softmax 聚合还原 43 的 max|Δ|={d:.3e} > 1e-6"


# --------------------------------------------------------------------------- #
# B. 生存风险 M0–M3
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("model", ["M0", "M1", "M2", "M3"])
def test_risk_matches_production(model):
    """M0–M3 风险对生产 ``_build_features`` + lock 的结果 max|Δ| ≤ 1e-5。"""
    import compass_os
    fx = load_fixtures()
    res = compass_os.predict(fx["expression"], fx["labels"]["cancer_type"],
                             clinical=fx["clinical"], model=model,
                             missing_gene_strategy="strict")
    d = max_abs_diff(res.risk[[model]], fx[f"risk_{model}"])
    assert d <= TOL_GOLDEN, f"{model} risk max|Δ|={d:.3e} > {TOL_GOLDEN:g}"


def test_default_model_is_m2():
    """默认模型由 manifest 决定（当前 M2），且四模型均可调用。"""
    from compass_os.api import available_models, default_model
    assert default_model() == "M2"
    assert set(available_models()) == {"M0", "M1", "M2", "M3"}


def test_all_models_in_one_call():
    """一次调用同时给出 M0–M3 + default，且与逐个调用一致。"""
    import compass_os
    fx = load_fixtures()
    res = compass_os.predict(fx["expression"], fx["labels"]["cancer_type"],
                             clinical=fx["clinical"], model="M0,M1,M2,M3")
    for mk in ("M0", "M1", "M2", "M3"):
        assert mk in res.risk.columns
        d = max_abs_diff(res.risk[[mk]], fx[f"risk_{mk}"])
        assert d <= TOL_GOLDEN, f"{mk} max|Δ|={d:.3e}"


# --------------------------------------------------------------------------- #
# C. 缺失基因策略
# --------------------------------------------------------------------------- #
def test_reference_mode_matches_production_rule():
    """``reference`` 用冻结的 log2 中位数填充（≡ TPM 7.76），与生产口径一致。"""
    import compass_os
    from compass_os.preprocessing import (align_expression, reference_fill_log2,
                                          reference_fill_tpm)
    fx = load_fixtures()
    assert abs(reference_fill_log2() - 3.130937933922) < 1e-12
    assert abs(reference_fill_tpm() - 7.76) < 5e-3
    # 去掉 20% 基因后仍能出结果，且 QC 记录缺失
    drop = list(fx["expression"].columns[::5])
    sub = fx["expression"].drop(columns=drop)
    res = compass_os.predict(sub, fx["labels"]["cancer_type"], clinical=fx["clinical"],
                             model="M2", missing_gene_strategy="reference")
    assert res.qc.n_missing_genes == len(drop)
    assert res.qc.missing_gene_strategy == "reference"
    assert abs(res.qc.gene_coverage - (1 - len(drop) / res.qc.n_required_genes)) < 1e-12
    # 填充值必须等于冻结参考值（而不是用户队列中位数）
    aligned = align_expression(sub, "reference")
    assert np.allclose(aligned.matrix[drop].to_numpy(), reference_fill_tpm(), rtol=0, atol=1e-9)


def test_zero_mode_scaled_space_equivalence():
    """``zero`` 的两种等价写法必须一致：A(scaled→0) vs B(TPM=2^data_min−1)，max|Δ| ≤ 1e-6。"""
    import compass_os
    import compass_os.representation as rep
    from compass_os.preprocessing import (align_expression, cancer_codes_for,
                                          zero_equivalent_tpm)
    fx = load_fixtures()
    drop = list(fx["expression"].columns[::5])
    sub = fx["expression"].drop(columns=drop)
    aligned = align_expression(sub, "zero")
    codes = cancer_codes_for(fx["labels"]["cancer_type"]).to_numpy()

    # A：标准化空间置 0（正式实现）
    a = rep.compute_representation(aligned.matrix, codes, missing_mask=aligned.missing_mask,
                                   strategy="zero")

    # B：输入侧等价式 TPM_g = 2^{data_min_g} − 1，再走 reference 之外的"不掩码"路径
    eq = zero_equivalent_tpm(aligned.missing_mask, rep.checkpoint_data_min(),
                             list(aligned.matrix.columns))
    tpm_b = aligned.matrix.copy()
    for g in drop:
        tpm_b[g] = eq[g].to_numpy()
    b = rep.compute_representation(tpm_b, codes, strategy="reference")

    for layer in ("signature_scores", "concept_scores"):
        d = max_abs_diff(a[layer], b[layer])
        assert d <= TOL_ZERO_EQ, f"zero 等价性失败（{layer}）max|Δ|={d:.3e} > {TOL_ZERO_EQ:g}"


def test_zero_differs_from_raw_tpm_zero():
    """证据留存：``zero``（训练空间）与把缺失基因写成 raw TPM=0 **不等价**。

    P0 审计：60.2% 的基因 ``data_min_ > 0``，raw TPM=0 会落到训练范围之外。
    """
    import compass_os
    import compass_os.representation as rep
    from compass_os.preprocessing import align_expression, cancer_codes_for
    fx = load_fixtures()
    drop = list(fx["expression"].columns[::5])
    sub = fx["expression"].drop(columns=drop)
    aligned = align_expression(sub, "zero")
    codes = cancer_codes_for(fx["labels"]["cancer_type"]).to_numpy()
    a = rep.compute_representation(aligned.matrix, codes, missing_mask=aligned.missing_mask,
                                   strategy="zero")
    raw0 = aligned.matrix.copy()
    raw0[drop] = 0.0
    c = rep.compute_representation(raw0, codes, strategy="reference")
    d = max_abs_diff(a["concept_scores"], c["concept_scores"])
    assert d > TOL_ZERO_EQ, ("期望 raw TPM=0 与训练空间 zero 不同，实测 max|Δ|="
                             f"{d:.3e}；若确为 0 则说明该夹具无法区分两种口径")


def test_strict_mode_raises():
    """``strict`` 缺任一必需基因即抛 MissingGenesError（用于 benchmark/严格复现）。"""
    import compass_os
    from compass_os.exceptions import MissingGenesError
    fx = load_fixtures()
    drop = list(fx["expression"].columns[:3])
    sub = fx["expression"].drop(columns=drop)
    try:
        compass_os.predict(sub, fx["labels"]["cancer_type"], clinical=fx["clinical"],
                           missing_gene_strategy="strict")
    except MissingGenesError as exc:
        assert exc.n_missing == 3
        assert set(exc.missing_genes) == set(drop)
        return
    raise AssertionError("strict 模式未抛出 MissingGenesError")


# --------------------------------------------------------------------------- #
# D. 官方 API 与生产（predict + capture）路径的逐位一致性
# --------------------------------------------------------------------------- #
def test_official_extract_matches_production_path():
    """正式包用官方 ``extract()``；其 132/43 必须与**生产路径产物**一致。

    夹具的 132/43 由原项目生产路径生成（``predict`` + ``enable_geneset_capture``，
    见 ``tools/build_fixtures.py``）。判据 1e-5；实测为 **0.0（逐位相同）**，
    因此改用官方 API 并移除 monkey patch 不改变任何数值结果。
    """
    import numpy as np
    fx = load_fixtures()
    r = _rep(missing_gene_strategy="strict")
    d_con = max_abs_diff(r.concept_scores, fx["concept_scores"])
    d_sig = max_abs_diff(r.signature_scores, fx["signature_scores"])
    assert d_con <= TOL_GOLDEN and d_sig <= TOL_GOLDEN, (
        f"官方 extract 与生产路径不一致：43 max|Δ|={d_con:.3e}，132 max|Δ|={d_sig:.3e}")
    # 实测残差量级 = 存档的 float32 十进制文本往返（43：~2.4e-7；132：~1e-16），
    # 远严于 1e-5 判据。注意：**包内** extract 与 predict+capture 直接对比为 0.0（逐位相同），
    # 这里比的是"包 vs 原项目存档产物"，残差来自存档精度而非算法差异。
    assert d_con <= 1e-6 and d_sig <= 1e-6, (
        f"残差超过 float32 存档精度（43={d_con:.3e}, 132={d_sig:.3e}）⇒ 需排查算法差异")
    assert np.isfinite(r.signature_scores.to_numpy()).all()
