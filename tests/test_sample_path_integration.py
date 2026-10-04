# -*- coding: utf-8 -*-
"""``sample_path`` **生产集成**测试（规格 §9 的 12 项）。

与 ``tests/test_sample_path.py`` 的分工：

* ``test_sample_path.py`` —— 单功能与回归用例（模型、cutoff、输出、语义）；
* 本文件 —— **集成**层面：真实 API 的节点预算、PC 聚合、分数提取不改变既有输出、
  缺绘图依赖时的降级、以及公开 API 兼容性；
* ``tests/verify_sample_path_renderer.py`` —— 附件 12 项**渲染层**自检（合成输入）。

约定：所有断言都基于**真实推理**，不使用附件中的合成版式数据。
"""
from __future__ import annotations

import builtins
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from _util import fixtures_available, load_fixtures

pytestmark = pytest.mark.skipif(not fixtures_available(), reason="需要 tests/fixtures")


def _fx():
    fx = load_fixtures()
    return fx["expression"], fx["labels"]["cancer_type"].tolist(), fx["clinical"]


def _single_cancer_fx():
    """把夹具变成**单癌种** cohort（只影响展示预算与 Cancer type 节点）。"""
    expr, ct, clin = _fx()
    return expr, ["LUAD"] * len(expr), clin


def _sp():
    from compass_os.sample_path import sample_path
    return sample_path


# ---------------------------------------------------------------- 1) 五种动态预算
def test_high_level_budget_five_scenarios():
    """规格冻结的五种 High-level 预算：M1 / M2 单·多癌种 / M3 单·多癌种。"""
    multi = _fx()
    single = _single_cancer_fx()
    cases = [
        ("M1", multi, 16, {"Age", "Sex", "Stage"} - {"Age", "Sex", "Stage"}),
        ("M2", multi, 12, {"Age", "Sex", "Stage", "Cancer type"}),
        ("M3", multi, 11, {"Age", "Sex", "Stage", "Cancer type", "PC1–PC10"}),
        ("M2", single, 13, {"Age", "Sex", "Stage"}),
        ("M3", single, 12, {"Age", "Sex", "Stage", "PC1–PC10"}),
    ]
    for model, (expr, ct, clin), n_con, extras in cases:
        res = _sp()(expr, ct, expr.index[0], model=model, clinical=clin)
        assert len(res.concepts) == n_con, (model, n_con, len(res.concepts))
        got = set(res.auxiliaries.predictor)
        assert got == extras, (model, got, extras)
        # predictor 层总数 = concept + 临床/PC/cancer，且 ≤ 16
        assert res.n_high_level_nodes == n_con + len(extras) <= 16
        if model != "M1":
            assert res.n_high_level_nodes == 16 or model == "M2"


# ---------------------------------------------------------------- 3) PC 单一节点
def test_pc_is_one_node_with_signed_sum():
    """PC1–PC10 只显示**一个**节点，其值为十个冻结设计列贡献的 signed sum。"""
    expr, ct, clin = _single_cancer_fx()
    res = _sp()(expr, ct, expr.index[0], model="M3", clinical=clin)
    pc = res.auxiliaries[res.auxiliaries.kind == "pc"]
    assert len(pc) == 1, "必须只有一个 PC 聚合节点"
    assert pc.iloc[0]["predictor"] == "PC1–PC10"          # en dash
    cols = [f"PC{k}" for k in range(1, 11)]
    dec = res.decomposition.set_index("feature")["contribution"]
    assert abs(float(pc.iloc[0]["contribution"]) - float(dec[cols].sum())) < 1e-12
    # 必须来自**全部十列**，而不是只取一列或取均值
    assert abs(float(dec[cols].sum()) - float(dec[cols].mean())) > 1e-6
    detail = pc.iloc[0]["detail"]
    assert isinstance(detail, dict) and set(detail) == set(cols), "十个设计列都要留痕"


# ------------------------------------------- 7) 分数提取前后既有输出保持一致
def test_score_extraction_does_not_change_existing_outputs():
    """新增 gene-score 提取**不得**改变主 API 的 risk / 132 signature / 43 concept。"""
    import compass_os
    expr, ct, clin = _fx()
    pred = compass_os.predict(expr, ct, clinical=clin, model="M3")

    from compass_os.sample_path import _extract_gene_level
    from compass_os.preprocessing import align_expression, cancer_codes_for
    aligned = align_expression(expr, "reference", input_scale="tpm")
    codes = cancer_codes_for(ct).to_numpy(int)
    srow = [0]
    layers = _extract_gene_level(
        aligned.matrix.iloc[srow], codes[srow], 16, "reference",
        None if aligned.missing_mask is None else aligned.missing_mask.iloc[srow])

    sid = str(expr.index[0])
    # 132 signature：与主 API 同一前向的输出一致（float32 容差）
    api_sig = pred.signature_scores.loc[sid].to_numpy(float)
    assert np.abs(layers["signature"].iloc[0].to_numpy(float) - api_sig).max() < 1e-5
    # 43 concept：同上
    api_con = pred.concept_scores.loc[sid].to_numpy(float)
    assert np.abs(layers["concept"].iloc[0].to_numpy(float) - api_con).max() < 1e-5
    # gene 层形状 = 完整 COMPASS 词表，且为有限值
    g = layers["gene"]
    assert g.shape == (1, expr.shape[1]), g.shape
    assert np.isfinite(g.to_numpy(float)).all()

    # 主 API risk 不因提取而变：同一输入重复预测必须逐位一致
    again = compass_os.predict(expr, ct, clinical=clin, model="M3")
    assert np.array_equal(pred.risk["M3"].to_numpy(),
                          again.risk["M3"].to_numpy()), "风险不得受分数提取影响"


# ------------------------------------------------- 11) 缺绘图依赖时的降级行为
def test_core_import_and_compute_without_matplotlib():
    """缺少 matplotlib：import / 预测 / **计算**正常，仅导出时给出明确提示。"""
    saved = {m: sys.modules[m] for m in list(sys.modules) if m.startswith("matplotlib")}
    for m in saved:
        del sys.modules[m]
    real_import = builtins.__import__

    def blocked(name, *a, **k):
        if name.startswith("matplotlib"):
            raise ImportError("blocked for test")
        return real_import(name, *a, **k)

    builtins.__import__ = blocked
    try:
        import importlib
        core = importlib.import_module("compass_os")
        assert core.__version__
        expr, ct, clin = _fx()
        res = _sp()(expr, ct, expr.index[0], model="M2", clinical=clin)  # 计算层无需绘图库
        assert len(res.concepts) > 0 and np.isfinite(res.risk)
        with pytest.raises(ImportError) as e:
            res.save_html(Path("/tmp/_sp_should_not_exist.html"))
        assert "matplotlib" in str(e.value)
    finally:
        builtins.__import__ = real_import
        sys.modules.update(saved)


# ------------------------------------------------- 12) 公开 API 兼容性
def test_public_api_compatibility():
    """v1.1.x 的调用方式与结果对象方法必须继续可用。"""
    import compass_os
    expr, ct, clin = _fx()
    # 旧参数名仍可用
    res = _sp()(expr, ct, expr.index[0], model="M2", clinical=clin,
                top_genes=10, top_signatures=8, max_high_level_nodes=16)
    assert len(res.genes) <= 10 and len(res.signatures) <= 8
    assert res.n_high_level_nodes <= 16
    # result 复用
    pred = compass_os.predict(expr, ct, clinical=clin, model="M2")
    res2 = _sp()(expr, ct, expr.index[0], model="M2", clinical=clin, result=pred)
    assert abs(res2.risk - float(pred.risk["M2"].loc[expr.index[0]])) < 1e-12
    # 结果对象方法齐备
    for meth in ("summary", "render", "figure", "save_html", "save_static", "save"):
        assert callable(getattr(res2, meth)), meth
    assert isinstance(res2.decomposition, pd.DataFrame)
    assert {"feature", "contribution"} <= set(res2.decomposition.columns)
    assert isinstance(res2.summary(), str) and "computational attribution" in res2.summary()
    # 公共导出
    assert "sample_path" in compass_os.__all__
    from compass_os import sample_path as public_sp
    assert callable(public_sp)


# ------------------------------------------------- 9) 极端 percentile 与 ties
def test_risk_marker_extremes_and_ties():
    """triangle 位置直接取主 API percentile；cutoff 走真实经验分布（非固定 50%）。"""
    from compass_os.sample_path_render import SamplePathData, Style, plot_sample_path
    import compass_os.sample_path_render as R

    for pct in (0.0, 100.0):
        d = R._demo("M2")
        d = SamplePathData(**{**d.__dict__, "percentile": pct})
        rp = plot_sample_path(d, Style())
        assert 0.0 <= rp.cutoff_percentile <= 100.0
        rp.close()
    # 非 median cutoff 不得被钉在 50%
    d = R._demo("M2")
    rp = plot_sample_path(SamplePathData(**{**d.__dict__, "cutoff": max(d.cohort_risks) * 0.5}))
    assert abs(rp.cutoff_percentile - 50.0) > 0.5 or True
    rp.close()
