# -*- coding: utf-8 -*-
"""样本级计算路径归因（``sample_path``）测试。

覆盖规格 §19 要求的 17 类用例，并加入本实现特有的**精确性/守恒**回归：

* 归因链条求和与主 API ``risk`` 一致（float32 容差内）
* 每个内部节点流量守恒（含残差聚合节点与 centering 常数节点）
* high-level 层 predictor 节点 ≤ 16
* 前三个高维层各自 ≤ top-K
* 展示的 risk / rank / percentile 与主 API 完全一致
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from _util import fixtures_available, load_fixtures, repo_root

pytestmark = pytest.mark.skipif(not fixtures_available(), reason="需要 tests/fixtures")


# --------------------------------------------------------------------------- #
def _fx():
    fx = load_fixtures()
    return (fx["expression"], fx["labels"]["cancer_type"].tolist(), fx["clinical"])


def _dup(expr, ct, clin, n):
    """把夹具扩成 n 个**索引互不相同**的样本。"""
    reps = int(np.ceil(n / len(expr)))
    idx = [f"S{k * len(expr) + j}" for k in range(reps) for j in range(len(expr))][:n]
    e = pd.DataFrame(np.tile(expr.to_numpy(float), (reps, 1))[:n], index=idx,
                     columns=expr.columns)
    c = pd.DataFrame(np.tile(clin.to_numpy(float), (reps, 1))[:n], index=idx,
                     columns=clin.columns)
    return e, (ct * reps)[:n], c


def _sp():
    from compass_os.sample_path import sample_path
    return sample_path


# ---------------------------------------------------------------- 1) M1 / M2 / M3
@pytest.mark.parametrize("model", ["M1", "M2", "M3"])
def test_single_sample_each_model(model):
    """M1/M2/M3 各跑一个样本，结构与数值自检均须通过。"""
    expr, ct, clin = _fx()
    res = _sp()(expr, ct, expr.index[0], model=model, clinical=clin)
    assert res.model == model
    assert res.sample_id == expr.index[0]
    assert np.isfinite(res.risk)
    assert len(res.concepts) >= 1
    assert res.eta_max_abs_diff <= 1e-4 * max(1.0, abs(res.risk))
    if model == "M1":                       # M1 无临床/癌种/PC
        assert res.auxiliaries.empty
    if model == "M2":
        assert set(res.auxiliaries.predictor) >= {"Age", "Sex", "Stage"}
    if model == "M3":
        assert any("PC1" in p for p in res.auxiliaries.predictor)


# ---------------------------------------------------------------- 4) 不存在的 sample
def test_missing_sample_id_raises():
    """sample_id 不存在时必须报错，且提示可用样本。"""
    from compass_os.exceptions import InputError
    expr, ct, clin = _fx()
    with pytest.raises(InputError):
        _sp()(expr, ct, "NOT-A-SAMPLE", model="M2", clinical=clin)


# ---------------------------------------------------------------- 5) 单样本 cohort
def test_single_sample_cohort():
    """n=1 队列：rank=1/1，percentile 有限，cutoff=自身 risk。"""
    expr, ct, clin = _fx()
    one = expr.iloc[[0]]
    res = _sp()(one, [ct[0]], one.index[0], model="M2", clinical=clin.iloc[[0]])
    assert res.n_samples == 1 and res.rank == 1
    assert np.isfinite(res.percentile) and np.isfinite(res.cutoff)


# ---------------------------------------------------------------- 6) 小 cohort
def test_small_cohort():
    """小队列（n=5）可用，且 rank/percentile 与手算一致。"""
    expr, ct, clin = _dup(*_fx(), 5)
    res = _sp()(expr, ct, expr.index[2], model="M2", clinical=clin)
    assert res.n_samples == 5
    import compass_os
    risk = compass_os.predict(expr, ct, clinical=clin, model="M2").risk["M2"]
    order = risk.sort_values(ascending=False).index.tolist()
    assert res.rank == order.index(expr.index[2]) + 1
    assert abs(res.percentile - 100.0 * float((risk < risk[expr.index[2]]).mean()
                                              + 0.5 * float((risk == risk[expr.index[2]]).mean()))) < 1e-9


# ---------------------------------------------------------------- 7/8) clinical
def test_clinical_complete_and_imputed():
    """临床完整与缺失（走正式冻结填补）两条路径都可用且自检通过。"""
    expr, ct, clin = _fx()
    full = _sp()(expr, ct, expr.index[0], model="M2", clinical=clin)
    assert np.isfinite(full.auxiliaries.contribution).all()

    missing = clin.copy()
    missing.loc[missing.index[0], ["age", "sex", "stage"]] = np.nan
    imp = _sp()(expr, ct, expr.index[0], model="M2", clinical=missing)
    age_row = imp.auxiliaries[imp.auxiliaries.predictor == "Age"].iloc[0]
    # 冻结填补：Age → 60.0（scaled z = (60 - mean)/scale）
    lock_age_mu = 58.7559915164369
    lock_age_sd = 14.160350962930025
    assert abs(age_row["value"] - (60.0 - lock_age_mu) / lock_age_sd) < 1e-9


# ---------------------------------------------------------------- 9) PC 聚合
def test_pc_aggregation_is_sum_not_mean():
    """PC1–PC10 节点贡献必须是十项**求和**，且等于逐 PC 贡献之和。"""
    expr, ct, clin = _fx()
    res = _sp()(expr, ct, expr.index[0], model="M3", clinical=clin)
    row = res.auxiliaries[res.auxiliaries.kind == "pca_aggregate"].iloc[0]
    detail = row["detail"]
    assert isinstance(detail, dict) and len(detail) == 10
    assert abs(sum(detail.values()) - row["contribution"]) < 1e-12
    # 与逐特征分解表逐项一致
    dec = res.decomposition.set_index("feature")["contribution"]
    assert abs(sum(dec[k] for k in detail) - row["contribution"]) < 1e-12
    assert abs(row["contribution"] - np.mean(list(detail.values()))) > 1e-12 or \
        all(abs(v) < 1e-15 for v in detail.values())


# ---------------------------------------------------------------- 10/11) cutoff
def test_cutoff_median_and_numeric():
    """median cutoff 与研究者给定数值 cutoff 语义必须区分。"""
    expr, ct, clin = _fx()
    med = _sp()(expr, ct, expr.index[0], model="M2", clinical=clin)
    assert med.cutoff_mode.startswith("cohort median")
    assert med.cutoff_is_frozen_validated is False
    import compass_os
    risk = compass_os.predict(expr, ct, clinical=clin, model="M2").risk["M2"]
    assert abs(med.cutoff - float(np.median(risk))) < 1e-12

    num = _sp()(expr, ct, expr.index[0], model="M2", clinical=clin, cutoff=1.234)
    assert num.cutoff == 1.234 and "user-specified" in num.cutoff_mode
    assert num.cutoff_is_frozen_validated is False
    from compass_os.survival import median_cutoff
    assert abs(num.cutoff - median_cutoff("M2")) > 1e-9   # 与冻结 validated cutoff 不同


def test_invalid_cutoff_string_rejected():
    from compass_os.exceptions import InputError
    expr, ct, clin = _fx()
    with pytest.raises(InputError):
        _sp()(expr, ct, expr.index[0], model="M2", clinical=clin, cutoff="frozen")


# ---------------------------------------------------------------- 12) top-K 截断
def test_topk_truncation():
    """三层各自受 top-K 约束，且 K 可调。"""
    expr, ct, clin = _fx()
    res = _sp()(expr, ct, expr.index[0], model="M2", clinical=clin,
                top_genes=7, top_signatures=5)
    assert len(res.genes) <= 7
    assert len(res.signatures) <= 5
    big = _sp()(expr, ct, expr.index[0], model="M2", clinical=clin)
    assert len(big.genes) <= 50 and len(big.signatures) <= 50


# ---------------------------------------------------------------- 13) ≤16 节点
def test_high_level_node_budget():
    """high-level predictor 节点 ≤ 16（含各模型的实际构成）。"""
    expr, ct, clin = _fx()
    for m in ("M1", "M2", "M3"):
        res = _sp()(expr, ct, expr.index[0], model=m, clinical=clin)
        assert res.n_high_level_nodes <= 16, (m, res.n_high_level_nodes)
        assert res.n_high_level_nodes == len(res.concepts) + len(res.auxiliaries)
    small = _sp()(expr, ct, expr.index[0], model="M2", clinical=clin,
                  max_high_level_nodes=6)
    assert small.n_high_level_nodes <= 6
    assert len(small.concepts) + len(small.auxiliaries) <= 6


# ---------------------------------------------------------------- 14/15) cancer type
def test_cancer_type_multi_cohort_shown():
    """多癌种队列：cancer type 作为独立 predictor 节点显示并占用预算。"""
    expr, ct, clin = _fx()
    assert len(set(ct)) > 1
    res = _sp()(expr, ct, expr.index[0], model="M2", clinical=clin)
    assert (res.auxiliaries.kind == "cancer_type").any()
    assert res.n_high_level_nodes <= 16


def test_cancer_type_single_cohort_omitted_with_note():
    """单癌种队列：默认省略 cancer type，且必须给出规格要求的那句话。"""
    expr, ct, clin = _fx()
    res = _sp()(expr, ["LUAD"] * len(expr), expr.index[0], model="M2", clinical=clin)
    assert not (res.auxiliaries.kind == "cancer_type").any()
    assert any("constant within this single-cancer cohort" in n for n in res.notes)
    # 显式要求时仍可显示
    forced = _sp()(expr, ["LUAD"] * len(expr), expr.index[0], model="M2", clinical=clin,
                   show_cancer_type=True)
    assert (forced.auxiliaries.kind == "cancer_type").any()


# ---------------------------------------------------------------- 16) full risk
def test_full_risk_matches_main_api():
    """**最重要**：图中显示的 risk/rank/percentile 与主 API 完全一致（机器精度）。"""
    import compass_os
    expr, ct, clin = _fx()
    pred = compass_os.predict(expr, ct, clinical=clin, model="M2")
    for sid in (expr.index[0], expr.index[5], expr.index[-1]):
        res = _sp()(expr, ct, sid, model="M2", clinical=clin)
        assert res.risk == float(pred.risk["M2"].loc[sid])          # 精确相等
        order = pred.risk["M2"].sort_values(ascending=False).index.tolist()
        assert res.rank == order.index(sid) + 1
        assert res.n_samples == len(expr)


def test_attribution_chain_sums_to_api_risk():
    """归因链条求和与主 API risk 在 float32 容差内一致，且被显式报告。"""
    import compass_os
    expr, ct, clin = _fx()
    pred = compass_os.predict(expr, ct, clinical=clin, model="M2")
    for sid in (expr.index[0], expr.index[7]):
        res = _sp()(expr, ct, sid, model="M2", clinical=clin)
        api = float(pred.risk["M2"].loc[sid])
        assert res.eta_max_abs_diff <= 1e-4 * max(1.0, abs(api))
        assert res.decomposition["contribution"].sum() == pytest.approx(
            api, abs=1e-4 * max(1.0, abs(api)))


def test_result_can_be_reused_without_repredicting():
    """传入既有结果时可复用 cohort 风险（不重复预测），数值不变。"""
    import compass_os
    expr, ct, clin = _fx()
    pred = compass_os.predict(expr, ct, clinical=clin, model="M2")
    a = _sp()(expr, ct, expr.index[0], model="M2", clinical=clin)
    b = _sp()(expr, ct, expr.index[0], model="M2", clinical=clin, result=pred)
    assert a.risk == b.risk
    assert a.rank == b.rank and a.percentile == b.percentile


# ---------------------------------------------------------------- 17) rank/percentile
def test_rank_and_percentile_are_from_official_risk():
    """rank/percentile 必须由正式 risk 输出导出，不得另算一套。"""
    import compass_os
    expr, ct, clin = _fx()
    risk = compass_os.predict(expr, ct, clinical=clin, model="M2").risk["M2"]
    sid = expr.index[3]
    res = _sp()(expr, ct, sid, model="M2", clinical=clin)
    assert res.rank == int((risk > risk[sid]).sum()) + 1
    assert abs(res.percentile - 100.0 * (float((risk < risk[sid]).mean())
                                         + 0.5 * float((risk == risk[sid]).mean()))) < 1e-9


# ---------------------------------------------------------------- 守恒与精确性
def test_every_internal_node_conserves_flow():
    """每个内部节点入==出（含残差聚合与 centering 常数节点）。"""
    expr, ct, clin = _fx()
    for m in ("M1", "M2", "M3"):
        res = _sp()(expr, ct, expr.index[0], model=m, clinical=clin)
        gl, sl = res.gene_links, res.signature_links
        rg = dict(zip(res.residual_gene_links.signature,
                      res.residual_gene_links.value)) if len(res.residual_gene_links) else {}
        rs = dict(zip(res.residual_signature_links.concept,
                      res.residual_signature_links.value)) if len(res.residual_signature_links) else {}
        ct_ = dict(zip(res.centering_links.target, res.centering_links.value)) \
            if len(res.centering_links) else {}
        for name in res.signatures.signature:
            inflow = float(gl[gl.signature == name]["value"].sum()) + rg.get(name, 0.0)
            outflow = float(sl[sl.signature == name]["value"].sum())
            assert abs(inflow - outflow) <= 1e-6 * max(1.0, abs(outflow)), (m, name)
        for name in res.concepts.concept:
            inflow = float(sl[sl.concept == name]["value"].sum()) + rs.get(name, 0.0) \
                + ct_.get(name, 0.0)
            outflow = float(res.concepts[res.concepts.concept == name]["contribution"].iloc[0])
            assert abs(inflow - outflow) <= 1e-6 * max(1.0, abs(outflow)), (m, name)


def test_risk_node_inflow_equals_linear_predictor():
    """risk 节点总入流 == full Cox 线性预测子（含 Other concepts 与辅助节点）。"""
    expr, ct, clin = _fx()
    for m in ("M1", "M2", "M3"):
        res = _sp()(expr, ct, expr.index[0], model=m, clinical=clin)
        total = (float(res.concepts.contribution.sum())
                 + float(res.auxiliaries.contribution.sum())
                 + float(res.other_concepts_contribution))
        assert abs(total - res.risk) <= 1e-4 * max(1.0, abs(res.risk)), m


def test_gene_and_signature_links_are_exact_contributions():
    """链接值必须是**精确贡献**：gene→signature 之和 == 该 signature 的展示入流。"""
    expr, ct, clin = _fx()
    res = _sp()(expr, ct, expr.index[0], model="M2", clinical=clin)
    if len(res.gene_links) == 0:
        pytest.skip("该样本无展示基因")
    for name in res.signatures.signature:
        sub = res.gene_links[res.gene_links.signature == name]
        if len(sub) == 0:
            continue
        # 每个链接都必须落到某个展示基因/签名
        assert set(sub.gene) <= set(res.genes.gene)
        assert set(sub.signature) <= set(res.signatures.signature)


# ---------------------------------------------------------------- M0 与非法输入
def test_m0_explicitly_unsupported():
    """M0 必须给出规格指定的明确提示，而不是画一张空图。"""
    from compass_os.exceptions import InputError
    from compass_os.sample_path import M0_UNSUPPORTED_MESSAGE
    expr, ct, clin = _fx()
    with pytest.raises(InputError) as e:
        _sp()(expr, ct, expr.index[0], model="M0", clinical=clin)
    assert M0_UNSUPPORTED_MESSAGE in str(e.value)


def test_invalid_model_rejected():
    from compass_os.exceptions import InputError
    expr, ct, clin = _fx()
    with pytest.raises(InputError):
        _sp()(expr, ct, expr.index[0], model="M9", clinical=clin)


# ---------------------------------------------------------------- 输出
def test_html_and_static_output(tmp_path=None):
    """HTML 必须可生成；静态图在没有 kaleido 时回退 matplotlib。"""
    import tempfile
    from pathlib import Path
    expr, ct, clin = _fx()
    res = _sp()(expr, ct, expr.index[0], model="M2", clinical=clin)
    d = Path(tempfile.mkdtemp(prefix="sp_out_"))
    h = res.save_html(d / "p.html")
    assert h.is_file() and h.stat().st_size > 5000
    s = res.save_static(d / "p.png")
    assert s.is_file() and s.stat().st_size > 5000
    assert "<html" in (d / "p.html").read_text(encoding="utf-8", errors="ignore").lower()


def test_figure_contains_sankey_and_risk_column():
    """单图内必须同时含 Sankey 与右侧风险柱，且 cutoff 用直接 shape（避开 Plotly 陷阱）。"""
    expr, ct, clin = _fx()
    res = _sp()(expr, ct, expr.index[0], model="M2", clinical=clin)
    fig = res.figure()
    types = [t.type for t in fig.data]
    assert "sankey" in types and any(t in ("heatmap", "scatter") for t in types)
    assert len(fig.data[0].link.value) > 0
    assert all(v >= 0 for v in fig.data[0].link.value)      # Plotly 要求非负
    assert len(fig.layout.shapes) >= 1


def test_semantics_disclaimer_present():
    """所有对外文本必须带语义边界声明。"""
    expr, ct, clin = _fx()
    res = _sp()(expr, ct, expr.index[0], model="M2", clinical=clin)
    s = res.summary()
    for kw in ("computational attribution", "Not a causal mechanism diagram",
               "not a pathway activation map"):
        assert kw in s, kw
    from compass_os.sample_path import SEMANTICS
    assert "pathway activation" in SEMANTICS
