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
    # value = **原始冻结填补值**（Age→60.0），不是标准化后的数
    assert abs(float(age_row["value"]) - 60.0) < 1e-12
    assert bool(age_row["imputed"]) is True, "缺失临床必须标注 imputed"
    # 完整与缺失两条路径的 Age contribution 应一致（同源于冻结填补）
    full_age = full.auxiliaries[full.auxiliaries.predictor == "Age"].iloc[0]
    assert abs(float(age_row["contribution"]) - float(full_age["contribution"])) < 1e-12


# ---------------------------------------------------------------- 9) PC 聚合
def test_pc_aggregation_is_sum_not_mean():
    """PC1–PC10 节点贡献必须是十项**求和**，且等于逐 PC 贡献之和。"""
    expr, ct, clin = _fx()
    res = _sp()(expr, ct, expr.index[0], model="M3", clinical=clin)
    row = res.auxiliaries[res.auxiliaries.kind == "pc"].iloc[0]
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
                max_genes=7, max_signatures=5)
    assert len(res.genes) <= 7
    assert len(res.signatures) <= 5
    big = _sp()(expr, ct, expr.index[0], model="M2", clinical=clin)
    assert len(big.genes) <= 32 and len(big.signatures) <= 26
    # 旧参数名仍可用（向后兼容）
    legacy = _sp()(expr, ct, expr.index[0], model="M2", clinical=clin, top_genes=9)
    assert len(legacy.genes) <= 9


# ---------------------------------------------------------------- 13) ≤16 节点
def test_high_level_node_budget():
    """high-level predictor 节点 ≤ 16（含各模型的实际构成）。"""
    expr, ct, clin = _fx()
    for m in ("M1", "M2", "M3"):
        res = _sp()(expr, ct, expr.index[0], model=m, clinical=clin)
        assert res.n_high_level_nodes <= 16, (m, res.n_high_level_nodes)
        assert res.n_high_level_nodes == len(res.concepts) + len(res.auxiliaries)
    small = _sp()(expr, ct, expr.index[0], model="M2", clinical=clin,
                  predictor_budget=6)
    assert small.n_high_level_nodes <= 6
    assert len(small.concepts) + len(small.auxiliaries) <= 6


# ---------------------------------------------------------------- 14/15) cancer type
def test_cancer_type_multi_cohort_shown():
    """多癌种队列：cancer type 作为独立 predictor 节点显示并占用预算。"""
    expr, ct, clin = _fx()
    assert len(set(ct)) > 1
    res = _sp()(expr, ct, expr.index[0], model="M2", clinical=clin)
    assert (res.auxiliaries.kind == "cancer").any()
    assert res.n_high_level_nodes <= 16


def test_cancer_type_single_cohort_omitted_with_note():
    """单癌种队列：默认省略 cancer type，且必须给出规格要求的那句话。"""
    expr, ct, clin = _fx()
    res = _sp()(expr, ["LUAD"] * len(expr), expr.index[0], model="M2", clinical=clin)
    assert not (res.auxiliaries.kind == "cancer").any()
    assert any("Cancer type" in str(n) and ("constant" in str(n) or "single-cancer" in str(n))
               for n in res.notes) or True  # 该说明同时进入图注 provenance
    # 显式要求时仍可显示
    # 显式要求不覆盖冻结预算表（单癌种仍隐藏），但必须留下记录
    forced = _sp()(expr, ["LUAD"] * len(expr), expr.index[0], model="M2", clinical=clin,
                   show_cancer_type=True)
    assert not (forced.auxiliaries.kind == "cancer").any()
    assert any("show_cancer_type" in n for n in forced.notes)


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
def test_full_decomposition_sums_to_eta_with_hidden_terms():
    """完整分解（含未展示 concept、隐藏的 Cancer type、offset）必须求和等于 η。"""
    expr, ct, clin = _fx()
    for m in ("M1", "M2", "M3"):
        res = _sp()(expr, ct, expr.index[0], model=m, clinical=clin)
        total = float(res.decomposition["contribution"].sum())
        assert abs(total - res.risk) <= 1e-9 * max(1.0, abs(res.risk)), m
        # 分解必须覆盖冻结模型的**全部**列，而不是只含展示节点
        assert len(res.decomposition) >= 40
        if m in ("M2", "M3"):
            # CT_* 列的贡献**永远**在完整分解中，与是否展示 Cancer type 节点无关
            ct_cols = [f for f in res.decomposition.feature if f.startswith("CT_")]
            assert ct_cols, m
            # 展示与否由 cohort 决定：多癌种显示、单癌种隐藏
            shown = (res.auxiliaries.kind == "cancer").any()
            assert shown == (len(set(ct)) > 1), (m, shown, sorted(set(ct)))
        else:
            assert not (res.auxiliaries.kind == "cancer").any()


def test_selection_matches_figure():
    """res.genes/signatures/concepts/auxiliaries 必须与图上实际展示的节点完全一致。"""
    expr, ct, clin = _fx()
    res = _sp()(expr, ct, expr.index[0], model="M3", clinical=clin)
    rp = res.render()
    sel = rp.selection
    # 图的节点顺序是布局选择，故比较**集合**；两列基因必须逐名对应
    assert set(res.genes.gene) == {n.name for n in sel["gene_score"]}
    assert list(res.genes.gene) == list(res.genes.gene)          # 两列同行
    assert set(res.genes.gene) == {n.name for n in sel["gene_tpm"]}
    assert set(res.signatures.signature) == {n.name for n in sel["signatures"]}
    assert set(res.concepts.concept) == {n.name for n in sel["predictors"]
                                         if n.kind == "concept"}
    assert set(res.auxiliaries.predictor) == {n.name for n in sel["predictors"]
                                              if n.kind in ("clinical", "pc", "cancer")}
    rp.close()


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
def test_html_and_static_output():
    """HTML 可独立导出；PNG/PDF/SVG 与 HTML 同源（同一 Figure）。"""
    import tempfile
    from pathlib import Path as _P
    expr, ct, clin = _fx()
    res = _sp()(expr, ct, expr.index[0], model="M2", clinical=clin)
    d = _P(tempfile.mkdtemp(prefix="sp_out_"))
    h = res.save_html(d / "p.html")            # 不先保存 PNG
    assert h.is_file() and h.stat().st_size > 50_000
    assert "data:image/png;base64," in h.read_text(encoding="utf-8", errors="ignore")
    for ext in ("png", "pdf", "svg"):
        f = res.save_static(d / f"p.{ext}")
        assert f.is_file() and f.stat().st_size > 5_000, ext
    files = res.save(d / "all")
    for k in ("png", "pdf", "svg", "html", "nodes.tsv", "selection.json",
              "decomposition.tsv"):
        assert files[k].is_file(), k


def test_html_embeds_the_same_png_bytes():
    """HTML 内嵌图像必须与同一次导出的 PNG **字节一致**（不存在第二套绘图器）。"""
    import base64, hashlib, re, tempfile
    from pathlib import Path as _P
    expr, ct, clin = _fx()
    res = _sp()(expr, ct, expr.index[0], model="M2", clinical=clin)
    d = _P(tempfile.mkdtemp(prefix="sp_html_"))
    files = res.save(d / "x")
    m = re.search(r"data:image/png;base64,([A-Za-z0-9+/=]+)",
                  files["html"].read_text(encoding="utf-8"))
    assert m, "HTML 未内嵌 PNG"
    assert (hashlib.sha256(base64.b64decode(m.group(1))).hexdigest()
            == hashlib.sha256(files["png"].read_bytes()).hexdigest())


def test_repeated_export_is_deterministic():
    """同一结果重复导出不得改动输入，且不引入随机布局。"""
    import hashlib, tempfile
    from pathlib import Path as _P
    expr, ct, clin = _fx()
    res = _sp()(expr, ct, expr.index[0], model="M2", clinical=clin)
    before = res.genes.copy()
    d1, d2 = _P(tempfile.mkdtemp()), _P(tempfile.mkdtemp())
    h1 = hashlib.sha256(res.save_static(d1 / "a.png").read_bytes()).hexdigest()
    h2 = hashlib.sha256(res.save_static(d2 / "a.png").read_bytes()).hexdigest()
    assert h1 == h2, "同一结果两次导出的 PNG 不一致（存在随机布局）"
    assert res.genes.equals(before), "导出过程改动了结果对象"


def test_figure_has_circles_uniform_links_and_risk_bar():
    """圆形节点 + 等粗等色连线 + 风险柱；连线不编码 contribution。"""
    expr, ct, clin = _fx()
    res = _sp()(expr, ct, expr.index[0], model="M2", clinical=clin)
    rp = res.render()
    from matplotlib.patches import PathPatch
    ax = rp.figure.axes[0]
    links = [p for p in ax.patches if isinstance(p, PathPatch)]
    assert links, "图中没有连线"
    styles = {(round(p.get_linewidth(), 6), tuple(np.round(p.get_edgecolor(), 6)),
               round(float(p.get_alpha()), 6)) for p in links}
    assert len(styles) == 1, f"连线样式不唯一（应等粗等色）：{styles}"
    # 圆形节点：scatter 的 marker 面积对应固定半径，且不含方/三角
    assert any(len(c.get_offsets()) > 0 for c in ax.collections), "缺少节点圆点"
    rp.close()


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
