# -*- coding: utf-8 -*-
"""QC 分级的合成覆盖度夹具测试（P3.5 §7）。

构造目标 global / Gsig 覆盖度的表达矩阵，逐一跑 ``compass_os.analyze``，
确认实际分级与**冻结的 qc_config** 完全一致，并覆盖组合规则
``overall = worse(global tier, Gsig tier)``。

同时检查**真实队列典型值**（global ≈ 0.967、Gsig ≈ 0.962）不得被无合理依据地标 warning。
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from _util import fixtures_available, load_fixtures

pytestmark = pytest.mark.skipif(not fixtures_available(), reason="需要 tests/fixtures")


def _tier(value: float, rec, warn) -> str:
    """与 compass_os.qc._tier_for 同语义的独立实现（测试端不复用被测代码）。"""
    if rec is not None and value >= rec:
        return "recommended"
    if warn is not None and value >= warn:
        return "warning"
    return "below_warning"


def _fixture_for(global_cov: float, gsig_cov: float) -> pd.DataFrame:
    """按目标覆盖度裁剪：先保 Gsig，再从 non-Gsig 里删到目标 global 覆盖度。"""
    from compass_os.preprocessing import feature_names
    from compass_os.qc import gene_sets

    fx = load_fixtures()
    expr = fx["expression"]
    genes = list(feature_names())
    sig = sorted({g for gl in gene_sets()["gene_list"] for g in gl})
    non_sig = [g for g in genes if g not in set(sig)]

    n_sig_drop = int(round((1.0 - gsig_cov) * len(sig)))
    dropped_sig = set(sig[:n_sig_drop])
    # 目标 global 覆盖度 = 保留的全部基因 / 15672
    keep_total = int(round(global_cov * len(genes)))
    n_non_drop = len(genes) - len(dropped_sig) - keep_total
    dropped_non = set(non_sig[:max(0, n_non_drop)])
    drop = sorted(dropped_sig | dropped_non)
    return expr.drop(columns=[g for g in drop if g in expr.columns])


def _analyze_cov(expr):
    import compass_os
    fx = load_fixtures()
    return compass_os.analyze(expr, fx["labels"]["cancer_type"].tolist(),
                              clinical=fx["clinical"], robustness=False)


def test_frozen_thresholds_are_on_clean_measured_levels():
    """冻结阈值必须落在实测/整数级覆盖度上，不出现虚假精度。"""
    from compass_os.qc import qc_config
    cfg = qc_config()
    assert cfg is not None and cfg.get("status") == "complete"
    for k in ("recommended_global_coverage", "warning_global_coverage",
              "recommended_Gsig_coverage", "warning_Gsig_coverage"):
        v = cfg.get(k)
        assert v is not None, f"{k} 缺失"
        assert abs(round(v, 2) - v) < 1e-9, f"{k}={v} 不是两位以内的实测级覆盖度"
    assert cfg["recommended_global_coverage"] > cfg["warning_global_coverage"]
    assert cfg["recommended_Gsig_coverage"] > cfg["warning_Gsig_coverage"]


@pytest.mark.parametrize("global_cov,gsig_cov", [
    (1.00, 1.00),
    (0.95, 0.995),
    (0.90, 0.95),
    (0.80, 0.90),
    (0.70, 0.85),
    (0.60, 0.80),
])
def test_analyze_tier_matches_frozen_config(global_cov, gsig_cov):
    """合成覆盖度下 analyze 的分级必须与 qc_config 的阈值逐一吻合。"""
    from compass_os.qc import qc_config
    cfg = qc_config()
    res = _analyze_cov(_fixture_for(global_cov, gsig_cov))
    got = res.qc.tiers
    g_exp = _tier(res.qc.gene_coverage, cfg["recommended_global_coverage"],
                  cfg["warning_global_coverage"])
    gsig = res.qc.gsig_coverage            # 与阈值同口径（916 基因中被观测到的比例）
    assert abs(gsig - gsig_cov) < 0.01, f"夹具未命中目标 Gsig 覆盖度：{gsig:.3f}"
    s_exp = _tier(gsig, cfg["recommended_Gsig_coverage"], cfg["warning_Gsig_coverage"])
    assert got["global"] == g_exp, f"global tier 不一致：{got['global']} vs {g_exp}"
    assert got["Gsig"] == s_exp, f"Gsig tier 不一致：{got['Gsig']} vs {s_exp}"
    worst = max([g_exp, s_exp], key=lambda t: ("recommended", "warning",
                                               "below_warning").index(t))
    assert got["overall"] == worst, f"overall 组合规则不一致：{got['overall']} vs {worst}"
    assert got["rule"] == "worse(global, Gsig)"


def test_analyze_summary_reports_overall_tier():
    """摘要中的 QC 一行必须反映 overall 分级（而不是只报 global）。"""
    res = _analyze_cov(_fixture_for(0.80, 0.90))
    txt = res.summary()
    assert "caution" in txt.lower() or "Recommended range" in txt
    assert "signature" in txt.lower()


def test_real_cohort_typical_coverage_is_recommended():
    """真实队列典型覆盖度（global≈0.967、Gsig≈0.962）必须判为 recommended。

    这是 P3.5 的核心回归点：阈值若报错轴（例如把 Gsig 阈值写成 0.994），
    绝大多数真实队列会被无合理依据地标 warning。
    """
    from compass_os.qc import qc_config
    cfg = qc_config()
    assert cfg["recommended_Gsig_coverage"] <= 0.962, (
        f"recommended_Gsig_coverage={cfg['recommended_Gsig_coverage']} 高于真实队列典型值 "
        "0.962 ⇒ 绝大多数真实队列会被误标 warning（阈值轴可能报错）")
    assert cfg["recommended_global_coverage"] <= 0.967
    res = _analyze_cov(_fixture_for(0.967, 0.962))
    assert res.qc.tiers["overall"] == "recommended", res.qc.tiers
    assert "Recommended range" in res.summary()


def test_below_warning_language_is_plain():
    """低于 warning 档时必须用人能读懂的语言，而不是只给 flag。"""
    res = _analyze_cov(_fixture_for(0.60, 0.80))
    txt = res.summary()
    assert res.qc.tiers["overall"] == "below_warning"
    assert "avoid quantitative interpretation" in txt
