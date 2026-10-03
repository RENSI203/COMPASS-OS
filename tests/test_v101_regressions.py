# -*- coding: utf-8 -*-
"""v1.0.1 regression tests —— 针对本轮 7 个 bug 的**可失败**回归测试。

每个测试都对应一个具体缺陷；在 v1.0.0 代码上应当 FAIL，在 v1.0.1 上 PASS：

* ``test_analyze_log2_scale_matches_tpm``            —— Bug 1：``analyze`` 未透传 ``input_scale``
* ``test_auto_robustness_uses_overall_qc``           —— Bug 2：auto 只看 global coverage
* ``test_auto_robustness_triggers_on_low_Gsig_high_global`` —— Bug 2：双轴 Case B/C
* ``test_summary_qc_and_robustness_are_consistent``  —— Bug 2b：``_coverage_ok`` 与 overall tier 不一致
* ``test_stratification_respects_min_cohort``        —— Bug 3：小队列仍给 risk_group
* ``test_robustness_documentation_is_current``       —— Bug 7：遗留 ``not_calibrated`` 文案
* ``test_quick_start_output_matches_readme``         —— Bug 5：README 展示结果与示例不一致
"""
from __future__ import annotations

import os
import re
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from _util import fixtures_available, load_fixtures, repo_root

pytestmark = pytest.mark.skipif(not fixtures_available(), reason="需要 tests/fixtures")

README = repo_root() / "README.md"


# ─────────────────────────── Bug 1：input_scale 透传 ───────────────────────────
def test_analyze_log2_scale_matches_tpm():
    """同一表达量的 TPM 与 log2(TPM+1) 两种写法，analyze() 结果必须一致。

    v1.0.0 缺陷：``analyze`` 未把 ``input_scale`` 传给 ``_api.predict``，
    于是 ``input_scale="log2_tpm1"`` 的数据被当作 TPM 处理 ⇒ 结果与 TPM 路径不一致。
    """
    import compass_os

    fx = load_fixtures()
    tpm = fx["expression"]
    log2_tpm1 = np.log2(tpm + 1.0)
    ct = fx["labels"]["cancer_type"].tolist()

    a = compass_os.analyze(tpm, ct, clinical=fx["clinical"], input_scale="tpm",
                           robustness=False)
    b = compass_os.analyze(log2_tpm1, ct, clinical=fx["clinical"], input_scale="log2_tpm1",
                           robustness=False)

    assert np.allclose(a.risk.to_numpy(float), b.risk.to_numpy(float), atol=1e-8), \
        "TPM 与 log2(TPM+1) 两条路径的 risk 不一致 ⇒ input_scale 未正确透传"
    assert np.allclose(a.concept_scores.to_numpy(float), b.concept_scores.to_numpy(float),
                       atol=1e-8)
    assert np.allclose(a.signature_scores.to_numpy(float), b.signature_scores.to_numpy(float),
                       atol=1e-8)
    assert a.qc.gene_coverage == b.qc.gene_coverage
    assert a.qc.tiers == b.qc.tiers


def test_advanced_paths_also_handle_log2_scale():
    """Level 2 的 predict / get_representation 在 log2 尺度下与 TPM 路径一致。"""
    import compass_os

    fx = load_fixtures()
    tpm, ct = fx["expression"], fx["labels"]["cancer_type"].tolist()
    log2_tpm1 = np.log2(tpm + 1.0)

    p1 = compass_os.predict(tpm, ct, clinical=fx["clinical"], model="M2", input_scale="tpm")
    p2 = compass_os.predict(log2_tpm1, ct, clinical=fx["clinical"], model="M2",
                            input_scale="log2_tpm1")
    assert np.allclose(p1.risk.to_numpy(float), p2.risk.to_numpy(float), atol=1e-8)

    r1 = compass_os.get_representation(tpm, ct, input_scale="tpm")
    r2 = compass_os.get_representation(log2_tpm1, ct, input_scale="log2_tpm1")
    assert np.allclose(r1.concept_scores.to_numpy(float), r2.concept_scores.to_numpy(float),
                       atol=1e-8)


def test_analyze_rejects_unknown_input_scale():
    """尺度必须显式且合法；不得静默猜测。"""
    import compass_os
    from compass_os.exceptions import InputError

    fx = load_fixtures()
    with pytest.raises((InputError, ValueError)):
        compass_os.analyze(fx["expression"], fx["labels"]["cancer_type"].tolist(),
                           input_scale="microarray", robustness=False)


# ─────────────────── Bug 2 / 2b：auto robustness 服从双轴 overall QC ───────────────────
def _fixture_for(global_cov: float, gsig_cov: float) -> pd.DataFrame:
    """构造目标 global / Gsig 覆盖度的表达矩阵（与 test_qc_tiers 同法）。"""
    from compass_os.preprocessing import feature_names
    from compass_os.qc import gene_sets

    expr = load_fixtures()["expression"]
    genes = list(feature_names())
    sig = sorted({g for gl in gene_sets()["gene_list"] for g in gl})
    non_sig = [g for g in genes if g not in set(sig)]
    n_sig_drop = int(round((1.0 - gsig_cov) * len(sig)))
    dropped_sig = set(sig[:n_sig_drop])
    keep_total = int(round(global_cov * len(genes)))
    n_non_drop = len(genes) - len(dropped_sig) - keep_total
    dropped_non = set(non_sig[:max(0, n_non_drop)])
    drop = sorted(dropped_sig | dropped_non)
    return expr.drop(columns=[g for g in drop if g in expr.columns])


def _analyze(expr, robustness="auto"):
    import compass_os
    fx = load_fixtures()
    return compass_os.analyze(expr, fx["labels"]["cancer_type"].tolist(),
                              clinical=fx["clinical"], robustness=robustness)


@pytest.mark.parametrize("g_cov,s_cov,expect_ran", [
    (1.00, 1.000, False),   # Case A：双轴 recommended → 不跑
    (0.95, 0.995, False),   # Case A'
    (0.95, 0.850, True),    # Case B：global 好、Gsig 差 → 必须跑（v1.0.0 会漏）
    (0.80, 0.950, True),    # Case C：global 差、Gsig 好 → 必须跑
    (0.60, 0.800, True),    # Case D：below_warning → 必须跑
])
def test_auto_robustness_uses_overall_qc(g_cov, s_cov, expect_ran):
    """``robustness="auto"`` 必须由 **overall = worse(global, Gsig)** 决定。"""
    res = _analyze(_fixture_for(g_cov, s_cov))
    tier = res.qc.tiers["overall"]
    assert (tier == "recommended") is (not expect_ran), \
        f"夹具未命中预期档位：global={g_cov} Gsig={s_cov} → overall={tier}"
    ran = res.robustness is not None
    assert ran is expect_ran, (
        f"auto robustness 行为与 overall QC 不一致：overall={tier} 但 "
        f"{'运行了' if ran else '未运行'}稳健性比较")


def test_auto_robustness_triggers_on_low_Gsig_high_global():
    """**本轮最重要的回归点**：global 达标但 Gsig 不达标时，auto 必须触发。

    v1.0.0 缺陷：auto 只比较 ``gene_coverage`` 与 recommended_global_coverage，
    因此 global=0.95 / Gsig=0.85（overall=warning）会被错误地判定为"无需稳健性检查"。
    """
    res = _analyze(_fixture_for(0.95, 0.85))
    assert res.qc.tiers["global"] == "recommended"
    assert res.qc.tiers["Gsig"] != "recommended"
    assert res.qc.tiers["overall"] != "recommended"
    assert res.robustness is not None, "Gsig 未达标时 auto robustness 必须运行"


def test_summary_qc_and_robustness_are_consistent():
    """摘要不得出现 ``QC = 需谨慎`` 却 ``Robustness check = Not required`` 的自相矛盾。"""
    res = _analyze(_fixture_for(0.95, 0.85))
    txt = res.summary()
    qc_line = next(l for l in txt.splitlines() if l.startswith("QC"))
    rob_line = next(l for l in txt.splitlines() if l.startswith("Robustness check"))
    assert "Recommended range" not in qc_line, f"QC 行不应为 recommended：{qc_line}"
    assert "Not required" not in rob_line, (
        f"QC 非 recommended 时不得声称无需稳健性检查：\n  {qc_line}\n  {rob_line}")


def test_robustness_false_never_runs_even_when_flagged():
    """显式 ``robustness=False`` 时即使 QC 触发也不运行。"""
    res = _analyze(_fixture_for(0.60, 0.80), robustness=False)
    assert res.robustness is None


def test_robustness_true_always_runs():
    """显式 ``robustness=True`` 时即使 QC 达标也运行。"""
    res = _analyze(_fixture_for(1.00, 1.00), robustness=True)
    assert res.robustness is not None


# ─────────────────────── Bug 3：min_cohort_for_stratification ───────────────────────
def _n_samples(n: int) -> tuple:
    """由夹具复制出 **n 个索引互不相同** 的样本（避免重复 index 干扰对齐）。"""
    fx = load_fixtures()
    base, clin0 = fx["expression"], fx["clinical"]
    b, m = base.to_numpy(float), clin0.to_numpy(float)
    idx = [f"S{k * len(base) + j}" for k in range(int(np.ceil(n / len(base)))) for j in
           range(len(base))][:n]
    expr = pd.DataFrame(np.tile(b, (int(np.ceil(n / len(b))), 1))[:n], index=idx,
                        columns=base.columns)
    clin = pd.DataFrame(np.tile(m, (int(np.ceil(n / len(m))), 1))[:n], index=idx,
                        columns=clin0.columns)
    ct = (fx["labels"]["cancer_type"].tolist() * int(np.ceil(n / len(base))))[:n]
    return expr, ct, clin


@pytest.mark.parametrize("n,expect_rank,expect_group", [
    (1, False, False),    # 单样本：两者皆无
    (5, True, False),     # 小队列：只有秩
    (29, True, False),    # 阈值下方：只有秩
    (30, True, True),     # 阈值：两者皆有
    (50, True, True),     # 阈值上方
])
def test_stratification_respects_min_cohort(n, expect_rank, expect_group):
    """分层与排序必须按文档分离：n < min_cohort 时不提供 risk_group。"""
    import compass_os

    expr, ct, clin = _n_samples(n)
    res = compass_os.predict(expr, ct, clinical=clin, model="M2")
    assert (res.risk_rank is not None) is expect_rank
    assert (res.risk_group is not None) is expect_group
    if not expect_group and n >= 2:
        assert any("risk_group is not assigned" in w for w in res.warnings)


def test_min_cohort_for_stratification_is_configurable():
    """Advanced 参数确实生效：min_cohort=10 时 n=15 应该有 risk_group。"""
    import compass_os

    expr, ct, clin = _n_samples(15)
    off = compass_os.predict(expr, ct, clinical=clin, model="M2",
                             min_cohort_for_stratification=30)
    on = compass_os.predict(expr, ct, clinical=clin, model="M2",
                            min_cohort_for_stratification=10)
    assert off.risk_group is None
    assert on.risk_group is not None and on.risk_group.shape[0] == 15


def test_level1_relative_split_is_separate_from_frozen_group():
    """Level 1 的 display-only 队列相对切分**不得**混同为冻结切点分层字段。"""
    import compass_os

    expr, ct, clin = _n_samples(25)
    res = compass_os.analyze(expr, ct, clinical=clin, robustness=False)
    # 25 < 30 ⇒ 冻结切点分组不可用；若给出相对切分，必须是独立字段且标注 display-only
    assert res.risk_group is None
    if res.risk_group_relative is not None:
        assert res.risk_group_relative.name == "risk_group_relative"
        assert any("display only" in note for note in res.notes)


# ─────────────────────────── Bug 7：文案已更新 ───────────────────────────
def test_robustness_documentation_is_current():
    """不得遗留 ``not_calibrated`` / "stress test 完成后校准" 这类过期表述。"""
    import compass_os

    fx = load_fixtures()
    r = compass_os.check_robustness(fx["expression"].iloc[:6],
                                    fx["labels"]["cancer_type"].iloc[:6],
                                    fx["clinical"].iloc[:6], model="M2")
    assert r.robustness_flag != "not_calibrated"
    assert "not_calibrated" not in (r.note or "")
    blob = (r.note or "") + (r.robustness_flag or "")
    assert "stress test 完成后" not in blob and "尚未校准" not in blob


def test_no_stale_calibration_placeholders_in_repo():
    """全仓（源码 + README + docs）不得残留过期的"待校准"占位表述。"""
    # "thresholds not frozen" 仅允许出现在**明确限定为本安装**的 fallback 里（异常安装场景）
    pat = re.compile(r"not_calibrated|calibration pending|stress test pending|"
                     r"signature-targeted calibration pending|"
                     r"thresholds not frozen(?! in this installation)")
    roots = [repo_root() / "src", repo_root() / "README.md", repo_root() / "docs"]
    bad = []
    for r in roots:
        files = [r] if r.is_file() else [p for p in r.rglob("*") if p.suffix in (".py", ".md")]
        for f in files:
            if "__pycache__" in str(f):
                continue
            for i, line in enumerate(f.read_text(encoding="utf-8").splitlines(), 1):
                if pat.search(line):
                    bad.append(f"{f.relative_to(repo_root())}:{i}: {line.strip()[:90]}")
    assert not bad, "发现过期校准占位表述：\n  " + "\n  ".join(bad)


# ─────────────────────── Bug 5：README 与示例输出一致 ───────────────────────
def test_quick_start_output_matches_readme():
    """README Quick Start 展示的结果必须与 ``examples/quick_start.py`` 真实输出一致。

    v1.0.0 缺陷：代码用 ``cancer_type="LUAD"``（全队列单一癌种），
    展示结果却写 "9 cancer types" —— 不是同一段代码的输出。
    """
    import compass_os

    data = repo_root() / "examples" / "example_data"
    expr = pd.read_csv(data / "example_expression.tsv.gz", sep="\t", index_col=0)
    clin = pd.read_csv(data / "example_clinical.tsv", sep="\t", index_col=0)
    lab = pd.read_csv(data / "example_labels.tsv", sep="\t", index_col=0)
    surv = lab.rename(columns={"os_time_days": "os_time_days", "os_event": "os_event"})

    res = compass_os.analyze(expr, "LUAD", clinical=clin, survival=surv)
    real = res.summary()

    txt = README.read_text(encoding="utf-8")
    block = re.search(r"## 2\. Quick start.*?```text\n(.*?)```", txt, re.S)
    assert block, "README 中找不到 Quick start 的结果代码块"
    shown = block.group(1)

    def grab(text, key):
        for line in text.splitlines():
            if line.startswith(key):
                return " ".join(line[len(key):].split())
        return None

    for key in ("Samples", "Cancer type", "Risk estimates", "Risk stratification",
                "Missing-gene strategy"):
        a, b = grab(shown, key), grab(real, key)
        assert a is not None, f"README 展示块缺少 {key!r} 行"
        assert a == b, f"{key} 不一致：README={a!r} 实际={b!r}"
    assert "cancer types (" not in shown, (
        "Quick Start 用单一 cancer_type，展示结果不得出现多癌种计数")
    assert "9 cancer types" not in shown

# ─────────── Bug 8：上游 loadcompass 会删除位于临时目录下的 checkpoint ───────────
def test_checkpoint_survives_loading_from_temp_directory(tmp_path=None):
    """模型文件若位于系统临时目录，加载后**不得被删除**。

    上游 ``compass.loadcompass`` 末尾有 ``if file.startswith(tempfile.gettempdir()):
    os.remove(file)``（本意是清理 URL 下载的临时文件），会把任何位于 /tmp 下的
    checkpoint 静默删除 ⇒ 装在 /tmp 下的 Python 环境（容器 / CI / pip --target /tmp）
    第一次预测后即失效。v1.0.1 改为直接 ``torch.load`` 本地冻结资产。
    """
    import shutil
    import tempfile

    import compass_os
    from compass_os import representation as R
    from compass_os._paths import asset

    tmp = Path(tempfile.mkdtemp(prefix="compass_ck_"))
    try:
        (tmp / "models").mkdir()
        for name in ("locked_M0.json", "locked_M1.json", "locked_M2.json", "locked_M3.json",
                     "locked_pca_M3.npz", "reference_quantiles.json",
                     "model_manifest.json", "qc_config.json"):
            shutil.copy2(asset(f"models/{name}"), tmp / "models" / name)
        ck = tmp / "models" / "pretrainer.pt"
        shutil.copy2(asset("models/pretrainer.pt"), ck)
        assert ck.is_file()

        saved = os.environ.get("COMPASS_OS_ROOT")
        os.environ["COMPASS_OS_ROOT"] = str(tmp)          # 让资产解析指向临时目录
        try:
            from compass_os import _paths
            _paths.assets_root.cache_clear()
            R.load_model.cache_clear()
            model = R.load_model()
            assert type(model).__name__ == "PreTrainer"
        finally:
            if saved is None:
                os.environ.pop("COMPASS_OS_ROOT", None)
            else:
                os.environ["COMPASS_OS_ROOT"] = saved
            from compass_os import _paths
            _paths.assets_root.cache_clear()
            R.load_model.cache_clear()

        assert ck.is_file(), (
            "位于临时目录的 checkpoint 在加载后被删除 —— 上游 loadcompass 的 "
            "os.remove 副作用未被规避")
        assert ck.stat().st_size == asset("models/pretrainer.pt").stat().st_size
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_model_loading_does_not_use_destructive_upstream_loader():
    """``load_model`` 不得调用上游 ``loadcompass``（其含 os.remove 副作用）。"""
    import inspect

    from compass_os import representation as R
    src = inspect.getsource(R.load_model)
    assert "loadcompass(" not in src, "load_model 仍调用上游 loadcompass（会删 /tmp 下的权重）"
    assert "torch.load" in src
