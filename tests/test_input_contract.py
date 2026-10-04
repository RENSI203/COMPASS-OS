# -*- coding: utf-8 -*-
"""输入契约回归测试（第 1 阶段数据限制审查发现的三处缺陷）。

覆盖：

1. ``±Inf`` 必须被**显式拒绝**为 :class:`InputError`，而不是静默流入 scaler
   后在 sklearn 里抛出与输入无关的报错；
2. ``cancer_type`` 长度不符必须是 :class:`InputError`，而不是 pandas 的裸 ``ValueError``；
3. 负 COMPASS 码值（``NORMAL`` = 正常组织）必须给出明确错误，
   而不是 embedding 层的裸 ``IndexError``。

同时锁定**不应改变**的既有语义：``NaN`` 仍是受支持的"缺失"标记。
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from _util import fixtures_available, load_fixtures
from compass_os.exceptions import InputError, UnknownCancerTypeError
from compass_os.preprocessing import align_expression, cancer_codes_for, feature_names

pytestmark = pytest.mark.skipif(not fixtures_available(), reason="需要 tests/fixtures")

G = list(feature_names())


def _small(n_genes=300, n_samples=4):
    rng = np.random.default_rng(0)
    return pd.DataFrame(rng.lognormal(2, 1.5, size=(n_samples, n_genes)),
                        columns=G[:n_genes])


# ---------------------------------------------------------------- 1) 非有限值
@pytest.mark.parametrize("bad", [np.inf, -np.inf])
def test_infinite_values_are_rejected(bad):
    """±Inf 必须是 InputError（明确指向输入），不得静默通过。"""
    e = _small()
    e.iloc[0, 0] = bad
    with pytest.raises(InputError) as exc:
        align_expression(e, "reference")
    msg = str(exc.value)
    assert "非有限值" in msg and ("+Inf" in msg or "-Inf" in msg)
    assert e.columns[0] in msg, "错误信息应指出出问题的基因"


def test_infinite_values_rejected_in_predict():
    """通过公共 predict 入口同样给出 InputError（而非 sklearn 的 ValueError）。"""
    import compass_os
    e = _small()
    e.iloc[0, 0] = np.inf
    with pytest.raises(InputError):
        compass_os.predict(e, ["LUAD"] * len(e), model="M2")


def test_nan_is_still_treated_as_missing():
    """NaN 是**受支持**的缺失标记：不得因为拒绝 Inf 而把 NaN 一并拒绝。"""
    e = _small()
    e.iloc[0, 0] = np.nan
    a = align_expression(e, "reference")
    # 列存在 ⇒ 队列层覆盖度仍按列计；逐样本统计另给
    assert 0.0 < a.qc.gene_coverage <= 1.0
    assert a.qc.per_sample is not None
    # NaN 位置被当作缺失（observed_mask 为 False）并按策略填补
    assert bool(a.observed_mask.iloc[0, 0]) is False
    assert np.isfinite(a.matrix.iloc[0, 0])


def test_negative_values_still_rejected():
    """负值检查不能被 Inf 检查取代。"""
    e = _small()
    e.iloc[0, 0] = -1.0
    with pytest.raises(InputError):
        align_expression(e, "reference")


# ---------------------------------------------------------------- 2) 长度
@pytest.mark.parametrize("n", [1, 2, 5])
def test_cancer_type_length_mismatch_is_input_error(n):
    """长度不符必须是 InputError，而不是 pandas 的裸 ValueError。"""
    import compass_os
    e = _small()
    with pytest.raises(InputError) as exc:
        compass_os.predict(e, ["LUAD"] * n, model="M2")
    assert "cancer_type 长度" in str(exc.value)


# ---------------------------------------------------------------- 3) NORMAL
def test_normal_is_rejected_with_clear_message():
    """NORMAL 是正常组织标记（compass_code = -1），必须明确拒绝。"""
    with pytest.raises(UnknownCancerTypeError) as exc:
        cancer_codes_for(["NORMAL", "LUAD"])
    assert "正常组织" in str(exc.value)


def test_normal_rejected_in_predict():
    """不能以 embedding 的裸 IndexError 形式暴露给用户。"""
    import compass_os
    e = _small()
    with pytest.raises(UnknownCancerTypeError):
        compass_os.predict(e, ["NORMAL"] * len(e), model="M2")


def test_all_table_tumour_types_still_resolve():
    """码表中**非负**码值的肿瘤类型必须全部仍可用（修复不得误伤）。"""
    from compass_os.qc import cancer_codes
    table = cancer_codes()
    tumour = [t for t, c in zip(table["cancer_type"], table["compass_code"]) if int(c) >= 0]
    assert len(tumour) >= 30
    codes = cancer_codes_for(tumour)
    assert (codes.to_numpy() >= 0).all()
    # NORMAL 是唯一的负码值
    assert [t for t, c in zip(table["cancer_type"], table["compass_code"])
            if int(c) < 0] == ["NORMAL"]


# ---------------------------------------------------------------- 输入尺度
def test_only_declared_scales_are_accepted():
    """只有 tpm / log2_tpm1 两个尺度；其他声明被拒绝。"""
    e = _small()
    for bad in ("counts", "microarray", "zscore", "TPM", ""):
        with pytest.raises(InputError):
            align_expression(e, "reference", input_scale=bad)


def test_counts_like_matrix_is_not_detectable():
    """**已知限制**：非负 counts/强度矩阵无法与 TPM 区分，会被静默当作 TPM。

    本测试锁定该事实（而非"修复"它）——因为库无法只凭数值判定尺度，
    尺度由调用方用 ``input_scale`` 显式声明。一旦未来加入启发式检测，
    本测试应当失败并需要相应更新文档。
    """
    rng = np.random.default_rng(1)
    counts = pd.DataFrame(rng.integers(0, 5000, size=(4, 300)).astype(float), columns=G[:300])
    a = align_expression(counts, "reference")           # 不报错
    assert a.matrix.shape == (4, len(G))
    assert np.isfinite(a.matrix.to_numpy()).all()
