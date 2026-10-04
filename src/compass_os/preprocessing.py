# -*- coding: utf-8 -*-
"""表达预处理：基因对齐 + 缺失基因策略 + QC。

输入契约
--------
* ``expression``：``samples × genes`` 的 DataFrame，**列名为基因符号**，索引为样本 ID；
* 尺度由 ``input_scale`` 显式声明（默认 ``"tpm"``）：
  ``"tpm"`` = 线性 TPM；``"log2_tpm1"`` = ``log2(TPM+1)``（内部转回 TPM）。
  模型内部（冻结 Datascaler）自行完成 ``log2(x+1)`` 与 MinMax 标准化，
  因此这里**不得**再做任何变换；
* ``cancer_type`` 单独提供，缺失即报错（**不得**静默推断）。

缺失基因策略（三种，语义严格区分）
----------------------------------
``reference``（默认）
    用**训练期冻结的 TCGA reference median** 填充：``log2(TPM+1) = 3.130937933922``
    ⟺ ``TPM = 2^3.130937933922 − 1 ≈ 7.76``。数值取自
    ``models/reference_quantiles.json`` 的 ``median(values)``。
    **绝不**用用户队列自身的 mean/median 估计——那会让同一样本的预测依赖同批其他样本。

``zero``
    训练空间掩码（training-space masking）：先正常经过冻结 Datascaler，
    再把缺失基因对应的**标准化值**置 0。这与 COMPASS 对比预训练中
    ``RandomMaskAugmentor`` 的 ``x[1:][mask] = 0`` 处在**同一归一化输入空间**。
    ⚠ 措辞边界：训练期掩码只作用于正样本视图、概率 0.41 且随机选择；推理期是
    **确定性地对真实缺失基因**做同一操作 ⇒ 只能说「同一归一化输入空间」，
    **不得**声称「与训练增强分布一致」。详见 docs/MISSING_GENES.md。
    **不得**实现成 ``raw TPM = 0``：实测 60.2% 的基因其训练 ``data_min_ > 0``，
    置 0 会落到训练分布之外（负标准化值）。

``strict``
    任一必需基因缺失即抛 :class:`MissingGenesError`（用于 benchmark / 严格复现）。

本模块只做**纯输入变换**；``zero`` 的标准化空间掩码由
:mod:`compass_os.representation` 在调用冻结 scaler 之后施加（见 ``_ZeroMaskedScaler``）。
"""
from __future__ import annotations

import functools
import json
from dataclasses import dataclass
from typing import Iterable, Sequence

import numpy as np
import pandas as pd

from ._paths import asset
from .exceptions import InputError, MissingGenesError, UnknownCancerTypeError
from .qc import CoverageQC, cancer_codes, gene_coverage

STRATEGIES = ("reference", "zero", "strict")
INPUT_SCALES = ("tpm", "log2_tpm1")


@dataclass
class AlignedExpression:
    """对齐后的表达与缺失信息。

    ``matrix``：samples × 15,672，TPM 尺度，列序 = 词表顺序（模型按列名取数，顺序不影响数值）；
    ``log2``：同形状的 ``log2(TPM+1)``，供 M3 的**锁定 PCA** 使用；
    ``observed_mask`` / ``missing_mask``：True 分别为「真实观测」/「被填补」。
    """

    matrix: pd.DataFrame
    log2: pd.DataFrame
    observed_mask: pd.DataFrame
    missing_mask: pd.DataFrame
    qc: CoverageQC


@functools.lru_cache(maxsize=1)
def feature_names() -> tuple:
    """模型词表（15,672 基因，顺序 = 冻结 checkpoint 的 ``feature_name``）。

    取自派生资产 ``src/compass_os/data/gene_vocabulary.txt``（由
    ``tools/build_assets.py`` 从 checkpoint 导出；``tests/test_assets.py`` 复核一致），
    因此 QC 路径无需加载 torch。
    """
    path = asset("src/compass_os/data/gene_vocabulary.txt")
    return tuple(path.read_text(encoding="utf-8").split())


@functools.lru_cache(maxsize=1)
def reference_fill_log2() -> float:
    """冻结的参考填充值（log2 尺度）= ``median(reference_quantiles['values'])``。"""
    ref = json.loads(asset("models/reference_quantiles.json").read_text(encoding="utf-8"))
    return float(np.median(np.asarray(ref["values"], dtype=float)))


def reference_fill_tpm() -> float:
    """冻结的参考填充值（TPM 尺度）= ``2^median − 1``。"""
    return float(2.0 ** reference_fill_log2() - 1.0)


def _to_tpm(expression: pd.DataFrame, input_scale: str) -> pd.DataFrame:
    if input_scale not in INPUT_SCALES:
        raise InputError(f"input_scale 必须是 {INPUT_SCALES} 之一，收到 {input_scale!r}")
    if not isinstance(expression, pd.DataFrame):
        raise InputError("expression 必须是 pandas DataFrame（samples × genes）")
    if expression.shape[0] == 0 or expression.shape[1] == 0:
        raise InputError(f"expression 为空：shape={expression.shape}")
    if expression.columns.duplicated().any():
        dup = expression.columns[expression.columns.duplicated()].unique()[:5].tolist()
        raise InputError(f"expression 存在重复基因列名（前几个：{dup}）；请先合并重复列")
    if expression.index.duplicated().any():
        dup = expression.index[expression.index.duplicated()].unique()[:5].tolist()
        raise InputError(f"expression 存在重复样本 ID（前几个：{dup}）")
    vals = expression.to_numpy(dtype=np.float64)
    # ±Inf 必须显式拒绝：np.isfinite 过滤会把它们排除在负值检查之外，
    # 从而静默流入 scaler（下游报出与输入无关的 sklearn 报错）。
    # NaN 不同——它是**受支持的"缺失"标记**，由 missing_gene_strategy 处理，故不在此拒绝。
    inf_mask = np.isinf(vals)
    if inf_mask.any():
        rows, cols = np.where(inf_mask)
        rc = sorted({str(expression.index[i]) for i in rows})[:3]
        cc = sorted({str(expression.columns[j]) for j in cols})[:3]
        kinds = sorted({"+Inf" if v > 0 else "-Inf" for v in vals[inf_mask]})
        raise InputError(
            f"表达矩阵存在非有限值 {kinds}（共 {int(inf_mask.sum())} 处；样本如 {rc}，"
            f"基因如 {cc}）。Inf/-Inf 不是有效表达值，请先修正为有限值；"
            "若某基因确实缺失请用 NaN（NaN 会按 missing_gene_strategy 处理）。")
    finite = vals
    if finite.size and finite.min() < -1e-9:
        raise InputError(
            f"表达矩阵存在负值（min={finite.min():.4g}）。TPM 不能为负；"
            "若输入是 log2(TPM+1) 请设 input_scale='log2_tpm1'。")
    if input_scale == "log2_tpm1":
        with np.errstate(over="ignore"):
            return pd.DataFrame(np.power(2.0, vals) - 1.0,
                                index=expression.index, columns=expression.columns)
    return expression.astype(np.float64)


def align_expression(expression: pd.DataFrame, missing_gene_strategy: str = "reference",
                     *, input_scale: str = "tpm",
                     genes: Sequence[str] | None = None) -> AlignedExpression:
    """把用户表达矩阵对齐到模型词表，并按策略填补缺失基因（含 QC）。

    填补只作用于**缺失基因**；输入中已有的 NaN 视同缺失，按同一策略处理。
    """
    if missing_gene_strategy not in STRATEGIES:
        raise InputError(f"missing_gene_strategy 必须是 {STRATEGIES} 之一，"
                         f"收到 {missing_gene_strategy!r}")
    genes = [str(g) for g in (genes if genes is not None else feature_names())]
    tpm = _to_tpm(expression, input_scale)
    tpm.columns = [str(c) for c in tpm.columns]

    qc = gene_coverage(tpm, genes, missing_gene_strategy)
    # 尺度自检（**提示，不阻断**）：与 analysis.py 既有启发式同一口径（max < 30 ⇒ 疑似 log2）。
    # 放在共享路径上，使 predict() 与 analyze() 行为一致（此前只有 analyze() 有该提示）。
    raw = expression.to_numpy(dtype=np.float64)
    if raw.size:
        mx = float(np.nanmax(raw))
        if np.isfinite(mx):
            if input_scale == "tpm" and mx < 30:
                qc.warnings.append(
                    f"input_scale='tpm' but the maximum input value is {mx:.3g} (< 30). If this "
                    "matrix is log2(TPM+1), re-run with input_scale='log2_tpm1'.")
            elif input_scale == "log2_tpm1" and mx > 30:
                qc.warnings.append(
                    f"input_scale='log2_tpm1' but the maximum input value is {mx:.3g} (> 30). If "
                    "this matrix is linear TPM, re-run with input_scale='tpm'.")
    if missing_gene_strategy == "strict" and qc.n_missing_genes:
        raise MissingGenesError(qc.missing_gene_names, qc.n_required_genes)

    aligned = tpm.reindex(columns=genes)
    observed = aligned.notna()
    missing_mask = ~observed

    if missing_gene_strategy == "reference":
        aligned = aligned.fillna(reference_fill_tpm())
    else:
        # zero：此处先填 NaN（真实缺失基因），标准化空间置 0 由 representation 层完成；
        # 填 0 只是占位，因为随后会被 scaler 结果覆盖为 0。
        aligned = aligned.fillna(0.0)

    log2 = pd.DataFrame(np.log2(aligned.to_numpy(dtype=np.float64) + 1.0),
                        index=aligned.index, columns=genes)
    return AlignedExpression(matrix=aligned, log2=log2, observed_mask=observed,
                             missing_mask=missing_mask, qc=qc)


def zero_equivalent_tpm(missing_mask: pd.DataFrame, data_min: Sequence[float],
                        genes: Sequence[str]) -> pd.DataFrame:
    """``zero`` 策略的**输入侧等价写法**：``TPM_g = 2^{data_min_g} − 1``（缺失处）。

    仅供 golden test 证明「标准化空间置 0」与「输入侧等价式」一致（判据 max|Δ| ≤ 1e-6）；
    正式实现走 representation 层的标准化空间掩码（与训练同一操作）。
    非缺失位置返回 NaN（不参与比较）。
    """
    genes = [str(g) for g in genes]
    dmin = np.asarray(data_min, dtype=np.float64)
    if dmin.size != len(genes):
        raise InputError(f"data_min 长度 {dmin.size} != 基因数 {len(genes)}")
    eq = np.power(2.0, dmin) - 1.0
    mask = missing_mask.reindex(columns=genes).to_numpy(bool)
    return pd.DataFrame(np.where(mask, eq[None, :], np.nan),
                        index=missing_mask.index, columns=genes)


def cancer_codes_for(cancer_type: Iterable) -> pd.Series:
    """用户 ``cancer_type`` → **真实 COMPASS 癌种码**（绝不使用占位码）。"""
    table = cancer_codes().set_index("cancer_type")
    ct = pd.Series(list(cancer_type), name="cancer_type")
    if ct.isna().any():
        raise UnknownCancerTypeError(
            f"cancer_type 存在缺失（{int(ct.isna().sum())} 处）；癌种必须逐个提供，不得推断")
    unknown = sorted(set(ct.astype(str)) - set(table.index))
    if unknown:
        raise UnknownCancerTypeError(
            f"未知癌种 {unknown}；可用取值见 src/compass_os/data/cancer_codes.tsv（TCGA 缩写）")
    codes = ct.astype(str).map(table["compass_code"]).astype(int)
    # 负码值是**正常组织**标记（NORMAL = -1），不是模型支持的肿瘤类型：
    # 直接送入 embedding 会抛 IndexError，故在此给出明确错误。
    non_tumour = sorted(set(ct.astype(str)[codes < 0]))
    if non_tumour:
        raise UnknownCancerTypeError(
            f"癌种 {non_tumour} 在 COMPASS 码表中是正常组织标记（compass_code < 0），"
            "不是本模型支持的肿瘤类型；请改用真正的 TCGA 肿瘤缩写"
            "（可用取值见 src/compass_os/data/cancer_codes.tsv）。")
    return codes


def ct_onehot_columns(cancer_type: Iterable) -> pd.DataFrame:
    """``CT_*`` one-hot（32 列，与锁定特征一致）。

    不在列中的癌种（训练队列的 drop-first 参照水平）⇒ 该样本 **32 列全 0**，
    这与生产外部验证口径（``12_validate_external.py`` 的 D38 分支）一致。
    """
    table = cancer_codes().set_index("cancer_type")
    ct_cols = [c for c in _lock_feature_names("M2") if c.startswith("CT_")]
    ct = pd.Series(list(cancer_type), name="cancer_type").astype(str)
    unknown = sorted(set(ct) - set(table.index))
    if unknown:
        raise UnknownCancerTypeError(f"未知癌种 {unknown}")
    non_tumour = sorted(set(ct[ct.map(table["compass_code"]).astype(int) < 0]))
    if non_tumour:
        raise UnknownCancerTypeError(
            f"癌种 {non_tumour} 是正常组织标记（compass_code < 0），不是支持的肿瘤类型")
    X = pd.DataFrame(0.0, index=pd.RangeIndex(len(ct)), columns=ct_cols)
    ct_col = table["ct_column"]
    for i, name in enumerate(ct):
        col = ct_col.get(name, "")
        if isinstance(col, str) and col in ct_cols:
            X.iat[i, ct_cols.index(col)] = 1.0
    return X


@functools.lru_cache(maxsize=4)
def _lock_feature_names(model: str) -> tuple:
    from .survival import lock_path
    lock = json.loads(lock_path(model).read_text(encoding="utf-8"))
    return tuple(lock["feature_names"])
