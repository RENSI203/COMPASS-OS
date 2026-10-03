# -*- coding: utf-8 -*-
"""COMPASS 表征计算：132 gene-signature 分数与 43 concept 分数。

数值链路（与生产 pipeline 逐位对齐）
------------------------------------
::

    输入 TPM (samples × 15,672)
      → 首列插入真实 cancer_code（int）           ← 模型要求首列为癌种码
      → 冻结 Datascaler：log2(x+1) + MinMax        ← 在 checkpoint 内，禁止重拟合
      → Compass encoder（gene tokens + PID + CANCER token）
      → DisentangledProjector.forward
             ├─ GeneSetProjector     → geneset scores (132)
             └─ CellPathwayProjector → concept  scores (43)   （由 132 经 softmax 注意力聚合）

取数路径：**官方 API**
---------------------
使用上游 ``PreTrainer.extract()``（一次返回 132 gene-set 层与 43 cell-pathway 层）；
**不再使用 monkey patch**。

依据（golden 实测，12 例）：官方 ``extract()`` 与生产期使用的
``predict()`` + ``enable_geneset_capture`` 路径结果**逐位相同**——
132 与 43 两层均为 ``max|Δ| = 0.000e+00``、Pearson = 1.0000000000。
因此改用官方 API **不改变任何数值结果**（见 ``tests/test_reproducibility.py``）。

⚠ 癌种码口径：一律使用**真实 cancer type → 真实 COMPASS code**。
项目内的机制分析存档（``analysis_out_concepts``，脚本 ``41``）曾用**常量占位码 5**，
两者产物可差到 ~1.4e-3；本包**不采用**占位码。
"""
from __future__ import annotations

import functools
from contextlib import contextmanager
from typing import Dict, Mapping, Sequence

import numpy as np
import pandas as pd

from ._paths import ensure_compass


@functools.lru_cache(maxsize=2)
def load_model(device: str = "cpu"):
    """加载冻结的 COMPASS checkpoint（带缓存；进程内只加载一次）。

    ⚠ v1.0.1 起**不再调用上游 ``compass.loadcompass``**。该函数末尾有：

    .. code-block:: python

        if file.startswith(tempfile.gettempdir()):
            os.remove(file)

    原意是清理"从 URL 下载得到的临时文件"，但判据是**路径前缀**：凡是位于系统临时目录
    （如 ``/tmp``）下的模型文件，加载后都会被 ``os.remove`` **静默删除**。因此任何装在
    ``/tmp`` 下的 Python 环境（容器、CI、``pip install --target /tmp/...`` 等）在第一次
    预测后就会丢失 checkpoint，后续调用全部失败。

    本函数用 ``torch.load`` 直接加载本地冻结资产（与上游同一个调用），并复刻上游那两处
    非破坏性调整，从而在**不修改 vendored 上游源码**的前提下消除该副作用。
    数值行为与上游完全一致（golden 复现验证）。
    """
    ensure_compass()
    import torch
    from ._paths import asset

    path = asset("models/pretrainer.pt")
    model = torch.load(str(path), weights_only=False, map_location=device)
    # —— 复刻上游 loadcompass 的非破坏性后处理 ——
    if hasattr(model, "with_wandb") and model.with_wandb:
        model.wandb._settings = ""
    if device == "cpu":
        model.device = "cpu"
    return model


def checkpoint_feature_names() -> tuple:
    """模型 ``feature_name``（15,672 基因）。"""
    return tuple(str(g) for g in load_model().feature_name)


def checkpoint_data_min() -> np.ndarray:
    """冻结 MinMaxScaler 的逐基因 ``data_min_``（log2 空间下界）。"""
    sc = load_model().scaler.scaler
    return np.asarray(sc.data_min_, dtype=np.float64)


def concept_weights() -> Dict[str, tuple]:
    """43 概念 → (成员基因集名, 冻结 softmax 注意力权重)。

    权重取自 checkpoint 的 ``CellPathwayAttentionAggregator``（非负、组内和为 1），
    用于 :func:`compass_os.qc.concept_input_coverage`。
    """
    import torch
    lp = load_model().model.latentprojector
    cp = lp.cellpathwayprojector
    agg = cp.cellpathway_aggregator.aggregator
    try:
        names = list(cp.cellpathway_names)
    except AttributeError:  # 上游属性名兼容
        names = list(lp.CELLPATHWAY.index)
    out: Dict[str, tuple] = {}
    for k, concept in enumerate(names):
        wt = agg.attention_weights[f"cellpathway_{k}"].detach()
        w = torch.softmax(wt, dim=0).numpy().ravel().astype(float)
        members = lp.CELLPATHWAY.iloc[k].tolist() if hasattr(
            lp.CELLPATHWAY.iloc[k], "tolist") else list(lp.CELLPATHWAY.iloc[k])
        member_names = [str(lp.GENESET.index[i]) for i in members]
        out[str(concept)] = (member_names, w)
    return out


class _ZeroMaskedScaler:
    """把冻结 Datascaler 的**标准化结果**中缺失基因列置 0（训练空间掩码）。

    只在 ``missing_gene_strategy="zero"`` 时临时包裹；``fit`` 被禁止，防止误拟合。
    """

    def __init__(self, inner, genes: Sequence[str]):
        self._inner = inner
        self._genes = [str(g) for g in genes]
        self.scale_method = getattr(inner, "scale_method", None)

    def transform(self, dfcx: pd.DataFrame) -> pd.DataFrame:
        out = self._inner.transform(dfcx)
        cols = [g for g in self._genes if g in out.columns]
        if cols:
            out = out.copy()
            out[cols] = 0.0
        return out

    def fit(self, *a, **k):  # pragma: no cover - 明确禁止
        raise RuntimeError("禁止在用户数据上拟合 scaler：模型资产必须只读使用")

    def fit_transform(self, *a, **k):  # pragma: no cover - 明确禁止
        raise RuntimeError("禁止在用户数据上拟合 scaler：模型资产必须只读使用")


@contextmanager
def _zero_masked(model, genes: Sequence[str]):
    inner = model.scaler
    model.scaler = _ZeroMaskedScaler(inner, genes)
    try:
        yield
    finally:
        model.scaler = inner


# --------------------------------------------------------------------------- #
# 主入口
# --------------------------------------------------------------------------- #
def compute_representation(tpm: pd.DataFrame, cancer_code: Sequence[int], *,
                           missing_mask: pd.DataFrame | None = None,
                           strategy: str = "reference",
                           batch_size: int = 16, device: str = "cpu",
                           model=None) -> Dict[str, pd.DataFrame]:
    """从 TPM（samples × 15,672，列序 = 词表）计算 132 与 43 两层分数。

    ``missing_mask`` 仅在 ``strategy="zero"`` 时需要（给出待掩码的基因位置）。
    返回 ``{"signature_scores": (n,132), "concept_scores": (n,43)}``。
    """
    model = model if model is not None else load_model(device)
    genes = [str(g) for g in model.feature_name]

    X = tpm.reindex(columns=genes).astype(np.float64).copy()
    X.insert(0, "cancer_code", np.asarray(list(cancer_code), dtype=int))
    if X.shape[0] != len(list(cancer_code)):
        raise InputError("cancer_code 长度与样本数不一致")

    masked_genes: list = []
    if strategy == "zero" and missing_mask is not None:
        mm = missing_mask.reindex(columns=genes).fillna(False).to_numpy(bool)
        masked_genes = [g for j, g in enumerate(genes) if mm[:, j].any()]
        # 输入侧统一填 0 占位（随后被 _ZeroMaskedScaler 在标准化空间覆盖为 0）
        X.loc[:, masked_genes] = 0.0

    with _zero_masked(model, masked_genes):
        # 上游 `PreTrainer.extract(with_gene_level=False)` 返回 **(dfgs, dfct)** 两个值
        dfgs, dfct = model.extract(X, batch_size=batch_size, num_workers=0,
                                   with_gene_level=False)
    sig = dfgs.drop(columns=[c for c in ("CANCER", "PID") if c in dfgs.columns])
    con = dfct.drop(columns=[c for c in ("CANCER", "PID") if c in dfct.columns])

    sig.index = [str(i) for i in sig.index]
    con.index = [str(i) for i in con.index]
    return {"signature_scores": sig, "concept_scores": con}
