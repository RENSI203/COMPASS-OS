# -*- coding: utf-8 -*-
"""生存风险：从冻结 lock 组装 Cox 设计矩阵并计算风险。

**所有参数只读冻结资产**：scaler、PCA、Cox 系数、临床填充值、癌种 one-hot 列、
中位切点、基线累积风险。本模块**不拟合任何东西**，并显式禁止在用户数据上 fit。

风险口径（与生产一致）::

    X = [CT_*(32) | Age, Sex, Stage | 43 concepts | (M3: PC1..PC10)]
    X[scaled_features] = (X[scaled_features] − scaler_mean) / scaler_scale
    risk = X @ beta                      # 线性预测子（相对风险）

绝对生存（**派生输出**，非核心卖点）::

    S_i(t) = exp( −H0(t) · exp(risk) )    # H0 取 lock 的 base_time / base_cumhazard

⚠ 校准限制：``base_time/base_cumhazard`` 来自训练域（TCGA）的 Breslow 基线，外部验证中
strict 口径的绝对风险校准明显弱于判别度。因此 ``survival_probability`` 属于
baseline-hazard-based derived output，**不得**宣传为稳定的个体绝对风险预测。
"""
from __future__ import annotations

import functools
import json
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd

from ._paths import asset, model_manifest
from .exceptions import InputError

CLINICAL_FIELDS = ("age", "sex", "stage")


def default_model() -> str:
    """当前默认模型（**由 ``models/model_manifest.json`` 决定**，改默认无需改代码）。"""
    return str(model_manifest().get("default_model", "M2"))


def available_models() -> tuple:
    return tuple(model_manifest()["models"].keys())


def lock_path(model: str) -> Path:
    mm = model_manifest()
    entry = mm["models"].get(str(model))
    if entry is None:
        raise InputError(f"未知模型 {model!r}；可用：{sorted(mm['models'])}")
    return asset(entry["lock"])


@functools.lru_cache(maxsize=8)
def load_lock(model: str) -> dict:
    return json.loads(lock_path(model).read_text(encoding="utf-8"))


@functools.lru_cache(maxsize=2)
def load_pca() -> dict:
    z = np.load(asset("models/locked_pca_M3.npz"), allow_pickle=True)
    return {"genes": [str(g) for g in z["genes"]], "mean": np.asarray(z["mean"], float),
            "scale": np.asarray(z["scale"], float),
            "components": np.asarray(z["components"], float)}


def build_design_matrix(concept_scores: pd.DataFrame, clinical: pd.DataFrame | None,
                        log2_matrix: pd.DataFrame, model: str,
                        cancer_type: Iterable) -> pd.DataFrame:
    """按 lock 的 ``feature_names`` 严格顺序组装设计矩阵（含 M3 的锁定 PCA）。"""
    from .preprocessing import ct_onehot_columns

    lock = load_lock(model)
    feats = list(lock["feature_names"])
    idx = [str(i) for i in concept_scores.index]
    X = pd.DataFrame(index=idx)

    for c in concept_scores.columns:
        if c in feats:
            X[c] = concept_scores[c].astype(np.float64).to_numpy()
    if any(f in feats for f in ("Age", "Sex", "Stage")):
        clin = _clinical_frame(clinical, idx)
        for nm, src in (("Age", "age"), ("Sex", "sex"), ("Stage", "stage")):
            if nm in feats:
                fill = float(lock["clinical_fill"][nm]["fill"])
                X[nm] = clin[src].fillna(fill).astype(np.float64).to_numpy()

    ct = ct_onehot_columns(list(cancer_type))
    ct.index = idx
    for c in [c for c in feats if c.startswith("CT_")]:
        X[c] = ct[c].to_numpy() if c in ct.columns else 0.0

    if any(f.startswith("PC") for f in feats):
        z = load_pca()
        mat = log2_matrix.reindex(index=idx, columns=z["genes"]).to_numpy(np.float64)
        std = (mat - z["mean"][None, :]) / z["scale"][None, :]
        pcs = std @ z["components"].T
        for i in range(pcs.shape[1]):
            X[f"PC{i + 1}"] = pcs[:, i]

    missing = [f for f in feats if f not in X.columns]
    if missing:
        raise InputError(f"设计矩阵缺少锁定特征 {missing[:8]}（共 {len(missing)} 个）")
    return X[feats]


def _clinical_frame(clinical: pd.DataFrame | None, idx: Sequence[str]) -> pd.DataFrame:
    """标准化用户临床表（列名 ``age`` / ``sex`` / ``stage``；缺失 → NaN → 用锁定常数填）。"""
    out = pd.DataFrame(index=list(idx))
    if clinical is None:
        for f in CLINICAL_FIELDS:
            out[f] = np.nan
        return out
    if not isinstance(clinical, pd.DataFrame):
        raise InputError("clinical 必须是 DataFrame（index = sample_id）")
    clin = clinical.copy()
    clin.index = [str(i) for i in clin.index]
    clin = clin.reindex([str(i) for i in idx])
    for f in CLINICAL_FIELDS:
        out[f] = pd.to_numeric(clin[f], errors="coerce") if f in clin.columns else np.nan
    # stage 只接受 1–4（生产口径：0/缺失 → NaN → 锁定常数）
    out["stage"] = out["stage"].where(out["stage"].between(1, 4), np.nan)
    return out


def risk_from_design(X: pd.DataFrame, model: str) -> np.ndarray:
    """标准化（仅 ``scaled_features``）后与冻结 ``beta`` 相乘。"""
    lock = load_lock(model)
    sc = list(lock.get("scaled_features") or [])
    Z = X.astype(np.float64).copy()
    if sc and lock.get("scaler_mean") is not None:
        mu = np.asarray(lock["scaler_mean"], dtype=np.float64)
        sd = np.asarray(lock["scaler_scale"], dtype=np.float64)
        if sd.size != len(sc):
            raise InputError(f"{model} 的 scaler_scale 长度 {sd.size} != scaled_features {len(sc)}")
        Z[sc] = (Z[sc].to_numpy(np.float64) - mu) / sd
    beta = np.asarray(lock["beta"], dtype=np.float64)
    if beta.size != Z.shape[1]:
        raise InputError(f"{model} 的 beta 长度 {beta.size} != 设计矩阵列数 {Z.shape[1]}")
    return Z.to_numpy(np.float64) @ beta


def median_cutoff(model: str) -> float:
    """冻结的中位风险切点（训练集 risk 中位数），用于队列内 high/low 分层。"""
    return float(load_lock(model)["median_cutoff"])


def survival_probability(risk: Sequence[float], times: Sequence[float], model: str,
                         sample_ids: Sequence[str] | None = None) -> pd.DataFrame:
    """``S_i(t) = exp(−H0(t)·exp(risk))``（派生输出；见模块 docstring 的校准限制）。"""
    lock = load_lock(model)
    bt = np.asarray(lock["base_time"], dtype=np.float64)
    h0 = np.asarray(lock["base_cumhazard"], dtype=np.float64)
    order = np.argsort(bt)
    bt, h0 = bt[order], h0[order]
    t = np.atleast_1d(np.asarray(times, dtype=np.float64))
    er = np.exp(np.asarray(risk, dtype=np.float64))
    idx = np.clip(np.searchsorted(bt, t, side="right") - 1, 0, len(bt) - 1)
    h = h0[idx].copy()
    h[t < bt[0]] = 0.0
    idx_out = [str(s) for s in sample_ids] if sample_ids is not None else \
        [f"sample_{i}" for i in range(len(er))]
    return pd.DataFrame(np.exp(-np.outer(h, er)).T, index=idx_out, columns=t)
