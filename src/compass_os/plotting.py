# -*- coding: utf-8 -*-
"""绘图：KM 曲线 + log-rank、风险分布、concept/signature 谱。

约定
----
* 只依赖 matplotlib（无 seaborn/lifelines）；log-rank 用 ``sksurv.compare``；
* 生存曲线用 KM 估计（非模型预测），风险分层默认按**冻结 median_cutoff**；
* **不提供**绝对生存概率的可视化默认项——绝对风险校准限制见
  :mod:`compass_os.survival`。
"""
from __future__ import annotations

from typing import Sequence

import numpy as np
import pandas as pd


def _plt():
    import matplotlib
    matplotlib.use("Agg", force=False)
    import matplotlib.pyplot as plt
    return plt


def km_estimate(time: Sequence[float], event: Sequence[bool]) -> tuple:
    """Kaplan–Meier 估计：返回 (times, survival)，含 t=0 起点。"""
    t = np.asarray(time, dtype=float)
    e = np.asarray(event, dtype=bool)
    ok = np.isfinite(t)
    t, e = t[ok], e[ok]
    order = np.argsort(t)
    t, e = t[order], e[order]
    n = len(t)
    times, surv = [0.0], [1.0]
    s, at_risk = 1.0, n
    i = 0
    while i < n:
        ti = t[i]
        d = int(e[i:][t[i:] == ti].sum())
        c = int((~e[i:][t[i:] == ti]).sum())
        if d:
            s *= (1.0 - d / at_risk)
            times.append(float(ti))
            surv.append(float(s))
        at_risk -= (d + c)
        i += d + c
    return np.asarray(times), np.asarray(surv)


def logrank(time: Sequence[float], event: Sequence[bool], group: Sequence) -> tuple:
    """log-rank 检验（``sksurv.compare.compare_survival``）：返回 (chi2, p)。"""
    from sksurv.compare import compare_survival
    from sksurv.util import Surv
    g = np.asarray(group)
    ok = pd.notna(g)
    y = Surv.from_arrays(np.asarray(event, dtype=bool)[ok], np.asarray(time, float)[ok])
    chisq, p = compare_survival(y, g[ok].astype(str), return_stats=False)
    return float(chisq), float(p)


def km_plot(risk: Sequence[float], time: Sequence[float], event: Sequence[bool], *,
            group: Sequence | None = None, cutoff: float | None = None,
            title: str = "Kaplan-Meier by risk group", ax=None, label_high: str = "high risk",
            label_low: str = "low risk"):
    """按风险分层的 KM 曲线（含 log-rank p）。

    ``group`` 未给时按 ``cutoff``（应传冻结的 ``median_cutoff``）二分为 high/low。
    """
    plt = _plt()
    if ax is None:
        _fig, ax = plt.subplots(figsize=(4.6, 4.0), dpi=160)
    r = np.asarray(risk, dtype=float)
    if group is None:
        if cutoff is None:
            cutoff = float(np.median(r))
        group = np.where(r >= cutoff, "high", "low")
    g = np.asarray(group).astype(str)
    p = None
    for name, lab, color in (("high", label_high, "#d62728"), ("low", label_low, "#1f77b4")):
        m = g == name
        if m.sum() == 0:
            continue
        tt, ss = km_estimate(np.asarray(time)[m], np.asarray(event)[m])
        ax.step(tt, ss, where="post", color=color, lw=1.4, label=f"{lab} (n={int(m.sum())})")
    if len(set(g)) > 1:
        try:
            _c, p = logrank(time, event, g)
        except Exception:  # 组内事件不足等
            p = None
    if p is not None:
        ax.set_title(f"{title}\nlog-rank p = {p:.3g}", fontsize=8)
    else:
        ax.set_title(title, fontsize=8)
    ax.set_xlabel("time (days)", fontsize=8)
    ax.set_ylabel("overall survival", fontsize=8)
    ax.set_ylim(0, 1.02)
    ax.tick_params(labelsize=8)
    ax.legend(fontsize=7, frameon=False)
    return ax


def risk_distribution(risk: Sequence[float], *, cutoff: float | None = None, ax=None,
                      title: str = "Risk score distribution"):
    """风险分布直方图（可叠加冻结切点）。"""
    plt = _plt()
    if ax is None:
        _fig, ax = plt.subplots(figsize=(4.6, 3.4), dpi=160)
    r = np.asarray(risk, dtype=float)
    ax.hist(r, bins=min(30, max(5, len(r) // 5)), color="#8c8c8c", alpha=0.85)
    if cutoff is not None:
        ax.axvline(cutoff, color="#d62728", ls="--", lw=1.2,
                   label=f"frozen cutoff = {cutoff:.3g}")
        ax.legend(fontsize=7, frameon=False)
    ax.set_xlabel("linear predictor (relative risk)", fontsize=8)
    ax.set_ylabel("samples", fontsize=8)
    ax.set_title(title, fontsize=8)
    ax.tick_params(labelsize=8)
    return ax


def representation_profile(scores: pd.DataFrame, *, top: int = 30, ax=None,
                           title: str = "Representation profile", by: str = "mean"):
    """逐样本的 concept/signature 谱（折线）或热图（``top`` 个方差最大的特征）。

    ``scores``：samples × features。样本数 > ``top`` 时画热图，否则画折线。
    """
    plt = _plt()
    if by == "variance" or scores.shape[0] > top:
        v = scores.var(axis=0).sort_values(ascending=False)
        cols = list(v.index[:top])
        mat = scores[cols].to_numpy(float)
        if ax is None:
            _fig, ax = plt.subplots(figsize=(6.0, 0.18 * scores.shape[0] + 1.6), dpi=160)
        im = ax.imshow(mat, aspect="auto", cmap="RdBu_r",
                       vmin=-np.nanmax(np.abs(mat)), vmax=np.nanmax(np.abs(mat)))
        ax.set_yticks(range(scores.shape[0]))
        ax.set_yticklabels([str(i) for i in scores.index], fontsize=5)
        ax.set_xticks(range(len(cols)))
        ax.set_xticklabels(cols, rotation=90, fontsize=5)
        ax.figure.colorbar(im, ax=ax, shrink=0.7, label="score")
        ax.set_title(title, fontsize=8)
        return ax
    if ax is None:
        _fig, ax = plt.subplots(figsize=(6.4, 3.0), dpi=160)
    for sid, row in scores.iterrows():
        ax.plot(range(scores.shape[1]), row.to_numpy(float), lw=0.9, label=str(sid))
    ax.set_xticks(range(scores.shape[1]))
    ax.set_xticklabels(list(scores.columns), rotation=90, fontsize=5)
    ax.set_ylabel("score", fontsize=8)
    ax.set_title(title, fontsize=8)
    ax.tick_params(labelsize=7)
    if scores.shape[0] <= 8:
        ax.legend(fontsize=6, frameon=False)
    return ax


def plot_robustness(curves: dict, *, xlabel: str = "gene coverage", ax=None,
                    title: str = "Missing-gene robustness (continuous)"):
    """稳健性曲线：``curves`` = {label: (coverage_array, stability_array)}。

    ⚠ 本图**不画任何阈值线**——阈值需由 stress test 数据确定后写入 ``qc_config.json``。
    """
    plt = _plt()
    if ax is None:
        _fig, ax = plt.subplots(figsize=(4.8, 3.6), dpi=160)
    for lab, (x, y) in curves.items():
        ax.plot(np.asarray(x, float), np.asarray(y, float), marker="o", ms=3, lw=1.2, label=lab)
    ax.set_xlabel(xlabel, fontsize=8)
    ax.set_ylabel("stability (sample-wise Pearson r)", fontsize=8)
    ax.set_title(title, fontsize=8)
    ax.tick_params(labelsize=8)
    ax.legend(fontsize=7, frameon=False)
    return ax
