# -*- coding: utf-8 -*-
"""compass_os —— COMPASS-OS：泛癌 OS 风险预测 + 机制导向表征的一体化工具包。

Public API（只有这些是给使用者用的）
------------------------------------
**Level 1 · Quick（一键式）**

* :func:`analyze` —— 表达（+ 癌种，可选临床与随访）→ 风险 + 43 concept / 132 signature
  表征 + 输入质控 + 人类可读报告；结果对象 :class:`AnalysisResult`
  （``summary()`` / ``save_report()``）。

**Level 2 · Advanced（研究型）**

* :func:`predict` —— 指定模型（M0/M1/M2/M3）、缺失策略、队列分层
* :func:`get_representation` —— 只要 132 / 43 表征
* :func:`check_robustness` —— 显式的缺失基因稳健性比较
  结果对象 :class:`PredictionResult` / :class:`RepresentationResult` / :class:`RobustnessResult`

异常：:class:`CompassOSError` 及其子类。

**内部实现**（``preprocessing`` / ``representation`` / ``survival`` / ``qc`` / ``plotting``）
不属公开契约；其函数签名可能变动，README 与 docs 不把它们作为使用文档。
模型定义（scaler / PCA / Cox 系数 / 词表 / 癌种码 / 参考中位数）**只能读冻结资产**，
不提供任何 refit / fit / 自定义码接口——需要改动模型定义请 fork 源码。

措辞纪律
--------
43 concepts 与 132 signatures 是 **mechanism-oriented biological representations /
hypothesis-generating features**，用于 **downstream prioritisation**；
不得表述为 pathway activation / suppression / causal mechanism。
"""
from __future__ import annotations

from .analysis import AnalysisResult, analyze
from .api import (PredictionResult, RepresentationResult, RobustnessResult,
                  check_robustness, get_representation, predict)
from .exceptions import (AssetNotFoundError, CompassOSError, InputError,
                         MissingGenesError, UnknownCancerTypeError)

__version__ = "1.0.0"

#: Level 1（Quick）与 Level 2（Advanced）正式导出；其余模块为内部实现。
__all__ = [
    # Level 1
    "analyze", "AnalysisResult",
    # Level 2
    "predict", "get_representation", "check_robustness",
    "PredictionResult", "RepresentationResult", "RobustnessResult",
    # exceptions
    "CompassOSError", "MissingGenesError", "UnknownCancerTypeError",
    "AssetNotFoundError", "InputError",
    "__version__",
]
