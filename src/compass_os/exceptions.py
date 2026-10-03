# -*- coding: utf-8 -*-
"""compass_os 的异常类型。"""
from __future__ import annotations

from typing import Sequence


class CompassOSError(Exception):
    """本包所有异常的基类。"""


class AssetNotFoundError(CompassOSError):
    """模型资产（checkpoint / lock / PCA / 数据表）缺失或无法解析。"""


class MissingGenesError(CompassOSError):
    """``missing_gene_strategy="strict"`` 下存在缺失的必需基因。"""

    def __init__(self, missing: Sequence[str], n_required: int, message: str | None = None):
        self.missing_genes = sorted(str(g) for g in missing)
        self.n_missing = len(self.missing_genes)
        self.n_required = int(n_required)
        head = ", ".join(self.missing_genes[:10])
        more = "" if self.n_missing <= 10 else f" …（共 {self.n_missing} 个）"
        super().__init__(message or (
            f"输入缺少 {self.n_missing}/{self.n_required} 个 COMPASS 必需基因"
            f"（missing_gene_strategy='strict' 不允许填补）：{head}{more}"))


class UnknownCancerTypeError(CompassOSError):
    """``cancer_type`` 缺失或不在冻结的 COMPASS 癌种码表内。"""


class InputError(CompassOSError):
    """输入表达矩阵的形状/尺度/取值不合法。"""
