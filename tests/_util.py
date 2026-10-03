# -*- coding: utf-8 -*-
"""测试公用工具：夹具加载 + 现场重跑生产实现（可选）。

夹具的期望值由原项目生产实现生成，见 ``tools/build_fixtures.py`` 与
``tests/fixtures/golden_provenance.json``。本模块只读夹具，不修改任何资产。
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
FIXTURES = REPO / "tests" / "fixtures"
SRC = REPO / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

#: 判据（P1 规格）：golden 复现 max|Δ| ≤ 1e-5；zero 两写法等价 max|Δ| ≤ 1e-6
TOL_GOLDEN = 1e-5
TOL_ZERO_EQ = 1e-6


def fixtures_available() -> bool:
    return (FIXTURES / "golden_expected.tsv.gz").is_file()


def load_fixtures() -> dict:
    """返回夹具 dict：``expression`` / ``clinical`` / ``labels`` / 各 ``expected`` 表。"""
    expr = pd.read_csv(FIXTURES / "golden_expression.tsv.gz", sep="\t", index_col=0)
    clin = pd.read_csv(FIXTURES / "golden_clinical.tsv", sep="\t", index_col=0)
    lab = pd.read_csv(FIXTURES / "golden_labels.tsv", sep="\t", index_col=0)
    wide = pd.read_csv(FIXTURES / "golden_expected.tsv.gz", sep="\t", index_col=0)
    groups: dict = {}
    for col in wide.columns:
        head, name = col.split("::", 1)
        groups.setdefault(head, {})[name] = wide[col]
    out = {"expression": expr, "clinical": clin, "labels": lab}
    for head, cols in groups.items():
        out[head] = pd.DataFrame(cols)
    prov = json.loads((FIXTURES / "golden_provenance.json").read_text(encoding="utf-8"))
    out["provenance"] = prov
    return out


def max_abs_diff(got: pd.DataFrame, expected: pd.DataFrame) -> float:
    """对齐列/行后取 max|Δ|（缺列/缺行视为失败）。"""
    cols = [c for c in expected.columns if c in got.columns]
    if len(cols) != expected.shape[1]:
        missing = sorted(set(expected.columns) - set(got.columns))[:5]
        raise AssertionError(f"缺少列（前几个）：{missing}")
    idx = [i for i in expected.index if i in set(got.index)]
    if len(idx) != expected.shape[0]:
        raise AssertionError("缺少样本行")
    a = got.loc[idx, cols].to_numpy(np.float64)
    b = expected.loc[idx, cols].to_numpy(np.float64)
    return float(np.abs(a - b).max())


def pearson(got: pd.DataFrame, expected: pd.DataFrame) -> float:
    idx = [i for i in expected.index if i in set(got.index)]
    cols = [c for c in expected.columns if c in got.columns]
    return float(np.corrcoef(got.loc[idx, cols].to_numpy(np.float64).ravel(),
                             expected.loc[idx, cols].to_numpy(np.float64).ravel())[0, 1])


def production_source_root() -> Path | None:
    """原项目根（仅当 ``COMPASS_OS_SOURCE_ROOT`` 环境变量给出且存在时返回）。"""
    root = os.environ.get("COMPASS_OS_SOURCE_ROOT")
    if root and (Path(root) / "scripts").is_dir():
        return Path(root)
    return None
