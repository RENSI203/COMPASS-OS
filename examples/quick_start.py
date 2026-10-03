#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Quick start —— 三行完成一次标准分析（Level 1 API）。

运行：
    MPLCONFIGDIR=/tmp/mplcfg python examples/quick_start.py
产物：examples/output/quick_start/（报告目录）
"""
from __future__ import annotations

import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

import compass_os  # noqa: E402

DATA = REPO / "examples" / "example_data"


def main() -> int:
    result = compass_os.analyze(
        DATA / "example_expression.tsv.gz",        # samples × genes（TPM）
        cancer_type="LUAD",                         # TCGA 缩写；单一字符串 = 全队列同一癌种
        clinical=DATA / "example_clinical.tsv",     # 可选：age / sex / stage
        survival=DATA / "example_labels.tsv",       # 可选：自动识别 time / event 列
    )

    result.summary()
    files = result.save_report(REPO / "examples" / "output" / "quick_start")
    print(f"\nsave_report 生成 {len(files)} 个文件 → "
          f"{REPO / 'examples' / 'output' / 'quick_start'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
