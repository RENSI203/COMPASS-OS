#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Sample-level computation-path Sankey —— M1 / M2 / M3 三个示例。

运行：
    MPLCONFIGDIR=/tmp/mplcfg python examples/sample_path_example.py

产物（examples/output/sample_path/）：
    sample_path_M1.html / .png
    sample_path_M2.html / .png
    sample_path_M3.html / .png
    <sample>.decomposition.tsv（完整逐特征 Cox 分解）

语义边界：这是 **sample-specific computational attribution / representation-flow
visualization**，不是因果机制图、不是通路激活图、也不是对最终风险的完整因果解释。
"""
from __future__ import annotations

import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

import pandas as pd  # noqa: E402

import compass_os  # noqa: E402
from compass_os import sample_path  # noqa: E402

DATA = REPO / "examples" / "example_data"
OUT = REPO / "examples" / "output" / "sample_path"


def main() -> int:
    expr = pd.read_csv(DATA / "example_expression.tsv.gz", sep="\t", index_col=0)
    clin = pd.read_csv(DATA / "example_clinical.tsv", sep="\t", index_col=0)
    lab = pd.read_csv(DATA / "example_labels.tsv", sep="\t", index_col=0)
    cancer = lab["cancer_type"].tolist()
    sample_id = str(expr.index[0])

    OUT.mkdir(parents=True, exist_ok=True)
    print(f"cohort: {len(expr)} samples | {len(set(cancer))} cancer type(s) | "
          f"sample of interest: {sample_id}\n")

    for model in ("M1", "M2", "M3"):
        res = sample_path(expr, cancer, sample_id, model=model, clinical=clin)
        html = res.save_html(OUT / f"sample_path_{model}.html")
        try:
            img = res.save_static(OUT / f"sample_path_{model}.png")
            static = f"{img.name} ({img.stat().st_size // 1024} KB)"
        except Exception as exc:  # noqa: BLE001
            static = f"静态导出不可用（{type(exc).__name__}）"
        res.decomposition.to_csv(OUT / f"{sample_id}.{model}.decomposition.tsv",
                                 sep="\t", index=False)
        print(f"=== {model} ===")
        print(res.summary())
        print(f"→ {html.name} ({html.stat().st_size // 1024} KB) | {static}\n")

    print(f"全部产物见 {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
