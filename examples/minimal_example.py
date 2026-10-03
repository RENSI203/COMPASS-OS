#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""最小示例：单样本预测（含 132 signature / 43 concept / M0–M3 风险 / QC）。

运行：
    MPLCONFIGDIR=/tmp/mplcfg python examples/minimal_example.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

import compass_os  # noqa: E402

DATA = REPO / "examples" / "example_data"


def main() -> int:
    expr = pd.read_csv(DATA / "example_expression.tsv.gz", sep="\t", index_col=0)
    clin = pd.read_csv(DATA / "example_clinical.tsv", sep="\t", index_col=0)
    lab = pd.read_csv(DATA / "example_labels.tsv", sep="\t", index_col=0)

    one = expr.iloc[[0]]
    sid = str(one.index[0])
    cancer = lab.loc[sid, "cancer_type"]
    print(f"样本 {sid}｜癌种 {cancer}｜输入 {one.shape[0]} 例 × {one.shape[1]} 基因")

    res = compass_os.predict(one, [cancer], clinical=clin.loc[[sid]])
    from compass_os.api import available_models
    print(f"\n默认模型：{res.default_model}（可用：{available_models()}）")
    print("相对风险（线性预测子）：")
    print(res.risk.T.to_string())
    print("\nQC：")
    for k, v in res.qc.to_dict().items():
        if k in ("missing_gene_names", "warnings"):
            continue
        print(f"  {k}: {v}")
    print(f"  signature 覆盖（132 维中位）：{res.qc.signature_gene_coverage.to_numpy().mean():.4f}")
    print(f"  concept 输入覆盖（43 维中位）：{res.qc.concept_input_coverage.to_numpy().mean():.4f}")

    top = res.concept_scores.iloc[0].sort_values(ascending=False).head(5)
    print("\n该样本 concept 分数最高/最低的 5 个（**仅作 association / prioritization**，"
          "不代表 activation/suppression）：")
    print("  最高：", ", ".join(f"{k}={v:.3f}" for k, v in top.items()))
    low = res.concept_scores.iloc[0].sort_values().head(5)
    print("  最低：", ", ".join(f"{k}={v:.3f}" for k, v in low.items()))

    print("\n⚠ 单样本不产生 high/low 分层（risk_group=None）："
          f"{res.risk_group is None or res.risk_group.shape[0] == 1}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
