#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""缺失基因策略示例：reference / zero / strict 的差别与 QC + 稳健性检查。

运行：
    MPLCONFIGDIR=/tmp/mplcfg python examples/missing_genes_example.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

import compass_os  # noqa: E402
from compass_os.exceptions import MissingGenesError  # noqa: E402
from compass_os.preprocessing import reference_fill_log2, reference_fill_tpm  # noqa: E402

DATA = REPO / "examples" / "example_data"


def main() -> int:
    expr = pd.read_csv(DATA / "example_expression_partial.tsv.gz", sep="\t", index_col=0)
    lab = pd.read_csv(DATA / "example_labels.tsv", sep="\t", index_col=0)
    ct = lab["cancer_type"].tolist()
    print(f"部分表达矩阵：{expr.shape[0]} 例 × {expr.shape[1]} 基因（词表需 15,672）")
    print(f"冻结参考填充值：log2(TPM+1) = {reference_fill_log2():.12f} ⟺ TPM = "
          f"{reference_fill_tpm():.4f}（**不使用**本队列自身的中位数）\n")

    for strategy in ("reference", "zero"):
        res = compass_os.predict(expr, ct, missing_gene_strategy=strategy)
        d = res.qc.to_dict()
        print(f"[{strategy:9s}] 覆盖 {d['n_observed_genes']}/{d['n_required_genes']}"
              f"（{d['gene_coverage']:.2%}）｜缺失 {d['n_missing_genes']}｜"
              f"默认风险 = {res.default_risk.round(4).tolist()}")

    try:
        compass_os.predict(expr, ct, missing_gene_strategy="strict")
    except MissingGenesError as exc:
        print(f"[strict   ] 按预期拒绝：{exc}")
        print(f"            缺失基因示例（前 5）：{exc.missing_genes[:5]}")

    print("\n--- 稳健性检查（连续指标；阈值尚未校准）---")
    r = compass_os.check_robustness(expr, ct)
    print(f"gene_coverage = {r.gene_coverage:.4f}｜robustness_flag = {r.robustness_flag}")
    print(f"风险差（reference − zero）：{r.risk_difference.round(4).tolist()}")
    print(f"concept 逐样本相关：{ {k: round(v, 5) for k, v in
                                 r.concept_correlation['reference_vs_zero'].items()} }")
    print(f"signature 逐样本相关：{ {k: round(v, 5) for k, v in
                                   r.signature_correlation['reference_vs_zero'].items()} }")
    print("\n注：本示例数据缺失 ~75% 基因，仅用于演示接口与 QC；"
          "真实可用区间需以 missing-gene stress test 结果为准。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
