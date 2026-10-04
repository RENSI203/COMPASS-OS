#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""队列示例：多样本预测 → 风险分层 → KM / log-rank → 机制导向表征谱。

运行：
    MPLCONFIGDIR=/tmp/mplcfg python examples/cohort_example.py
产物（写入 ``examples/output/``，示例目录内，可安全删除）：
    cohort_risk.tsv / km_risk_group.png / concept_heatmap.png
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

import compass_os  # noqa: E402
from compass_os import plotting  # noqa: E402
from compass_os.survival import median_cutoff  # noqa: E402

DATA = REPO / "examples" / "example_data"
OUT = REPO / "examples" / "output"


def main() -> int:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    OUT.mkdir(exist_ok=True)
    expr = pd.read_csv(DATA / "example_expression.tsv.gz", sep="\t", index_col=0)
    clin = pd.read_csv(DATA / "example_clinical.tsv", sep="\t", index_col=0)
    lab = pd.read_csv(DATA / "example_labels.tsv", sep="\t", index_col=0)

    res = compass_os.predict(expr, lab["cancer_type"].tolist(), clinical=clin)
    print(f"队列 {expr.shape[0]} 例｜默认模型 {res.default_model}")
    print("风险（含 default）：")
    print(res.risk.round(4).to_string())

    tab = res.risk.copy()
    # risk_group 只在队列 >= min_cohort_for_stratification 时给出（默认 30）；
    # 本示例仅 5 例 ⇒ 必为 None。这里显式处理，避免在示例里抛
    # `TypeError: 'NoneType' object is not subscriptable`。
    if res.risk_group is None:
        tab["risk_group(默认模型)"] = "(未分层：队列 < 30 例)"
        print(f"\n⚠ 队列仅 {len(expr)} 例（< 30），不产生冻结切点 high/low 分层 ⇒ "
              "risk_group=None；已用占位列写入，数值请勿用于分层结论。")
    else:
        tab["risk_group(默认模型)"] = res.risk_group[res.default_model]
    tab["os_time_days"] = lab["os_time_days"]
    tab["os_event"] = lab["os_event"]
    tab.to_csv(OUT / "cohort_risk.tsv", sep="\t")
    print(f"\n已写出 {OUT / 'cohort_risk.tsv'}")

    if len(expr) >= 4 and lab["os_event"].sum() >= 2:
        fig, ax = plt.subplots(figsize=(4.6, 4.0), dpi=160)
        plotting.km_plot(res.default_risk.to_numpy(), lab["os_time_days"].to_numpy(),
                         lab["os_event"].to_numpy(bool), cutoff=median_cutoff(res.default_model),
                         ax=ax, title=f"KM by risk group (frozen cutoff, {res.default_model})")
        fig.tight_layout()
        fig.savefig(OUT / "km_risk_group.png")
        plt.close(fig)
        print(f"已写出 {OUT / 'km_risk_group.png'}")
    else:
        print("样本/事件过少，跳过 KM（仅示例数据）")

    fig, ax = plt.subplots(figsize=(6.0, 2.4), dpi=160)
    plotting.representation_profile(res.concept_scores, top=43, ax=ax,
                                    title="43 concept scores (per sample; latent coordinates)")
    fig.tight_layout()
    fig.savefig(OUT / "concept_heatmap.png")
    plt.close(fig)
    print(f"已写出 {OUT / 'concept_heatmap.png'}")

    print("\n风险秩（队列内）：")
    print(res.risk_rank.round(1).to_string())
    print(f"\n⚠ 分层切点 = 冻结的 median_cutoff（训练集中位风险），不是本队列中位数。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
