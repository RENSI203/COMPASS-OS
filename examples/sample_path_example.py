#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""样本级计算路径图 —— 真实队列 M2/M3 端到端示例。

用法：
    MPLCONFIGDIR=/tmp/mplcfg python examples/sample_path_example.py \
        --cohort GSE39582 --out examples/output/sample_path

对**同一批** 3 个样本（按 **M3** percentile 最接近 10 / 50 / 90 程序化选择）
分别绘制 M2 与 M3，共用同一份 ``Style``，输出 6 组 PNG/PDF/SVG/HTML +
完整 decomposition + selection 审计。
"""
from __future__ import annotations

import argparse
import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

import compass_os  # noqa: E402
from compass_os import sample_path  # noqa: E402
from compass_os.sample_path_render import Style  # noqa: E402

#: 表达矩阵目录。默认取仓库上一级的和声化缓存；不存在时用 --expr-dir 指定。
DEFAULT_EXPR_DIR = REPO.parent / "scripts_dev" / "analysis_out_harmonize" / "expr_log2"

#: M2/M3 **共用**的视觉规格：画布、节点尺寸、列位置、字体规则与颜色尺度一致。
#: 颜色范围按真实分布设定（非每图 min-max），详见 docs/SAMPLE_PATH_INTEGRATION_VALIDATION.md
SHARED_STYLE = Style(
    gene_threshold=0.25,            # 作用于真实 gene-token score 尺度（仅显示参数）
    max_genes=32, max_signatures=26, predictor_budget=16,
    expression_limits=(0.0, 1704.3),        # 队列 TPM 的 p99.5
    gene_score_limits=(-2.0, 2.0),          # 实测 |gene score| max ≈1.98
    signature_limits=(-1.0, 1.0),
    concept_limits=(-1.0, 1.0),
    contribution_limits=(-1.5, 1.5),        # M2 ±0.93 / M3 ±1.50（PC 聚合）
)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cohort", default="GSE39582")
    ap.add_argument("--expr-dir", default=None,
                    help="含 <cohort>.tsv.gz 的表达矩阵目录（默认：仓库同级 "
                         "scripts_dev/analysis_out_harmonize/expr_log2）")
    ap.add_argument("--out", default=str(REPO / "examples" / "output" / "sample_path"))
    a = ap.parse_args()

    expr_dir = Path(a.expr_dir) if a.expr_dir else DEFAULT_EXPR_DIR
    if not (expr_dir / f"{a.cohort}.tsv.gz").is_file():
        raise SystemExit(f"找不到表达矩阵 {expr_dir / (a.cohort + '.tsv.gz')}；"
                         "请用 --expr-dir 指定目录")

    ph = pd.read_csv(REPO / "validation" / "inputs" / "pheno.tsv", sep="\t")
    sub = ph[ph.cohort == a.cohort]
    expr = pd.read_csv(expr_dir / f"{a.cohort}.tsv.gz", sep="\t", index_col=0)
    expr = expr.loc[[i for i in expr.index if i in set(sub.sample_id)]]
    ct = sub.set_index("sample_id").loc[expr.index, "cancer_type"].tolist()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    print(f"cohort={a.cohort}  n={len(expr)}  cancer_types={sorted(set(ct))}")

    # 一次预测同时得到 M2/M3（同一 cohort、同一批样本）
    pred = compass_os.predict(expr, ct, model="M2,M3", input_scale="log2_tpm1")
    r3 = pred.risk["M3"].to_numpy(float)
    pct3 = np.array([100.0 * ((r3 < v).mean() + 0.5 * (r3 == v).mean()) for v in r3])
    idx = list(map(str, expr.index))
    levels = {"low": 10.0, "mid": 50.0, "high": 90.0}
    picked = {k: idx[int(np.argmin(np.abs(pct3 - t)))] for k, t in levels.items()}
    print("selection (by M3 percentile):",
          {k: (v, f"{pct3[idx.index(v)]:.1f}%") for k, v in picked.items()})

    rows = []
    for level, sid in picked.items():
        for model in ("M2", "M3"):
            res = sample_path(expr, ct, sid, model=model, style=SHARED_STYLE,
                              input_scale="log2_tpm1", result=pred)
            files = res.save(out / level / model)
            d = json.loads(files["selection.json"].read_text(encoding="utf-8"))
            rows.append({"level": level, "model": model, "sample_id": sid,
                         "risk": d["risk"], "rank": d["rank"], "n": d["n"],
                         "percentile": d["percentile"], "cutoff": d["cutoff"],
                         "cutoff_percentile": d["cutoff_percentile"],
                         **{f"n_{k}": v for k, v in d["counts"].items()},
                         "decomposition_rows": len(res.decomposition),
                         "eta_diff": res.eta_max_abs_diff,
                         "imputed": d["metadata"]["clinical_imputed"]})
            print(f"  {level}/{model}: {sid} η={d['risk']:.4f} rank={d['rank']}/{d['n']} "
                  f"pct={d['percentile']:.1f}% genes={d['counts']['gene_score']} "
                  f"pred={d['counts']['predictors']}")

    audit = pd.DataFrame(rows)
    audit.to_csv(out / "REAL_M2_M3_selection_audit.tsv", sep="\t", index=False)
    (out / "shared_style.json").write_text(
        json.dumps({k: (list(v) if isinstance(v, tuple) else v)
                    for k, v in SHARED_STYLE.__dict__.items()}, indent=1), encoding="utf-8")
    print(f"\naudit → {out / 'REAL_M2_M3_selection_audit.tsv'}")
    print(f"figures → {out}")
    return 0


if __name__ == "__main__":
    warnings.filterwarnings("ignore")
    sys.exit(main())
