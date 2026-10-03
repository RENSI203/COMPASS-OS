#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""missing_gene_stress_test.py —— 缺失基因稳健性压力测试（**阈值由本实验确定**）。

设计（对应 P0 报告 §5）
----------------------
* **主分析用 39 个外部队列**：所有 TCGA 样本都已参与最终 Cox 系数拟合，
  因此 TCGA 内部**没有**"未参与拟合"的样本；外部队列从未参与任何 Cox 参数拟合。
* TCGA val/test 折（若提供）只能作为 **in-domain perturbation stability reference**，
  **不得**表述为 independent model validation。
* mask 层级：0 / 5 / 10 / 20 / 30 / 40 / 50 %；策略：``reference`` / ``zero``；
  重复 ≥20 次；mask 作用于两类基因集合：
  ``all``（全部 15,672）与 ``signature``（132 基因集涉及的 916 个基因）。
* 每个重复内**所有样本使用同一 mask 基因集合**（保证可比）。

输入契约（两个参数，无绝对路径默认）
------------------------------------
``--pheno``        TSV：``sample_id  cohort  cancer_type  os_time_days  os_event``
``--expr-dir``     目录，内含每个队列一个 ``<cohort>.tsv.gz``（样本 × 15,672，
                   **log2(TPM+1)** 尺度，索引 = sample_id），或单个 ``.tsv.gz`` 文件。

用法示例::

    python validation/missing_gene_stress_test.py \
        --pheno validation/inputs/pheno.tsv \
        --expr-dir validation/inputs/harmonized \
        --out validation/output --repeats 20

产物：``stress_raw.tsv``（逐次重复）、``stress_curves.tsv``（coverage × 指标中位）、
``stress_report.md``、``fig_stress_*.png``。

⚠ 本脚本**不设定**任何"80%/90% 安全"阈值；阈值须由本实验产出后写入 ``qc_config.json``。
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

import compass_os  # noqa: E402
from compass_os.survival import median_cutoff  # noqa: E402

MASK_LEVELS = (0.0, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50)
SCOPES = ("all", "signature")
STRATEGIES = ("reference", "zero")


def _signature_genes() -> list:
    from compass_os.qc import gene_sets
    return sorted({g for gl in gene_sets()["gene_list"] for g in gl})


def load_inputs(pheno_path: Path, expr_dir: Path) -> tuple:
    pheno = pd.read_csv(pheno_path, sep="\t")
    need = {"sample_id", "cohort", "cancer_type", "os_time_days", "os_event"}
    miss = need - set(pheno.columns)
    if miss:
        raise SystemExit(f"--pheno 缺列 {sorted(miss)}；需要 {sorted(need)}")
    pheno = pheno.set_index("sample_id")

    if expr_dir.is_file():
        expr = pd.read_csv(expr_dir, sep="\t", index_col=0)
    else:
        parts = []
        for coh, sub in pheno.groupby("cohort"):
            f = expr_dir / f"{coh}.tsv.gz"
            if not f.is_file():
                print(f"  跳过 {coh}：缺 {f}")
                continue
            d = pd.read_csv(f, sep="\t", index_col=0)
            d = d.reindex([i for i in sub.index if i in d.index])
            parts.append(d)
        if not parts:
            raise SystemExit(f"--expr-dir {expr_dir} 下没有可用的 <cohort>.tsv.gz")
        expr = pd.concat(parts)
    common = [i for i in pheno.index if i in set(expr.index)]
    if not common:
        raise SystemExit("表型与表达矩阵没有交集样本（检查 sample_id 命名）")
    return pheno.loc[common], expr.loc[common]


def uno_c(risk: np.ndarray, time: np.ndarray, event: np.ndarray) -> float:
    try:
        from sksurv.metrics import concordance_index_ipcw
        from sksurv.util import Surv
        y = Surv.from_arrays(event.astype(bool), time)
        return float(concordance_index_ipcw(y, y, risk)[0])
    except Exception:
        return float("nan")


def run(args) -> int:
    t0 = time.time()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    pheno, expr = load_inputs(Path(args.pheno), Path(args.expr_dir))
    print(f"输入：{pheno.shape[0]} 例 / {pheno['cohort'].nunique()} 队列 / "
          f"{expr.shape[1]} 基因（log2 尺度）")
    sig_genes = _signature_genes()
    gene_pool = {"all": list(expr.columns), "signature": [g for g in sig_genes
                                                          if g in expr.columns]}

    by_cohort = {c: list(sub.index) for c, sub in pheno.groupby("cohort")}
    t_all = pheno["os_time_days"].to_numpy(float)
    e_all = pheno["os_event"].to_numpy(float) >= 0.5

    # ---- 完整输入基线 ----
    base = compass_os.predict(expr, pheno["cancer_type"].tolist(),
                              missing_gene_strategy="reference", input_scale="log2_tpm1")
    risk0 = base.default_risk.to_numpy(float)
    con0, sig0 = base.concept_scores, base.signature_scores
    c0_tot = uno_c(risk0, t_all, e_all)
    c0_by = {c: uno_c(risk0[[i for i, s in enumerate(pheno.index) if s in set(by_cohort[c])]],
                      pheno.loc[by_cohort[c], "os_time_days"].to_numpy(float),
                      pheno.loc[by_cohort[c], "os_event"].to_numpy(float) >= 0.5)
             for c in by_cohort}
    cut = median_cutoff(base.default_model)
    grp0 = (risk0 >= cut).astype(int)
    print(f"基线：Uno C = {c0_tot:.4f}（分子队列中位 "
          f"{np.nanmedian(list(c0_by.values())):.4f}）")

    rows = []
    for level in MASK_LEVELS:
        for scope in SCOPES:
            pool = gene_pool[scope]
            n_mask = int(round(level * len(pool)))
            for rep in range(args.repeats if level > 0 else 1):
                idx = (rng.choice(len(pool), size=n_mask, replace=False)
                       if n_mask else np.array([], dtype=int))
                masked = [pool[i] for i in idx]
                sub = expr.drop(columns=masked) if masked else expr
                cov = 1.0 - len(masked) / expr.shape[1]
                for strategy in STRATEGIES:
                    try:
                        res = compass_os.predict(sub, pheno["cancer_type"].tolist(),
                                                 missing_gene_strategy=strategy,
                                                 input_scale="log2_tpm1")
                    except Exception as exc:  # noqa: BLE001
                        rows.append({"mask_level": level, "scope": scope, "repeat": rep,
                                     "strategy": strategy, "gene_coverage": cov,
                                     "error": f"{type(exc).__name__}: {exc}"})
                        continue
                    r = res.default_risk.to_numpy(float)
                    cc = _row_stats(con0, res.concept_scores)
                    ss = _row_stats(sig0, res.signature_scores)
                    c_by = {c: uno_c(
                        r[[i for i, s in enumerate(pheno.index) if s in set(by_cohort[c])]],
                        pheno.loc[by_cohort[c], "os_time_days"].to_numpy(float),
                        pheno.loc[by_cohort[c], "os_event"].to_numpy(float) >= 0.5)
                        for c in by_cohort}
                    rows.append({
                        "mask_level": level, "scope": scope, "repeat": rep,
                        "strategy": strategy, "gene_coverage": cov,
                        "risk_pearson": _corr(risk0, r, "pearson"),
                        "risk_spearman": _corr(risk0, r, "spearman"),
                        "risk_mad": float(np.mean(np.abs(r - risk0))),
                        "risk_rank_spearman": _corr(pd.Series(risk0).rank().to_numpy(),
                                                    pd.Series(r).rank().to_numpy(), "spearman"),
                        "group_agreement": float(np.mean((r >= cut).astype(int) == grp0)),
                        "uno_c": uno_c(r, t_all, e_all),
                        "uno_c_delta": uno_c(r, t_all, e_all) - c0_tot,
                        "uno_c_cohort_median": float(np.nanmedian(list(c_by.values()))),
                        "uno_c_cohort_delta": float(np.nanmedian(list(c_by.values()))
                                                    - np.nanmedian(list(c0_by.values()))),
                        "concept_sample_r": cc["sample_r"], "concept_feature_r": cc["feature_r"],
                        "concept_mad": cc["mad"],
                        "signature_sample_r": ss["sample_r"],
                        "signature_feature_r": ss["feature_r"],
                        "signature_mad": ss["mad"],
                        "n_masked": len(masked), "q": res.qc.gene_coverage,
                    })
            print(f"  level={level:.0%} scope={scope:9s} 完成（{time.time() - t0:.0f}s）")

    raw = pd.DataFrame(rows)
    raw.to_csv(out / "stress_raw.tsv", sep="\t", index=False)

    ok = raw[raw.get("error").isna()] if "error" in raw.columns else raw
    keys = ["risk_pearson", "risk_spearman", "risk_rank_spearman", "group_agreement",
            "uno_c", "uno_c_delta", "concept_sample_r", "concept_feature_r",
            "signature_sample_r", "signature_feature_r", "risk_mad", "concept_mad",
            "signature_mad"]
    curves = (ok.groupby(["scope", "strategy", "gene_coverage"])[keys]
              .median().reset_index())
    curves.to_csv(out / "stress_curves.tsv", sep="\t", index=False)

    _figures(out, curves)
    _report(out, curves, raw, base, c0_tot, pheno, args)
    print(f"\n完成（{time.time() - t0:.0f}s）⇒ {out}")
    return 0


def _corr(a, b, kind: str) -> float:
    from scipy.stats import pearsonr, spearmanr
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 3 or np.std(a[ok]) == 0 or np.std(b[ok]) == 0:
        return float("nan")
    return float((pearsonr if kind == "pearson" else spearmanr)(a[ok], b[ok])[0])


def _row_stats(base: pd.DataFrame, cur: pd.DataFrame) -> dict:
    """逐样本相关 / 逐特征相关 / MAD（先按共同行与列对齐）。"""
    idx = [i for i in base.index if i in set(cur.index)]
    cols = [c for c in base.columns if c in set(cur.columns)]
    A, B = base.loc[idx, cols].to_numpy(float), cur.loc[idx, cols].to_numpy(float)
    sr = [np.corrcoef(A[i], B[i])[0, 1] if np.std(A[i]) > 0 and np.std(B[i]) > 0 else np.nan
          for i in range(A.shape[0])]
    fr = [np.corrcoef(A[:, j], B[:, j])[0, 1] if np.std(A[:, j]) > 0 and np.std(B[:, j]) > 0
          else np.nan for j in range(A.shape[1])]
    return {"sample_r": float(np.nanmedian(sr)), "feature_r": float(np.nanmedian(fr)),
            "mad": float(np.mean(np.abs(A - B)))}


def _figures(out: Path, curves: pd.DataFrame) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    panels = [("risk_pearson", "risk stability (Pearson r)"),
              ("concept_sample_r", "43-concept stability (sample-wise r)"),
              ("signature_sample_r", "132-signature stability (sample-wise r)")]
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.4), dpi=160)
    for ax, (col, lab) in zip(axes, panels):
        for (scope, strategy), sub in curves.groupby(["scope", "strategy"]):
            s = sub.sort_values("gene_coverage")
            ax.plot(s["gene_coverage"], s[col], marker="o", ms=3, lw=1.2,
                    label=f"{scope}/{strategy}")
        ax.set_xlabel("gene coverage", fontsize=8)
        ax.set_ylabel(lab, fontsize=8)
        ax.tick_params(labelsize=8)
        ax.legend(fontsize=6, frameon=False)
    fig.suptitle("Missing-gene stress test (no thresholds drawn: pending calibration)",
                 fontsize=9)
    fig.tight_layout()
    fig.savefig(out / "fig_stress_curves.png")
    fig, ax = plt.subplots(figsize=(5.2, 3.4), dpi=160)
    for (scope, strategy), sub in curves.groupby(["scope", "strategy"]):
        s = sub.sort_values("gene_coverage")
        ax.plot(s["gene_coverage"], s["uno_c"], marker="o", ms=3, lw=1.2,
                label=f"{scope}/{strategy}")
    ax.set_xlabel("gene coverage", fontsize=8)
    ax.set_ylabel("Uno C (pooled external)", fontsize=8)
    ax.tick_params(labelsize=8)
    ax.legend(fontsize=6, frameon=False)
    fig.tight_layout()
    fig.savefig(out / "fig_stress_uno_c.png")


def _report(out: Path, curves: pd.DataFrame, raw: pd.DataFrame, base, c0: float,
            pheno: pd.DataFrame, args) -> None:
    L = ["# Missing-gene stress test 报告（**阈值待定**）", "",
         f"- 生成时间：{time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}",
         f"- 样本：{pheno.shape[0]} 例 / {pheno['cohort'].nunique()} 队列",
         f"- 默认模型：{base.default_model}；mask 层级：{list(MASK_LEVELS)}；"
         f"重复：{args.repeats}；mask 集合：{list(SCOPES)}；策略：{list(STRATEGIES)}",
         f"- 基线 Uno C（池化）= **{c0:.4f}**", "",
         "> ⚠ 本报告**不给出** recommended/warning/unreliable 阈值。",
         "> 阈值须由本实验数据确定后写入 `qc_config.json`（带版本号）。", "",
         "## 曲线（中位，逐 coverage）", "",
         "| scope | strategy | coverage | risk_r | risk_rank_r | group_agree | "
         "UnoC | ΔUnoC | concept_r | signature_r |", "| --- | --- | --- | --- | --- | "
         "--- | --- | --- | --- | --- |"]
    for r in curves.sort_values(["scope", "strategy", "gene_coverage"]).itertuples():
        L.append(f"| {r.scope} | {r.strategy} | {r.gene_coverage:.3f} | "
                 f"{r.risk_pearson:.4f} | {r.risk_rank_spearman:.4f} | "
                 f"{r.group_agreement:.4f} | {r.uno_c:.4f} | {r.uno_c_delta:+.4f} | "
                 f"{r.concept_sample_r:.4f} | {r.signature_sample_r:.4f} |")
    L += ["", "## 读数提示（写讨论时必须遵守）", "",
          "1. 本实验只支持**在测试过的覆盖区间内**描述稳健性；不得写成"
          "“COMPASS 支持任意不完整转录组”或“对缺失基因不变”。",
          "2. `zero` 是**训练空间掩码**：与训练增强**同一归一化输入空间**，"
          "但训练期只作用于正样本视图且概率 0.41，**不是**同一分布。",
          "3. 若某 coverage 下 `reference` 与 `zero` 曲线分离明显，"
          "说明结论对缺失填补策略敏感，必须在 README/论文中同时报告两者。", "",
          "## 产物", "", "| 文件 | 内容 |", "| --- | --- |",
          "| `stress_raw.tsv` | 逐次重复的全部指标 |",
          "| `stress_curves.tsv` | coverage × 指标中位（画图数据） |",
          "| `fig_stress_curves.png` | risk / 43-concept / 132-signature 稳健性曲线 |",
          "| `fig_stress_uno_c.png` | 外部队列 Uno C vs coverage |", ""]
    (out / "stress_report.md").write_text("\n".join(L) + "\n", encoding="utf-8")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="缺失基因稳健性压力测试")
    ap.add_argument("--pheno", required=True, help="TSV：sample_id/cohort/cancer_type/"
                                                  "os_time_days/os_event")
    ap.add_argument("--expr-dir", required=True,
                    help="<cohort>.tsv.gz 目录（log2(TPM+1)，样本 × 15,672）")
    ap.add_argument("--out", default="validation/output")
    ap.add_argument("--repeats", type=int, default=20)
    ap.add_argument("--seed", type=int, default=20261003)
    return run(ap.parse_args(argv))


if __name__ == "__main__":
    sys.exit(main())
