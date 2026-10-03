#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""summarize_validation.py —— 汇总 stress test 结果：cohort-level / pan-cohort / 曲线 / 图 / 阈值。

输入：``validation/results/raw*.tsv``（各分片逐条件结果）+ ``platform_masks.json``
输出：
    ``stress_cohort_summary.tsv``     cohort × mask_type × level × strategy（中位）
    ``stress_pancohort_summary.tsv``  pan-cohort：以 **cohort 为独立单位** 的 median/IQR/bootstrap CI
    ``coverage_qc_validation.tsv``    signature/concept 覆盖度 vs 不稳定性
    ``qc_config.json``                由数据确定的 QC 阈值（**不预设**）
    ``figures/fig_{A..E}*.{pdf,png}`` + 对应 source TSV

原则
----
* pan-cohort **不把患者混在一起**：先算每个 cohort 的指标，再对 cohort 做汇总
  （median / IQR / bootstrap 95% CI）。
* 若曲线无断点，不制造"科学阈值"，改用 *empirically calibrated descriptive tiers*，
  并在 qc_config 中标注 ``threshold_kind``。
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

KEYS = ["risk_pearson", "risk_spearman", "risk_rank_spearman", "group_agreement",
        "risk_mad", "risk_delta_median", "risk_delta_p95", "delta_uno",
        "concept_sample_r_median", "concept_sample_r_q05", "concept_feature_r_median",
        "concept_mad", "concept_rmse",
        "signature_sample_r_median", "signature_sample_r_q05", "signature_feature_r_median",
        "signature_mad", "signature_rmse"]


def load_raw(out: Path) -> pd.DataFrame:
    files = sorted(out.glob("raw*.tsv"))
    if not files:
        raise SystemExit(f"{out} 下没有 raw*.tsv")
    df = pd.concat([pd.read_csv(f, sep="\t") for f in files], ignore_index=True)
    if "error" in df.columns:
        bad = df[df["error"].notna()]
        if len(bad):
            print(f"  ⚠ {len(bad)} 条条件失败（前 3 条）：")
            print(bad[["mask_type", "level", "strategy", "cohort", "error"]].head(3).to_string(index=False))
        df = df[df["error"].isna()].copy()
    df = df.drop_duplicates(subset=["mask_type", "level", "repeat", "strategy", "cohort"],
                            keep="last")
    return df


def boot_ci(x, n=2000, seed=7) -> tuple:
    x = np.asarray([v for v in x if np.isfinite(v)], float)
    if x.size < 3:
        return (np.nan, np.nan)
    rng = np.random.default_rng(seed)
    med = np.median(rng.choice(x, size=(n, x.size), replace=True), axis=1)
    return (float(np.percentile(med, 2.5)), float(np.percentile(med, 97.5)))


def summarize(df: pd.DataFrame, out: Path) -> tuple:
    # ---- cohort × condition（对 repeats 取中位）----
    g = (df.groupby(["mask_type", "level", "strategy", "cohort"])[KEYS]
         .median().reset_index())
    cov = (df.groupby(["mask_type", "level", "strategy", "cohort"])["gene_coverage"]
           .median().reset_index())
    coh = g.merge(cov, on=["mask_type", "level", "strategy", "cohort"])
    coh.to_csv(out / "stress_cohort_summary.tsv", sep="\t", index=False)

    # ---- pan-cohort：以 cohort 为独立单位 ----
    rows = []
    for (mt, lv, st), sub in coh.groupby(["mask_type", "level", "strategy"]):
        rec = {"mask_type": mt, "level": lv, "strategy": st,
               "n_cohorts": int(sub["cohort"].nunique()),
               "gene_coverage": float(sub["gene_coverage"].median())}
        for k in KEYS:
            v = sub[k].to_numpy(float)
            lo, hi = boot_ci(v)
            rec[f"{k}_median"] = float(np.nanmedian(v))
            rec[f"{k}_q25"] = float(np.nanpercentile(v, 25))
            rec[f"{k}_q75"] = float(np.nanpercentile(v, 75))
            rec[f"{k}_ci_lo"] = lo
            rec[f"{k}_ci_hi"] = hi
        rows.append(rec)
    pan = pd.DataFrame(rows).sort_values(["mask_type", "strategy", "level"])
    pan.to_csv(out / "stress_pancohort_summary.tsv", sep="\t", index=False)
    return coh, pan


def coverage_qc_validation(df: pd.DataFrame, out: Path) -> pd.DataFrame:
    """覆盖度 QC 是否有预测力：逐 (cohort, condition) 的**逐特征覆盖度 vs 逐特征误差**。"""
    rows = []
    for f in sorted(out.glob("_tmp_sigdelta_s*.npy")) + sorted(out.glob("_tmp_sigcov_s*.npy")):
        pass
    # 逐 signature / concept：把 (覆盖度, |Δ|) 在各 cohort × 条件下汇总
    sig_cov, sig_del = _stack(out, "_tmp_sigcov_s", "_tmp_sigdelta_s")
    con_cov, con_del = _stack(out, "_tmp_concov_s", "_tmp_condelta_s")
    from scipy.stats import spearmanr
    if sig_cov is not None:
        for j in range(sig_cov.shape[1]):
            c, d = sig_cov[:, j], sig_del[:, j]
            ok = np.isfinite(c) & np.isfinite(d)
            rows.append({"level": "signature", "index": j, "n": int(ok.sum()),
                         "coverage_mean": float(np.nanmean(c)),
                         "delta_mean": float(np.nanmean(d)),
                         "spearman_coverage_vs_delta":
                             float(spearmanr(c[ok], d[ok])[0]) if ok.sum() > 10 else np.nan})
    if con_cov is not None:
        for j in range(con_cov.shape[1]):
            c, d = con_cov[:, j], con_del[:, j]
            ok = np.isfinite(c) & np.isfinite(d)
            rows.append({"level": "concept", "index": j, "n": int(ok.sum()),
                         "coverage_mean": float(np.nanmean(c)),
                         "delta_mean": float(np.nanmean(d)),
                         "spearman_coverage_vs_delta":
                             float(spearmanr(c[ok], d[ok])[0]) if ok.sum() > 10 else np.nan})
    tab = pd.DataFrame(rows)
    if not len(tab):
        prev = out / "coverage_qc_validation.tsv"
        if prev.is_file():
            print("  （无 _tmp 中间件：沿用既有 coverage_qc_validation.tsv）")
            return pd.read_csv(prev, sep="\t")
    tab.to_csv(out / "coverage_qc_validation.tsv", sep="\t", index=False)
    # 汇总：跨 feature 的 coverage–instability 关系
    if len(tab):
        from compass_os.qc import gene_sets
        names = gene_sets()["gene_set"].tolist()
        tab["name"] = [names[i] if i < len(names) else str(i) for i in tab["index"]]
        tab.to_csv(out / "coverage_qc_validation.tsv", sep="\t", index=False)
    return tab


def _stack(out: Path, cov_pat: str, del_pat: str):
    """读回逐条件的 (coverage, |delta|) 中间件。

    ⚠ 中间件是**每个 (shard, cohort) 最后一次条件**的覆盖；若已被清理（发布前清理中间产物），
    则 ``_stack`` 返回 (None, None)，此时 ``coverage_qc_validation`` 会**沿用既有结果表**
    而不是把它清空。
    """
    import glob
    covs = sorted(glob.glob(str(out / f"{cov_pat}*.npy")))
    dels = sorted(glob.glob(str(out / f"{del_pat}*.npy")))
    if not covs or not dels:
        return None, None
    C = np.vstack([np.load(f)[None, :] for f in covs])
    D = np.vstack([np.load(f)[None, :] for f in dels])
    n = min(C.shape[0], D.shape[0])
    return C[:n], D[:n]


#: 每种 masking 的**阈值轴**定义（阈值必须报在这个轴上，否则就是"口径混用"）：
#:   global  → 保留的全部基因比例（gene_coverage）
#:   signature → Gsig 覆盖度 = 1 − level（只 mask 916 基因）
#:   non_sig → non-Gsig 覆盖度 = 1 − level（916 基因恒为 100%）
AXIS = {
    "global": ("gene_coverage", "global gene coverage"),
    "signature": ("axis_coverage", "signature (Gsig) coverage"),
    "non_sig": ("axis_coverage", "non-signature gene coverage"),
}


def add_axis(pan: pd.DataFrame) -> pd.DataFrame:
    """给汇总表加统一的 ``axis_coverage``（各 mask 类型自己的覆盖度轴）。"""
    pan = pan.copy()
    pan["axis_coverage"] = np.where(pan["mask_type"] == "global",
                                    pan["gene_coverage"], 1.0 - pan["level"])
    pan["axis_name"] = pan["mask_type"].map(lambda m: AXIS.get(m, ("", ""))[1])
    return pan


def threshold_audit(pan: pd.DataFrame, out: Path) -> pd.DataFrame:
    """逐 coverage level 输出判据满足情况（**不得只给最终数字**）。"""
    rows = []
    for mt in ("global", "signature", "non_sig"):
        sub = pan[(pan["mask_type"] == mt) & (pan["strategy"] == "reference")]
        if not len(sub):
            continue
        axis_col, axis_name = AXIS[mt]
        for r in sub.sort_values("axis_coverage", ascending=False).itertuples():
            sp, ga = float(r.risk_spearman_median), float(r.group_agreement_median)
            rows.append({
                "mask_type": mt,
                "coverage_axis": axis_name,
                "coverage": round(float(getattr(r, axis_col)), 6),
                "mask_fraction_of_pool": round(float(r.level), 4),
                "gene_coverage": round(float(r.gene_coverage), 6),
                "risk_spearman": round(sp, 4),
                "group_agreement": round(ga, 4),
                "recommended_condition": bool(sp >= 0.95 and ga >= 0.90),
                "warning_condition": bool(sp >= 0.80 and ga >= 0.70),
                "n_cohorts": int(r.n_cohorts),
            })
    tab = pd.DataFrame(rows)
    tab.to_csv(out / "qc_threshold_audit.tsv", sep="\t", index=False)
    return tab


def figures(pan: pd.DataFrame, coh: pd.DataFrame, out: Path, cov_tab) -> None:
    """Figure A/B/C（global / Gsig-only / non-Gsig-only 同图对比）、D、E。"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    figdir = out / "figures"
    figdir.mkdir(parents=True, exist_ok=True)

    def save(fig, name, source: pd.DataFrame):
        for ext in ("pdf", "png"):
            fig.savefig(figdir / f"{name}.{ext}", bbox_inches="tight")
        source.to_csv(figdir / f"{name}_source.tsv", sep="\t", index=False)
        plt.close(fig)

    sty = {("global", "reference"): ("#1f77b4", "o", "-"),
           ("global", "zero"): ("#1f77b4", "^", "--"),
           ("signature", "reference"): ("#d62728", "s", "-"),
           ("signature", "zero"): ("#d62728", "v", "--"),
           ("non_sig", "reference"): ("#2ca02c", "D", "-"),
           ("non_sig", "zero"): ("#2ca02c", "d", "--")}
    axis_lab = {"global": "global gene coverage",
                "signature": "signature (Gsig) coverage",
                "non_sig": "non-signature gene coverage"}

    def _panel(ax, col, ylab, title):
        for (mt, st), sub in pan.groupby(["mask_type", "strategy"]):
            s = sub.sort_values("axis_coverage")
            c, m, ls = sty.get((mt, st), ("grey", "o", "-"))
            ax.plot(s["axis_coverage"], s[col], marker=m, ms=3.5, lw=1.2, color=c, ls=ls,
                    label=f"{mt}/{st}")
        ax.set_xlabel("coverage on the masking-specific axis", fontsize=8)
        ax.set_ylabel(ylab, fontsize=8)
        ax.set_title(title, fontsize=8)
        ax.tick_params(labelsize=8)
        ax.legend(fontsize=6, frameon=False)

    # Figure A：三臂风险稳定性
    fig, ax = plt.subplots(figsize=(5.6, 3.8), dpi=160)
    _panel(ax, "risk_spearman_median", "risk Spearman (vs full input)",
           "A. Coverage → risk stability (three masking arms)")
    save(fig, "figA_coverage_vs_risk", pan)

    # Figure B：ΔUno C
    fig, ax = plt.subplots(figsize=(5.6, 3.8), dpi=160)
    _panel(ax, "delta_uno_median", "ΔUno C (perturbed − full)",
           "B. Coverage → ΔUno C")
    ax.axhline(0, color="grey", lw=0.6)
    save(fig, "figB_coverage_vs_unoC", pan)

    # Figure C：表征稳定性（两 panel）
    fig, axes = plt.subplots(1, 2, figsize=(10.0, 3.6), dpi=160)
    _panel(axes[0], "concept_sample_r_median_median",
           "43-concept stability (sample r)", "43 concepts")
    _panel(axes[1], "signature_sample_r_median_median",
           "132-signature stability (sample r)", "132 signatures")
    fig.suptitle("C. Coverage → representation stability (three masking arms)", fontsize=9)
    save(fig, "figC_coverage_vs_representation", pan)

    # Figure D：真实平台覆盖分布
    cov = out / "external_coverage_summary.tsv"
    if cov.is_file():
        c = pd.read_csv(cov, sep="\t")
        ok = c[c["status"] == "ok"]
        fig, ax = plt.subplots(figsize=(5.4, 3.8), dpi=160)
        ax.scatter(ok["global_gene_coverage"], ok["Gsig_coverage"], s=22, color="#2ca02c")
        for r in ok.itertuples():
            ax.annotate(r.cohort, (r.global_gene_coverage, r.Gsig_coverage), fontsize=4,
                        xytext=(2, 2), textcoords="offset points")
        ax.set_xlabel("global COMPASS gene coverage", fontsize=8)
        ax.set_ylabel("Gsig (916 signature genes) coverage", fontsize=8)
        ax.set_title("D. Real platform coverage (pre-imputation)", fontsize=8)
        ax.tick_params(labelsize=8)
        save(fig, "figD_external_coverage", ok)

    # Figure E：signature 覆盖度 vs 扰动误差
    if cov_tab is not None and len(cov_tab):
        s = cov_tab[cov_tab["level"] == "signature"]
        if len(s):
            fig, ax = plt.subplots(figsize=(5.2, 3.8), dpi=160)
            ax.scatter(s["coverage_mean"], s["delta_mean"], s=10, color="#9467bd")
            ax.set_xlabel("signature_gene_coverage", fontsize=8)
            ax.set_ylabel("mean |Δ signature score|", fontsize=8)
            ax.set_title("E. Signature coverage vs perturbation error", fontsize=8)
            ax.tick_params(labelsize=8)
            save(fig, "figE_signature_coverage_vs_error", s)


def write_qc_config(out: Path, pan: pd.DataFrame, cov_tab, design: dict) -> dict:
    """由**数据**确定 QC 阈值并冻结（不预设）。无断点时退化为描述性分层。"""
    import datetime as _dt
    cfg = {
        "version": "1.0.0",
        "validation_date": _dt.datetime.now(_dt.timezone.utc).strftime("%Y-%m-%d"),
        "validation_dataset_count": int(design.get("n_loaded", 0)),
        "validation_sample_count": int(sum(design.get("samples", {}).values())),
        "metrics_used": ["risk_spearman", "group_agreement", "risk_pearson",
                         "concept_sample_r", "signature_sample_r", "delta_uno"],
        "cutoffs_used": {"recommended": {"risk_spearman": 0.95, "group_agreement": 0.90},
                         "warning": {"risk_spearman": 0.80, "group_agreement": 0.70}},
        "threshold_kind": "empirically_calibrated_descriptive_tiers",
        "status": "complete",
        "note": ("阈值由 missing-gene stress test 的实测曲线确定，属**工程 QC 阈值**，"
                 "不是生物学边界；曲线若无明显断点，则它们是描述性分层而非‘安全线’。"),
    }

    # 语义：**满足判据的最低 coverage**（安全区下边界）。要求 ≥4 个不同 coverage 点，
    # 否则该轴样本不足 ⇒ 返回 None 并注明（绝不给出由单点外推的"阈值"）。
    MIN_POINTS = 4

    def lowest(mt, st, sp, ga):
        """满足判据的**最低实测覆盖度**（在**该 mask 类型自己的轴**上；不做插值）。"""
        sub = pan[(pan["mask_type"] == mt) & (pan["strategy"] == st)]
        pts = sub["axis_coverage"].nunique()
        if pts < MIN_POINTS:
            return None, (f"{mt} 轴只有 {pts} 个实测 coverage 点（<{MIN_POINTS}），"
                          "不足以确定阈值")
        hit = sub[(sub["risk_spearman_median"] >= sp) & (sub["group_agreement_median"] >= ga)]
        if not len(hit):
            return None, f"{mt} 轴没有任何实测 coverage 满足判据（Spearman≥{sp} 且 agreement≥{ga}）"
        # 只用**实测点**；取满足判据的最低实测覆盖度（不插值、不外推）
        return round(float(hit["axis_coverage"].min()), 2), ""

    notes = []
    for key, mt, sp, ga in (("recommended_global_coverage", "global", 0.95, 0.90),
                            ("warning_global_coverage", "global", 0.80, 0.70),
                            ("recommended_Gsig_coverage", "signature", 0.95, 0.90),
                            ("warning_Gsig_coverage", "signature", 0.80, 0.70)):
        val, why = lowest(mt, "reference", sp, ga)
        cfg[key] = val
        if why:
            notes.append(f"{key}: {why}")
    if notes:
        cfg["notes_missing"] = notes

    sig_pts = pan[pan["mask_type"] == "signature"]["axis_coverage"].nunique()
    if cov_tab is not None and len(cov_tab) and sig_pts >= MIN_POINTS:
        s = cov_tab[cov_tab["level"] == "signature"]
        if len(s) and s["delta_mean"].notna().any():
            floor = float(np.nanpercentile(s["delta_mean"], 10))
            bad = s[s["delta_mean"] > 2 * floor] if floor > 0 else s.iloc[0:0]
            thr = float(bad["coverage_mean"].max()) if len(bad) else None
            if thr is None:
                cfg["signature_coverage_warning_threshold"] = None
                cfg.setdefault("notes_missing", []).append(
                    "signature_coverage_warning_threshold 未校准且**不予发布**：可复核的"
                    "覆盖度–误差表要求 (signature × coverage) 逐条件累积，当前保留的结果表"
                    "中 coverage_mean 对 132 个 signature 恒为 1.0（无覆盖度变异）⇒ 无法从"
                    "中估计阈值。per-signature coverage 继续以连续值报告"
                    "（signature_gene_coverage，n×132）；**不使用**任何不可复现的数值阈值。")
            else:
                cfg["signature_coverage_warning_threshold"] = {
                    "units": "signature_gene_coverage (per signature, 0-1)",
                    "value": round(float(thr), 2),
                    "rule": ("warn for an individual signature when its "
                             "signature_gene_coverage is below this value"),
                    "basis": ("mean perturbation error of that signature exceeded 2x the "
                              "low-error floor (10th percentile of mean |delta| across the "
                              "132 signatures); the threshold is the highest average "
                              "coverage among those high-error signatures"),
                    "precision_note": "data-derived; rounded to 2 decimals",
                }
    def _n(mt):
        s = pan[pan["mask_type"] == mt]
        return int(s["n_cohorts"].max()) if len(s) else 0

    cfg["curve_cohorts"] = _n("global")
    cfg["empirical_cohorts"] = _n("empirical")
    cfg["signature_axis_points"] = int(sig_pts)
    if sig_pts < MIN_POINTS:
        cfg.setdefault("notes_missing", []).append(
            f"signature_coverage_warning_threshold 未定：signature 轴仅 {sig_pts} 个 "
            f"coverage 点（<{MIN_POINTS}）")
    (out / "qc_config.json").write_text(json.dumps(cfg, ensure_ascii=False, indent=2),
                                        encoding="utf-8")
    return cfg


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="validation/results")
    a = ap.parse_args(argv)
    out = Path(a.out)
    df = load_raw(out)
    print(f"载入 {len(df)} 条结果（{df['cohort'].nunique()} 队列，"
          f"{df['mask_type'].nunique()} 类 mask）")
    coh, pan = summarize(df, out)
    pan = add_axis(pan)
    pan.to_csv(out / "stress_pancohort_summary.tsv", sep="\t", index=False)
    audit = threshold_audit(pan, out)
    print(f"阈值审计表：{len(audit)} 行 ⇒ qc_threshold_audit.tsv")
    print(f"cohort 级 {len(coh)} 行；pan-cohort {len(pan)} 行")
    cov_tab = coverage_qc_validation(df, out)
    design = json.loads((out / "design.json").read_text()) if (out / "design.json").is_file() else {}
    cfg = write_qc_config(out, pan, cov_tab, design)
    print("qc_config.json 已写出：")
    print(json.dumps({k: v for k, v in cfg.items()
                      if "coverage" in k or k == "threshold_kind"},
                     ensure_ascii=False, indent=1))
    figures(pan, coh, out, cov_tab)
    print("图与汇总已写出（validation/results/figures/）")
    return 0


if __name__ == "__main__":
    sys.exit(main())
