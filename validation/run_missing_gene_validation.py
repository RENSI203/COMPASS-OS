#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""run_missing_gene_validation.py —— 正式缺失基因稳健性验证（三类 missingness）。

三类实验
--------
``global``    在全部 15,672 冻结 feature 中按**每队列同一集合**随机 mask（模拟平台缺失）。
``signature`` 只在 916 个 signature 相关基因中 mask（控制性缺失，暴露小 gene set 的敏感性）。
``empirical`` 直接 **replay 每个队列真实平台的 pre-imputation 观测集合**（不回推、不猜测：
              由 series matrix × platform annotation 现算，见 ``extract_platform_masks.py``）。

比较基准（对每个样本）
----------------------
``FULL``     完整输入预测（perturbation reference）
``REFERENCE`` 上述 mask 下用冻结参考中位数填补
``ZERO``      上述 mask 下用训练空间掩码（标准化空间置 0）
⇒ 所有指标都是 **FULL vs perturbed** 的配对量。

成本口径（**必须如实报告**）
---------------------------
COMPASS 前向实测 ≈ 33.5 ms/样本（CPU, batch=128, 24 线程；加大线程无效）。
总成本 = 前向次数 × 每批样本数 × 33.5 ms，因此本脚本采用：
每队列 **至多 --max-per-cohort 例**、每 level **--repeats** 次重复；
``empirical`` 不受此限（每个队列只需 2 次前向）。全量 9,359 例 × 全网格约需 40+ 小时，
故默认走"计算受限但完整"的设计；扩展只需调大 ``--max-per-cohort`` / ``--cohorts``。

可续跑：结果逐条 append 到 ``raw.tsv``，重跑自动跳过已完成的 (mask_type, level, repeat, strategy)。
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
from compass_os import survival as _csu  # noqa: E402
from compass_os.survival import median_cutoff  # noqa: E402

GLOBAL_LEVELS = (0.05, 0.10, 0.20, 0.30, 0.40, 0.50)
SIG_LEVELS = (0.01, 0.02, 0.05, 0.10, 0.20, 0.30)
#: non-Gsig 对照：只在 G_all − G_sig 中 mask（Gsig 覆盖度恒为 1.0）
NONSIG_LEVELS = (0.05, 0.10, 0.20, 0.30, 0.40, 0.50)
STRATEGIES = ("reference", "zero")
CANCER_ALIAS = {"CRC": "COAD", "COADREAD": "COAD", "COLORECTAL": "COAD"}


# --------------------------------------------------------------------------- #
# 输入
# --------------------------------------------------------------------------- #
def load_pheno(path: Path) -> pd.DataFrame:
    ph = pd.read_csv(path, sep="\t")
    ph["cancer_type"] = ph["cancer_type"].replace(CANCER_ALIAS)
    return ph


def signature_genes() -> list:
    from compass_os.qc import gene_sets
    return sorted({g for gl in gene_sets()["gene_list"] for g in gl})


def pick_cohorts(pheno: pd.DataFrame, names: str | None, max_n: int) -> list:
    if names:
        want = [c.strip() for c in names.split(",") if c.strip()]
        return [c for c in want if c in set(pheno["cohort"])]
    # 每队列最多 max_n 例，按队列样本量降序取前 --n-cohorts
    return sorted(pheno["cohort"].unique())


# --------------------------------------------------------------------------- #
# 指标
# --------------------------------------------------------------------------- #
def _safe_corr(a, b, kind="spearman") -> float:
    from scipy.stats import pearsonr, spearmanr
    a, b = np.asarray(a, float), np.asarray(b, float)
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 3 or np.std(a[ok]) == 0 or np.std(b[ok]) == 0:
        return float("nan")
    f = pearsonr if kind == "pearson" else spearmanr
    return float(f(a[ok], b[ok])[0])


def _row_stats(full: pd.DataFrame, pert: pd.DataFrame) -> dict:
    """逐样本相关 / 逐特征相关 / MAD / RMSE（零方差 → NaN，不填 0）。"""
    idx = [i for i in full.index if i in set(pert.index)]
    cols = [c for c in full.columns if c in set(pert.columns)]
    A, B = full.loc[idx, cols].to_numpy(float), pert.loc[idx, cols].to_numpy(float)
    sr = [_safe_corr(A[i], B[i], "pearson") for i in range(A.shape[0])]
    fr = [_safe_corr(A[:, j], B[:, j], "pearson") for j in range(A.shape[1])]
    d = B - A
    return {"sample_r_median": float(np.nanmedian(sr)),
            "sample_r_q05": float(np.nanpercentile(sr, 5)) if np.isfinite(sr).any() else np.nan,
            "sample_r_q25": float(np.nanpercentile(sr, 25)) if np.isfinite(sr).any() else np.nan,
            "sample_r_q75": float(np.nanpercentile(sr, 75)) if np.isfinite(sr).any() else np.nan,
            "feature_r_median": float(np.nanmedian(fr)),
            "mad": float(np.nanmean(np.abs(d))),
            "rmse": float(np.sqrt(np.nanmean(d ** 2)))}


def _uno(risk, time, event) -> float:
    if int(np.sum(event)) < 10 or len(risk) < 30:
        return float("nan")
    try:
        from sksurv.metrics import concordance_index_ipcw
        from sksurv.util import Surv
        y = Surv.from_arrays(np.asarray(event, bool), np.asarray(time, float))
        return float(concordance_index_ipcw(y, y, np.asarray(risk, float))[0])
    except Exception:
        return float("nan")


def cohort_metrics(full_risk, pert_risk, full_con, pert_con, full_sig, pert_sig,
                   time, event, cut) -> dict:
    """一个 cohort × 一个条件的全部配对指标。"""
    fr, pr = np.asarray(full_risk, float), np.asarray(pert_risk, float)
    d = pr - fr
    out = {
        "risk_pearson": _safe_corr(fr, pr, "pearson"),
        "risk_spearman": _safe_corr(fr, pr, "spearman"),
        "risk_mad": float(np.mean(np.abs(d))),
        "risk_delta_median": float(np.median(np.abs(d))),
        "risk_delta_p95": float(np.percentile(np.abs(d), 95)),
        "risk_rank_spearman": _safe_corr(pd.Series(fr).rank(), pd.Series(pr).rank(), "spearman"),
        "group_agreement": float(np.mean((fr >= cut).astype(int) == (pr >= cut).astype(int))),
        "uno_full": _uno(fr, time, event),
        "uno_perturbed": _uno(pr, time, event),
    }
    out["delta_uno"] = (out["uno_perturbed"] - out["uno_full"]
                        if np.isfinite(out["uno_full"]) and np.isfinite(out["uno_perturbed"])
                        else float("nan"))
    for pre, a, b in (("concept", full_con, pert_con), ("signature", full_sig, pert_sig)):
        st = _row_stats(a, b)
        for k, v in st.items():
            out[f"{pre}_{k}"] = v
    # coverage–instability 关联（逐 signature / 逐 concept）
    out["_sig_delta"] = np.abs(pert_sig.loc[full_sig.index].to_numpy(float)
                               - full_sig.to_numpy(float)).mean(axis=0)
    out["_con_delta"] = np.abs(pert_con.loc[full_con.index].to_numpy(float)
                               - full_con.to_numpy(float)).mean(axis=0)
    return out


# --------------------------------------------------------------------------- #
# 主流程
# --------------------------------------------------------------------------- #
def run(args) -> int:
    t0 = time.time()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    tag = "" if args.nshards == 1 else f"_shard{args.shard}of{args.nshards}"
    raw_path = out / f"raw{tag}.tsv"
    # 溯源侧车文件与 raw **同 tag 配对**：`design.json` 是单文件、每次运行被覆盖，
    # 无法同时描述 12 队列曲线轴与 38 队列经验轴两次运行（见 docs/audit/04）。
    design_path = out / f"raw{tag}.design.json"

    def design_signature() -> dict:
        """**影响数值**的运行参数指纹；变化时必须拒绝续跑。"""
        return {"cohorts": list(args.cohorts.split(",")) if args.cohorts else None,
                "max_per_cohort": int(args.max_per_cohort),
                "repeats": int(args.repeats),
                "seed": int(args.seed),
                "mask_types": str(args.mask_types),
                "nonsig_levels": [float(x) for x in args.nonsig_levels]}

    done = set()
    prev_design = None
    if design_path.is_file():
        prev_design = json.loads(design_path.read_text(encoding="utf-8")).get("params")
    if raw_path.is_file() and not args.restart:
        if prev_design is not None and prev_design != design_signature():
            raise SystemExit(
                f"{raw_path} 已存在，但其溯源参数与本次调用不一致：\n"
                f"  已有: {prev_design}\n  本次: {design_signature()}\n"
                "续跑会把两套定义**混进同一个 raw 文件**。请改用 --restart 重新开始，"
                "或使参数一致（含 --seed / --max-per-cohort / --repeats / --cohorts / "
                "--mask-types / --nonsig-levels）。")
        prev = pd.read_csv(raw_path, sep="\t")
        done = {(r.mask_type, float(r.level), int(r.repeat), r.strategy, str(r.cohort))
                for r in prev.itertuples()}
        print(f"续跑：已有 {len(prev)} 行结果（溯源参数一致）")

    pheno = load_pheno(Path(args.pheno))
    cohorts = pick_cohorts(pheno, args.cohorts, args.max_per_cohort)
    sig_genes = signature_genes()

    # ---- 逐队列载入表达 + FULL 基线 ----
    data = {}
    for coh in cohorts:
        f = Path(args.expr_dir) / f"{coh}.tsv.gz"
        if not f.is_file():
            print(f"  跳过 {coh}：无 {f.name}")
            continue
        sub = pheno[pheno["cohort"] == coh].sort_values("sample_id")
        ids = list(sub["sample_id"])
        mat = pd.read_csv(f, sep="\t", index_col=0)
        ids = [i for i in ids if i in set(mat.index)][: args.max_per_cohort]
        if len(ids) < 20:
            print(f"  跳过 {coh}：可用样本 {len(ids)} < 20")
            continue
        X = mat.loc[ids]
        ct = sub.set_index("sample_id").loc[ids, "cancer_type"].tolist()
        data[coh] = {"X": X, "ct": ct, "ids": ids,
                     "time": sub.set_index("sample_id").loc[ids, "os_time_days"].to_numpy(float),
                     "event": sub.set_index("sample_id").loc[ids, "os_event"].to_numpy(float) >= 0.5}
    if args.torch_threads > 0:
        import torch
        torch.set_num_threads(int(args.torch_threads))
    print(f"载入 {len(data)} 个队列 / {sum(len(d['ids']) for d in data.values())} 例")

    full = {}
    for coh, d in data.items():
        r = compass_os.predict(d["X"], d["ct"], input_scale="log2_tpm1",
                               batch_size=args.batch_size)
        full[coh] = {"risk": r.default_risk.to_numpy(float),
                     "con": r.concept_scores, "sig": r.signature_scores,
                     "cov": r.qc.signature_gene_coverage,
                     "ccov": r.qc.concept_input_coverage}
    print(f"FULL 基线完成（{time.time() - t0:.0f}s）")
    cut = median_cutoff(_csu.default_model())
    design = {"cohorts": cohorts, "n_loaded": len(data),
              "samples": {c: len(d["ids"]) for c, d in data.items()},
              "cut": cut, "max_per_cohort": args.max_per_cohort,
              "repeats": args.repeats,
              # 与 raw 配对的参数指纹（守卫用）与本次运行标识
              "params": design_signature(), "tag": tag,
              "raw_file": raw_path.name}
    json.dump(design, open(design_path, "w"), ensure_ascii=False, indent=1)
    # 兼容既有读者（validation/summarize_validation.py 读 out/"design.json"）：
    # 显式标注这是"最后一次运行"，完整溯源请看配对侧车文件。
    design_last = dict(design); design_last["_note"] = (
        "LAST RUN ONLY — 该文件每次运行被覆盖。多次运行（如 12 队列曲线轴 + 38 队列"
        "经验轴）请查阅与各自 raw 配对的 raw<tag>.design.json。")
    json.dump(design_last, open(out / "design.json", "w"), ensure_ascii=False, indent=1)

    rows = []
    fh = open(raw_path, "a", encoding="utf-8")
    header_written = raw_path.stat().st_size > 0

    def emit(rec: dict):
        nonlocal header_written
        if not header_written:
            fh.write("\t".join(rec.keys()) + "\n")
            header_written = True
        fh.write("\t".join("" if v is None or (isinstance(v, float) and not np.isfinite(v))
                           else str(v) for v in rec.values()) + "\n")
        fh.flush()

    def do_condition(mask_type: str, level: float, repeat: int, strategy: str,
                     masks: dict):
        """masks: cohort → 待 mask 基因列表（队列内统一）。"""
        for coh, d in data.items():
            key = (mask_type, float(level), int(repeat), strategy, str(coh))
            if key in done:
                continue
            drop = masks.get(coh, [])
            Xp = d["X"].drop(columns=drop) if drop else d["X"]
            try:
                r = compass_os.predict(Xp, d["ct"], input_scale="log2_tpm1",
                                       missing_gene_strategy=strategy,
                                       batch_size=args.batch_size)
            except Exception as exc:  # noqa: BLE001
                emit({"mask_type": mask_type, "level": level, "repeat": repeat,
                      "strategy": strategy, "cohort": coh, "n_masked": len(drop),
                      "error": f"{type(exc).__name__}: {exc}"})
                continue
            m = cohort_metrics(full[coh]["risk"], r.default_risk.to_numpy(float),
                               full[coh]["con"], r.concept_scores,
                               full[coh]["sig"], r.signature_scores,
                               d["time"], d["event"], cut)
            sig_cov = r.qc.signature_gene_coverage.to_numpy(float).mean(axis=0)
            np.save(out / f"_tmp_sigdelta_s{args.shard}_{coh}.npy", m.pop("_sig_delta"))
            np.save(out / f"_tmp_condelta_s{args.shard}_{coh}.npy", m.pop("_con_delta"))
            np.save(out / f"_tmp_sigcov_s{args.shard}_{coh}.npy", sig_cov)
            ccov = r.qc.concept_input_coverage.to_numpy(float).mean(axis=0)
            np.save(out / f"_tmp_concov_s{args.shard}_{coh}.npy", ccov)
            emit({"mask_type": mask_type, "level": level, "repeat": repeat,
                  "strategy": strategy, "cohort": coh,
                  "n_samples": len(d["ids"]), "n_masked": len(drop),
                  "gene_coverage": 1.0 - len(drop) / d["X"].shape[1],
                  **{k: v for k, v in m.items()}, "error": ""})

    # ---------------- 条件表（可跨进程分片） ----------------
    pool = list(next(iter(data.values()))["X"].columns)
    #: non-Gsig 对照的基因池（**程序计算**，不硬编码）：G_all − G_sig
    _sigset = set(sig_genes)
    non_sig_genes = [g for g in pool if g not in _sigset]
    print(f"基因池：G_all={len(pool)}｜G_sig={sum(1 for g in pool if g in _sigset)}"
          f"｜G_non_sig={len(non_sig_genes)}", flush=True)
    kinds = [k.strip() for k in str(args.mask_types).split(",") if k.strip()]
    conditions = []
    if "global" in kinds:
        for level in args.global_levels:
            for rep in range(args.repeats if level > 0 else 1):
                conditions.append(("global", float(level), rep))
    if "signature" in kinds:
        for level in args.sig_levels:
            for rep in range(args.repeats):
                conditions.append(("signature", float(level), rep))
    if "non_sig" in kinds:
        for level in args.nonsig_levels:
            for rep in range(args.repeats):
                conditions.append(("non_sig", float(level), rep))
    if "empirical" in kinds and (out / "platform_masks.json").is_file():
        conditions.append(("empirical", 0.0, 0))

    mine = [c for i, c in enumerate(conditions) if i % args.nshards == args.shard]
    print(f"分片 {args.shard}/{args.nshards}：{len(mine)}/{len(conditions)} 个条件", flush=True)

    def cond_rng(mask_type, level, rep):
        """按条件独立播种 ⇒ 同一条件在任何分片下 mask 完全一致（可复现、可分片）。"""
        return np.random.default_rng([args.seed, int(round(level * 10000)), rep,
                                      {"global": 1, "signature": 2, "non_sig": 3, "empirical": 4}[mask_type]])

    def make_masks(mask_type, level, rep):
        rng = cond_rng(mask_type, level, rep)
        if mask_type == "global":
            n = int(round(level * len(pool)))
            idx = rng.choice(len(pool), size=n, replace=False) if n else []
            return {c: [pool[i] for i in idx] for c in data}
        if mask_type == "signature":
            n = max(1, int(round(level * len(sig_genes))))
            idx = rng.choice(len(sig_genes), size=n, replace=False)
            return {c: [sig_genes[i] for i in idx] for c in data}
        if mask_type == "non_sig":
            # 对照：**只动 non-Gsig 基因**，916 个 signature 基因保持 100% 覆盖
            n = max(1, int(round(level * len(non_sig_genes))))
            idx = rng.choice(len(non_sig_genes), size=n, replace=False)
            return {c: [non_sig_genes[i] for i in idx] for c in data}
        pm = json.loads((out / "platform_masks.json").read_text(encoding="utf-8"))["cohorts"]
        return {c: pm.get(c, {}).get("missing_genes", []) for c in data}

    for mask_type, level, rep in mine:
        masks = make_masks(mask_type, level, rep)
        for strategy in STRATEGIES:
            do_condition(mask_type, level, rep, strategy, masks)
        print(f"  {mask_type} {level:.0%} rep{rep} 完成（{time.time() - t0:.0f}s）", flush=True)

    fh.close()
    print(f"\n完成（{time.time() - t0:.0f}s）⇒ {raw_path}")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="正式缺失基因稳健性验证")
    ap.add_argument("--pheno", required=True)
    ap.add_argument("--expr-dir", required=True)
    ap.add_argument("--out", default="validation/results")
    ap.add_argument("--cohorts", default=None, help="逗号分隔；默认全部")
    ap.add_argument("--max-per-cohort", type=int, default=60)
    ap.add_argument("--repeats", type=int, default=20)
    ap.add_argument("--batch-size", type=int, default=128)
    ap.add_argument("--seed", type=int, default=20261003)
    ap.add_argument("--global-levels", type=float, nargs="+", default=list(GLOBAL_LEVELS))
    ap.add_argument("--sig-levels", type=float, nargs="+", default=list(SIG_LEVELS))
    ap.add_argument("--restart", action="store_true")
    ap.add_argument("--torch-threads", type=int, default=0)
    ap.add_argument("--nshards", type=int, default=1)
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nonsig-levels", type=float, nargs="+", default=list(NONSIG_LEVELS))
    ap.add_argument("--mask-types", default="global,signature,non_sig,empirical",
                    help="逗号分隔：global/signature/empirical")
    return run(ap.parse_args(argv))


if __name__ == "__main__":
    sys.exit(main())
