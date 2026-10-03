#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""extract_platform_masks.py —— 从原始 GEO series matrix 还原**pre-imputation 观测基因集合**。

原则
----
**只现算，不回推。** 不从已 imputed 的和声化缓存反推缺失模式；而是复刻生产口径
（``32``/``47``）：series matrix 的探针 × 平台注释 ``probe2symbol`` → 基因符号
（同符号取平均表达最大探针）→ 与 15,672 词表求交 ⇒ 得到

* ``observed_genes``：该队列真实观测到的 COMPASS 基因
* ``missing_genes``：需要填补的基因（= 词表 − observed）

同时记录 ``n_observed_Gsig`` / ``Gsig_coverage``（916 个 signature 相关基因的覆盖），
因为全局覆盖相同不代表 Gsig 覆盖相同。

输出：``validation/results/platform_masks.json``（含逐队列掩码与覆盖统计）
      ``validation/results/external_coverage_summary.tsv``
"""
from __future__ import annotations

import argparse
import gzip
import json
import re
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
from compass_os.preprocessing import feature_names  # noqa: E402
from compass_os.qc import gene_sets  # noqa: E402


def _norm_index(x) -> str:
    """探针 ID 归一：去引号/空白；浮点化的整数探针（7892501.0）还原为 7892501。

    ⚠ 踩过的坑：用 ``dtype=np.float64`` 读入时，Affymetrix Gene ST 的**数字探针 ID**
    会被解析成 float，``astype(str)`` 得到 ``"7892501.0"``，与注释文件里的 ``"7892501"``
    对不上 ⇒ 整个平台覆盖率算成 0。GPL6244/GPL5175/GPL17692 三个平台都踩到这里。
    """
    s = str(x).strip().strip('"')
    if s.endswith(".0") and s[:-2].isdigit():
        return s[:-2]
    return s


def read_series_matrix(path: Path) -> pd.DataFrame:
    """series matrix 的表格段 → DataFrame(探针 × GSM)。

    用 pandas C 解析器 + ``comment='!'`` + ``quotechar='"'``：元信息行都以 ``!`` 开头
    （含结束标记），表内单元格带引号，一次 ``read_csv`` 即可。
    数值列交给 pandas 自行推断（数值表直接得到 float64，无需逐列 to_numeric）。
    """
    df = pd.read_csv(path, sep="\t", comment="!", quotechar='"', index_col=0,
                     low_memory=False)
    df.index = [_norm_index(i) for i in df.index]
    if df.dtypes.eq(object).any():
        df = df.apply(pd.to_numeric, errors="coerce")
    else:
        df = df.astype(np.float64, copy=False)
    return df[df.index != ""]


def map_probes(mat: pd.DataFrame, ann: pd.DataFrame) -> pd.DataFrame:
    """探针 → 符号；同一符号多个探针取**平均表达最大**者（与生产同口径）。"""
    # 两侧探针 ID 都过同一个归一器：注释文件若被读成 int64（Affymetrix 数字探针），
    # 直接用会与矩阵的字符串索引完全对不上（实测 GPL6244 → 0 个符号）。
    a = pd.Series(ann[ann.columns[-1]].astype(str).to_numpy(),
                  index=[_norm_index(p) for p in ann[ann.columns[0]]])
    a = a[~a.index.duplicated(keep="first")]
    probes = np.asarray([_norm_index(p) for p in mat.index], dtype=object)
    keep_mask = np.isin(probes, a.index.to_numpy())
    m = mat.loc[keep_mask]
    sym = a.reindex(probes[keep_mask])
    means = m.mean(axis=1).to_numpy(float)
    # 向量化选"同符号中平均表达最大"的探针
    d = pd.DataFrame({"sym": sym.to_numpy(), "mean": means,
                      "probe": np.asarray(probes[keep_mask], dtype=object)})
    d = d.sort_values("mean", ascending=False).drop_duplicates("sym", keep="first")
    out = m.loc[d["probe"].to_numpy()]
    out.index = d["sym"].to_numpy()
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="还原 pre-imputation 平台观测集合")
    ap.add_argument("--source-root", required=True,
                    help="只读源项目根（含 data_geo/、scripts_dev/）；本脚本仅用于 validation")
    ap.add_argument("--pheno", required=True)
    ap.add_argument("--out", default="validation/results")
    a = ap.parse_args(argv)
    src = Path(a.source_root)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    vocab = list(feature_names())
    vs = set(vocab)
    sig = set(g for gl in gene_sets()["gene_list"] for g in gl)
    pheno = pd.read_csv(a.pheno, sep="\t")

    status = pd.read_csv(src / "scripts_dev/analysis_out_geo_validate/cohort_platform_status.tsv",
                         sep="\t")
    plat_of = {r["GEO数据集"]: str(r["平台"]) for _, r in status.iterrows() if pd.notna(r["平台"])}
    ann_dirs = [src / "data_geo/platform_annot_soft", src / "data_geo/platform_annot"]

    masks, cov_rows = {}, []
    for coh in sorted(pheno["cohort"].unique()):
        if coh == "METABRIC":
            cov_rows.append({"cohort": coh, "platform": "METABRIC", "status": "unavailable",
                             "reason": "非 GEO 平台；和声化缓存无 pre-imputation 记录"})
            continue
        cands = sorted((src / "data_geo/GEO_cohorts" / coh).glob("*series_matrix*"))
        if not cands:
            cov_rows.append({"cohort": coh, "platform": plat_of.get(coh, ""),
                             "status": "unavailable", "reason": "未找到 series matrix"})
            continue
        f = cands[0]
        plat = plat_of.get(coh, "")
        mm = re.search(r"(GPL\d+)", f.name)
        if mm:
            plat = mm.group(1)
        ann = None
        for d in ann_dirs:
            p = d / f"{plat}.probe2symbol.tsv"
            if p.is_file():
                ann = pd.read_csv(p, sep="\t")
                break
        if ann is None:
            cov_rows.append({"cohort": coh, "platform": plat, "status": "unavailable",
                             "reason": f"无平台注释 {plat}.probe2symbol.tsv"})
            continue
        try:
            mat = read_series_matrix(f)
            mapped = map_probes(mat, ann)
        except Exception as exc:  # noqa: BLE001
            cov_rows.append({"cohort": coh, "platform": plat, "status": "unavailable",
                             "reason": f"{type(exc).__name__}: {exc}"})
            continue
        obs = [g for g in mapped.index.astype(str) if g in vs]
        obs_sig = [g for g in obs if g in sig]
        miss = [g for g in vocab if g not in set(obs)]
        masks[coh] = {"platform": plat,
                      "n_probes_in_matrix": int(mat.shape[0]),
                      "n_symbols_mapped": int(mapped.shape[0]),
                      "n_observed_COMPASS_genes": len(obs),
                      "n_missing_COMPASS_genes": len(miss),
                      "global_gene_coverage": len(obs) / len(vocab),
                      "n_observed_Gsig": len(obs_sig),
                      "Gsig_coverage": len(obs_sig) / len(sig),
                      "observed_genes": obs,
                      "missing_genes": miss}
        cov_rows.append({"cohort": coh, "platform": plat, "status": "ok",
                         "n_observed_COMPASS_genes": len(obs),
                         "n_missing_COMPASS_genes": len(miss),
                         "global_gene_coverage": len(obs) / len(vocab),
                         "n_observed_Gsig": len(obs_sig),
                         "Gsig_coverage": len(obs_sig) / len(sig),
                         "reason": ""})
        print(f"  {coh:12s} {plat:9s} observed {len(obs):5d}/{len(vocab)} "
              f"({len(obs)/len(vocab):.3f})｜Gsig {len(obs_sig)}/{len(sig)} "
              f"({len(obs_sig)/len(sig):.3f})", flush=True)

    (out / "platform_masks.json").write_text(
        json.dumps({"n_cohorts": len(masks), "cohorts": masks}, ensure_ascii=False),
        encoding="utf-8")
    cov = pd.DataFrame(cov_rows)
    cov.to_csv(out / "external_coverage_summary.tsv", sep="\t", index=False)
    ok = cov[cov["status"] == "ok"]
    print(f"\n完成（{time.time() - t0:.0f}s）：可还原 {len(ok)}/{len(cov)} 个队列")
    if len(ok):
        print(f"  global coverage 中位 {ok['global_gene_coverage'].median():.3f} "
              f"（范围 {ok['global_gene_coverage'].min():.3f}–{ok['global_gene_coverage'].max():.3f}）")
        print(f"  Gsig   coverage 中位 {ok['Gsig_coverage'].median():.3f} "
              f"（范围 {ok['Gsig_coverage'].min():.3f}–{ok['Gsig_coverage'].max():.3f}）")
    return 0


if __name__ == "__main__":
    sys.exit(main())
