#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""build_fixtures.py —— 生成 golden 测试夹具（**开发/打包工具，不属于运行时**）。

夹具的"期望值"来自**原项目的生产实现**，不是本包自己算的：

* 43 concepts：直接取原项目 `04_extract_embedding.py` 的产物
  `data_tcga/processed/04_embedding/embeddings_43.tsv.gz`（真实癌种码）；
* 132 gene sets：
  1) 原项目 `41_extract_all_concepts.py` 的产物
     `analysis_out_concepts/tcga_genesets.tsv.gz`（**占位癌种码 5**）；
  2) 用**生产同一调用**（`common.enable_geneset_capture` + `predict()`）在**真实癌种码**下
     重跑一次（本包采用的正是真实码口径），作为默认路径的期望值；
* M0–M3 风险：用原项目**生产函数** `scripts/12_validate_external.py::_build_features`
  组装设计矩阵，再按 `lock` 的 `scaled_features`/`scaler_mean`/`scaler_scale`/`beta` 计算
  （与 `12` 第 283–290 行逐行相同）。

输出（写入 `tests/fixtures/`）：
    golden_expression.tsv.gz   12 例 × 15,672 基因（原始 TPM，用户输入格式）
    golden_clinical.tsv        age / sex / stage（部分缺失去触发冻结参考值填充）
    golden_labels.tsv          cancer_type + os_time_days + os_event
    golden_expected.tsv.gz     期望的 132 / 43 / M0–M3 风险
    README.md                  夹具来源、哈希与重建命令
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
FIX = REPO / "tests" / "fixtures"
N_SAMPLES = 12


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def _load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="生成 golden 夹具（需原项目只读访问）")
    ap.add_argument("--source-root", default="/home/rensi/projects/202608读书汇报")
    ap.add_argument("--n", type=int, default=N_SAMPLES)
    a = ap.parse_args(argv)
    src = Path(a.source_root)
    FIX.mkdir(parents=True, exist_ok=True)

    sys.path.insert(0, str(src / "scripts"))
    sys.path.insert(0, str(src / "COMPASS-main"))
    common = _load_module("prod_common", src / "scripts" / "common.py")

    # ---- 选样：同时存在于 04 概念档案与 41 基因集档案，且标签齐全 ----
    emb = pd.read_csv(src / "data_tcga/processed/04_embedding/embeddings_43.tsv.gz",
                      sep="\t", index_col=0)
    genesets = pd.read_csv(src / "scripts_dev/analysis_out_concepts/tcga_genesets.tsv.gz",
                           sep="\t", index_col=0)
    labels = pd.read_csv(src / "data_tcga/processed/02_labels/labels_clinical.tsv.gz",
                         sep="\t", compression="gzip").set_index("patient")
    cand = [p for p in emb.index if p in set(genesets.index) and p in set(labels.index)]
    cand = cand[: a.n]
    print(f"夹具样本 {len(cand)} 例：{cand[:3]} …")

    # ---- 表达（原始 TPM，15,672 词表列）----
    model = common.load_checkpoint(str(src / "COMPASS-main/paper/checkpoint/latest/pretrainer.pt"),
                                   device="cpu")
    genes = [str(g) for g in model.feature_name]
    tpm = pd.read_csv(src / "data_tcga/processed/01_expression/expression_tcga_tpm.tsv.gz",
                      sep="\t", index_col=0, usecols=["Unnamed: 0"] + genes)
    expr = tpm.loc[cand].copy()
    expr.index.name = "sample_id"
    expr.to_csv(FIX / "golden_expression.tsv.gz", sep="\t", compression="gzip")

    lab = labels.loc[cand, ["cancer_type", "os_time_days", "os_event"]].copy()
    lab.index.name = "sample_id"
    lab.to_csv(FIX / "golden_labels.tsv", sep="\t")

    # ---- 临床：保留真实缺失以触发"冻结参考值填充"分支（每 4 例去掉 1 个字段）----
    clin = labels.loc[cand, ["age", "sex", "stage"]].copy()
    clin.index.name = "sample_id"
    for i in range(0, len(clin), 4):
        clin.iloc[i, i % 3] = np.nan
    if len(clin) > 5:
        clin.iloc[5, :] = np.nan           # 整行缺失 ⇒ 三个字段都走参考值
    clin.to_csv(FIX / "golden_clinical.tsv", sep="\t")

    # ---- 43 concepts（原项目 04 产物，真实癌种码）----
    con_expected = emb.loc[cand].copy()

    # ---- 132 gene sets：真实码（生产调用重跑）----
    codes = common.load_cancer_code(str(src / "COMPASS-main"))
    X = expr.copy()
    X.insert(0, "cancer_code", [int(codes[labels.loc[p, "cancer_type"]]) for p in X.index])
    common.enable_geneset_capture(True)
    dfe, _ = model.predict(X, batch_size=16, num_workers=0)
    gs_real = common.take_geneset_capture(index=list(X.index)).drop(columns=["CANCER"],
                                                                   errors="ignore")
    con_check = dfe.drop(columns=["CANCER"])

    # ---- M0–M3 风险：用生产函数 _build_features ----
    m12 = _load_module("prod_12", src / "scripts" / "12_validate_external.py")
    # ⚠ 生产函数 `_build_features` 期望 **基因 × 样本** 方向（外部 GEO pipeline 的口径），
    # 而本包内部统一是 样本 × 基因；此处按生产口径转置。
    log2_mat = pd.DataFrame(np.log2(expr.to_numpy(np.float64) + 1.0).T,
                            index=expr.columns, columns=expr.index)
    clin_std = pd.DataFrame({"age": clin["age"], "sex": clin["sex"], "stage": clin["stage"]})
    # 原项目 `41` 档案（**占位癌种码 5**）⇒ 用于独立复核 132 层
    gs_code5 = genesets.loc[cand].copy()

    expected = {"concept_scores": con_expected, "signature_scores": gs_real,
                "signature_scores_code5": gs_code5}
    # ⚠ 生产函数 `_build_features` 一次只接受**一个 cohort 癌种**（外部验证口径是"一个队列一个癌种"），
    # 而本夹具是混合癌种 ⇒ 必须**逐患者**用各自的癌种码调用，才是正确的逐样本期望值。
    for mk in ("M0", "M1", "M2", "M3"):
        lock = json.loads((src / f"data_tcga/processed/locked_model/locked_{mk}.json").read_text()
                          ) if mk in ("M2", "M3") else json.loads(
            (src / f"scripts_dev/analysis_out_cluster/ref_lock/locked_{mk}.json").read_text())
        sc = list(lock["scaled_features"])
        risks = {}
        for pid in cand:
            Xd = m12._build_features(con_expected.loc[[pid]], clin_std.loc[[pid]], lock,
                                     log2_mat[[pid]],
                                     src / "data_tcga/processed/locked_model",
                                     labels.loc[pid, "cancer_type"], None)
            Z = Xd.astype(np.float64).copy()
            if sc and lock.get("scaler_mean") is not None:
                Z[sc] = (Z[sc].to_numpy(np.float64) - np.asarray(lock["scaler_mean"], float)) / \
                    np.asarray(lock["scaler_scale"], float)
            risks[pid] = float((Z.to_numpy(np.float64) @ np.asarray(lock["beta"], float)).ravel()[0])
        expected[f"risk_{mk}"] = pd.Series(risks, name=mk).loc[cand]
        if mk == "M0":
            print(f"  {mk} 逐患者设计矩阵 {Xd.shape}（lock 特征 {len(lock['feature_names'])}）")

    # ---- 单调一致性：包内 43 应与 04 一致（真实码），132 与聚合一致 ----
    d = np.abs(con_check[con_expected.columns].to_numpy() - con_expected.to_numpy()).max()
    print(f"  自检：生产重跑的 43 与 04 档案 max|Δ| = {d:.3e}")

    frames = []
    for name, df in expected.items():
        part = pd.DataFrame(df).copy()
        part.columns = [f"{name}::{c}" for c in part.columns]
        frames.append(part)
    wide = pd.concat(frames, axis=1)
    wide.index.name = "sample_id"
    wide.to_csv(FIX / "golden_expected.tsv.gz", sep="\t", compression="gzip")

    # ---- 示例数据（前 5 例；含"完整"与"部分基因"两版，便于演示缺失基因策略）----
    ex = REPO / "examples" / "example_data"
    ex.mkdir(parents=True, exist_ok=True)
    n_ex = min(5, len(cand))
    ex_rows = cand[:n_ex]
    sample = expr.loc[ex_rows]
    sample.to_csv(ex / "example_expression.tsv.gz", sep="\t", compression="gzip")
    partial_cols = list(sample.columns[::4])          # 保留 1/4 基因 ⇒ 演示缺失处理
    sample[partial_cols].to_csv(ex / "example_expression_partial.tsv.gz", sep="\t",
                                compression="gzip")
    lab.loc[ex_rows].to_csv(ex / "example_labels.tsv", sep="\t")
    clin.loc[ex_rows].to_csv(ex / "example_clinical.tsv", sep="\t")
    print(f"示例数据已写出：{n_ex} 例（完整 {sample.shape[1]} 基因；"
          f"部分 {len(partial_cols)} 基因）")

    prov = {
        "n_samples": len(cand),
        "samples": cand,
        "provenance": {
            "expression": "data_tcga/processed/01_expression/expression_tcga_tpm.tsv.gz",
            "concepts_43": "data_tcga/processed/04_embedding/embeddings_43.tsv.gz (真实癌种码)",
            "genesets_132": "生产调用（predict + enable_geneset_capture，真实癌种码）",
            "risk": "scripts/12_validate_external.py::_build_features + lock(β/scaler)",
            "genesets_132_placeholder_code5":
                "scripts_dev/analysis_out_concepts/tcga_genesets.tsv.gz（占位码 5，另测）",
        },
        "files": {},
    }
    for f in sorted(FIX.glob("golden_*")):
        prov["files"][f.name] = {"size": f.stat().st_size, "sha256": sha256(f)}
    (FIX / "golden_provenance.json").write_text(
        json.dumps(prov, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print("夹具已写出：", ", ".join(sorted(prov["files"])))
    return 0


if __name__ == "__main__":
    sys.exit(main())
