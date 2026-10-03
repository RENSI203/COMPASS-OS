#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""build_assets.py —— 从**只读源项目**派生小体积数据资产并生成哈希清单。

**这是开发/打包工具，不是运行时依赖。**它记录的原项目绝对路径仅作为**溯源信息**
写进 `ASSET_MANIFEST.tsv` 的 `source_path` 列；正式 package 运行时**不读取**任何原项目路径。

用法：
    python tools/build_assets.py --source-root /path/to/project      # 重新派生 + 写清单
    python tools/build_assets.py --verify                            # 只校验现有资产哈希
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
PKG = REPO / "src" / "compass_os"
DATA = PKG / "data"

#: (资产名, 包内相对路径, 原项目相对路径, 类别, 用途, 推理必需, 测试必需, 示例必需, 是否可改)
ITEMS = [
    ("pretrainer_ckpt", "src/compass_os/assets/models/pretrainer.pt",
     "COMPASS-main/paper/checkpoint/latest/pretrainer.pt",
     "frozen_model", "COMPASS 预训练 checkpoint：encoder + 132 基因集 / 43 概念投影 + Datascaler",
     True, True, True, False),
    ("lock_M0", "src/compass_os/assets/models/locked_M0.json",
     "scripts_dev/analysis_out_cluster/ref_lock/locked_M0.json",
     "frozen_model", "M0 锁定 Cox 头（仅临床基线）", True, True, True, False),
    ("lock_M1", "src/compass_os/assets/models/locked_M1.json",
     "scripts_dev/analysis_out_cluster/ref_lock/locked_M1.json",
     "frozen_model", "M1 锁定 Cox 头（仅 43 概念）", True, True, True, False),
    ("lock_M2", "src/compass_os/assets/models/locked_M2.json",
     "data_tcga/processed/locked_model/locked_M2.json",
     "frozen_model", "M2 锁定 Cox 头（主模型：临床 + 43 概念）", True, True, True, False),
    ("lock_M3", "src/compass_os/assets/models/locked_M3.json",
     "data_tcga/processed/locked_model/locked_M3.json",
     "frozen_model", "M3 锁定 Cox 头（M2 + PC1–10，敏感性模型）", True, True, True, False),
    ("pca_M3", "src/compass_os/assets/models/locked_pca_M3.npz",
     "data_tcga/processed/locked_model/locked_pca_M3.npz",
     "frozen_model", "M3 的锁定 PCA（components/mean/scale/genes），禁止在用户队列重拟合",
     True, True, True, False),
    ("ref_quantiles", "src/compass_os/assets/models/reference_quantiles.json",
     "data_tcga/processed/locked_model/reference_quantiles.json",
     "frozen_input", "TCGA log2(TPM+1) 参考分位数（1024 点）：跨平台映射与缺失基因填充中位数来源",
     False, True, False, False),
    ("qc_config", "src/compass_os/assets/models/qc_config.json",
     "validation/results/qc_config.json",
     "frozen_config", "由 missing-gene stress test 实测曲线确定的 QC 阈值（global 轴已定；Gsig 待补）",
     True, True, True, False),
    ("upstream_license", "src/compass_os/assets/third_party/COMPASS_LICENSE",
     "COMPASS-main/LICENSE",
     "license", "上游 COMPASS 的 MIT 许可证（必须随附）", False, False, False, False),
    ("upstream_compass_pkg", "src/compass_os/assets/third_party/compass",
     "COMPASS-main/compass",
     "upstream_code", "上游 COMPASS 包源码 + tokenizer 资产（conception_processed.tsv 等）",
     True, True, True, False),
    ("concept_gene_sets", "src/compass_os/data/concept_gene_sets.tsv",
     "COMPASS-main/compass/tokenizer/conception_processed.tsv",
     "derived_data", "从冻结 CONCEPT 派生的 132 基因集→基因符号表（供 signature coverage 独立使用）",
     False, True, True, False),
    ("gene_vocabulary", "src/compass_os/data/gene_vocabulary.txt",
     "COMPASS-main/paper/checkpoint/latest/pretrainer.pt",
     "derived_data", "从冻结 checkpoint 的 feature_name 导出的 15,672 基因词表（顺序权威）",
     True, True, True, False),
    ("cancer_codes", "src/compass_os/data/cancer_codes.tsv",
     "COMPASS-main/compass/tokenizer/cancer_code.json",
     "derived_data", "从冻结 cancer_code.json 派生的癌种码表（真实 COMPASS code + CT_* 列名）",
     True, True, True, False),
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def sha256_tree(path: Path) -> str:
    """目录哈希：按相对路径排序后逐个更新（文件内容 + 相对路径）。

    跳过 ``__pycache__`` / ``*.pyc``——字节码缓存不是资产内容，导入后必然变化。
    """
    h = hashlib.sha256()
    files = [p for p in path.rglob("*") if p.is_file()
             and "__pycache__" not in p.parts and p.suffix != ".pyc"]
    for f in sorted(files):
        h.update(str(f.relative_to(path)).encode())
        h.update(sha256(f).encode())
    return h.hexdigest()


def build_derived(source_root: Path) -> None:
    """派生两个小数据资产（必须与上游冻结文件逐值一致，由 tests/test_assets.py 复核）。"""
    import pandas as pd
    DATA.mkdir(parents=True, exist_ok=True)

    concept = pd.read_csv(source_root / "COMPASS-main/compass/tokenizer/conception_processed.tsv",
                          sep="\t", index_col=0)
    out = pd.DataFrame({
        "gene_set": concept.index.astype(str),
        "broad_pathway": concept["BroadCelltypePathway"].astype(str),
        "n_genes": concept["Genes"].astype(str).str.split(":").apply(len),
        "genes": concept["Genes"].astype(str),
    }).reset_index(drop=True)
    out.to_csv(DATA / "concept_gene_sets.tsv", sep="\t", index=False)

    import torch
    ckpt = torch.load(source_root / "COMPASS-main/paper/checkpoint/latest/pretrainer.pt",
                      map_location="cpu", weights_only=False)
    genes = [str(g) for g in ckpt.feature_name]
    (DATA / "gene_vocabulary.txt").write_text("\n".join(genes) + "\n", encoding="utf-8")
    del ckpt

    codes = json.loads((source_root / "COMPASS-main/compass/tokenizer/cancer_code.json")
                       .read_text(encoding="utf-8"))
    ct = [c for c in json.loads((REPO / "src/compass_os/assets/models/locked_M2.json").read_text(encoding="utf-8"))
          ["feature_names"] if c.startswith("CT_")]
    rows = [{"cancer_type": k, "compass_code": int(v),
             "ct_column": (f"CT_{k}" if f"CT_{k}" in ct else ""),
             "is_reference_level": bool(v >= 0 and f"CT_{k}" not in ct)}
            for k, v in codes.items()]
    pd.DataFrame(rows).sort_values("compass_code").to_csv(
        DATA / "cancer_codes.tsv", sep="\t", index=False)
    print(f"派生完成：gene_vocabulary.txt（{len(genes)} 基因）、"
          f"concept_gene_sets.tsv（{len(out)} 行）、cancer_codes.tsv（{len(rows)} 行）")


def write_manifest(source_root: Path) -> None:
    rows = []
    for name, rel, src_rel, kind, purpose, inf, test, ex, mod in ITEMS:
        pkg_path = REPO / rel
        if not pkg_path.exists():
            print(f"!! 缺失：{rel}", file=sys.stderr)
            continue
        src_path = source_root / src_rel
        is_dir = pkg_path.is_dir()
        rows.append({
            "asset_name": name,
            "asset_type": kind,
            "source_path": str(src_path),
            "github_copy_path": rel,
            "copied_or_reference_only": "copied",
            "is_directory": "yes" if is_dir else "no",
            "file_size": "" if is_dir else pkg_path.stat().st_size,
            "sha256": sha256_tree(pkg_path) if is_dir else sha256(pkg_path),
            "purpose": purpose,
            "required_for_inference": "yes" if inf else "no",
            "required_for_testing": "yes" if test else "no",
            "required_for_example": "yes" if ex else "no",
            "modifiable": "yes" if mod else "no",
            "notes": "只复制自只读源项目；运行时通过包内相对路径解析" if not is_dir
                     else "目录哈希 = 各文件内容哈希 + 相对路径 的有序拼接",
        })
    df = __import__("pandas").DataFrame(rows)
    df.to_csv(REPO / "ASSET_MANIFEST.tsv", sep="\t", index=False)
    print(f"ASSET_MANIFEST.tsv 已写出：{len(df)} 项")

    model_meta = {
        "schema_version": 1,
        "default_model": "M2",
        "models": {
            "M0": {"label": "clinical baseline",
                   "lock": "src/compass_os/assets/models/locked_M0.json",
                   "features": "CancerType(one-hot) + Age + Sex + Stage",
                   "description": "仅临床变量的基线模型。"},
            "M1": {"label": "COMPASS concept model",
                   "lock": "src/compass_os/assets/models/locked_M1.json",
                   "features": "43 COMPASS concepts",
                   "description": "仅 COMPASS 43 概念，无临床变量。"},
            "M2": {"label": "COMPASS concepts + clinical (primary)",
                   "lock": "src/compass_os/assets/models/locked_M2.json",
                   "features": "CancerType(one-hot) + Age + Sex + Stage + 43 COMPASS concepts",
                   "description": "当前冻结的主模型（frozen primary model）。"},
            "M3": {"label": "COMPASS concepts + clinical + expression PCs",
                   "lock": "src/compass_os/assets/models/locked_M3.json",
                   "pca": "src/compass_os/assets/models/locked_pca_M3.npz",
                   "features": "M2 + PC1-PC10（锁定 PCA，作用于 log2(TPM+1)）",
                   "description": "表达主成分增强模型，作为附加模型提供。"},
        },
        "checkpoint": "src/compass_os/assets/models/pretrainer.pt",
        "reference_quantiles": "src/compass_os/assets/models/reference_quantiles.json",
        "missing_gene": {
            "strategies": ["reference", "zero", "strict"],
            "default": "reference",
            "reference_value_log2": 3.130937933922,
            "reference_value_tpm": 7.76,
            "reference_value_source": "median(reference_quantiles['values'])",
            "zero_space": "MinMax-scaled input space (post-Datascaler)",
        },
    }
    (REPO / "src/compass_os/assets/models" / "model_manifest.json").write_text(
        json.dumps(model_meta, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print("src/compass_os/assets/models/model_manifest.json 已写出")


def verify() -> int:
    import pandas as pd
    df = pd.read_csv(REPO / "ASSET_MANIFEST.tsv", sep="\t")
    bad = 0
    for r in df.itertuples(index=False):
        p = REPO / r.github_copy_path
        if not p.exists():
            print(f"MISSING  {r.github_copy_path}")
            bad += 1
            continue
        got = sha256_tree(p) if p.is_dir() else sha256(p)
        if got != r.sha256:
            print(f"CHANGED  {r.github_copy_path}\n  manifest={r.sha256}\n  actual  ={got}")
            bad += 1
    print(f"校验：{len(df)} 项 ⇒ 异常 {bad}")
    return 1 if bad else 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="派生资产 + 生成哈希清单")
    ap.add_argument("--source-root", default="/home/rensi/projects/202608读书汇报",
                    help="只读源项目根（仅用于 build/verify，运行时不使用）")
    ap.add_argument("--verify", action="store_true")
    a = ap.parse_args(argv)
    if a.verify:
        return verify()
    root = Path(a.source_root)
    build_derived(root)
    write_manifest(root)
    return verify()


if __name__ == "__main__":
    sys.exit(main())
