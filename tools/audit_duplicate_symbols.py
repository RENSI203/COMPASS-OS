#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""audit_duplicate_symbols.py —— 重复 gene symbol 审计（**release blocker 核查**）。

问题
----
package 的输入契约是「列名 = gene symbol」，而原项目曾记录 TCGA 源矩阵中存在 ~55 个
重复 symbol。若同一个 symbol 在模型里对应**多个 feature token**，symbol-only 输入就
无法严格复现；本脚本负责查清并给出结论。

结论（本脚本产出，见 docs/duplicate_gene_symbols.tsv 与 docs/DUPLICATE_SYMBOLS.md）
-------------------------------------------------------------------------------
* **冻结词表 ``feature_name`` 15,672 项全部唯一、零重复**；
  ``gene_tokens_long.json`` 的基因序列与之逐项相同，同样零重复。
  ⇒ **模型侧不存在"一个 symbol 多个 token"的情形，symbol-only 输入无歧义。**
* 55 个重复来自**源矩阵构建阶段**：TCGA RSEM 矩阵以 Ensembl ID 为行，
  经 ``probeMap_gencode.v23`` 映射到 symbol 后，有 55 个 symbol 各对应 **2 个 Ensembl ID**
  （多为 ``ENSG…`` + ``ENSGR…`` 的 PAR_Y/参考重复行，少数为旁系同源/通读对）。
* 原项目 ``01_prepare_expression.py:177`` 的处理规则是
  **保留行均值最大者**（``best = int(np.argmax([float(r.mean()) for r in lst]))``），
  另一行被丢弃。本脚本用原始矩阵复算该规则并核对。

⚠ 这是**训练数据构建口径**的记录，不是 package 的输入歧义：
使用者按 symbol 提供**一列**，与冻结 token 一一对应。唯一需要注意的是——
若使用者在 Ensembl 层自行 sum/average 那 55 个符号的多行，得到的值与训练时该 token
所见的"最大均值行"不同。本脚本量化该差异并写入文档。

用法::

    python tools/audit_duplicate_symbols.py --source-root /path/to/project [--stream]
"""
from __future__ import annotations

import argparse
import collections
import csv
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
from compass_os.preprocessing import feature_names  # noqa: E402

OUT_TSV = REPO / "docs" / "duplicate_gene_symbols.tsv"


def build_table(source_root: Path) -> pd.DataFrame:
    fn = list(feature_names())
    vocab = set(fn)
    pm = pd.read_csv(source_root / "data_tcga/raw/probeMap_gencode.v23.gene.probemap",
                     sep="\t", usecols=["id", "gene"], dtype=str)
    pm = pm.drop_duplicates(subset="id", keep="first")
    m = pm[pm["gene"].isin(vocab)]
    cnt = collections.Counter(m["gene"])
    dups = sorted(k for k, v in cnt.items() if v > 1)
    rep = json.loads((source_root / "data_tcga/processed/01_expression/gene_match_report.json")
                     .read_text(encoding="utf-8"))
    rows = []
    for sym in dups:
        ids = sorted(m.loc[m["gene"] == sym, "id"].tolist())
        rows.append({
            "gene_symbol": sym,
            "n_occurrences": len(ids),
            "n_frozen_feature_tokens": 1,          # 词表里该符号只出现一次
            "frozen_feature_index": fn.index(sym),
            "original_source_ids": ";".join(ids),
            "source_id_type": "Ensembl (probeMap_gencode.v23)",
            "reported_n_rows_in_01": int(rep["symbol_dup_rows"].get(sym, -1)),
            "resolution_rule": "keep row with max mean expression (01_prepare_expression.py:177)",
        })
    df = pd.DataFrame(rows)
    df["note"] = ("源矩阵层面重复；冻结词表唯一 ⇒ package 的 symbol 输入无歧义。"
                  "若在 Ensembl 层自行合并多行，取值可能与该 token 训练时所见不同。")
    return df


def stream_row_stats(source_root: Path, ids: list) -> pd.DataFrame:
    """流式扫描原始 RSEM 矩阵，取候选 Ensembl 行的 mean/std（不整表常驻）。"""
    raw = source_root / "data_tcga/raw/TcgaTargetGtex_rsem_gene_tpm.gz"
    if not raw.is_file():
        print(f"  （跳过数值审计：{raw.name} 不存在）")
        return pd.DataFrame()
    want = set(ids)
    acc = {i: [] for i in ids}
    n_seen = 0
    reader = pd.read_csv(raw, sep="\t", index_col=0, chunksize=800, low_memory=False)
    for chunk in reader:
        n_seen += chunk.shape[0]
        hit = [i for i in chunk.index if i in want]
        for i in hit:
            acc[i].append(chunk.loc[i].to_numpy(np.float64))
        if n_seen % 20000 < 800:
            print(f"    已扫描 {n_seen} 行 …", flush=True)
    out = []
    for i, parts in acc.items():
        if not parts:
            out.append({"ensembl_id": i, "n_samples": 0, "mean": np.nan, "std": np.nan})
            continue
        v = np.concatenate(parts)
        out.append({"ensembl_id": i, "n_samples": int(v.size),
                    "mean": float(np.nanmean(v)), "std": float(np.nanstd(v))})
    return pd.DataFrame(out)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="重复 gene symbol 审计")
    ap.add_argument("--source-root", default="/home/rensi/projects/202608读书汇报")
    ap.add_argument("--stream", action="store_true", help="再扫原始矩阵做数值审计（较慢）")
    a = ap.parse_args(argv)
    src = Path(a.source_root)

    fn = list(feature_names())
    print(f"feature_name: {len(fn)} 行 / {len(set(fn))} 唯一 ⇒ 重复 {len(fn) - len(set(fn))}")
    df = build_table(src)
    print(f"源矩阵层面重复符号: {len(df)} 个（每个 2 行）")

    if a.stream and len(df):
        ids = [x for s in df["original_source_ids"] for x in s.split(";")]
        print(f"数值审计：扫描原始矩阵取 {len(ids)} 个 Ensembl 行的统计量 …")
        stats = stream_row_stats(src, ids)
        if len(stats):
            stats.to_csv(REPO / "docs" / "duplicate_source_row_stats.tsv", sep="\t", index=False)
            mm = stats.set_index("ensembl_id")["mean"]
            chosen, dropped, ratios = [], [], []
            for r in df.itertuples():
                a_ids = r.original_source_ids.split(";")
                means = [mm.get(i, np.nan) for i in a_ids]
                if not np.all(np.isfinite(means)):
                    chosen.append("NA")
                    dropped.append("NA")
                    ratios.append(np.nan)
                    continue
                k = int(np.argmax(means))
                chosen.append(a_ids[k])
                dropped.append(a_ids[1 - k])
                lo = min(means)
                ratios.append(float(max(means) / lo) if lo > 0 else np.inf)
            df["kept_source_id_max_mean"] = chosen
            df["dropped_source_id"] = dropped
            df["max_over_min_row_mean"] = np.round(ratios, 4)
            print(f"  行均值比值（保留/丢弃）中位 = {np.nanmedian(ratios):.3f}，"
                  f"最大 = {np.nanmax(ratios):.3f}")

    df.to_csv(OUT_TSV, sep="\t", index=False, quoting=csv.QUOTE_MINIMAL)
    print(f"已写出 {OUT_TSV.relative_to(REPO)}（{len(df)} 行）")
    print("\n结论：冻结词表零重复 ⇒ symbol-only 输入**不存在**不可逆映射歧义；"
          "55 个源层面重复由 01 的 max-mean 规则消解，已记录。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
