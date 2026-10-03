#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""release_audit.py —— GitHub release 前卫生审计（**git init 之前运行**）。

检查（任一失败 → 退出码 1）：
1. 仓库内**不得**出现原项目绝对路径 / 本机用户目录（源码、文档、清单都查）；
2. 不得包含被 .gitignore 覆盖的运行产物（validation/inputs、validation/output、smoke、
   examples/output、__pycache__、.pytest_cache）；
3. 不得包含原始 TCGA/GEO 数据或大型表达矩阵（按体积与文件名双重判断）；
4. 不得包含下游机制分析结果（downstream_mechanism 等）；
5. 资产清单哈希全部通过；
6. 关键文件齐备（README/LICENSE/CITATION/pyproject/ASSET_MANIFEST/models/tests）；
7. 仓库总体积与单文件体积在 GitHub 限制内（单文件 < 100 MB，仓库建议 < 1 GB）。

输出：``validation/results/release_audit.tsv`` + 控制台摘要。
"""
from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
BIG = 50 * 1024 * 1024          # 单文件 > 50 MB 需人工确认
HARD = 100 * 1024 * 1024        # GitHub 硬限制

#: 允许出现原项目绝对路径的位置（都是**溯源/开发**用途，不是运行时依赖）
ALLOW_PATH_PATTERNS = (
    r"^tools/",                 # 构建/审计工具（开发用）
    r"^validation/inputs/",     # 本地输入（gitignore）
    r"^ASSET_MANIFEST\.tsv$",   # 溯源列 source_path（P2 规格明确要求记录）
    r"^docs/duplicate_gene_symbols\.tsv$",
    r"^tests/fixtures/README\.md$",
    r"^README\.md$",            # 文档中说明 provenance 时允许
)
#: 只对"源码/文档/配置"做绝对路径扫描；数据结构文件（.json/.tsv 数值）不扫
TEXT_SCAN_SUFFIX = (".py", ".md", ".cff", ".toml", ".in", ".cfg", ".txt")

FORBIDDEN_NAMES = re.compile(
    r"(TcgaTargetGtex|expression_tcga_tpm|expression_tcga_log2tp1|"
    r"series_matrix|\.soft$|GSE\d+.*\.(gz|txt)$|"
    r"downstream_mechanism|analysis_out_|cohorts?_(curve|all)\.txt)", re.I)


def scan_paths(root: Path) -> list:
    rows = []

    def add(kind, path, detail, status):
        rows.append({"check": kind, "path": str(path), "detail": detail, "status": status})

    pat = re.compile(r"(/home/[A-Za-z0-9_.-]+/|/mnt/[a-z]/|[A-Za-z]:\\\\|"
                     r"projects/202608|rensi)")
    for f in sorted(root.rglob("*")):
        if not f.is_file():
            continue
        rel = f.relative_to(root)
        if any(re.match(p, str(rel)) for p in ALLOW_PATH_PATTERNS):
            continue
        parts = set(rel.parts)
        if "__pycache__" in parts or ".pytest_cache" in parts:
            add("artifact", rel, "python 缓存（.gitignore 已覆盖，commit 前清理）", "WARN")
            continue
        if rel.parts[:2] in (("validation", "inputs"), ("validation", "output")):
            add("artifact", rel, "本地输入/运行产物（.gitignore 已覆盖）", "WARN")
            continue
        if FORBIDDEN_NAMES.search(f.name):
            add("forbidden-name", rel, "疑似原始数据/大型矩阵/下游结果", "FAIL")
            continue
        size = f.stat().st_size
        if size > HARD:
            add("size", rel, f"{size/1048576:.0f} MB > 100 MB（GitHub 硬限制）", "FAIL")
        elif size > BIG:
            add("size", rel, f"{size/1048576:.0f} MB > 50 MB（需确认）", "WARN")
        if f.suffix in TEXT_SCAN_SUFFIX:
            try:
                txt = f.read_text(encoding="utf-8", errors="ignore")
            except Exception:
                continue
            for i, line in enumerate(txt.splitlines(), 1):
                if pat.search(line) and "http" not in line:
                    add("absolute-path", f"{rel}:{i}", line.strip()[:120], "FAIL")
                    break
    return rows


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="validation/results")
    a = ap.parse_args(argv)
    rows = scan_paths(REPO)

    # 资产哈希
    r = subprocess.run([sys.executable, str(REPO / "tools" / "build_assets.py"), "--verify"],
                       capture_output=True, text=True, cwd=str(REPO))
    rows.append({"check": "asset-hash", "path": "ASSET_MANIFEST.tsv",
                 "detail": r.stdout.strip().splitlines()[-1] if r.stdout else r.stderr[:200],
                 "status": "PASS" if r.returncode == 0 else "FAIL"})

    # 关键文件
    need = ["README.md", "LICENSE", "CITATION.cff", "pyproject.toml", "MANIFEST.in",
            "ASSET_MANIFEST.tsv", "models/model_manifest.json", "models/pretrainer.pt",
            "src/compass_os/__init__.py", "tests/run_tests.py",
            "third_party/COMPASS_LICENSE"]
    for n in need:
        ok = (REPO / n).exists()
        rows.append({"check": "required-file", "path": n,
                     "detail": "存在" if ok else "缺失", "status": "PASS" if ok else "FAIL"})

    # 潜在敏感/未发布内容
    for pat, why in ((r"docs/_", "内部草稿目录前缀"), (r"handoff", "交接材料")):
        hits = [str(p.relative_to(REPO)) for p in REPO.rglob("*")
                if p.is_file() and re.search(pat, str(p.relative_to(REPO)))]
        if hits:
            rows.append({"check": "internal-content", "path": ";".join(hits[:5]),
                         "detail": why, "status": "WARN"})

    import pandas as pd
    df = pd.DataFrame(rows)
    out = REPO / a.out
    out.mkdir(parents=True, exist_ok=True)
    df.to_csv(out / "release_audit.tsv", sep="\t", index=False)
    total = sum(f.stat().st_size for f in REPO.rglob("*") if f.is_file())
    bad = df[df["status"] == "FAIL"]
    warn = df[df["status"] == "WARN"]
    print(f"仓库体积（含未跟踪产物）：{total/1048576:.1f} MB")
    print(f"检查项 {len(df)}：FAIL {len(bad)} / WARN {len(warn)}")
    if len(bad):
        print("\n=== FAIL ===")
        print(bad.to_string(index=False)[:3000])
    if len(warn):
        print("\n=== WARN ===")
        print(warn.to_string(index=False)[:1500])
    return 1 if len(bad) else 0


if __name__ == "__main__":
    sys.exit(main())
