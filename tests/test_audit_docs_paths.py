# -*- coding: utf-8 -*-
"""审查文档不得包含开发机绝对路径（第 5 阶段新增的机械化防护）。

**为什么需要这个测试**：`release_audit.py` 的 `absolute-path` 检查会扫描
`docs/**`。在审查过程中，"引用证据"时复制机读路径**已经发生过三次**
（Phase 0 §6.5 / Phase 4 §7.1 / Phase 5 §3.1），每次都只补了文字教训，
因此第三次仍然发生。本测试把该纪律变成可执行断言：

* 若某份审查文档引用了一条含开发机路径的证据，**写入时应脱敏**；
* 一旦漏脱敏，本测试立刻失败，而不是等到发布审计才暴露。

原始字符串仍可由 `python tools/release_audit.py` 复现，证据价值不损失。
"""
from __future__ import annotations

import re

from _util import repo_root

#: 与 tools/release_audit.py 的检测口径保持一致（见其 ALLOW_PATH_PATTERNS）
PATTERN = re.compile(r"(/home/[A-Za-z0-9_.-]+/|/mnt/[a-z]/|[A-Za-z]:\\\\|"
                     r"projects/202608|rensi)")


def _audit_files():
    d = repo_root() / "docs" / "audit"
    return sorted(list(d.glob("*.md")) + list(d.glob("*.tsv")))


def test_audit_docs_contain_no_machine_paths():
    """审查文档（含 manifest）不得出现开发机绝对路径或用户名。"""
    files = _audit_files()
    assert files, "docs/audit 下没有审查文档（路径解析可能出错）"
    bad = []
    for f in files:
        for i, line in enumerate(f.read_text(encoding="utf-8", errors="ignore").splitlines(), 1):
            if PATTERN.search(line) and "http" not in line:
                bad.append(f"{f.name}:{i}: {line.strip()[:80]}")
                break
    assert not bad, (
        "docs/audit 下出现开发机绝对路径/用户名，须按 docs/audit/05_MAIN_RESULTS_CLAIMS.md "
        "§3.1 的口径脱敏：\n  " + "\n  ".join(bad))


def test_audit_bundle_lists_every_evidence_file():
    """汇总包 README 必须索引 docs/audit 下的全部证据文件（避免漏交）。"""
    d = repo_root() / "docs" / "audit"
    idx = (d / "README.md").read_text(encoding="utf-8")
    missing = [f.name for f in _audit_files()
               if f.name != "README.md" and f.name not in idx]
    assert not missing, f"汇总包未索引：{missing}"
