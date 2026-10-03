# -*- coding: utf-8 -*-
"""资产完整性测试：哈希、派生数据与冻结源逐值一致、无绝对路径依赖。"""
from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from _util import REPO, fixtures_available, load_fixtures

pytestmark = pytest.mark.skipif(not fixtures_available(), reason="需要 tests/fixtures")


def test_asset_manifest_hashes():
    """ASSET_MANIFEST.tsv 的每一项 sha256 与实际文件一致。"""
    import subprocess
    import sys
    r = subprocess.run([sys.executable, str(REPO / "tools" / "build_assets.py"), "--verify"],
                       capture_output=True, text=True, cwd=str(REPO))
    assert r.returncode == 0, f"资产哈希校验失败：\n{r.stdout}\n{r.stderr}"
    assert "异常 0" in r.stdout


def test_gene_vocabulary_matches_checkpoint():
    """派生词表与冻结 checkpoint 的 ``feature_name`` 逐一相同（顺序也相同）。"""
    from compass_os.preprocessing import feature_names
    from compass_os.representation import checkpoint_feature_names
    assert feature_names() == checkpoint_feature_names()
    assert len(feature_names()) == 15672


def test_concept_gene_sets_match_vendored_concept():
    """派生基因集表与 vendored ``conception_processed.tsv`` 逐值一致。"""
    from compass_os.qc import gene_sets
    from compass_os._paths import asset
    ref = pd.read_csv(asset("third_party/compass/tokenizer/conception_processed.tsv"),
                      sep="\t", index_col=0)
    gs = gene_sets()
    assert list(gs["gene_set"]) == [str(i) for i in ref.index]
    assert (gs["genes"].tolist() == ref["Genes"].astype(str).tolist())
    assert (gs["broad_pathway"].tolist() == ref["BroadCelltypePathway"].astype(str).tolist())
    assert gs["n_genes"].sum() == ref["Genes"].astype(str).str.split(":").apply(len).sum()


def test_cancer_codes_match_frozen_json():
    """派生癌种码表与 vendored ``cancer_code.json`` 完全一致。"""
    from compass_os._paths import asset
    from compass_os.qc import cancer_codes
    ref = json.loads(asset("third_party/compass/tokenizer/cancer_code.json").read_text())
    cc = cancer_codes().set_index("cancer_type")
    assert set(cc.index) == set(ref)
    for k, v in ref.items():
        assert int(cc.loc[k, "compass_code"]) == int(v)
    # CT_* 列中恰有 1 个是 drop-first 参照水平
    assert int(cc["is_reference_level"].sum()) == 1


def test_all_signature_genes_in_vocabulary():
    """132 基因集引用的 916 个基因必须全部在 15,672 词表内（P0 已核，此处固化为测试）。"""
    from compass_os.preprocessing import feature_names
    from compass_os.qc import gene_sets
    vocab = set(feature_names())
    used = {g for gl in gene_sets()["gene_list"] for g in gl}
    missing = sorted(used - vocab)
    assert not missing, f"词表外基因：{missing[:10]}"
    assert len(used) == 916


def test_no_absolute_paths_in_package():
    """正式包源码不得出现原项目绝对路径 / 本机用户目录。"""
    # 模式由片段拼接：避免"检测绝对路径的代码"本身被 release 审计误报
    pat = re.compile("(" + "|".join([
        "/" + "home/[A-Za-z0-9_.-]+/", "/" + "mnt/[a-z]/",
        "[A-Za-z]:" + chr(92) * 2, "pro" + "jects/202608", "ren" + "si"]) + ")")
    hits = []
    for f in sorted((REPO / "src").rglob("*.py")):
        for i, line in enumerate(f.read_text(encoding="utf-8").splitlines(), 1):
            if pat.search(line) and "http" not in line:
                hits.append(f"{f.relative_to(REPO)}:{i}: {line.strip()[:100]}")
    assert not hits, "源码中出现绝对路径：\n" + "\n".join(hits)


def test_locks_are_read_only_assets():
    """lock 的关键字段齐备（β / scaler / 填充 / 切点 / 基线），且模型不可在用户数据上 fit。"""
    from compass_os.survival import load_lock
    for mk in ("M0", "M1", "M2", "M3"):
        lock = load_lock(mk)
        for key in ("feature_names", "beta", "clinical_fill", "median_cutoff"):
            assert key in lock, f"{mk} 缺 {key}"
        assert len(lock["beta"]) == len(lock["feature_names"])
        if lock.get("scaled_features"):
            assert len(lock["scaler_mean"]) == len(lock["scaled_features"])
            assert len(lock["scaler_scale"]) == len(lock["scaled_features"])
    # M3 必须有 PCA，且只在 lock 内
    from compass_os.survival import load_pca
    pca = load_pca()
    assert pca["components"].shape == (10, 15672)
    assert len(pca["genes"]) == 15672


def test_scaler_fit_is_forbidden():
    """显式拦截：禁止在用户数据上拟合 scaler（PCA/缩放泄漏防线）。"""
    from compass_os.representation import _ZeroMaskedScaler
    from compass_os.representation import load_model
    s = _ZeroMaskedScaler(load_model().scaler, ["A1BG"])
    for fn in (s.fit, s.fit_transform):
        try:
            fn(pd.DataFrame({"a": [1.0]}))
        except RuntimeError as exc:
            assert "禁止" in str(exc)
        else:
            raise AssertionError("fit 未被拦截")
