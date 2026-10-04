# -*- coding: utf-8 -*-
"""v1.0.1 packaging regression tests —— 资产必须真正随包分发。

覆盖 v1.0.0 的 packaging 缺陷：``models/`` 与 ``third_party/`` 位于**仓库根**，
``pip install .`` / wheel 安装后无法定位 ⇒ 推理不可用。

v1.0.1 起资产位于 ``src/compass_os/assets/``，以 package data 进入 wheel / sdist。
"""
from __future__ import annotations

import os
import zipfile
from pathlib import Path

import pytest

from _util import repo_root

PKG = repo_root() / "src" / "compass_os"
ASSETS = PKG / "assets"

#: wheel / sdist 必须携带的 runtime 资产（相对包内路径）
REQUIRED = [
    "assets/models/pretrainer.pt",
    "assets/models/locked_M0.json",
    "assets/models/locked_M1.json",
    "assets/models/locked_M2.json",
    "assets/models/locked_M3.json",
    "assets/models/locked_pca_M3.npz",
    "assets/models/reference_quantiles.json",
    "assets/models/model_manifest.json",
    "assets/models/qc_config.json",
    "assets/third_party/COMPASS_LICENSE",
    "assets/third_party/compass/__init__.py",
    "assets/third_party/compass/model/_distilled.pth",
    "assets/third_party/compass/tokenizer/gene_tokens_long.json",
    "assets/third_party/compass/tokenizer/conceptor_gene_map.csv",
    "data/gene_vocabulary.txt",
    "data/concept_gene_sets.tsv",
    "data/cancer_codes.tsv",
]


def test_assets_live_inside_the_package():
    """runtime 资产必须在包内（不得继续放在仓库根，否则 wheel 无法携带）。"""
    assert ASSETS.is_dir(), "包内缺少 assets/ 目录"
    assert not (repo_root() / "models").exists(), \
        "仓库根仍有 models/：wheel 安装将无法定位资产"
    assert not (repo_root() / "third_party").exists(), \
        "仓库根仍有 third_party/：wheel 安装将无法定位 vendored COMPASS"
    for rel in REQUIRED:
        assert (PKG / rel).is_file(), f"包内缺少 {rel}"


def test_assets_resolve_without_environment_variable():
    """不设 COMPASS_OS_ROOT 也能解析全部资产（wheel/源码树通用）。"""
    from compass_os import _paths
    saved = os.environ.pop("COMPASS_OS_ROOT", None)
    _paths.assets_root.cache_clear()
    try:
        _check_env_free_resolution(_paths)
    finally:
        if saved is not None:
            os.environ["COMPASS_OS_ROOT"] = saved
        _paths.assets_root.cache_clear()


def _check_env_free_resolution(_paths):
    root = _paths.assets_root()
    assert root.name == "assets" and (root / "models" / "model_manifest.json").is_file()
    for rel in REQUIRED:
        assert _paths.asset(rel).is_file(), f"解析失败：{rel}"
    assert _paths.qc_config_path() is not None


def test_compass_resolves_to_vendored_copy():
    """默认执行路径必须是**包内 vendored** COMPASS（保证与验证时一致）。"""
    from compass_os import _paths
    saved = os.environ.pop("COMPASS_OS_ROOT", None)
    _paths.assets_root.cache_clear()
    try:
        origin = Path(_paths.compass_origin())
        assert origin.is_relative_to(_paths.assets_root() / "third_party"), \
            f"compass 未从包内 vendored 副本导入：{origin}"
    finally:
        if saved is not None:
            os.environ["COMPASS_OS_ROOT"] = saved
        _paths.assets_root.cache_clear()


def test_compass_os_root_override_still_supported():
    """COMPASS_OS_ROOT 仍可作为**开发者覆盖**（指向 assets/ 目录）。"""
    from compass_os import _paths
    saved = os.environ.get("COMPASS_OS_ROOT")
    os.environ["COMPASS_OS_ROOT"] = str(ASSETS)
    _paths.assets_root.cache_clear()
    try:
        assert _paths.assets_root() == ASSETS.resolve()
    finally:
        if saved is None:
            os.environ.pop("COMPASS_OS_ROOT", None)
        else:
            os.environ["COMPASS_OS_ROOT"] = saved
        _paths.assets_root.cache_clear()


def test_setuptools_configuration_declares_package_data():
    """pyproject 必须声明 assets 为 package data，并排除其被当作子包发现。"""
    txt = (repo_root() / "pyproject.toml").read_text(encoding="utf-8")
    assert "assets/models/*.pt" in txt and "assets/models/*.json" in txt
    assert "assets/third_party/compass/**/*" in txt.replace(" ", "") or \
           "assets/third_party/compass/**/*.py" in txt
    assert 'exclude = ["compass_os.assets*"]' in txt, \
        "vendored compass 不得被 setuptools 当作 compass_os.assets.* 子包发现"
    import re as _re
    assert _re.search(r'version = "\d+\.\d+\.\d+"', txt), "缺少版本号"


def test_manifest_in_excludes_bytecode():
    """MANIFEST.in 必须排除字节码（vendored 包导入时会生成 __pycache__）。"""
    txt = (repo_root() / "MANIFEST.in").read_text(encoding="utf-8")
    assert "global-exclude" in txt and "__pycache__" in txt


@pytest.mark.skipif(not (repo_root() / "dist").is_dir(),
                    reason="需要先 python -m build（dist/ 不存在）")
def test_built_artifacts_contain_runtime_assets():
    """已构建的 wheel / sdist 必须包含全部 runtime 资产，且不含字节码。

    这是 v1.0.1 的核心 packaging 断言：仅"build 成功"不算通过。
    """
    wheels = sorted((repo_root() / "dist").glob("compass_os-*.whl"))
    sdists = sorted((repo_root() / "dist").glob("compass_os-*.tar.gz"))
    assert wheels, "dist/ 下没有 wheel"
    assert sdists, "dist/ 下没有 sdist"

    with zipfile.ZipFile(wheels[-1]) as z:
        names = z.namelist()
    missing = [r for r in REQUIRED if not any(n.endswith("compass_os/" + r) for n in names)]
    assert not missing, f"wheel 缺少资产：{missing}"
    assert not [n for n in names if "__pycache__" in n or n.endswith((".pyc", ".pyo"))], \
        "wheel 含字节码"
    assert any("dist-info/licenses/LICENSE" in n for n in names), "wheel 缺少 LICENSE"
    assert any("dist-info/licenses/NOTICE" in n for n in names), "wheel 缺少 NOTICE"

    import tarfile
    with tarfile.open(sdists[-1]) as tf:
        tn = tf.getnames()
    missing_s = [r for r in REQUIRED
                 if not any(n.endswith("src/compass_os/" + r) for n in tn)]
    assert not missing_s, f"sdist 缺少资产：{missing_s}"
    assert not [n for n in tn if "__pycache__" in n or n.endswith((".pyc", ".pyo"))], \
        "sdist 含字节码"
