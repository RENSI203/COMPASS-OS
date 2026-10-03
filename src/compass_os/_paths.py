# -*- coding: utf-8 -*-
"""包内资源解析：**安装后不依赖 Git 仓库、不依赖 ``COMPASS_OS_ROOT``**。

自 v1.0.1 起，全部 runtime 资产以 **package data** 形式随 wheel 分发::

    src/compass_os/assets/
    ├── models/            # pretrainer.pt, locked_M0-M3.json, locked_pca_M3.npz,
    │                      # reference_quantiles.json, model_manifest.json, qc_config.json
    └── third_party/       # vendored upstream COMPASS + COMPASS_LICENSE

解析顺序：

1. ``COMPASS_OS_ROOT``（**开发者覆盖**，可选）：指向包含 ``models/`` 的目录、
   直接指向 ``assets/``，或指向 ``src/compass_os``；
2. **package-internal assets**（wheel / sdist / editable / 源码树 通用），
   经 :mod:`importlib.resources` 解析 —— 这是 v1.0.1 的规范路径；
3. 均失败 ⇒ 抛出带明确指引的 :class:`AssetNotFoundError`。

上游 COMPASS 包解析顺序（**优先 vendored 副本**）：

1. ``<assets>/third_party/compass``（随包发布、版本固定、已通过 golden test）；
2. 已安装的 ``compass`` / ``immuno-compass``（仅在 vendored 副本缺失时）。

优先 vendored 是为了**可复现**：本包的全部数值验证都是在 vendored 副本上做的，
即使使用者环境里另有一份 editable/不同来源的 ``compass`` 也不会改变结果。
"""
from __future__ import annotations

import os
import sys
from functools import lru_cache
from pathlib import Path

from .exceptions import AssetNotFoundError

_ENV_ROOT = "COMPASS_OS_ROOT"


def _looks_like_assets(p: Path) -> bool:
    """该目录是否直接包含 ``models/model_manifest.json``。"""
    return (p / "models" / "model_manifest.json").is_file()


@lru_cache(maxsize=1)
def assets_root() -> Path:
    """返回包含 ``models/`` 与 ``third_party/`` 的资产根目录（**绝对路径**）。

    wheel 安装 ⇒ ``<site-packages>/compass_os/assets``；
    源码树 / editable ⇒ ``<repo>/src/compass_os/assets``。
    """
    env = os.environ.get(_ENV_ROOT)
    if env:
        p = Path(env).expanduser().resolve()
        for cand in (p, p / "assets", p / "src" / "compass_os" / "assets"):
            if _looks_like_assets(cand):
                return cand
        raise AssetNotFoundError(
            f"{_ENV_ROOT}={p} 下找不到 models/model_manifest.json；请指向包含 models/ 的"
            "目录、assets/ 目录，或取消该环境变量以使用包内资产。")

    # 规范路径：package-internal assets（wheel / editable / 源码树通用）
    try:
        from importlib.resources import files as _res_files
        cand = Path(str(_res_files(__package__))) / "assets"
        if _looks_like_assets(cand):
            return cand.resolve()
    except Exception:  # noqa: BLE001 —— 退回 __file__ 推断
        pass

    cand = Path(__file__).resolve().parent / "assets"
    if _looks_like_assets(cand):
        return cand.resolve()

    raise AssetNotFoundError(
        "无法定位模型资产目录：包内未找到 assets/models/model_manifest.json。"
        "该安装可能不完整（wheel 应包含 compass_os/assets/）；"
        f"也可用 {_ENV_ROOT} 显式指向资产目录。")


def package_root() -> Path:
    """包目录（``.../compass_os``），其中含 ``assets/`` 与 ``data/``。"""
    return assets_root().parent


def asset(rel: str) -> Path:
    """解析资产相对路径为绝对路径。

    支持两种书写习惯，便于源码树与 wheel 共用同一份调用代码：

    * ``models/locked_M2.json`` / ``third_party/...`` ⇒ ``<assets>/models/...``
    * ``data/gene_vocabulary.txt`` ⇒ ``<package>/data/...``
    * 历史写法 ``src/compass_os/data/...`` 会被自动归一化（去掉前缀）。
    """
    rel = str(rel).replace("\\", "/").lstrip("/")
    if rel.startswith("src/compass_os/"):
        rel = rel[len("src/compass_os/"):]
    for cand in (assets_root() / rel, package_root() / rel):
        if cand.exists():
            return cand
    raise AssetNotFoundError(f"资产不存在：{assets_root() / rel}")


def model_manifest() -> dict:
    import json
    return json.loads(asset("models/model_manifest.json").read_text(encoding="utf-8"))


def qc_config_path() -> Path | None:
    """返回 ``qc_config.json`` 路径；**缺失时返回 None**（安装不完整场景，不抛异常）。"""
    p = assets_root() / "models" / "qc_config.json"
    return p if p.is_file() else None


_compass_ready = False


def ensure_compass() -> None:
    """确保 ``import compass`` 指向 **vendored 副本**（不存在时才用已安装的包）。"""
    global _compass_ready
    if _compass_ready:
        return
    vendor = assets_root() / "third_party"
    if (vendor / "compass").is_dir():
        # 已安装的 compass（可能是 editable 指向别的源码树）优先度**低于**随包副本：
        # 先把它从 sys.modules 清掉，再让 assets/third_party 排在 sys.path 最前。
        for mod in [m for m in sys.modules if m == "compass" or m.startswith("compass.")]:
            del sys.modules[mod]
        vstr = str(vendor)
        if vstr in sys.path:
            sys.path.remove(vstr)
        sys.path.insert(0, vstr)
        import compass
        if not Path(compass.__file__).resolve().is_relative_to(vendor.resolve()):
            raise AssetNotFoundError(
                f"compass 未从 vendored 副本导入（实际：{compass.__file__}）")
    else:
        try:
            import compass  # noqa: F401
        except ImportError:
            raise AssetNotFoundError(
                f"vendored 副本不存在（{vendor / 'compass'}）且未安装 compass 包；"
                "请安装 immuno-compass==2.5.3 或使用完整 wheel") from None
    _compass_ready = True


def compass_origin() -> str:
    """当前 ``compass`` 的实际来源路径（用于自查/报告）。"""
    ensure_compass()
    import compass
    return str(Path(compass.__file__).resolve())


def repo_root() -> Path:
    """**已废弃**（v1.0.1 起资产随包分发）。返回资产根目录，保留仅为兼容旧调用。"""
    return assets_root()
