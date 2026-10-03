# -*- coding: utf-8 -*-
"""包内资源解析：**正式 package 不依赖原项目绝对路径**。

解析顺序：
1. 环境变量 ``COMPASS_OS_ROOT``（指向包含 ``models/`` 的仓库根）；
2. 由 ``__file__`` 推断仓库根（``src/compass_os/_paths.py`` → ``parents[2]``），
   适用于源码目录直接使用与 editable 安装；
3. 均失败 ⇒ 抛出带明确指引的 :class:`AssetNotFoundError`。

上游 COMPASS 包解析顺序（**优先 vendored 副本**）：

1. ``<repo>/third_party/compass``（随包发布、版本固定、已通过 golden test）；
2. 已安装的 ``compass`` / ``immuno-compass``（例如 wheel 安装且无 ``third_party/`` 时）。

优先 vendored 是为了**可复现**：本包的全部数值验证都是在 vendored 副本上做的，
即使使用者环境里另有一份 editable/不同来源的 ``compass`` 也不会改变结果。
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

from .exceptions import AssetNotFoundError

_ENV_ROOT = "COMPASS_OS_ROOT"


def repo_root() -> Path:
    """返回包含 ``models/`` 的仓库根（可由 ``COMPASS_OS_ROOT`` 覆盖）。"""
    env = os.environ.get(_ENV_ROOT)
    if env:
        p = Path(env).expanduser().resolve()
        if not (p / "models").is_dir():
            raise AssetNotFoundError(
                f"{_ENV_ROOT}={p} 下找不到 models/ 目录；请指向包含 models/ 的仓库根")
        return p
    p = Path(__file__).resolve().parents[2]
    if (p / "models").is_dir():
        return p
    raise AssetNotFoundError(
        "无法定位模型资产目录：既未设置环境变量 COMPASS_OS_ROOT，也未能在 "
        f"{p} 找到 models/。请设置 COMPASS_OS_ROOT=<包含 models/ 的仓库根>。")


def asset(rel: str) -> Path:
    """解析包内资产相对路径（如 ``models/locked_M2.json``）为绝对路径。"""
    p = repo_root() / rel
    if not p.exists():
        raise AssetNotFoundError(f"资产不存在：{p}")
    return p


def model_manifest() -> dict:
    import json
    return json.loads(asset("models/model_manifest.json").read_text(encoding="utf-8"))


_compass_ready = False


def ensure_compass() -> None:
    """确保 ``import compass`` 指向 **vendored 副本**（不存在时才用已安装的包）。"""
    global _compass_ready
    if _compass_ready:
        return
    vendor = repo_root() / "third_party"
    if (vendor / "compass").is_dir():
        # 已安装的 compass（可能是 editable 指向别的源码树）优先度**低于**随包副本：
        # 先把它从 sys.modules 清掉，再让 third_party 排在 sys.path 最前。
        for mod in [m for m in sys.modules if m == "compass" or m.startswith("compass.")]:
            del sys.modules[mod]
        if str(vendor) in sys.path:
            sys.path.remove(str(vendor))
        sys.path.insert(0, str(vendor))
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
                "请安装 immuno-compass==2.5.3 或使用完整仓库") from None
    _compass_ready = True


def compass_origin() -> str:
    """当前 ``compass`` 的实际来源路径（用于自查/报告）。"""
    ensure_compass()
    import compass
    return str(Path(compass.__file__).resolve())
