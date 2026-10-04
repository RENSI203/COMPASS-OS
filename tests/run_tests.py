#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""无第三方依赖的测试运行器（环境无 pytest 时使用）。

用法：
    python tests/run_tests.py            # 全部
    python tests/run_tests.py test_api   # 指定模块（可多个）

测试本身是 pytest 兼容的（模块级 ``test_*`` 函数 + 普通 assert），
有 pytest 时也可直接 ``pytest tests/``。
"""
from __future__ import annotations

import importlib
import inspect
import sys
import time
import traceback
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "src"))

MODULES = ["test_assets", "test_duplicate_genes", "test_reproducibility",
           "test_api", "test_analysis", "test_qc_tiers", "test_v101_regressions",
            "test_packaging", "test_sample_path"]


class _Skip(Exception):
    pass


class _Mark:
    """最小化 pytest.mark 替身（支持 skipif 与单参数 parametrize）。"""

    def __init__(self, condition: bool = False, reason: str = "", params=None):
        self.condition, self.reason, self.params = condition, reason, params

    def __call__(self, fn):
        if self.params is not None:
            fn.__params__ = self.params
        else:
            fn.__skip__ = self
        return fn


def _parametrize(argnames, argvalues):
    names = [n.strip() for n in str(argnames).split(",")]
    cases = [dict(zip(names, v if isinstance(v, (tuple, list)) else (v,)))
             for v in argvalues]
    return _Mark(params=cases)


class _Approx:
    def __init__(self, v, rel=1e-6, abs_=1e-12):
        self.v, self.rel, self.abs = v, rel, abs_

    def __eq__(self, other):
        return abs(float(other) - float(self.v)) <= max(self.abs, self.rel * abs(float(self.v)))


class _PytestShim:
    approx = staticmethod(lambda v, rel=1e-6, abs=1e-12: _Approx(v, rel, abs))
    mark = type("mark", (), {
        "skipif": staticmethod(lambda cond, reason="": _Mark(cond, reason)),
        "parametrize": staticmethod(_parametrize)})()

    class _Raises:
        def __init__(self, exc):
            self.exc = exc
            self.value = None          # 与 pytest 一致：`with raises(X) as e: e.value`

        def __enter__(self):
            return self

        def __exit__(self, et, ev, tb):
            if et is None:
                raise AssertionError(f"期望抛出 {self.exc.__name__}，但没有异常")
            self.value = ev
            return issubclass(et, self.exc)

    @staticmethod
    def raises(exc):
        return _PytestShim._Raises(exc)


def run_module(name: str) -> tuple:
    sys.modules.setdefault("pytest", _PytestShim())
    mod = importlib.import_module(name)
    passed, failed, skipped = [], [], []
    tests = [(n, f) for n, f in vars(mod).items()
             if n.startswith("test_") and inspect.isfunction(f)]
    for tname, fn in tests:
        mark = getattr(fn, "__skip__", None)
        if mark is not None and mark.condition:
            skipped.append((tname, mark.reason))
            continue
        cases = getattr(fn, "__params__", None) or [{}]
        for case in cases:
            label = tname + ("" if not case else
                             "[" + ",".join(f"{k}={v}" for k, v in case.items()) + "]")
            t0 = time.time()
            try:
                fn(**case)
            except Exception:  # noqa: BLE001
                failed.append((label, None, traceback.format_exc(), time.time() - t0))
            else:
                passed.append((label, time.time() - t0))
    return passed, failed, skipped


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    mods = argv or MODULES
    total_p = total_f = total_s = 0
    for m in mods:
        print(f"\n=== {m} ===")
        try:
            p, f, s = run_module(m)
        except Exception:  # noqa: BLE001 - 模块导入失败
            print(f"  !! 模块导入失败\n{traceback.format_exc()}")
            total_f += 1
            continue
        for n, dt in p:
            print(f"  PASS  {n}  ({dt:.1f}s)")
        for n, reason in s:
            print(f"  SKIP  {n}  ({reason})")
        for n, _exc, tb, dt in f:
            print(f"  FAIL  {n}  ({dt:.1f}s)\n{_indent(tb)}")
        total_p += len(p)
        total_f += len(f)
        total_s += len(s)
    print(f"\n合计：PASS {total_p} / FAIL {total_f} / SKIP {total_s}")
    return 1 if total_f else 0


def _indent(tb: str) -> str:
    return "\n".join("        " + ln for ln in tb.strip().splitlines()[-12:])


if __name__ == "__main__":
    sys.exit(main())
