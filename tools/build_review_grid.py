#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""把 redesign 渲染器的 6 张单图拼成 2×3 人工验收总览图（只拼贴，不重绘）。

用法：
    python tools/redesign_render/build_redesign_grid.py docs/figures/sample_path_redesign
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt

ORDER = [("M2","low"),("M2","mid"),("M2","high"),("M3","low"),("M3","mid"),("M3","high")]
COLS = ["Low (pct≈10)", "Mid (pct≈50)", "High (pct≈90)"]

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("root", help="含 <level>/<model>.{png,selection.json} 的目录")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    D = Path(a.root)
    fig, axes = plt.subplots(2, 3, figsize=(30, 12.4), dpi=100)
    for ax, (model, lvl) in zip(axes.ravel(), ORDER):
        img = D / lvl / f"{model}.png"
        if not img.is_file():
            ax.axis("off"); ax.set_title(f"{model} {lvl}: 缺少 {img.name}", fontsize=12); continue
        ax.imshow(mpimg.imread(img)); ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values(): s.set_edgecolor("#888")
        sel = json.loads((D / lvl / f"{model}.selection.json").read_text(encoding="utf-8"))
        c = sel.get("counts", {})
        ax.set_title(f"{model}  |  {sel.get('sample','?')}  |  risk={sel.get('risk',float('nan')):.4f}  "
                     f"rank={sel.get('rank','?')}/{sel.get('n','?')}  pct={sel.get('percentile',float('nan')):.1f}%  "
                     f"|  genes={c.get('gene_score','?')} sigs={c.get('signatures','?')} "
                     f"predictors={c.get('predictors','?')}", fontsize=12, pad=8)
    for c, lab in enumerate(COLS):
        axes[0, c].annotate(lab, xy=(0.5, 1.10), xycoords="axes fraction",
                            ha="center", fontsize=15, fontweight="bold")
    fig.suptitle("COMPASS-OS sample computation path — redesign renderer, human review sheet\n"
                 "GSE39582 (COAD, n=573) · REAL COMPASS-OS inference · M2 vs M3 on the same three "
                 "samples · one shared Style · composited from the rendered figures (no re-plotting)",
                 fontsize=15)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    out = Path(a.out) if a.out else D / "sample_path_redesign_review_grid.png"
    fig.savefig(out, dpi=100, bbox_inches="tight"); plt.close(fig)
    print(f"✅ {out} ({out.stat().st_size:,} bytes)")
    return 0

if __name__ == "__main__":
    sys.exit(main())
