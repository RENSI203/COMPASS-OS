# Sample computation path — visual review (round 1)

**Purpose.** Show what the **current, unmodified** implementation actually draws so the layout
can be reviewed before any further change.

**This round changed no code.** No merge, no tag, no push.

> ⚠️ **The figures in this round are visually poor, and that is the intended finding.**
> The list of defects in §4 is the deliverable; §5 records what I did *not* fix on purpose.

---

## 1. Cohort and programmatic sample selection

| | |
|---|---|
| Cohort | **GSE39582** |
| Cancer type | **COAD** (single-cancer cohort, 1 unique type) |
| Samples in cohort | **n = 573** |
| Expression scale | `log2_tpm1` (harmonised cache, read-only) |
| Genes | 15,672 (full COMPASS vocabulary) |
| Selection rule | M2 risk percentile closest to **10 / 50 / 90**, computed programmatically |
| M2 & M3 on the same cohort | ✅ both run on the identical matrix; **same 3 samples** used for both models |

No sample was hand-picked, no data was altered, no node was removed for appearance.

| level | sample_id | M2 risk | M3 risk | M2 rank | M3 rank | M2 pct | M3 pct |
|---|---|---|---|---|---|---|---|
| low | `GSE39582::GSM1681368` | 0.955266 | 0.057925 | 516 / 573 | 464 / 573 | 10.03 | 19.11 |
| mid | `GSE39582::GSM972085` | 1.246794 | 0.262943 | 287 / 573 | 286 / 573 | 50.00 | 50.17 |
| high | `GSE39582::GSM972003` | 1.530004 | 0.514420 | 58 / 573 | 77 / 573 | 89.97 | 86.65 |

**M2 → M3 rank shift on the same samples:** low `−52`, mid `−1`, high `+19`.
M3 pulls both extremes toward the middle; the PC block contribution to η is
**0.4685 / 0.9205 / 0.7141** for low / mid / high respectively.

---

## 2. Per-figure acceptance table

| Model | Sample | Risk | Percentile | Cutoff (cohort median) | Visible concepts | Clinical nodes | PC node | Imputed? | HTML | PNG |
|---|---|---|---|---|---|---|---|---|---|---|
| M2 | GSM1681368 | 0.9553 | 10.0 | 1.2468 | 13 | Age, Sex, Stage | — | **yes (all 3)** | ✅ 41 KB | ✅ 703 KB |
| M2 | GSM972085 | 1.2468 | 50.0 | 1.2468 | 13 | Age, Sex, Stage | — | **yes (all 3)** | ✅ 43 KB | ✅ 769 KB |
| M2 | GSM972003 | 1.5300 | 90.0 | 1.2468 | 13 | Age, Sex, Stage | — | **yes (all 3)** | ✅ 42 KB | ✅ 703 KB |
| M3 | GSM1681368 | 0.0579 | 19.1 | 0.2629 | 12 | Age, Sex, Stage | **PC1–PC10 (1 node)** | **yes (all 3)** | ✅ 42 KB | ✅ 740 KB |
| M3 | GSM972085 | 0.2629 | 50.2 | 0.2629 | 12 | Age, Sex, Stage | **PC1–PC10 (1 node)** | **yes (all 3)** | ✅ 43 KB | ✅ 763 KB |
| M3 | GSM972003 | 0.5144 | 86.6 | 0.2629 | 12 | Age, Sex, Stage | **PC1–PC10 (1 node)** | **yes (all 3)** | ✅ 41 KB | ✅ 742 KB |

*Predictor-node budget = 16 in every figure* (M2: 13 concepts + 3 clinical; M3: 12 concepts + 3
clinical + 1 PC aggregate). Rendering also includes 2 bookkeeping nodes (see §3), so the
high-level column actually contains 18 drawn nodes.

**Cancer type** is omitted in all six (single-cancer cohort) and the mandated caption is present
in `res.notes` and in the figure subtitle.

---

## 3. Files

```
docs/figures/sample_path_review/
├── M2_low.html / .png / .pdf        M2_mid.*        M2_high.*
├── M3_low.html / .png / .pdf        M3_mid.*        M3_high.*
├── sample_path_M2_M3_review_grid.png      ← 2×3 review sheet (composited, not re-plotted)
├── plotly.min.js                          ← shared, so each HTML works offline
├── mpl_fallback/*.png                     ← the package's own save_static() output (different rendering)
├── _selection.tsv                         ← programmatic selection record
└── _meta.json                             ← per-figure metadata
```

**PNG = a real raster of the same figure as the HTML.** They were produced by rendering the
self-contained HTML in headless Edge (Chromium) at a **fixed 1500 × 860 window** for all six,
so M2 and M3 share figure width, column positions, font size, node spacing, palette, margins
and legend position. `kaleido` is unavailable and PyPI is unreachable here, and the plotly CDN
is also unreachable, so:

* the committed HTML is **self-contained** (loads `plotly.min.js` from the same directory) —
  it opens and hovers offline;
* the PNGs are Edge screenshots rather than `fig.write_image()` output.

`mpl_fallback/` keeps what the package's own `save_static()` produces without `kaleido` — a
matplotlib contribution-bar rendering, **not** a Sankey. That difference is itself a finding.

---

## 4. Answers to the ten questions

### 1. Do the figures have separate “Gene TPM” and “COMPASS gene score” columns?

**No.** The current figure has **one** gene column:

`gene (symbol) → gene-set/signature → concept → [aux] → risk`

Gene **TPM** is encoded only as node colour, and `gene_score` appears only in the hover text and
in `res.genes` — it is **not** a column. So the intended five-stage chain
(Gene TPM → gene score → signature → high-level → risk) is currently rendered as four.

This is a genuine structural gap, and it is **implementable**: the model genuinely produces a
distinct gene-level score — `extract(with_gene_level=True)` returns `dfg` =
`geneset_scorer(encoder_output)[:, 2:]`, a learned 32→1 projection per gene — which is
numerically different from TPM. Nothing was fixed this round, as instructed.

### 2. Are all three high-dimensional layers really truncated to ≤ 50?

**Yes — by real node/path removal plus explicit residual aggregation, not by blanking labels.**

| layer | cap | actual (six figures) |
|---|---|---|
| gene | 50 | **50** in all six |
| signature | 50 | 45, 50, 49 (M2) / 47, 47, 46 (M3) |
| concept | 16-predictor budget | 13 (M2) / 12 (M3) |

Signature counts fall below 50 only because ranked entries whose contribution is exactly zero
are dropped, not because of a hidden cap. Every truncation is complemented by a residual node
(§3 of the feasibility audit) so no node silently loses flow.

### 3. Is the high-level column strictly ≤ 16 predictor nodes?

**Yes — exactly 16 in all six figures.** M2 (single-cancer): 13 concepts + Age + Sex + Stage.
M3 (single-cancer): 12 concepts + Age + Sex + Stage + PC1–PC10.
Two further nodes are *bookkeeping* and are reported separately
(`n_high_level_nodes_total = 18`): “Other concepts (aggregated)” and
“Frozen model centering (constant)”.

### 4. Does M3 use exactly one visible PC1–PC10 node?

**Yes** — a single node labelled `PC1–PC10`; the ten per-PC contributions are kept in
`res.auxiliaries['detail']` and appear in the node's hover text.

### 5. Is the PC aggregate contribution Σ β_PC × PC?

**Yes.** The node value equals the **sum** of the ten per-PC Cox contributions
(not a mean), verified against the per-feature decomposition table:
`0.4685` (low), `0.9205` (mid), `0.7141` (high). Guarded by
`test_pc_aggregation_is_sum_not_mean`.

### 6. Are risk / rank / percentile taken directly from the main API?

**Yes.** `res.risk` is exactly `float(predict(...).risk[model].loc[sample_id])`; rank and
percentile derive from the official cohort risk ordering. The attribution chain is a separate
float32 forward pass and its deviation from the API value is reported
(`eta_max_abs_diff` = `2.1e-07 … 1.0e-06` in these six figures).

### 7. How much does the ranking change between M2 and M3 on the same samples?

| level | M2 rank (pct) | M3 rank (pct) | Δrank | Δpct |
|---|---|---|---|---|
| low | 516 (10.0) | 464 (19.1) | **−52** | **+9.1** |
| mid | 287 (50.0) | 286 (50.2) | −1 | +0.2 |
| high | 58 (90.0) | 77 (86.6) | **+19** | **−3.3** |

Adding PC1–PC10 compresses the tails: the low-risk sample rises out of the bottom decile and the
high-risk sample falls back from the top decile.

### 8. Do the PNG and the HTML express the same information?

**In this round, yes** — the PNGs are Chromium screenshots of the very same HTML, so hover
content, node names and contributions are identical.

**But the package's own `save_static()` does not.** Without `kaleido` it falls back to a
matplotlib contribution-bar figure (`mpl_fallback/`) that shows the same *numbers* but a
completely different *layout* — no Sankey, no flow, no risk column. Relying on that fallback as
“the static deliverable” would misrepresent the feature, so it is quarantined in
`mpl_fallback/` rather than used as the review PNG.

### 9. Most obvious visual problems

1. **Vertical overcrowding.** 50 genes + ~50 signatures + 18 high-level nodes are compressed
   into ~800 px; node labels overlap into illegible blocks in every panel. This is the single
   biggest usability problem.
2. **Residual/bookkeeping nodes dominate.** “Other genes (not shown)”, “Other signatures (not
   shown)” and “Frozen model centering (constant)” carry far more mass than any individual gene
   or signature, producing a few overwhelming bands that drown the actual per-node detail.
3. **The risk column is effectively invisible.** It occupies only 13 % of the width, the
   red→green gradient is squeezed against the right edge, and the heatmap's y-axis ticks
   (`1, 0.5, 0, −0.5, −1`) leak through despite `showticklabels=False`.
4. **The sample marker does not encode cohort position.** The scale is computed as
   `lo = min(Σcontribution, risk)`, `hi = max(lo, risk)` — both are the sample's own values, so
   `frac` collapses to ≈ 0 or 1 and the triangle is pinned to the top/bottom of the bar. This is
   a **functional defect**, not an aesthetic one: the marker currently tells the reader nothing.
5. **The cutoff line is not drawn.** The implementation writes a shape with
   `color="rgba(0,0,0,0)"` and `width=0`; the cutoff exists only in the title annotation.
6. **Title/annotation collision.** The in-figure annotation row (“Risk = … Rank = …”, the cutoff
   caption and the link-legend line) overlaps the top node labels.
7. **No separate gene-score column** (see question 1).
8. **The caption still uses flow-conservation wording.** The brief asks for
   *“Link width represents contribution magnitude; link color represents contribution direction.
   Exact signed Cox accounting is reported separately.”* The current legend instead implies the
   visual flow equals η.
9. **Imputed clinical values are not marked.** This cohort has no clinical data, so Age, Sex and
   Stage are **all** frozen-imputed, yet the nodes read `Age`, `Sex`, `Stage` with no
   `(imputed)` marker.
10. **Weak layer colour separation.** Genes (grey), signatures (pale blue) and concepts (purple)
    are not distinct enough at this density to read the four stages fluently.
11. **Risk-bar scale is not comparable across figures.** Because the scale uses the sample's own
    values, the six bars are not on a shared axis; M2 and M3 bars cannot be compared by looking.

### 10. Things I believe need fixing, but did **not** change this round

| # | change | why it matters |
|---|---|---|
| 1 | Fix the risk-bar scale to the **cohort** risk range (min…max of the cohort's own risk output) | the marker is currently meaningless |
| 2 | Draw the cutoff as a **visible dashed line** + `Cutoff = …` label | required element of the design |
| 3 | Add a dedicated **“COMPASS gene score”** column between TPM and signature | matches the frozen target structure; the data already exists |
| 4 | Reduce default density (fewer shown nodes, or taller figure / larger canvas) | readability |
| 5 | De-emphasise residual/bookkeeping nodes (thinner links, distinct style, off to the side) | they currently dominate |
| 6 | Re-label imputed clinical nodes as `Age (imputed)` etc., incl. hover | avoid implying real inputs |
| 7 | Replace the caption with the required magnitude/direction wording and drop the “inflow = η” framing | scientific wording discipline |
| 8 | Give the risk column more width and suppress the leaking axis ticks | currently invisible |
| 9 | Define the risk bar on a **shared** scale across M2/M3 (e.g. both on the M2 or a common range) | required for direct comparison |
| 10 | Fix annotation/title overlap | readability |

---

## 5. Explicit statement of what was not done

* No source file was modified in this round.
* No merge, tag or push.
* No sample was swapped, no risk value, contribution, concept score or PC value was altered.
* No node was removed to improve appearance, and no PC contribution was scaled up to make M3
  look different.
* The poor appearance above **is** the measured state of the current implementation.
