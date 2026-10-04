# Sample-level computation path (`sample_path`)

> **What this is.** A *sample-specific computational attribution / representation-flow
> visualization*: for one sample in an already-analysed cohort, it shows how the frozen model
> graph produced that sample's Cox linear predictor, layer by layer, and where the sample sits
> in the cohort risk distribution.
>
> **What this is not.** Not a causal mechanism diagram, not a pathway-activation map, not a
> biological causal network, and not a complete causal explanation of the predicted risk.
> Link colour encodes *contribution direction* (positive → higher model-estimated risk,
> negative → lower); representation colour encodes *high / low model representation value*.
> Neither may be read as "pathway activated / inhibited". The 43 concepts and 132 signatures
> remain mechanism-oriented **hypothesis-generating** representations.

---

## 1. Quick start

```python
import pandas as pd
import compass_os
from compass_os import sample_path

expr = pd.read_csv("expression.tsv", sep="\t", index_col=0)   # samples × genes, TPM
cancer = pd.read_csv("labels.tsv", sep="\t", index_col=0)["cancer_type"].tolist()
clinical = pd.read_csv("clinical.tsv", sep="\t", index_col=0)

res = sample_path(
    expr, cancer, "TCGA-19-1787",      # the sample of interest
    model="M2",                        # M1 / M2 / M3
    clinical=clinical,
)

res.summary()                          # plain-text overview
res.save_html("sample_path.html")      # interactive Plotly figure (required output)
res.save_static("sample_path.png")     # static: kaleido, else matplotlib fallback
```

## 2. What the figure contains

Left to right:

| column | content | how nodes are selected |
|---|---|---|
| gene expression | ≤ `top_genes` (default 50) genes | **by exact contribution** to the displayed concepts, not by TPM |
| gene-set / granular representation | ≤ `top_signatures` (50) of the 132 | by the sample's own contribution |
| high-level concepts | top concepts by `\|x_c β_c\|` | **not** by raw concept score |
| auxiliary predictors | Age, Sex, Stage, PC1–PC10 (M3), cancer type | direct Cox predictors of the frozen head |
| cohort-relative risk | continuous vertical column | red (higher) → green (lower) within the cohort |

Node colour and node *selection* are deliberately separate rules: colour reflects the
sample's expression / representation value, while selection reflects contribution.

### Exact flow conservation

Every internal node conserves flow. Where a layer is truncated for readability, an explicit
residual node carries the remainder:

* **Other genes (not shown)** → the displayed gene-sets
* **Other signatures (not shown)** → the displayed concepts
* **Other concepts (aggregated)** → risk
* **Frozen model centering (constant)** → the displayed concepts

Consequently the inflow at the risk node sums to the **full Cox linear predictor** — the figure
never implies a partial sum is the whole.

## 3. The decomposition is exact, and the risk is the main API's risk

Nothing is re-predicted for display. All quantities come from frozen artefacts and the official
extraction API:

| quantity | source | measured agreement |
|---|---|---|
| gene / 132 / 43 layers | `PreTrainer.extract(..., with_gene_level=True)` | — |
| `gene-set score = Σ_g a_{g,j}·gene_score_g` | frozen attention + shared scorer | `max\|Δ\| = 1.8e-07` |
| `concept score = Σ_j b_{j,c}·gene-set score_j` | frozen attention | `max\|Δ\| = 3.3e-08` |
| per-feature Cox contributions | `survival.build_design_matrix` + frozen `beta` | sums to η at `2.2e-16 … 8.9e-16` (M1/M2/M3) |
| displayed **risk / rank / percentile** | `compass_os.predict()` output | **identical** |

The attribution chain and the cohort prediction run as two independent float32 forward passes;
their linear predictors agree to ≈`1e-06`–`1e-05`. The difference is reported as
`SamplePathResult.eta_max_abs_diff`, and the figure always displays the **main API** risk.

## 4. Parameters

| parameter | default | meaning |
|---|---|---|
| `model` | `"M2"` | `M1` / `M2` / `M3`. **`M0` is not supported** — it has no COMPASS concept branch. |
| `top_genes` / `top_signatures` | `50` | per-layer display caps |
| `top_concepts` | `None` | explicit cap; otherwise the remaining high-level budget is used |
| `max_high_level_nodes` | `16` | **predictor nodes** in the high-level column (concepts + clinical/PC/cancer) |
| `show_cancer_type` | `None` (auto) | `None`: omitted for single-cancer cohorts (with a caption), shown for multi-cancer cohorts |
| `cutoff` | `"median"` | `"median"` = **cohort** median of the current risk output (visualization only) or a numeric risk threshold |
| `result` | `None` | pass an existing `PredictionResult` / `AnalysisResult` to reuse the cohort risk |
| `input_scale` / `missing_gene_strategy` | as in `predict` | identical alignment and imputation path |

### High-level node budget

The 16-node budget counts **Cox predictors** only (spec: concepts + clinical + PC aggregate +
other required Cox predictors). The two bookkeeping nodes (*Other concepts*, *centering
constant*) are reported separately as `n_high_level_nodes_total` and do not consume the budget.

| model | auxiliary predictors shown | concepts shown | total |
|---|---|---|---|
| M1 | — | 16 | 16 |
| M2 | Age, Sex, Stage, Cancer type | 12 | 16 |
| M3 | Age, Sex, Stage, Cancer type, PC1–PC10 | 11 | 16 |

(Single-cancer cohort: cancer type is omitted → one more concept slot.)

## 5. Reading the result object

```python
res.risk                       # Cox linear predictor for this sample (= main API value)
res.rank, res.n_samples        # 1 = highest risk in the cohort
res.percentile                 # 0–100, from the official risk ordering
res.cutoff, res.cutoff_mode    # and res.cutoff_is_frozen_validated (always False here)
res.genes / .signatures / .concepts / .auxiliaries     # displayed layers + contributions
res.decomposition              # full per-feature Cox decomposition (feature/value/beta/contribution)
res.gene_links / .signature_links                       # exact Sankey topology
res.other_concepts_contribution                         # aggregated hidden concepts
res.eta_max_abs_diff           # attribution-chain vs main-API deviation (float32)
res.notes                      # truncation, single-cancer, float32 caveats
```

## 6. Cutoffs: two different things

* `cutoff="median"` (default) — the **cohort median of the current risk output**. A
  visualization / cohort-relative grouping convenience. It is **not** a validated prognostic
  threshold, and it moves with the cohort.
* `cutoff=<number>` — a researcher-supplied risk threshold for display.
* The **frozen validated cutoff** (`median_cutoff` in the model lock, e.g. `1.0783` for M2) is a
  separate concept and is **not** selectable through this parameter;
  `cutoff_is_frozen_validated` is therefore always `False` here.

## 7. Dependencies

`plotly` is required for figures and is declared as the optional extra
`compass-os[plot]` (`plotly` + `kaleido`). **Core inference does not depend on it.**
Static export uses `kaleido` when present and otherwise falls back to a
dependency-free matplotlib rendering of the same attribution data (layered contribution bars).

```bash
pip install "compass-os[plot]"
```

## 8. Bundled examples

Reference outputs produced by `examples/sample_path_example.py` on the 5-sample example cohort
(sample `TCGA-19-1787`):

| file | content |
|---|---|
| [`figures/sample_path/sample_path_M1.html`](figures/sample_path/sample_path_M1.html) | interactive, M1 (16 concepts, no clinical/PCA) |
| [`figures/sample_path/sample_path_M2.html`](figures/sample_path/sample_path_M2.html) | interactive, M2 (12 concepts + Age/Sex/Stage/Cancer type) |
| [`figures/sample_path/sample_path_M3.html`](figures/sample_path/sample_path_M3.html) | interactive, M3 (11 concepts + those + aggregated PC1–PC10) |
| `figures/sample_path/sample_path_M{1,2,3}.png` | static (matplotlib fallback; kaleido not installed here) |
| [`figures/sample_path/example_M2.decomposition.tsv`](figures/sample_path/example_M2.decomposition.tsv) | full per-feature Cox decomposition for the M2 example |

## 9. Interpreting responsibly

* A gene's contribution may be **negative**; width encodes magnitude, colour encodes direction.
  A large width is not "biological importance".
* Concepts are **learned latent coordinates** (20/43 are reverse-encoded at the model level),
  and the default model uses them as Cox predictors — so any association with risk is a
  *model-linked representation association*, not independent evidence.
* For a single-cancer cohort the cancer-type term is constant by construction; the figure says
  so rather than hiding it.
* The right-hand column is *cohort-relative model risk*, **not** an absolute death probability.
