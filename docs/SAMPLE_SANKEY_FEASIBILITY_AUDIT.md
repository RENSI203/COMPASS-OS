# Sample-level computation-path Sankey — feasibility audit

**Scope.** This audit answers whether a *sample-specific computational attribution /
representation-flow visualization* can be built for COMPASS-OS **without** touching frozen
scientific assets, without changing the risk algorithm, and without duplicating prediction
logic.

**Semantic boundary (binding for all code, figures and docs).** What is visualized is a
*sample-specific computational attribution* through the frozen model graph. It is **not** a
causal mechanism diagram, **not** a pathway-activation map, **not** a biological causal
network, and **not** a complete causal explanation of the risk. Node and link wording must
stay inside that boundary (see §11).

---

## 0. Reference script — resolved paths

| | |
|---|---|
| Requested path | `build/lib/compass_os/assets/third_party/compass/utils/sankey.py` |
| **Authoritative resolved path** | `<repo>/src/compass_os/assets/third_party/compass/utils/sankey.py` (resolved, then reported relative to the repository root so the document stays portable) |
| Build-tree copy | `<repo>/build/lib/compass_os/assets/third_party/compass/utils/sankey.py` (gitignored build artefact; byte-identical to the above) |
| Size / lines | 9,652 bytes / 326 lines |
| SHA-256 (both copies, identical) | `20ed31d9eb38f438…` |

The `build/` tree is a transient copy produced by `python -m build`; the `src/…` path is the
vendored authoritative copy that ships inside the wheel. Both were read and are byte-identical.

Audited entry points: `plot_sankey_diagram()`, `get_projector_weights()`.

---

## 1. What layers does the original Sankey actually have?

Reading the code (**not** the example images), `plot_sankey_diagram()` emits exactly **two link
groups** and therefore **three node columns**:

| column | object in code | identity | count |
|---|---|---|---|
| 1 | `gene_name[geneset_idx[j]]` where `gene_name = model.feature_name` | gene symbols | 15,672 (subset shown) |
| 2 | `model.model.latentprojector.genesetprojector.genesets_names` | **gene-set / granular representation** | **132** |
| 3 | `model.model.latentprojector.cellpathwayprojector.cellpathway_names` | **high-level concept** | **43** |

Link groups: `"gene->geneset"` (lines 118–141) and `"geneset->celltype"` (lines 108–116).

**There is no Cox layer, no clinical layer, no PCA layer and no risk node in the original
script.** It stops at the concept column.

Naming note: the code calls column 3 `cellpathway` / `celltype`; the frozen COMPASS-OS data
dictionary and the package public API call the same 43 objects **concepts**. This document uses
*concept* and quotes the code identifier where it matters.

### 1.1 Stability of the accessed objects in the current frozen checkpoint

All objects required by the original script were verified present and structurally unchanged on
the shipped `pretrainer.pt` (introspection, read-only):

| attribute | status | shape / type |
|---|---|---|
| `model.feature_name` | ✅ | `Index`, 15,672 symbols |
| `…genesetprojector.genesets_names` | ✅ | `Index`, 132 (`Bcell_l_Danaher17`, …) |
| `…genesetprojector.genesets_indices` | ✅ | `list[132]` of gene-index lists (gene space, 0-based) |
| `…genesetprojector.geneset_aggregator.aggregator.attention_weights` | ✅ | `ParameterDict`, 132 keys, e.g. `geneset_0` → `(8, 1)` |
| `…genesetprojector.geneset_scorer.genesetscorer.fc` | ✅ | `nn.Linear(32 → 1)`, **shared across all 132 sets** |
| `…cellpathwayprojector.cellpathway_names` | ✅ | `Index`, 43 (`Bcell_general`, …) |
| `…cellpathwayprojector.cellpathway_indices` | ✅ | `list[43]` of gene-set index lists |
| `…cellpathwayprojector.cellpathway_aggregator.aggregator.attention_weights` | ✅ | `ParameterDict`, 43 keys, e.g. `cellpathway_0` → `(2, 1)` |

Both aggregators are `mode="attention"` with `softmax_mean=True`.

**Correction to the original script's vocabulary:** the code uses `genesets_names` /
`genesets_indices` (plural *genesets*) on the gene-set projector. The singular spellings
`geneset_names` / `geneset_indices` do **not** exist. The original script uses the correct
plural names; new code must too.

---

## 2. What does the original script really display? (top-K behaviour)

**The original default is `topK_vis=10`, not 50.** More importantly, `topK_vis` does **not**
remove nodes:

```python
tail_nodes = (df.groupby(["group","source"]).weights.sum().sort_values(ascending=False)
                .loc["gene->geneset"].iloc[topK_vis:].index)
filtered_labels = ["" if label in tail_nodes else label for label in labels]
```

* The slicing is applied **only** to the `gene->geneset` group, i.e. to **genes**.
* It only blanks the **label**; the node and all its links are still drawn.

Findings, summarised:

| question | answer |
|---|---|
| separate top-K for genes? | Only label-visibility, and only for genes. Nodes/links remain. |
| separate top-K for gene-sets? | **None.** |
| separate top-K for concepts? | **None** — but the real filter is `concept2plot` (default `["NKcell"]`), which selects which concepts (and hence which gene-sets and genes) enter the diagram at all. |
| default `topK_vis` | `10` |
| default `concept2plot` | `["NKcell"]` (one concept) |
| duplicate `(source, target)` pairs | **averaged**, not summed: `df.groupby(["source","target"])["weights"].mean()` |

### 2.1 The original is *not* a sample-specific path

This is the single most important finding for the new feature. All link values come from
**static model parameters**:

```python
geneset_weights = …cellpathway_aggregator.aggregator.attention_weights   # frozen softmax params
gene_weights    = …geneset_aggregator.aggregator.attention_weights       # frozen softmax params
```

re-softmaxed with `F.softmax(..., dim=0)` and optionally scaled by `concept_imp_dict`
(default: all ones). **No expression value, no per-sample quantity and no Cox coefficient
enters the original computation.** The original renders the model's static architecture for a
chosen concept — it cannot be reused as-is for a sample-specific view.

The new feature therefore keeps the original's *layered flow idea* and its **layer semantics**,
but derives every value from an exact per-sample decomposition (§4).

---

## 3. Is an exact per-sample decomposition available?

Yes. The frozen graph is:

```
expression → [frozen scaler] → inputencoder → gene encoding  e_{i,g} ∈ R^32   (15,672 × 32)
             → gene-level score        gene_score_g(i) = w · e_{i,g} + b     ← official API
             → gene-set attention      geneset_score_j(i) = Σ_{g∈G_j} a_{g,j} · gene_score_g(i)
             → concept attention       concept_c(i)     = Σ_{j∈C_c} b_{j,c} · geneset_score_j(i)
             → frozen Cox head         η_i = Σ_k β_k · z_k(i)
```

`gene_score` is an **official** output: `PreTrainer.extract(..., with_gene_level=True)` returns
`(dfg, dfgs, dfct)`, where `dfg` is exactly
`genesetprojector.geneset_scorer(encoding)[:, 2:]` (`compass/model/tune.py:289-291`). No hooks,
no monkey-patching and no re-implementation of the encoder are required.

### 3.1 Measured exactness (12 golden samples, reference strategy)

| identity | measured `max|Δ|` | verdict |
|---|---|---|
| re-implemented `scorer(aggregator(enc[:,2:]))` vs official `dfgs` | `1.8e-07` | exact (float32) |
| manual attention aggregation vs model aggregator | `3.6e-07` | exact (float32) |
| manual gene-level vs official `dfg` | `4.8e-07` | exact (float32) |
| **`Σ_g a_{g,j}·gene_score_g` vs official `dfgs` (gene→gene-set)** | **`1.8e-07`** | **exact** |
| **`Σ_j b_{j,c}·geneset_score_j` vs official `dfct` (gene-set→concept)** | **`3.3e-08`** | **exact** |
| gene→concept end-to-end | `1.8e-07` | exact |

The identity `geneset_score_j = Σ_g a_{g,j} · gene_score_g` holds because the gene-set score is
an attention-weighted sum of gene embeddings passed through **one shared linear scorer**, and
the attention weights sum to 1 (measured: `0.99999990 … 1.00000009`).

**Trap found and handled:** 2 of the 132 sets list a gene index **twice**
(e.g. `geneset_21`: 15 entries, 14 unique). A naive `A[j, idxs] = w` assignment silently drops
the duplicate and produced a `6.7e-02` error. The correct construction is a **scatter-add**
(`np.add.at`). With scatter-add the identity is exact. Any implementation must do the same.

### 3.2 Measured exactness of the Cox decomposition

Using the **package's own** `build_design_matrix()` + `loader` beta, with
`z_scaled = (x − scaler_mean) / scaler_scale` applied exactly to `lock["scaled_features"]`:

| model | design columns | `Σ(per-feature contribution)` vs main-API `risk` | groups |
|---|---|---|---|
| M1 | 43 | **`2.2e-16`** | 43 concepts |
| M2 | 78 | **`4.4e-16`** | 32 `CT_*` + Age + Sex + Stage + 43 concepts |
| M3 | 88 | **`8.9e-16`** | 32 `CT_*` + Age + Sex + Stage + **PC1–PC10** + 43 concepts |

The decomposition of `η_i` is therefore **exact at machine precision**, not an approximation.

Frozen-head facts that the implementation must respect:

* `scaled_features` = Age + 43 concepts (M2: 44 entries; M3: 54 = 44 + PC1–PC10).
  **Sex and Stage are NOT scaled** — they enter the linear predictor raw.
* `cancer_onehot_cols` = 32 `CT_*` columns; the reference level is dropped, so a sample whose
  cancer type is the reference (ACC) has **all** `CT_* = 0`.
* M1 has **no** clinical, cancer or PC columns — concepts only.

---

## 4. Answers to the twelve audit questions

**Q1 — What does each original layer correspond to?**
Gene (15,672 symbols) → gene-set / granular representation (`genesets_names`, 132) → high-level
concept (`cellpathway_names`, 43). No risk layer. See §1.

**Q2 — What is the original top-K behaviour?**
`topK_vis=10`, applies to genes only, and only blanks labels; nodes and links stay.
`concept2plot` (default one concept) is the real filter. No per-layer top-K for gene-sets or
concepts. See §2.

**Q3 — Can top-50 per layer for the first three columns be implemented stably?**
Yes. The new implementation builds nodes and links itself from the exact decomposition, so
truncation is explicit and testable (`top_genes`, `top_signatures`, `top_concepts`). The
original's label-only truncation is **not** reused, because it cannot bound the node count.

**Q4 — Is the ≤16-node high-level layer implementable?**
Yes. The budget is enforced on the *rendered* high-level column
(`concepts_shown + auxiliary_nodes ≤ 16`), computed from the model's own predictor list.

**Q5 — How are exact Age/Sex/Stage Cox contributions obtained?**
Directly from the frozen lock:
`contrib_Age = β_Age · (age − mean_Age)/scale_Age`;
`contrib_Sex = β_Sex · sex`; `contrib_Stage = β_Stage · stage` (Sex/Stage are unscaled).
Imputation of missing clinical values is taken from the main API's frozen-reference fill, so
the visualized value equals the value used by the risk computation. Verified: the sum of all
contributions reproduces the main-API risk to `4.4e-16`.

**Q6 — How is cancer type handled?**
The 32 `CT_*` columns genuinely enter the M2/M3 Cox head and may not be pretended away.
Within a **single-cancer cohort** the term is a per-sample constant
(`β_CTk` for every sample, or `0` for the reference level ACC) and therefore cannot change
cohort-internal ranking → omitted from the diagram by default, with a mandatory caption
stating that cancer type is part of the frozen model but constant in this cohort.
For a **multi-cancer cohort** it is shown as one predictor node connected directly to risk,
occupying one unit of the 16-node budget. (Note the term is *not* strictly proportional to
"number of cancer types"; it is the sample's own one-hot contribution.)

**Q7 — How are PC1–PC10 obtained and aggregated exactly?**
From the locked M3 PCA (`locked_pca_M3.npz`, already applied by `build_design_matrix`) and the
frozen β:
`C_PC(i) = Σ_{k=1..10} β_PCk · z_PCk(i)`, where `z_PCk` is the **scaled** PC
(PC1–PC10 are inside `scaled_features`). This is a sum of ten real contributions, **not** a
mean. Rendering merges the ten into one node; the ten per-PC contributions are retained in the
returned data frame and exposed in hover text. Verified exact as part of the M3 row above.

**Q8 — How are signed Cox contributions mapped to Sankey width/colour?**
Plotly requires `link.value ≥ 0`, so:
* **width** = `|contribution|`
* **colour** = sign — one hue family for *positive contribution (higher model-estimated risk)*,
  another for *negative contribution (lower)*, with the direction stated in the legend.
Negative values are never encoded as negative width. The same rule is applied to the
representation layers (`a_{g,j}·gene_score_g`, `b_{j,c}·geneset_score_j`), whose values are also
signed, because `gene_score` is an unbounded linear output.

**Q9 — Can the full risk reuse the main result instead of re-predicting?**
Yes, and it must. The implementation calls the **same** frozen functions the main API uses
(`survival.build_design_matrix` + `risk_from_design`), and both:
1. accepts an existing `PredictionResult` / `AnalysisResult` and reuses its risk when supplied;
2. asserts the reconstructed `η` equals the main-API risk (`max|Δ| ≤ 1e-9`) and surfaces the
   measured deviation in the returned result.
A regression test asserts equality against `compass_os.predict()` at machine precision.

**Q10 — Can Plotly render the Sankey and a continuous risk bar in one figure?**
Yes, with one important trap. A Sankey trace has **no x/y axes**, so axis-spanning helpers
(`add_hline`, `add_vline`) raise `PlotlyKeyError: Invalid property … 'xaxis'`. Verified
workaround: build the risk bar in a second subplot via
`make_subplots(specs=[[{"type":"sankey"},{"type":"xy"}]])` and add the cutoff line through
`fig.add_shape(..., xref="x2", yref="y2")` with explicit axis references.
This was spiked successfully (Sankey + gradient risk column + sample triangle marker + dashed
cutoff, single HTML figure).

**Q11 — Does the main API need modification?**
No. The feature is purely additive: a new module plus an optional public entry point. No change
to `api.py` risk logic, `analysis.py`, the design matrix, the frozen thresholds or the models.
`plotly` becomes an optional plotting dependency; the core inference dependency set is
unchanged.

**Q12 — Are frozen scientific assets touched?**
No. The checkpoint and the locked Cox JSONs are opened **read-only**; attention weights, betas
and scalers are read, never written. `ASSET_MANIFEST` hashes are unaffected. This feature adds
no asset.

---

## 5. Design decisions fixed by this audit

1. **Two extra layers are added relative to the original**: a clinical/auxiliary predictor layer
   and a terminal cohort-relative risk column. The first three columns keep the original
   semantics (gene → gene-set/granular representation → concept).
2. **Column 1 is gene expression (TPM by default)**, with node colour = the sample's value on
   the declared input scale. Selection of *which* genes to show is a **separate** rule from
   colour (per the brief): genes are ranked by their exact aggregate contribution to the
   concepts displayed, `|Σ_c β_c/scale_c · Σ_{j∈C_c} b_{j,c} · a_{g,j} · gene_score_g|`, and
   not by TPM.
3. **Column 2** is ranked by the sample's own gene-set contribution to the displayed concepts.
4. **Column 3** ranks the 132 granular representations by `|contribution|` the same way.
5. **High-level column budget**: `n_concepts_shown + n_auxiliary_shown ≤ 16`, with auxiliary
   nodes (Age/Sex/Stage/PC-aggregate/cancer type) placed first in the budget when present.
6. **Concept ranking uses `|x_j β_j|`**, i.e. the sample's actual Cox contribution — never the
   raw concept score.
7. **Completeness option chosen: "Other concepts" aggregate node** (option A of §14 in the
   brief) *and* an explicit note. An aggregate node makes the visible inflow to the risk node
   sum to the full linear predictor, which is less misleading than showing a partial set with a
   footnote alone. Its contribution is `Σ_{c ∉ shown} β_c z_c(i)`, computed exactly.
8. **Two cutoffs are distinguished**: `cutoff="median"` (cohort median of the current risk
   output, default) is a visualization/grouping convenience only; the frozen validated cutoff
   (`median_cutoff` in the lock) is a separate explicit mode and is **never** conflated with it.
9. **Cancer type**: omitted by default for single-cancer cohorts (with mandatory caption),
   shown for multi-cancer cohorts.
10. **M0 is rejected** with the exact message required by the brief.
11. **Static export** needs `kaleido`, which is not installed and cannot be installed here
    (no PyPI access). It is therefore an **optional extra** (`compass-os[plot]`), and the
    implementation additionally offers a dependency-free matplotlib renderer so a static figure
    is still producible offline.

---

## 6. Verdict

```text
GO
```

Every required quantity is available **exactly** from frozen, read-only artefacts and the
official extraction API:

* gene → gene-set → concept: exact (`1.8e-07` / `3.3e-08`, float32 precision);
* concept/clinical/cancer/PC → risk: exact (`2.2e-16` … `8.9e-16`);
* full risk reusable from the main API without duplicating prediction logic;
* Sankey + continuous risk column renderable in one Plotly figure (trap identified and
  workaround verified);
* no frozen asset modified, no risk algorithm change, no retraining, no approximation
  presented as a Cox contribution.

None of the STOP conditions apply.

### Residual limitations to be carried into the validation report

* Static PNG/SVG export requires the optional `kaleido` package (unavailable offline here).
* The reference level of the cancer one-hot (ACC) contributes exactly 0; for a single-cohort
  study all samples share that same term, so the visualization cannot show cancer-type
  variation within one cancer type — by construction, not by omission.
* Gene-level attribution is exact with respect to the frozen graph, but a gene's contribution
  may be negative (the shared scorer is an unbounded linear map); the diagram encodes magnitude
  as width and sign as colour and must not be read as "importance" in a biological sense.
