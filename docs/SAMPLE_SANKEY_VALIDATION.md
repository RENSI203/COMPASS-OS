# Sample-level computation path — validation report

Feature: `compass_os.sample_path` (v1.1.0).
Companion documents: [`SAMPLE_SANKEY_FEASIBILITY_AUDIT.md`](SAMPLE_SANKEY_FEASIBILITY_AUDIT.md)
(design + GO decision), [`SAMPLE_PATH.md`](SAMPLE_PATH.md) (usage).

**Verdict: READY WITH CONDITIONS** — see §7.

---

## 1. What was delivered

| item | path |
|---|---|
| Feasibility audit | `docs/SAMPLE_SANKEY_FEASIBILITY_AUDIT.md` |
| Implementation | `src/compass_os/sample_path.py` |
| Public export | `compass_os.sample_path`, `compass_os.SamplePathResult` |
| Tests | `tests/test_sample_path.py` — 26 cases |
| Usage guide | `docs/SAMPLE_PATH.md` |
| Examples | `examples/sample_path_example.py` (M1/M2/M3) |
| Interactive HTML | `examples/output/sample_path/sample_path_{M1,M2,M3}.html` |
| Static figures | `examples/output/sample_path/sample_path_{M1,M2,M3}.png` |
| Full decomposition tables | `examples/output/sample_path/<sample>.{M1,M2,M3}.decomposition.tsv` |

No frozen asset was modified; no risk algorithm was changed; the main API is untouched.

---

## 2. Exactness of the attribution chain

Measured on the 12-sample golden fixture (reference strategy), comparing the
re-implementation against the **official** extraction outputs:

| identity | measured `max|Δ|` |
|---|---|
| `Σ_g a_{g,j}·gene_score_g` vs official gene-set score | `1.8e-07` |
| `Σ_j b_{j,c}·gene-set score_j` vs official concept score | `3.3e-08` |
| gene → concept end-to-end | `1.8e-07` |

| model | design columns | `Σ(per-feature contribution)` vs main-API `risk` |
|---|---|---|
| M1 | 43 | `2.2e-16` |
| M2 | 78 | `4.4e-16` |
| M3 | 88 | `8.9e-16` |

All differences are float32 round-off, not modelling approximations.

### 2.1 Duplicate gene indices (trap)

2 of the 132 gene sets list a gene index twice (`geneset_21`: 15 entries / 14 unique).
Naive assignment produced `6.7e-02` error; **scatter-add** (`np.add.at`) makes the identity
exact. Guarded by `test_gene_and_signature_links_are_exact_contributions` and by the
per-node conservation test.

---

## 3. Flow conservation

Every internal node is checked at runtime; a violation raises `CompassOSError` instead of
drawing a misleading figure.

| node class | inflow | outflow |
|---|---|---|
| displayed gene | — (source) | Σ contributions to kept gene-sets |
| **Other genes (not shown)** | — (source) | residual = `L_j·gene-set_j − shown` |
| gene-set *j* | shown genes + residual | to displayed concepts |
| **Other signatures (not shown)** | — (source) | residual to displayed concepts |
| concept *c* | kept signatures + residual + centering | centred Cox contribution |
| **Frozen model centering (constant)** | — (source) | `−β_c·mean_c/scale_c` per concept |
| **Other concepts (aggregated)** | — (source) | `Σ_{c∉shown} β_c z_c` |
| auxiliary predictors | — (sources, by design: not connected upstream) | their exact Cox contribution |
| **risk** | all of the above | **= full Cox linear predictor** |

Design rationale for the sources with no inflow: the brief requires that clinical / PCA / cancer
nodes are **not** connected to the COMPASS representation layers (they are additional Cox
predictors, not products of the representation). The same convention is applied to the residual
aggregates so that no node silently "gains" flow.

Verified by `test_every_internal_node_conserves_flow` and
`test_risk_node_inflow_equals_linear_predictor`.

---

## 4. Risk, rank and percentile are the official ones

* `res.risk` is `float(predict(...).risk[model].loc[sample_id])` — **exact equality**
  (`test_full_risk_matches_main_api`).
* `res.rank`, `res.percentile` derive from the official cohort risk ordering
  (`test_rank_and_percentile_are_from_official_risk`).
* The attribution chain independently sums to η; the deviation from the API value is reported as
  `eta_max_abs_diff` and must stay within a relative `1e-4` or the call fails.
  Observed: `9.2e-06` (M1), `1.3e-06` (M2), `1.2e-06` (M3) — float32 forward-pass round-off
  between two independent passes, not a decomposition error.

---

## 5. Required test matrix (§19 of the brief)

| # | case | test |
|---|---|---|
| 1–3 | M1 / M2 / M3 single sample | `test_single_sample_each_model` |
| 4 | sample_id not present | `test_missing_sample_id_raises` |
| 5 | single-sample cohort | `test_single_sample_cohort` |
| 6 | small cohort | `test_small_cohort` |
| 7 | M2 complete clinical | `test_clinical_complete_and_imputed` |
| 8 | M2 missing clinical → frozen imputation | `test_clinical_complete_and_imputed` |
| 9 | M3 PC contribution aggregation | `test_pc_aggregation_is_sum_not_mean` |
| 10 | median cutoff | `test_cutoff_median_and_numeric` |
| 11 | numeric cutoff | `test_cutoff_median_and_numeric` |
| 12 | top-50 truncation | `test_topk_truncation` |
| 13 | high-level ≤ 16 | `test_high_level_node_budget` |
| 14 | multi-cancer cancer-type behaviour | `test_cancer_type_multi_cohort_shown` |
| 15 | single-cancer omission | `test_cancer_type_single_cohort_omitted_with_note` |
| 16 | **full risk identical to main API** | `test_full_risk_matches_main_api`, `test_attribution_chain_sums_to_api_risk` |
| 17 | rank / percentile identical | `test_rank_and_percentile_are_from_official_risk` |

Extras: M0 rejection with the mandated message, invalid model, result reuse without
re-prediction, HTML/static output, figure contains Sankey **and** risk column with non-negative
link values, semantics disclaimer present, flow conservation, risk-node inflow.

**Result: 26/26 pass.** Full suite: **125/125 pass** (no regression in the pre-existing 99).

---

## 6. Figures produced in this environment

| model | HTML | static PNG |
|---|---|---|
| M1 | 44 KB | 333 KB |
| M2 | 45 KB | 331 KB |
| M3 | 45 KB | 322 KB |

The HTML is the native Plotly rendering (Sankey + risk column + sample marker). The PNG was
produced by the **matplotlib fallback**, because `kaleido` is unavailable here and PyPI is
unreachable from this machine (see §7).

---

## 7. Remaining limitations (conditions)

1. **Static export depends on an optional package.** Plotly's `write_image` needs `kaleido`,
   which is not installed in this environment and cannot be installed offline. It is declared as
   the optional extra `compass-os[plot]` (`plotly>=5.18`, `kaleido>=0.2.1`).
   A dependency-free matplotlib fallback ships, so a static figure is always producible, but it
   is a *different rendering* (layered contribution bars), not a raster of the Sankey.
2. **Two float32 forward passes.** The attribution chain and the cohort prediction are separate
   passes, so η differs by ≈`1e-06`–`1e-05`. The figure displays the main-API risk; the deviation
   is reported and bounded (`ETA_TOL_REL = 1e-4`). If bit-identical attribution is ever required,
   the cohort would have to be extracted once with `with_gene_level=True`.
3. **Truncated layers carry residuals.** With `top_genes/top_signatures` at 50, the displayed
   nodes represent ≈63–66 % of the gene-set inflow in absolute terms; the rest is carried by
   explicit "not shown" nodes. The figure is conservative, but individual gene links are a
   *subset* of the true gene→gene-set flow.
4. **Cancer-type term in single-cancer cohorts** is constant by construction and therefore
   omitted (with a caption). The visualization cannot show within-cohort cancer-type variation
   for a single-cancer cohort — that is a property of the model, not a display choice.
5. **No causal reading.** Contributions are properties of the frozen computation, and a negative
   contribution means "pushes the linear predictor down", not "protective biology".
6. **`M0` unsupported** by design (no COMPASS branch); the mandated message is raised.

None of these conditions affects the correctness of the reported risk, rank, percentile or
decomposition.

---

## 8. Verdict

```text
READY WITH CONDITIONS
```

The feature is numerically exact, flow-conserving, fully tested (26 new tests, 125/125 total),
and correctly bounded in its scientific claims. The conditions in §7 are documentation- and
dependency-level only; none of them silently misrepresents the model.
