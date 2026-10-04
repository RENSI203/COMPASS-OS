# Sample path redesign — real-data wiring report

Adopted the uploaded `COMPASS_sample_path_redesign` renderer and wired it to **real
COMPASS-OS inference**. No merge, no tag, no push. **No `src/` or `tests/` change.**

---

## 1. What was adopted

| item | location | integrity |
|---|---|---|
| Renderer | `tools/redesign_render/sankey.py` | sha256 `daa89531…cfae6b5` — **byte-identical** to the package |
| Verifier | `tools/redesign_render/verify_renderer.py` | sha256 `72f7a244…f47a13` — byte-identical |
| Example driver | `tools/redesign_render/sample_path_example.py` | verbatim |
| Adapter (mine) | `tools/redesign_render/build_real_inputs.py` | real inference → `SamplePathData` |
| Grid tool (mine) | `tools/redesign_render/build_redesign_grid.py` | 2×3 review sheet |

Both package files match `CHECKSUMS.sha256` from the zip exactly.

`verify_renderer.py`: **12/12 PASS**.

---

## 2. Real-data wiring (field by field)

Cohort **GSE39582**, COAD, **n = 573**, harmonised cache, `input_scale="log2_tpm1"`.
Same three samples as the previous review round, selected programmatically on the **M2**
percentile nearest 10 / 50 / 90; M2 and M3 share samples, cohort and one Style.

| `SamplePathData` field | real source |
|---|---|
| `risk` | **main API** `compass_os.predict(...).risk[model]` — Cox linear predictor η |
| `rank`, `percentile` | main API cohort ordering; definition written into `rank_definition` |
| `cohort_risks` | all 573 η values from the same prediction call |
| **`gene_score`** | **`PreTrainer.extract(..., with_gene_level=True)` → `dfg`** = `genesetprojector.geneset_scorer(encoder_output)[:, 2:]`, the real gene-token score (shared `nn.Linear(32→1)`). **Not** attention weight, **not** TPM, **not** risk attribution |
| `gene_tpm` | `align_expression(...).matrix` — the TPM-scale matrix the model consumes; units labelled `TPM (harmonised cache, model input scale)` |
| `signatures` | real 132 gene-set/signature scores (`dfgs`) |
| concept predictors | real 43 concept scores for `score`; exact signed Cox contribution for `contribution`. Concept `score` and `contribution` both come from the **main-API** concept scores so colour and contribution are mutually consistent |
| clinical predictors | `Age` / `Sex` / `Stage`, `score == contribution` = exact signed contribution from the frozen design matrix (`Age` standardised, `Sex`/`Stage` raw) |
| `PC1–PC10` | `aggregate_cox_predictor("PC1–PC10", [PC1..PC10], {col: contribution})` — exact signed sum of the ten **frozen** design-column contributions |
| `imputed` | per field. This cohort supplied **no clinical data**, so Age/Sex/Stage are **all** frozen-imputed (`value` = 60.0 / 0.0 / 2.0) and are marked `imputed=True`, rendering as **`Age (imputed)`** etc. |
| edges | frozen topology: 1,281 gene→signature and 132 signature→concept edges, weights = frozen attention (used **only** for connection filtering) |
| `decomposition` | full untruncated signed decomposition: **78 rows (M2) / 88 rows (M3)** |
| `cancer_types` | cohort list `["COAD"]` → single-cancer ⇒ Cancer type hidden, its slot returned to concepts |

### One wiring decision worth recording

The decomposition is built from the **main API's** `concept_scores`, not from the adapter's own
single-sample extraction. Rationale: the main API computes η from exactly those concept scores,
so the decomposition is then **exactly** self-consistent with `risk`. Using the adapter's
separate float32 forward pass made `sum(decomposition)` differ from η by ~3.8e-07, which the
renderer's `sum == η` guard correctly rejected (it failed loudly on `M3_low`, η = 0.058, where
the relative tolerance is tightest). After the change the residual is **0.00e+00 … 2.22e-16**.

---

## 3. Verification results

| check | result |
|---|---|
| `verify_renderer.py` | **12/12 PASS** |
| decomposition `Σ contribution` vs η | **Δ = 0.00e+00 … 2.22e-16** (≤ 2.2e-16) |
| figure `risk` vs main API, all 6 | **exactly equal** |
| figure `rank` vs main API, all 6 | **exactly equal** |
| HTML embeds the same PNG bytes | **True** for all 6 |
| package regression | **125 / 125 PASS** |
| frozen asset hashes | **13/13 unchanged** |
| `src/` / `tests/` modified | **none** |

### Node budgets actually rendered

| figure | gene_tpm | gene_score | signatures | predictors |
|---|---|---|---|---|
| low / M2 | 32 | **32** | 26 | 16 |
| low / M3 | 28 | **28** | 26 | 16 |
| mid / M2 | 26 | **26** | 26 | 16 |
| mid / M3 | 29 | **29** | 26 | 16 |
| high / M2 | 32 | **32** | 26 | 16 |
| high / M3 | 30 | **30** | 26 | 16 |

The two gene columns are always equal in length **and** name-for-name identical, as required.

---

## 4. Shared Style (real distributions, one file for M2/M3)

`docs/figures/sample_path_redesign/shared_style.json`

| parameter | value | basis |
|---|---|---|
| `gene_threshold` | 0.25 | real gene-token score scale (median \|score\| = 0.40) — a **display** threshold, not a validated scientific one |
| `max_genes` / `max_signatures` | 32 / 26 | package defaults, ≤ 50 |
| `predictor_budget` | 16 | unchanged |
| `expression_limits` | (0, 1704.3) | cohort TPM p99.5 — the package default (0, 20) would have saturated everything |
| `gene_score_limits` | (−2, 2) | measured \|gene score\| max ≈ 1.98 |
| `signature_limits` / `concept_limits` | (−1, 1) | measured ±0.82 |
| `contribution_limits` | (−1.5, 1.5) | M2 ±0.93 / M3 ±1.50 (PC aggregate) |

---

## 5. Deliverables

```
docs/figures/sample_path_redesign/
├── low/{M2,M3}.{png,pdf,svg,html,nodes.tsv,selection.json,decomposition.tsv}
├── mid/{M2,M3}.{…}
├── high/{M2,M3}.{…}
├── M2_{low,mid,high}.json  M3_{low,mid,high}.json   ← real SamplePathData payloads
├── shared_style.json  _selection.json
└── sample_path_redesign_review_grid.png             ← 2×3 sheet, composited (no re-plotting)
```

`docs/figures/sample_path_redesign/` replaces the earlier hand-rolled Sankey for review
purposes; the old implementation and its review material are untouched under
`docs/figures/sample_path_review/`.

---

## 6. What the redesign fixes (vs the previous round)

* Separate **Gene expression** and **COMPASS gene score** columns — the structural gap reported
  last round is closed (question 1 of the previous review).
* Genuine readability: circular nodes, uniform grey links, layer colour bars with explicit
  limits, and legible per-layer node labels.
* Risk column is **visible** and correct: cohort-percentile axis, red→yellow gradient, triangle
  marker positioned from the **main API percentile**, dashed cutoff at its true empirical
  percentile, and `Risk / Rank / Percentile` labels.
* Residual/bookkeeping nodes are gone from the main figure; hidden contributions remain in the
  full decomposition table.
* No Plotly / kaleido / CDN dependency: one Matplotlib figure renders PNG, PDF and SVG, and the
  HTML embeds the identical PNG bytes.

## 7. Open items (for the user to decide, not changed here)

1. `docs/figures/sample_path_redesign/` is **~15 MB** and currently untracked; decide whether it
   ships in git.
2. The **old** `src/compass_os/sample_path.py` renderer is still in place and still exported.
   Adopting the redesign as the package's renderer (and retiring `_build_figure`) is a code
   change that has **not** been made.
3. `gene_threshold = 0.25` is a display parameter; if this becomes a published figure it should
   be justified on the real gene-score distribution rather than inherited.
4. Clinical values here are **100 % imputed** (the cohort has no clinical data), so the
   Age/Sex/Stage contributions are identical for every sample and cannot illustrate variation.
