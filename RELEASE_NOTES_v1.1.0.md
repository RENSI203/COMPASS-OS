# RELEASE_NOTES_v1.1.0.md

**Version**: `1.1.0` ｜ **Status**: **prepared, not tagged** ｜ **Date**: to be set at release

> This file is written in advance of tagging. The tag, the release commit and the release
> date are **not** set yet; nothing here should be cited as an existing release.

---

## What this release adds

### New: sample-level computation-path visualization (`compass_os.sample_path`)

For one sample of an analysed cohort, draw its computational attribution through the frozen
graph and its position in the cohort risk distribution:

```
Gene expression → COMPASS gene score → granular signature score
    → high-level concepts + clinical / PC → cohort-relative risk
```

* **Real scores, not weights.** The gene-score column is the official
  `extract(with_gene_level=True)` gene-token scalar
  (`geneset_scorer(encoder_output)[:, 2:]`, a shared `nn.Linear(32→1)`); the 132 signature and
  43 concept columns come from the **same forward pass as η**. Attention weights decide only
  *which connections* are drawn, never node values.
* **Clinical / cancer type / PC1–PC10** are aggregated from the frozen Cox design matrix
  (multi-column terms use the signed sum; `PC1–PC10` is one node). Fields that were actually
  imputed are labelled `(imputed)`.
* **Complete accounting.** The full untruncated decomposition — including undisplayed
  concepts and the hidden cancer term — satisfies `sum(contribution) == η`; the renderer
  enforces it.
* **One call, unified export.** `sample_path(...)` performs input checks, inference, score
  extraction, decomposition, selection and export. PNG / PDF / SVG come from a single
  Matplotlib figure; the HTML embeds the **identical PNG bytes** plus an expandable exact-node
  table. There is no second renderer and no silent fallback.

Display budgets (frozen): M1 16 concepts; M2 single-cancer 13 + Age/Sex/Stage; M2
multi-cancer 12 + clinical + Cancer type; M3 single-cancer 12 + clinical + one `PC1–PC10`;
M3 multi-cancer 11 + clinical + Cancer type + `PC1–PC10`. Each high-dimensional column ≤ 50.

### Changed

* The sample-path drawing layer is now the package's unified renderer. The earlier Plotly
  Sankey and its Matplotlib fallback are **removed**, so there are no longer two divergent
  appearances.
* `SamplePathResult.figure()` now returns a `RenderedPath`; `save_html()` / `save_static()` /
  `save()` are the export entry points. `save_html()` works standalone.
* Legacy keyword aliases `top_genes` / `top_signatures` / `max_high_level_nodes` are still
  accepted; `top_concepts` is deprecated.

### Fixed

Seven defects were found by the pre-release audit and fixed in this release line:

| # | Defect | Effect |
|---|---|---|
| 1 | `±Inf` passed input validation silently | surfaced later as an unrelated sklearn error |
| 2 | `cancer_type` length mismatch raised a bare pandas `ValueError` | the package's own length check never ran |
| 3 | `cancer_type="NORMAL"` raised a bare `IndexError` | and on the one-hot path silently equalled the ACC reference level |
| 4 | the input-scale sanity note ran only inside `analyze()` | `predict()` accepted a wrong `input_scale` silently |
| 5 | non-numeric clinical values (e.g. `Female`/`Male`) were reported only as "missing" | an encoding error looked like absent data |
| 6 | the robustness runner resumed from a cache keyed without `--seed` / `--max-per-cohort` | two definitions could be mixed into one result file |
| 7 | `examples/cohort_example.py` crashed on its own demo data | `risk_group` is `None` below 30 samples |
| 8 | **`torch` and `seaborn` were undeclared runtime dependencies** | `pip install compass-os` succeeded but the first `predict()` raised `ModuleNotFoundError` |

Defect 8 is the one that matters most to users: it was invisible to `--no-deps` installs and is
why this release requires a dependency-resolving install.

---

## Installation

```bash
git clone https://github.com/RENSI203/COMPASS-OS.git
cd COMPASS-OS
pip install .          # resolves PyTorch, seaborn, scipy and the rest automatically
```

Requires Python ≥ 3.10. **PyTorch is a core runtime dependency**, not an extra.

---

## Verification performed for this release

| Check | Result |
|---|---|
| Package test suite | 152 / 152 pass (build first, then test) |
| Renderer self-checks | 12 / 12 pass |
| Release hygiene audit | FAIL 0 / WARN 0 |
| Frozen model assets | 13 / 13 unchanged; byte-identical to v1.0.1 |
| Wheel reproducibility | byte-identical across rebuilds at the pinned epoch |
| sdist acceptance | 162 entries, all assets present, zero bytecode; builds, installs and runs |
| Installed-package end-to-end | core prediction, `sample_path`, asset resolution and HTML/PNG/PDF/SVG export verified outside the repository |
| Reproducibility of published numbers | 132 signatures `1.1e-16`, 43 concepts `2.4e-07`, M0–M3 risks `≤1.06e-06` (criterion 1e-5) |

---

## Known limitations (please read before citing)

1. **No significant increment over the clinical baseline was detected.** In the
   clinical-baseline-evaluable primary analysis set (k = 24), paired DL ΔUnoC(M2−M0) =
   **+0.008075**, 95 % CI **[−0.003864, +0.020015]**, p = 0.185. No equivalence margin was
   pre-specified and no equivalence test was performed, so this is **not** evidence of
   equivalence — only of no detected increment. A cohort-set sensitivity exists and is
   disclosed: over all 39 cohorts the same contrast is +0.015619, CI [+0.004217, +0.027022].
   The 24-cohort set is used because M0 is degenerate (constant risk) in the excluded cohorts,
   which makes the contrast uninterpretable there — **not** because it is more conservative.
2. **M1 (concepts only) is significantly worse than the clinical baseline** (−0.0495,
   CI [−0.0829, −0.0160]).
3. **Risk scores order samples; they are not calibrated absolute risks.** `S(t)` is a
   baseline-hazard-derived output; no calibration evidence (calibration curves, Brier, ICI) is
   included. Discrimination must not be presented as absolute-risk calibration.
4. **Microarray input requires external harmonisation.** There is no
   `input_scale="microarray"`. Only `"tpm"` and `"log2_tpm1"` are accepted.
5. **A non-negative counts matrix is indistinguishable from TPM** and will be accepted
   silently. Declaring the scale correctly is the caller's responsibility.
6. **Robustness to missing genes holds only within the tested range** — do not describe the
   model as invariant to missing genes.
7. **QC coverage tiers are empirically calibrated descriptive engineering tiers**, not
   biological boundaries or safety lines. The per-signature threshold is **not calibrated**
   (`null` in `qc_config.json`).
8. **Representation-layer sensitivity differs by strategy**: the risk ranking is nearly
   insensitive to `reference` vs `zero` (Δ ≤ 0.0003), but the representation layer is not —
   `zero` shows 52 %–76 % larger mean absolute deviation. The two strategies are **not**
   interchangeable for representation readouts.
9. **43 concepts / 132 signatures / gene scores are representation coordinates and internal
   model scores.** They do not establish pathway activation, suppression, or causal mechanism.
10. **Display thresholds and colour ranges are display parameters**, not validated scientific
    or clinical thresholds. Figure truncation does not change predictions, but the figure is
    not the complete model.
11. **No claim is made about immunotherapy response.** The upstream COMPASS publication
    concerns immunotherapy; COMPASS-OS is a frozen-asset inference and engineering package and
    makes no such claim.
12. **The frozen median cutoff need not split a new cohort 50:50.** An unequal split is
    expected when the risk distribution shifts; it is not by itself evidence that the cutoff
    failed. What requires dedicated verification is the post-stratification prognostic
    separation and the applicable range.

---

## Attribution

COMPASS-OS is a frozen-asset inference package. The pretrained COMPASS model, the 132 gene-set
and 43 concept attention topology and the gene encoder come from the upstream COMPASS work
(Shen, W., Moon, I., Nguyen, T.H. et al. *Nat Med* **32**, 3010–3022 (2026),
<https://doi.org/10.1038/s41591-026-04502-7>), vendored under its MIT license
(`assets/third_party/COMPASS_LICENSE`). This package adds the frozen inference wrapper, the
Cox heads M0–M3 and locked PCA, the missing-gene strategies, coverage QC, the sample-path
attribution and renderer, and the packaging and release auditing. No training, re-fitting or
tuning happens at prediction time.
