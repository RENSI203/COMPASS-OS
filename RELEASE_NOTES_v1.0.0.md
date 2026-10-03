# COMPASS-OS v1.0.0

**Pan-cancer survival prediction and mechanism-oriented transcriptomic representation
framework.**

This is the first public release: a frozen-asset inference package that turns a bulk
expression matrix into prognostic risk, cohort stratification, a 43-concept COMPASS
representation, a 132-gene-signature representation and input-quality QC.

## Core capabilities

* **One-click analysis** — `compass_os.analyze(expression, cancer_type, clinical=None,
  survival=None)` returns a result object with `summary()` and `save_report()`.
* **Frozen prognostic models** — `M0` (clinical baseline), `M1` (43 concepts),
  **`M2` (primary / default)**, `M3` (M2 + expression PCs). No training, re-fitting or
  tuning at prediction time.
* **43 concept representations** and **132 gene-signature representations**, reported as
  mechanism-oriented, hypothesis-generating features.
* **Missing-gene QC** — gene coverage, per-signature coverage (*n* × 132) and per-concept
  input coverage (*n* × 43), with calibrated tiers.
* **Three missing-gene strategies** — `reference` (default), `zero`
  (training-space masking), `strict`.
* **Robustness assessment** — `check_robustness()`, or `robustness="auto"` which only runs
  the extra comparison when coverage falls below the recommended tier.
* **Controlled advanced API** — `predict`, `get_representation`, `check_robustness` with a
  limited, documented parameter set (model, strategy, input scale, batch/device, cohort-size
  threshold). All model internals are read-only frozen assets with no public setter.

## Validation

* **Golden reproducibility against the original model-development pipeline:** 132
  signatures `max|Δ| = 1.1e-16`, 43 concepts `max|Δ| = 2.4e-07`, M0–M3 risks
  `max|Δ| ≤ 1.06e-06` (criterion 1e-5).
* **External gene coverage** (38 restorable public cohorts, pre-imputation): global
  coverage median 96.7 % (63.5–98.0 %), signature-gene coverage median 96.2 %
  (80.3–97.8 %).
* **Missing-gene perturbation validation** (global / signature-only / non-signature-only
  masking, 12 cohorts × 20 repeats each, plus replay of real platform masks in 38
  cohorts): across the tested masking range, changes in Uno's C remained within
  approximately ±0.02, with no consistent deterioration across the evaluated cohorts.
  Preservation of signature-related genes accounted for a substantial component of
  prognostic robustness (50 % of the remaining 14,756 genes masked → risk Spearman
  ≈ 0.9998).
* **QC tiers** derived from the measured curves: recommended 0.90 / warning 0.70 for both
  the global and the signature-gene axis; overall grade = worse of the two.

## Input

* Bulk expression, `samples × genes`, gene-symbol columns, in **TPM** or **log2(TPM+1)**
  (declared by the caller; never auto-detected).
* Cancer type required (TCGA acronym).
* Optional clinical variables (`age`, `sex`, `stage`) and survival follow-up
  (time + event) for KM / log-rank / Cox summaries.

## Important limitations

* Cancer type is **required** and is never inferred.
* **Raw microarray intensity input is not supported** in v1; such data must be harmonised
  externally first.
* The 43 concepts and 132 signatures are **hypothesis-generating representations**, not
  causal or mechanistic conclusions; concept scores are latent coordinates and should not
  be read as pathway activation or suppression.
* **Absolute survival calibration is limited**; the primary output is relative risk.
* Robustness to missing genes is characterised only within the tested range; the model
  should not be described as invariant to missing genes.
