# Methods snapshot — COMPASS-OS v1.0.0

> Final v1.0.0 methods only. Superseded or abandoned designs are not recorded here.

## Software architecture

`compass-os` 1.0.0 is a frozen-asset inference package with two public layers:

* **Level 1 (`compass_os.analyze`)** — one call takes expression, cancer type, optional
  clinical variables and optional survival follow-up, and returns an `AnalysisResult` with
  `summary()` (plain-text report) and `save_report()` (TSV + PDF/PNG). Model, missing-gene
  strategy and QC thresholds are fixed internally.
* **Level 2 (`predict` / `get_representation` / `check_robustness`)** — exposes model
  choice (M0–M3), missing-gene strategy, input scale, batching/device and the
  cohort-size threshold for stratification.

All model internals (COMPASS scaler, gene vocabulary, PCA loadings, Cox coefficients,
cancer code table, reference median, signature→concept attention weights, baseline hazard)
are read-only frozen assets with no public setter; `fit` on the frozen scaler is blocked at
runtime. The upstream COMPASS package (v2.5.3) is vendored in `third_party/compass/` and is
the default execution path.

## Model definitions

| model | features |
|---|---|
| `M0` | cancer-type one-hot (32 columns) + age + sex + stage |
| `M1` | 43 COMPASS concepts |
| **`M2`** | cancer-type one-hot + age + sex + stage + 43 concepts — **frozen primary / default** |
| `M3` | M2 + PC1–PC10 of log2(TPM+1) (locked PCA) |

All four are ridge-penalised Cox models with Efron ties, fitted previously on TCGA
pan-cancer overall survival; the package only applies the locked coefficients
(`risk = X β`, with `X[scaled_features]` standardised by the locked mean/scale).
Missing clinical values are filled with frozen training reference values
(age 60, sex 0, stage 2). Stratification uses the frozen pan-cancer median cut-off, not the
submitted cohort's median.

## COMPASS representation extraction

Expression (TPM) is aligned to the frozen 15,672-gene symbol vocabulary; the cancer type is
supplied to the model, from the same canonical value, both as the encoder's cancer token
and as the Cox cancer-type indicator. The frozen `Datascaler` (log2(x+1) + MinMax) is
applied inside the checkpoint. The representation is read through the upstream
`PreTrainer.extract()` call, which returns both projection layers:
**132 gene-set-level scores** and **43 cell-pathway-level scores**, the latter being a
softmax-attention convex combination of the former. Using the official API is equivalent to
the production `predict()` + capture path at `max|Δ| = 0.0`.

## Missing-gene strategies

* **`reference` (default)** — missing genes are filled with the frozen TCGA reference
  median, `log2(TPM+1) = 3.130937933922` (TPM ≈ 7.76), taken from
  `models/reference_quantiles.json`. No statistic is estimated from the submitted cohort.
* **`zero`** — training-space masking: the frozen scaler is applied normally and the
  *normalised* value of each missing gene is set to 0, i.e. the same operation used by the
  COMPASS contrastive-pretraining mask augmentor. It is **not** `TPM = 0`: 60.2 % of the
  vocabulary has `data_min_ > 0`, so a raw zero would fall outside the training range.
  Training masking was applied to the positive view only, with probability 0.41, at random;
  inference applies the same operation deterministically to the genes that are missing, so
  the two share the input space but not the distribution.
* **`strict`** — raises `MissingGenesError` listing the missing genes.

## QC definition

Coverage is reported at three levels: global gene coverage, per-signature coverage
(*n* × 132) and per-concept input coverage (*n* × 43, the frozen softmax-weighted average of
its member signatures). The graded axis for signature coverage is the proportion of the 916
signature-associated genes that were observed. Calibrated tiers
(`models/qc_config.json`, v1.0.0): recommended 0.90 and warning 0.70 on both the global and
the signature-gene axis; the overall grade is `worse(global tier, signature tier)`. These
are empirically calibrated descriptive engineering tiers — the underlying curves are smooth
with no breakpoint, and the values are measured coverage levels without interpolation.

## Robustness validation

Three controlled masking arms were run on 12 external cohorts (≤ 60 samples each,
20 repeats per level, all samples of a cohort sharing one mask set):
**global** (random over all 15,672 genes), **signature-only** (random within the 916
signature genes) and **non-signature-only** (random within the remaining 14,756 genes, with
signature coverage pinned at 1.0); plus **empirical replay** of the real pre-imputation
platform masks of 38 cohorts. Every perturbed input was compared against the full-input
prediction on the same samples (`FULL` vs `reference`-filled vs `zero`-masked), with risk
Pearson/Spearman, mean/median/p95 |Δ risk|, risk-rank Spearman, risk-group agreement under
the frozen cut-off, ΔUno's C where event counts allowed, and sample-wise / feature-wise
correlations plus MAD/RMSE for both representation layers. Real platform masks were
recovered from GEO series matrices and platform annotations (probe → symbol, keeping the
highest-mean probe), never back-inferred from imputed matrices; one cohort (METABRIC, not a
GEO platform) was marked unavailable. Chromosome-level compute used four processes with six
torch threads each; the full-sample grid over all 9,359 samples was not run.

## Golden reproducibility

`tests/test_reproducibility.py` compares the package against artefacts of the original
model-development pipeline on 12 TCGA patients: 132 signatures `max|Δ| = 1.1e-16`,
43 concepts `max|Δ| = 2.4e-07` (independent re-check against a second archive with a
different cancer-code convention: 3.0e-07), M0–M3 risks `max|Δ| ≤ 1.06e-06`, against an
acceptance criterion of 1e-5. The alternative `zero` implementation
(input-side `TPM = 2^{data_min} − 1`) reproduces the normalised-space masking at
`max|Δ| = 0.0`.

## Interpretation boundaries

43 concepts and 132 signatures are **mechanism-oriented biological representations /
hypothesis-generating features** for downstream prioritisation. Concept scores are learned
latent coordinates: every gene-set score is one shared `nn.Linear(32→1)` projection of the
encoder embedding with `ReLU` disabled, and a dedicated audit found 20/43 concepts
reverse-encoded at the model level. The package therefore reports association and
prioritisation only, never pathway activation, suppression or causal mechanism. Any
association between concept scores and the risk score is a *model-linked representation
association*, because the default model uses those concepts as Cox predictors.
