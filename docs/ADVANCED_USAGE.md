# Advanced usage (Level 2 API)

Level 1 (`compass_os.analyze`) fixes the model, the missing-gene strategy and the QC
thresholds. Level 2 exposes the choices a computational researcher legitimately needs —
and nothing below them.

```python
import compass_os

pred = compass_os.predict(expr, cancer_type, clinical=clin, model="M2")
rep  = compass_os.get_representation(expr, cancer_type)
rob  = compass_os.check_robustness(expr, cancer_type, strategies=("reference", "zero"))
```

---

## 1. Models

| model | features | when to use it |
|---|---|---|
| `M0` | cancer type + age/sex/stage | clinical baseline; the comparator for "does the molecular block add anything" |
| `M1` | 43 COMPASS concepts | molecular-only model; use it to separate representation effects from clinical ones |
| **`M2`** | M0 + 43 concepts | **the frozen primary model**; default in both Level 1 and Level 2 |
| `M3` | M2 + PC1–PC10 (locked PCA of log2(TPM+1)) | expression-PC-augmented model, provided additionally |

```python
pred = compass_os.predict(expr, cancer_type, clinical=clin, model="M0,M1,M2,M3")
pred.risk          # DataFrame with one column per requested model
pred.risk["M2"]
```

**M2 is the frozen primary model.** The default is defined in
`models/model_manifest.json`; changing it requires editing that file, not the code.

---

## 2. Missing-gene strategies

| strategy | behaviour | data space | use when |
|---|---|---|---|
| `reference` (default) | fill with the frozen TCGA reference median (`log2(TPM+1) = 3.130937933922` ⟺ TPM ≈ 7.76) | input space | standard inference; the production convention |
| `zero` | the frozen scaler is applied normally, then missing genes' **normalised** values are set to 0 | the same normalised space used by COMPASS contrastive pretraining | you want a masking-based alternative, or a sensitivity comparison |
| `strict` | raise `MissingGenesError` listing the missing genes | — | benchmarks and exact reproduction |

```python
rep = compass_os.get_representation(expr, cancer_type, missing_gene_strategy="strict")
```

`zero` is **not** `TPM = 0`: 60.2 % of the frozen vocabulary has `data_min_ > 0`, so a raw
zero would fall outside the training range. And note the wording boundary — training
masking was applied to the positive view only, with probability 0.41, at random, whereas
inference applies the same operation deterministically to the genes that are actually
missing. The correct description is therefore "the masking operation is applied in the
same normalized input space used during COMPASS contrastive pretraining", not "identical
to the training augmentation distribution". See [`MISSING_GENES.md`](MISSING_GENES.md).

---

## 3. Expression scale

```python
compass_os.predict(expr_tpm,  ct)                       # input_scale="tpm" (default)
compass_os.predict(expr_log2, ct, input_scale="log2_tpm1")
```

Scales are declared, never guessed, and `analyze()` forwards `input_scale` to the same
resolution path as Level 2 (fixed in v1.0.1 — in v1.0.0 the one-click entry point ignored it).
Raw microarray intensities need external harmonisation first; there is no
`input_scale="microarray"` in v1.

---

## 4. Robustness analysis

```python
rob = compass_os.check_robustness(expr, cancer_type,
                                  strategies=("reference", "zero"),
                                  model="M2")
rob.risk_difference                        # Series: reference − zero
rob.concept_correlation["reference_vs_zero"]
rob.signature_correlation["reference_vs_zero"]
rob.gene_coverage
rob.concept_input_coverage["reference"]    # n × 43
rob.robustness_flag                        # "continuous_only": no binary pass/fail in v1.x
```

Each strategy costs one forward pass. `analyze(..., robustness="auto")` is the cheap version:
it performs the comparison only when the combined QC grade
(`overall = worse(global tier, signature tier)`) is not `recommended`.

---

## 5. Accessing the representations

```python
rep = compass_os.get_representation(expr, cancer_type, signature_detail=True)
rep.signature_scores          # n × 132
rep.concept_scores            # n × 43
rep.qc.signature_gene_coverage    # n × 132 per-signature input coverage
rep.qc.concept_input_coverage     # n × 43
```

Interpretation discipline (see [`CONCEPT_SEMANTICS.md`](CONCEPT_SEMANTICS.md)):

* the 132 gene-signature layer is the more directly interpretable one;
* the 43 concepts are learned latent coordinates — 20/43 are reverse-encoded at the model
  level, so do not read them as abundances or activities;
* any association between concept scores and risk is a **model-linked representation
  association**, because M2 uses those concepts as predictors.

---

## 6. Cohort stratification

```python
pred.risk_group        # high/low by the frozen pan-cancer median cutoff
pred.risk_rank         # rank within the submitted cohort
```

Ranking and stratification are separate guarantees, controlled by
`min_cohort_for_stratification` (default 30):

| cohort size | `risk_rank` | `risk_group` |
|---|---|---|
| `n == 1` | not defined | not defined |
| `2 ≤ n < min_cohort_for_stratification` | available | **not assigned** |
| `n ≥ min_cohort_for_stratification` | available | frozen-cutoff `high`/`low` |

Two things to keep in mind:

* the cutoff is a **single pan-cancer value** from the training set, not your cohort median,
  so a sample's group does not depend on the other samples submitted;
* risk levels shift by cancer type, so with external cohorts the frozen cutoff frequently
  places every sample in one group. At Level 1, `analyze` then reports a *cohort-relative*
  split (`result.risk_group_relative`, labelled `cohort_relative_median_display_only`) for
  display only and says so explicitly — it is never merged into `risk_group`. At Level 2 you
  get the frozen grouping unchanged — check that both groups are non-empty before using it.

---

## 7. Absolute survival probabilities (limitations)

```python
pred.survival_probability([365, 1095, 1825])
```

`S_i(t) = exp(−H0(t)·exp(η))` where `η = risk` is the Cox linear predictor, with the Breslow
baseline shipped in the model lock.
This is a **baseline-hazard-derived output**: the baseline comes from the training domain,
and in external validation the strict-transport absolute-risk calibration was substantially
weaker than the discrimination. It is not a validated individual absolute-risk predictor;
prefer relative risk for ranking and comparison.

---

## 8. Performance and devices

```python
compass_os.predict(expr, ct, batch_size=128, device="cpu")
compass_os.predict(expr, ct, min_cohort_for_stratification=30)
```

The COMPASS forward pass costs ≈28–34 ms per sample on CPU and does not benefit from more
torch threads; for large screens, run several processes over disjoint sample sets rather
than raising thread counts (this is what the validation runs did).

---

## 9. What is deliberately not exposed

The following are read-only frozen assets with **no public setter**, and there is no
`refit_pca` / `fit_scaler` / `fit_model` / `custom_cancer_code` entry point:

COMPASS scaler · gene vocabulary · PCA components / mean / scale · Cox coefficients ·
Cox-feature scaler · reference median · cancer numeric codes · signature→concept attention
weights · baseline hazard internals · masking probability · checkpoint internals.

If you need to change the model definition, fork the repository and modify the source; that
scenario is not covered by user documentation. The internal modules
(`preprocessing`, `representation`, `survival`, `qc`, `plotting`) are implementation
details and not part of the public contract.
