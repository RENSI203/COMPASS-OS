# COMPASS-OS

**Pan-cancer survival prediction and mechanism-oriented transcriptomic representation
framework.**

* **Version:** 1.0.0 · **Python package:** `compass_os` · **License:** MIT
* **Repository:** `https://github.com/RENSI203/COMPASS-OS`

Built on bulk transcriptomic expression + cancer type (+ optional clinical variables), it
returns:

| | |
|---|---|
| **Prognostic risk** | a relative risk score per sample |
| **Cohort risk stratification** | high/low grouping within a submitted cohort |
| **43 concept representations** | COMPASS latent coordinates |
| **132 gene-signature representations** | gene-set-level profiles |
| **Input / robustness QC** | gene, signature and concept-level coverage |

It is a frozen-asset inference package: no training, no re-fitting, no tuning at
prediction time. The representations are hypothesis-generating features for downstream
prioritisation — this is not automatic mechanism discovery.

---

## 1. What is COMPASS-OS?

| | |
|---|---|
| **Input** | bulk expression (`samples x genes`, TPM or log2(TPM+1)), a cancer type per sample, optionally clinical variables and survival follow-up |
| **Output** | relative risk score + risk ranking/stratification · 43 concept profiles · 132 gene-signature profiles · coverage QC · a plain-text report and TSV/PDF files |
| **Model** | a frozen ridge Cox head on the COMPASS representation (default: cancer type + age/sex/stage + 43 concepts) |
| **Trained on** | TCGA pan-cancer overall survival |

Three kinds of output, deliberately kept separate:

1. **Prognostic outputs** — risk score, ranking, stratification.
2. **Mechanism-oriented representations** — 43 concepts, 132 gene signatures.
   These are **hypothesis-generating features** for downstream prioritisation; they are
   *not* evidence of pathway activation, suppression, or causal mechanism.
3. **Input / robustness QC** — gene coverage, signature coverage, concept input coverage.

---

## 2. Quick start

```bash
pip install -e .
export COMPASS_OS_ROOT=$(pwd)     # only if models/ is not found automatically
```

```python
import compass_os

result = compass_os.analyze(
    "expression.tsv",              # samples x genes, TPM
    cancer_type="LUAD",            # TCGA acronym
    clinical="clinical.tsv",       # optional: age / sex / stage
    survival="survival.tsv",       # optional: time + event columns
)

result.summary()
result.save_report("compass_results/")
```

That is the whole workflow. Real output on a small cohort:

```text
COMPASS-OS analysis
────────────────────────────────────────────────────────────────
Samples                         48
Cancer type                     9 cancer types (BRCA, CESC, ESCA, GBM, …)
Model                           Default prognostic model (M2)

INPUT QUALITY
Gene coverage                   100.0% (15672/15672 genes)
Signature-related coverage      100.0% (median across 132 sets)
Concept input coverage          100.0% (median across 43 concepts)
QC                              Recommended range (≥ 90% of required genes observed)

SURVIVAL
Risk estimates                  Available (relative risk, linear predictor)
Risk stratification             Available (frozen cutoff)
High risk                       28
Low risk                        20
Log-rank p (risk groups)        0.0353

BIOLOGICAL REPRESENTATIONS
43 concept profiles             Generated
132 gene-signature profiles     Generated
Interpretation                  hypothesis-generating features

ROBUSTNESS
Missing-gene strategy           Reference imputation (frozen median)
Robustness check                Not required (coverage in recommended range)
```

You do **not** need to know about model variants, scaling, PCA, masking spaces, cancer
tokens or design matrices. Everything that must be fixed is fixed internally.

Runnable scripts: [`examples/quick_start.py`](examples/quick_start.py) (Level 1),
[`examples/cohort_example.py`](examples/cohort_example.py),
[`examples/missing_genes_example.py`](examples/missing_genes_example.py),
[`examples/minimal_example.py`](examples/minimal_example.py) (Level 2).

---

## 3. Understanding the output

| kind | object | what it is | how to use it |
|---|---|---|---|
| **Prognostic** | `result.risk` | relative risk (linear predictor); higher = worse prognosis | ranking samples *within* a cohort |
| | `result.risk_rank` | rank within the analysed cohort | ordering |
| | `result.risk_group` | `high`/`low` by the **frozen** pan-cancer median cutoff | group comparison inside one cohort |
| | `result.risk_group_relative` | cohort-relative median split, **display only** | used when the frozen cutoff fails to separate a cohort |
| **Representation** | `result.concept_scores` | 43 COMPASS latent coordinates | downstream prioritisation, clustering, association testing |
| | `result.signature_scores` | 132 gene-set-level scores | the more directly interpretable layer |
| **QC** | `result.qc` | coverage at gene / signature / concept level | decide how far to trust the readouts |

Two cautions that matter for interpretation:

* **Risk scores are relative, not absolute.** They order samples; they are not calibrated
  probabilities. `survival_probability(t)` exists but is a baseline-hazard-derived output
  (see [ADVANCED_USAGE](docs/ADVANCED_USAGE.md)).
* **Risk comparisons across cancer types are not meaningful.** The model contains
  cancer-type terms, so risk distributions differ by cancer type. Compare within a cohort.

---

## 4. Input requirements

* `expression`: `samples x genes`, **columns are gene symbols**, index is a sample ID.
  Declare the scale with `input_scale` — `"tpm"` (default) or `"log2_tpm1"`.
  Scales are never auto-converted; if the values look inconsistent with the declared
  scale you get a note asking you to check `input_scale`.
* `cancer_type`: required, one TCGA acronym per sample (a single string applies to the
  whole cohort). Accepted values: `src/compass_os/data/cancer_codes.tsv`.
  It is used, from the same value, both as the COMPASS encoder token and as the model's
  cancer-type indicator.
* `clinical` (optional): columns `age`, `sex`, `stage` (stage 1–4). Missing values are
  filled with **frozen training reference values** and listed in the report.
* `survival` (optional): a time column and an event column; common names
  (`time`, `os_time_days`, … / `event`, `os_event`, `status`, …) are auto-detected.

**Direct public inference currently accepts TPM or log2(TPM+1)-scale expression.**
Raw microarray intensity matrices are not accepted directly — they first need
cross-platform harmonisation (see [Validation](#8-validation)).

Missing genes are handled with the frozen TCGA reference median; this never uses
statistics estimated from your own cohort, so a sample's prediction does not depend on
which other samples were submitted alongside it.

---

## 5. Advanced usage

For model selection (`M0`–`M3`), explicit missing-gene strategies, batch/device control,
and direct access to the representation:

```python
rep = compass_os.get_representation(expr, cancer_type, missing_gene_strategy="strict")
pred = compass_os.predict(expr, cancer_type, clinical=clin, model="M0,M1,M2,M3")
rob = compass_os.check_robustness(expr, cancer_type, strategies=("reference", "zero"))
```

Full details, including when to use each model and what each strategy means:
[`docs/ADVANCED_USAGE.md`](docs/ADVANCED_USAGE.md).

Everything below Level 2 — the COMPASS scaler, gene vocabulary, PCA loadings, Cox
coefficients, cancer code table, reference median, signature→concept weights and baseline
hazard — is **read-only frozen assets with no public setter**. To change the model
definition, fork the repository and modify the source.

---

## 6. Missing genes and robustness

Real public cohorts rarely contain all 15,672 COMPASS genes. Measured across 38 external
cohorts (before any imputation) the global COMPASS gene coverage was **0.635–0.980
(median 0.967)** and the 916-gene signature coverage **0.803–0.978 (median 0.962)**.

Three strategies exist (`reference` is the default used by `analyze`):

| strategy | behaviour |
|---|---|
| `reference` | fill missing genes with the frozen TCGA reference median |
| `zero` | training-space masking: the frozen scaler is applied normally, then missing genes' normalised values are set to 0 |
| `strict` | raise `MissingGenesError` if any required gene is missing |

`analyze(..., robustness="auto")` runs a second forward pass **only** when coverage falls
below the frozen recommended threshold, and reports the comparison in plain language.

Coverage QC is always returned, at three levels: gene coverage, per-signature coverage
(`n x 132`) and per-concept input coverage (`n x 43`). Calibrated thresholds live in
`models/qc_config.json`, derived from a stress test — never hard-coded:

| threshold | value | basis |
|---|---|---|
| recommended global coverage | **0.90** | lowest measured coverage with risk Spearman ≥ 0.95 and group agreement ≥ 0.90 |
| warning global coverage | **0.70** | lowest measured coverage with risk Spearman ≥ 0.80 and group agreement ≥ 0.70 |
| recommended signature-gene coverage | **0.90** | same criterion, on the signature-gene axis |
| warning signature-gene coverage | **0.70** | same criterion, on the signature-gene axis |

These are **empirically calibrated descriptive engineering tiers**, not biological or
"safe" thresholds: the underlying curves are smooth with no breakpoint, and the values are
measured coverage levels (no interpolation). The overall QC grade is
**`overall = worse(global tier, signature tier)`**. Per-signature coverage
(`signature_gene_coverage`, *n* × 132) is reported as a **continuous** QC value; no
per-signature numeric threshold is shipped.

Three masking arms were compared (global / signature-only / non-signature-only,
12 cohorts × 20 repeats each):

| coverage on the arm's own axis | 0.95 | 0.90 | 0.80 | 0.70 | 0.60 | 0.50 |
|---|---|---|---|---|---|---|
| **global** — risk Spearman | 0.986 | 0.959 | 0.923 | 0.845 | 0.793 | 0.741 |
| **signature-only** — risk Spearman | 0.983 | 0.964 | 0.924 | 0.849 | — | — |
| **non-signature-only** — risk Spearman | 1.000 | 1.000 | 0.9999 | 0.9999 | 0.9998 | 0.9998 |

Preservation of the signature-related genes accounted for a substantial component of
prognostic robustness: with the 916 signature genes fully retained, removing half of the
remaining 14,756 genes left risk Spearman at 0.9998, group agreement at 1.00 and both
representation layers at a correlation of 1.000, whereas removing the same *fraction* of
the signature genes degraded risk Spearman to 0.849. The non-signature genes are not
without effect — the mean absolute risk shift grew from 0.002 to 0.022 across that range,
and ΔUno C rose to +0.0008 — so the supported statement is that preserving the signature
genes accounts for a substantial component of the robustness, not that the remaining genes
contribute nothing.

`concept_input_coverage` is *input coverage*, **not** a confidence in any biological
interpretation, and must not be renamed to "confidence".

Details and the underlying evidence: [`docs/MISSING_GENES.md`](docs/MISSING_GENES.md).

---

## 7. Biological interpretation

The 43 concepts are **learned latent coordinates** of the COMPASS model, not expression
abundance or pathway-activity scores: every gene-set score is one shared
`nn.Linear(32→1)` projection of the encoder embedding (with `ReLU` disabled, so signs are
unconstrained), and the 43 concepts are a softmax-attention convex combination of the 132
gene-set scores. A dedicated audit found **20/43 concepts reverse-encoded at the model
level**. Concept names are **semantic anchors**.

Practical guidance:

* Use the **132 gene-signature layer** as the primary route to biological interpretation;
  treat the 43 concepts as model-internal coordinates.
* Reported vocabulary: *mechanism-oriented representation*, *candidate biological
  program*, *hypothesis-generating feature*, *risk-associated representation*.
* Do **not** report "activated mechanism", "suppressed pathway" or "causal mechanism".
* Any association between concept scores and the risk score is a **model-linked
  representation association** — the default model uses those very concepts as Cox
  predictors, so this is not independent evidence. Report files label it as such.

See [`docs/CONCEPT_SEMANTICS.md`](docs/CONCEPT_SEMANTICS.md).

---

## 8. Validation

### Reproducibility

The frozen package reproduces production outputs from the original model-development
pipeline: 132 signatures `max|Δ| = 1.1e-16`, 43 concepts `max|Δ| = 2.4e-07`, and M0–M3
risks `max|Δ| ≤ 1.06e-06` (acceptance criterion 1e-5). The official upstream
`PreTrainer.extract()` and the production `predict()` + capture path agree to
`max|Δ| = 0.0`.

### External input availability

Across 38 restorable public cohorts (measured before any imputation):

| | median | range |
|---|---|---|
| global COMPASS gene coverage | **96.7 %** | 63.5 % – 98.0 % |
| signature-gene (916) coverage | **96.2 %** | 80.3 % – 97.8 % |

### Missing-gene robustness

Across the tested masking range, changes in Uno's C remained within approximately ±0.02,
with no consistent deterioration across the evaluated cohorts. Risk *ranking* did degrade
smoothly with coverage (Spearman 0.986 / 0.959 / 0.923 / 0.845 / 0.793 / 0.741 at 95 / 90 /
80 / 70 / 60 / 50 % global coverage), so coverage and ranking stability are reported
together rather than described as "unaffected".

### Signature-related genes

**Preservation of signature-related genes accounted for a substantial component of
prognostic robustness.** With all 916 signature-related genes retained, masking 50 % of the
remaining 14,756 genes preserved risk rankings (Spearman ≈ 0.9998, group agreement 1.00,
43/132 representation correlations 1.000). Masking the same *fraction* of the signature
genes themselves reduced risk Spearman to 0.849. The remaining genes are not irrelevant —
the mean absolute risk shift grew from 0.002 to 0.022 over that range — so the supported
statement is the one above, not that non-signature genes have no contribution.

### Real platform missing patterns and QC thresholds

Replaying the real pre-imputation platform masks of 38 cohorts (training-space masking)
left risk ranking unchanged (Spearman 1.000, group agreement 1.000) while representation
readouts shifted by ≈6 %. The corresponding `reference` arm is a **production-pipeline
self-consistency control** — it returns the input unchanged by construction, because the
harmonised cache is itself reference-imputed — and is not independent robustness evidence.
Coverage QC was predictive: the correlation between per-feature input coverage and
perturbation error was negative for 43/43 concepts and 132/132 signatures (median
Spearman −0.68 and −0.60).

Calibrated QC tiers (see section 6) and the per-level audit are in
`validation/results/qc_threshold_audit.tsv`; design and cost accounting in
[`docs/VALIDATION_DESIGN.md`](docs/VALIDATION_DESIGN.md).

## 9. Limitations

* Trained on TCGA; external performance is dataset-dependent. Discrimination transfers
  better than absolute-risk calibration.
* The frozen pan-cancer median cutoff often fails to separate an external cohort (risk
  levels shift by cancer type). `analyze` then reports a cohort-relative split for
  display only — never for cross-cohort comparison.
* Single samples are not stratified.
* Microarray input requires external harmonisation; an `input_scale="microarray"` mode does
  not exist in v1.
* Robustness to missing genes is characterised only within the tested range; do not
  describe the model as invariant to missing genes.
* Signature-coverage thresholds are not yet calibrated.

---

## 10. Installation

```bash
pip install -e .                    # from the repository root
export COMPASS_OS_ROOT=$(pwd)       # if models/ is not auto-discovered
```

Requires Python ≥ 3.10, PyTorch, pandas, numpy, scikit-learn, scikit-survival, matplotlib.
Model assets (~12 MB) ship in `models/`; the upstream COMPASS package is **vendored** in
`third_party/compass/` (MIT, see `third_party/COMPASS_LICENSE`) and is the default
execution path. An installed `immuno-compass==2.5.3` is only a fallback when `third_party/`
is absent (optional extra `compass-os[upstream]`).

Reproduce the checks:

```bash
python tests/run_tests.py      # no pytest needed; or: pytest tests/
```

`ASSET_MANIFEST.tsv` records every shipped asset (source path, size, SHA-256);
`models/model_manifest.json` records model composition and the default model.

---

## 11. Citation

See [`CITATION.cff`](CITATION.cff). If you use this software, please cite it together with
the upstream COMPASS publication.

---

## 12. License and attribution

Code: [`LICENSE`](LICENSE) (MIT).

**COMPASS-OS incorporates vendored components from COMPASS under its upstream license.**
The upstream COMPASS source is redistributed unmodified in `third_party/compass/` together
with its MIT licence text at [`third_party/COMPASS_LICENSE`](third_party/COMPASS_LICENSE);
the pretrained checkpoint (`models/pretrainer.pt`) is likewise distributed under that
upstream licence. Those components are **not** original work of the COMPASS-OS authors —
please cite the upstream COMPASS publication as well.

> Copyright (c) 2026 RENSI203 (MIT). Software citation authors: Jiahao Ren, Junyi Xin
> (see [`CITATION.cff`](CITATION.cff)).
