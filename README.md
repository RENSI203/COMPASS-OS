<div align="center">

# COMPASS-OS

**Pan-cancer survival risk prediction with COMPASS representations**

Frozen inference · 43 concepts · 132 gene signatures · Coverage QC · Sample computation paths

**English** · [简体中文](README.zh-CN.md)

[![Release](https://img.shields.io/badge/release-v1.1.0-2463eb)](https://github.com/RENSI203/COMPASS-OS/releases/tag/v1.1.0)
[![Python](https://img.shields.io/badge/python-%E2%89%A53.10-3776ab)](pyproject.toml)
[![License](https://img.shields.io/badge/license-MIT-16a085)](LICENSE)

[Quick start](#quick-start) · [Sample path](#sample-path) · [Documentation](#documentation) · [Citation](#citation)

</div>

---

COMPASS-OS is a Python package for **overall-survival risk prediction from bulk transcriptomic expression**. It combines a frozen COMPASS encoder with TCGA-trained Cox models and returns risk scores, molecular representations, and input-quality information in one workflow.

Model assets are included in the package. Inference uses frozen parameters without training, re-fitting, or tuning on the submitted cohort.

| What you get | What it provides |
|---|---|
| **Prognostic risk** | Per-sample Cox linear predictor, η = Xβ |
| **Cohort stratification** | Risk ranking and, when the cohort is large enough, high/low grouping |
| **43 concept representations** | Learned COMPASS coordinates for downstream association and exploration |
| **132 gene-signature representations** | Gene-set-level model readouts |
| **Coverage QC** | Global, signature-related, per-signature, and concept input coverage |
| **Reports and sample paths** | Tables, report figures, and a layered view of a sample's computation path |

The molecular readouts support hypothesis generation. Their numerical direction does not, by itself, establish pathway activation, inhibition, or causality.

## Installation

Download the wheel from [GitHub Releases](https://github.com/RENSI203/COMPASS-OS/releases/tag/v1.1.0), then install it:

```bash
python -m pip install compass_os-1.1.0-py3-none-any.whl
```

Alternatively, install from the repository:

```bash
git clone https://github.com/RENSI203/COMPASS-OS.git
cd COMPASS-OS
python -m pip install .
```

Python ≥3.10 is declared in the package metadata. Runtime dependencies, including PyTorch and torchvision, are installed by pip. The package currently requires `matplotlib>=3.7,<3.9`; the full dependency list is in [pyproject.toml](pyproject.toml).

The COMPASS checkpoint, frozen model assets, and vendored upstream code ship with the package. **No repository path or environment-variable setup is required after installation.**

## Quick start

Prepare `expression.tsv` with **samples as rows**, **gene symbols as columns**, and sample IDs in the first column. For a single-cancer cohort:

```python
import compass_os

result = compass_os.analyze(
    "expression.tsv",
    cancer_type="LUAD",
    input_scale="tpm",
)

print(result.summary())
result.save_report("compass_results/")
```

Clinical variables and survival follow-up are optional:

```python
result = compass_os.analyze(
    "expression.tsv",
    cancer_type="LUAD",
    clinical="clinical.tsv",    # age, sex, stage; indexed by sample ID
    survival="survival.tsv",    # time and event; indexed by sample ID
    input_scale="tpm",
)
```

`analyze()` selects the model specified by the frozen manifest: **M2 in v1.1.0**. M3 is also available through the parameterized API below. For mixed-cancer cohorts, provide one cancer acronym per sample in expression-row order.

To run the included example from a cloned repository:

```bash
python examples/quick_start.py
```

## Sample path

<p align="center">
  <img src="docs/figures/sample_path_redesign/high/M2.png" alt="COMPASS-OS circular sample computation path: gene expression, gene score, signature score, high-level predictors, and cohort-relative risk" width="100%">
</p>

*Example: a high-risk sample from GSE39582 (COAD), using M2. The figure shows a selected computation path, not a causal mechanism or pathway-activation map.*

The circular diagram separates gene expression, COMPASS gene scores, signature scores, and high-level predictors. Node colors encode the displayed values; thin, uniform links show selected model connections. The risk bar shows the sample's cohort-relative position, a triangle marker, and a cutoff line.

```python
import pandas as pd
import compass_os

expr = pd.read_csv("expression.tsv", sep="\t", index_col=0)
expr.index = expr.index.map(str)
sample_id = expr.index[0]

path = compass_os.sample_path(
    expr,
    cancer_type="LUAD",
    sample_id=sample_id,
    model="M3",
    input_scale="tpm",
)

print(path.summary())
path.save("compass_results/sample_path")
```

Use a cohort containing the sample of interest to give its rank and percentile context. `save()` exports HTML, PNG, PDF, SVG, visible-node values, selection metadata, and the full Cox decomposition. Individual exports are also available:

```python
path.save_html("sample_path.html")
path.save_static("sample_path.pdf")    # also .png or .svg
```

- HTML embeds the same PNG image and includes an expandable exact-node table; the main figure is static.
- Clinical fields filled from frozen references are marked **imputed**. PC1–PC10 is shown as one aggregated node.
- The default diagram cutoff is the **submitted cohort's median**, a visualization boundary distinct from the frozen prediction cutoff.
- M1, M2, and M3 are supported. M0 has no COMPASS concept branch and is not supported by this visualization.

Display thresholds, node limits, and colors can be configured. See the [sample-path guide](docs/SAMPLE_PATH.md) and [runnable example](examples/sample_path_example.py).

## Inputs and outputs

### Input contract

| Input | Requirement |
|---|---|
| Expression | Samples × genes; gene-symbol columns; sample-ID index |
| Expression scale | `input_scale="tpm"` or `input_scale="log2_tpm1"` |
| Cancer type | Required; supported TCGA tumor acronym per sample. A single string is supported by `analyze()` and `sample_path()` for a single-cancer cohort |
| Clinical data | Optional numeric `age`, `sex`, `stage`; use the documented frozen-model encoding and matching sample IDs |
| Survival follow-up | Optional time/event columns with consistent units and matching sample IDs |

Missing clinical values are filled from frozen training references and reported. Missing genes follow the selected strategy. Declare the expression scale correctly: **nonnegative counts or raw intensities can resemble TPM numerically**, and the package cannot reliably infer their units. Raw microarray intensities are not a directly supported input scale.

### Reading an analysis result

| Field | Meaning |
|---|---|
| `result.risk` | Cox linear predictor η; higher values indicate higher model-estimated hazard |
| `result.risk_rank` | Within-cohort ranking, when available |
| `result.risk_group` | High/low grouping using the frozen cutoff, when available |
| `result.risk_group_relative` | Separately labeled cohort-median grouping for display, when returned |
| `result.concept_scores` | Samples × 43 concept representations |
| `result.signature_scores` | Samples × 132 gene-signature representations |
| `result.qc` | Input-coverage measures and engineering QC grades |

`risk` is **not a probability**: η is on the log-relative-hazard scale; `exp(η)` is the relative-hazard multiplier. Absolute survival outputs derived from the training baseline hazard have separate calibration limitations.

Single samples do not receive a prediction rank or high/low group. By default, prediction ranking is available for at least two samples, while grouping requires at least 30. A frozen cutoff does not guarantee a balanced split in an external cohort. Interpret rankings in a defined cohort and account for cancer composition when comparing samples.

## Parameterized usage

Use the one-call `analyze()` workflow for routine analysis. For explicit model selection and missing-gene handling, use `predict()`, `get_representation()`, or `check_robustness()`:

```python
import pandas as pd
import compass_os

expr = pd.read_csv("expression.tsv", sep="\t", index_col=0)
expr.index = expr.index.map(str)
ct = ["LUAD"] * len(expr)       # low-level API: one label per sample

pred = compass_os.predict(
    expr,
    ct,
    model="M3",
    input_scale="tpm",
    missing_gene_strategy="reference",
)

print(pred.risk)
print(pred.concept_scores.shape)
print(pred.signature_scores.shape)
```

| Model | Frozen predictors |
|---|---|
| M0 | Cancer type + age/sex/stage |
| M1 | 43 COMPASS concepts |
| **M2 — default** | Cancer type + age/sex/stage + 43 concepts |
| M3 | M2 + PC1–PC10 from the locked expression PCA |

| Missing-gene strategy | Behavior |
|---|---|
| **`reference` — default** | Fill from the frozen reference median |
| `zero` | Apply the frozen scaler, then set missing genes' normalized input values to zero |
| `strict` | Raise an error for missing required genes |

`zero` is a normalized-input masking operation, not a claim that unmeasured genes have zero TPM. Global and signature-related QC grades are engineering tiers; per-signature coverage remains a continuous measure. Model internals are frozen and have no public re-fitting interface.

See [advanced usage](docs/ADVANCED_USAGE.md) for clinical coding, batching, robustness comparison, and derived survival outputs.

## Evidence and interpretation

The release documents numerical reproduction within a declared tolerance, external input coverage, and missing-gene stress tests. In 38 restorable external cohorts, pre-imputation median coverage was **96.7% globally** and **96.2% for the release's 916 signature-related genes**. See [missing-gene evidence](docs/MISSING_GENES.md) and the [validation design](docs/VALIDATION_DESIGN.md) for definitions and tested conditions.

Important boundaries:

- Relative to the clinical baseline, the audited 24-cohort M2 comparison **did not detect a significant incremental gain**. This is not an equivalence claim or proof of superiority.
- Concept and signature scores are model representations, not independently validated pathway-activity measurements. Their signs alone do not establish biological activation or suppression.
- Associations between risk and features used by the Cox model are model-linked associations, not independent mechanistic evidence.
- These readouts can complement transcriptome-wide downstream analysis; signature membership should not restrict the genes eligible for mechanism discovery.
- Coverage QC and missing-gene robustness apply to the tested settings. They do not establish universal cross-platform performance or robustness to every missingness mechanism.
- Individual absolute-survival calibration and transfer of stratification to new cohorts require separate validation.

## Documentation

| Guide | Contents |
|---|---|
| [Advanced usage](docs/ADVANCED_USAGE.md) | Model selection and the parameterized API |
| [API reference](docs/API.md) | Public interfaces and result objects |
| [Sample paths](docs/SAMPLE_PATH.md) | Figure semantics, display parameters, and exports |
| [Missing genes](docs/MISSING_GENES.md) | Imputation, masking, coverage QC, and robustness evidence |
| [Concept semantics](docs/CONCEPT_SEMANTICS.md) | Interpretation of the 43/132 representation layers |
| [Validation design](docs/VALIDATION_DESIGN.md) | Validation setup and tested conditions |
| [Release notes](RELEASE_NOTES_v1.1.0.md) | Changes in v1.1.0 |

Examples: [quick start](examples/quick_start.py) · [cohort analysis](examples/cohort_example.py) · [missing genes](examples/missing_genes_example.py) · [sample paths](examples/sample_path_example.py).

## Citation

Please cite **COMPASS-OS and the upstream COMPASS publication** when using this software. Software citation metadata are maintained in [CITATION.cff](CITATION.cff); GitHub's **Cite this repository** option provides citation formats.

COMPASS-OS software authors: **Jiahao Ren and Junyi Xin**. The pretrained COMPASS model and vendored upstream components retain their original attribution.

## License and contact

COMPASS-OS is released under the [MIT License](LICENSE). Upstream attribution is recorded in [NOTICE](NOTICE) and the [vendored COMPASS license](src/compass_os/assets/third_party/COMPASS_LICENSE).

Questions and reproducible bug reports: [GitHub Issues](https://github.com/RENSI203/COMPASS-OS/issues).

Contact: **Jiahao Ren** — [rjh2623826975@stu.njmu.edu.cn](mailto:rjh2623826975@stu.njmu.edu.cn).

---

[简体中文](README.zh-CN.md) · [Release downloads](https://github.com/RENSI203/COMPASS-OS/releases/tag/v1.1.0) · [Back to top](#compass-os)
