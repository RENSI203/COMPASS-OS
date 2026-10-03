# COMPASS-OS v1.0.1

Bug-fix and packaging release. **No scientific model changes.**

## Bug fixes

* **`analyze()` now propagates `input_scale` correctly.** In v1.0.0 the one-click entry point
  accepted `input_scale="log2_tpm1"` but did not forward it to the prediction layer, so
  log2(TPM+1) input was processed as if it were linear TPM. Both scales now give identical
  results. `batch_size` and `device` remain fixed internally and are still not exposed at
  Level 1.
* **`robustness="auto"` now follows the combined two-axis QC grade.** It previously keyed off
  global gene coverage alone, so a cohort with adequate global coverage but inadequate
  signature-gene coverage (or vice versa) skipped the robustness comparison despite a
  non-`recommended` overall grade. Auto now runs whenever
  `overall = worse(global tier, signature tier)` is not `recommended`, using the same frozen
  thresholds — no new threshold was introduced.
* **QC verdict and robustness text are now consistent.** The summary can no longer report
  `QC: usable with caution` together with `Robustness check: Not required`.
* **`min_cohort_for_stratification` now matches its documented behaviour.** Ranking and
  stratification are separate guarantees:
  `n == 1` → neither; `2 ≤ n < min_cohort` → `risk_rank` only; `n ≥ min_cohort` → `risk_rank`
  plus frozen-cutoff `high`/`low`. Previously `risk_group` was emitted for any cohort with
  `n ≥ 2`, with only a warning. The Level 1 display-only cohort-relative split
  (`risk_group_relative`, labelled `cohort_relative_median_display_only`) remains a separate
  field and is unchanged.
* **Checkpoint no longer deleted when the environment lives under the temp directory.**
  The package used to load the frozen model through upstream `compass.loadcompass()`, whose
  final lines are

  ```python
  if file.startswith(tempfile.gettempdir()):
      os.remove(file)
  ```

  That was intended to clean up a model downloaded from a URL, but the test is a *path
  prefix*: any checkpoint located under the system temp directory was **silently deleted right
  after loading**. Any environment installed under `/tmp` — containers, CI runners,
  `pip install --target /tmp/...`, some conda/Docker layouts — therefore worked exactly once
  and then failed on every later call. `compass_os` now loads the local frozen asset with
  `torch.load` directly and reproduces the two non-destructive post-processing steps of the
  upstream loader, so the vendored upstream source stays byte-identical while the side effect
  is gone. Numerical behaviour is unchanged.
* **Wheel / standard-install packaging fixed.** Model assets and the vendored upstream COMPASS
  package now ship as Python package data under `compass_os/assets/`, so `pip install .` and
  `pip install compass_os-1.0.1-py3-none-any.whl` resolve everything from `site-packages`.
  No Git checkout and no `COMPASS_OS_ROOT` are required any more; the environment variable is
  retained only as an optional developer override. Bytecode is excluded from both wheel and
  sdist.
* **Quick Start example and documented output now agree.** The README previously showed an
  output block ("48 samples / 9 cancer types") that the shown code could not produce, since it
  passed a single `cancer_type="LUAD"`. The README now shows the real output of
  `examples/quick_start.py`, and a regression test compares the two.
* **Cox risk terminology clarified.** `risk` (a.k.a. `linear_predictor`) is the Cox linear
  predictor η = Xβ on the log-relative-hazard scale, described as the *prognostic risk score*;
  `exp(η)` — not `η` — is the relative hazard. The numerical API is unchanged.
* **Obsolete robustness wording updated.** `robustness_flag` is now `"continuous_only"` instead
  of the stale `"not_calibrated"`, and the accompanying note states that coverage QC tiers are
  calibrated and reported separately, while the reference-versus-zero comparison has no binary
  pass/fail threshold in v1.x.

## Packaging changes

* Runtime assets moved from the repository root into the package:
  `models/` → `src/compass_os/assets/models/`,
  `third_party/` → `src/compass_os/assets/third_party/`.
* Asset resolution rewritten around `importlib.resources` with a `COMPASS_OS_ROOT` developer
  override. Repository-relative resolution is no longer used by the installed package.
* `LICENSE` stays verbatim MIT; third-party attribution lives in `NOTICE` and
  `assets/third_party/COMPASS_LICENSE` (both shipped in wheel and sdist).

## No scientific model changes

* M0–M3 definitions, coefficients and the default model (M2) are unchanged.
* The COMPASS checkpoint, the locked PCA and the reference quantiles are byte-identical
  (SHA-256 unchanged; see `ASSET_MANIFEST.tsv`).
* The 132 gene-signature and 43 concept definitions are unchanged.
* The `reference` / `zero` / `strict` missing-gene rules are unchanged.
* QC thresholds (recommended 0.90 / warning 0.70, both axes) are unchanged.
* Validation results, manuscript numbers and paper figures are unchanged.

## Upgrading

```bash
git pull            # or download the v1.0.1 wheel
pip install .
```

`v1.0.0` remains an immutable earlier release; use `v1.0.1` for the fixes above.
