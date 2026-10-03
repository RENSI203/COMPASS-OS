# Figure 1 — workflow / software architecture (diagram specification)

> **Status: specification only.** The final publication figure is to be drawn separately
> in a vector editor; no auto-generated placeholder is shipped (P4 §17).

## Panel content

```
                 Expression matrix (samples x genes, TPM or log2(TPM+1))
                 + cancer type (TCGA acronym)
                 + optional clinical (age / sex / stage)
                 + optional survival (time, event)
                                  |
                                  v
                        COMPASS-OS  analyze()
                                  |
        +-------------------------+-------------------------+
        |                                                   |
   gene alignment                                    coverage QC
   (frozen 15,672 vocabulary)                        (gene / per-signature /
        |                                             per-concept)
        v                                                   |
   missing-gene handling                                     |
   reference | zero | strict                                 |
        |                                                   |
        v                                                   |
   frozen COMPASS representation  ----------------------------+
        |                                                   |
   +----+-----------------------------+                      |
   |                                  |                      |
   v                                  v                      |
 132 gene-signature scores        43 concept scores          |
   |                                  |                      |
   +---------------+------------------+                      |
                   v                                         |
        frozen Cox model (M0 / M1 / M2* / M3)                 |
                   |                                         |
                   v                                         v
     relative risk + cohort stratification          QC tier / robustness
                   |                                (worse of global & Gsig)
                   v
     report: summary.txt, prediction.tsv, concept_scores.tsv,
             signature_scores.tsv, qc.tsv, KM / risk plots
```

`*` default model.

## Drawing guidance

* Two visual layers must be distinguishable: the **one-click path** (solid arrows,
  left-to-right spine) and the **advanced / QC branches** (lighter arrows).
* The three *output families* must be visually separated, because they carry different
  interpretation status:
  1. **Prognostic outputs** (risk, ranking, stratification) — solid,
  2. **Mechanism-oriented representations** (43 concepts, 132 signatures) — annotated
     "hypothesis-generating",
  3. **Input / robustness QC** — annotated "engineering tiers".
* Frozen components (COMPASS representation, Cox heads, vocabulary, scalers, PCA,
  baseline hazard) should be marked with a lock/freeze glyph and the note
  "read-only frozen assets — no public setter".
* Do not draw the 43 concepts as *derived from* the risk model; they are parallel
  representation readouts, and the default model additionally uses them as predictors
  (mark that edge as "model-linked").
* Colour palette: keep the three output families colour-coded consistently with Figures 2
  and S1–S3.
