# FAST Track 2.2

This repository is organized around two related user surfaces:

- [training](training/README.md): fit ENET/XGB treatment-response models and
  emit reproducible model artifacts.
- [inference](inference/README.md): replay training preprocessing and score
  saved models on labeled future datasets.

The shared artifact contract is:

- `manifest.json`
- `models/reports/preprocessing.csv`
- follow-up-specific model artifacts under `models/FU*/`

`models/` contains the artifacts needed for inference and high-level review.
`data/` contains subject-level training artifacts such as model matrices,
predictions, fold assignments, and row maps.

## Current Entry Point

Training currently loads from `training/main.R`:

```r
source(file.path("training", "main.R"))
manifest <- FAST_treatment_ML(...)
```

Inference loads from `inference/main.R`:

```r
source(file.path("inference", "main.R"))
pred <- FAST_treatment_predict(pheno, omics, manifest_path = "runs/run_a/manifest.json")
```

## Documentation

- [Training README](training/README.md)
- [Inference README](inference/README.md)
- [Input and output schemas](training/INPUTS_OUTPUTS.md)
- [Code walkthrough](training/CODE_WALKTHROUGH.md)

## Tests

```bash
Rscript tests/run_tests.R
```
