# FAST Track 2.2

This repository is organized around two related user surfaces:

- [training](training/README.md): fit ENET/XGB treatment-response models and
  emit reproducible model artifacts.
- [inference](inference/README.md): replay training preprocessing and score
  saved models on future datasets. This surface is planned but not implemented
  yet.

The shared artifact contract is:

- `manifest.json`
- `reports/preprocessing.csv`
- follow-up-specific model artifacts under `FU*/`

`reports/preprocessing.csv` is both an audit table and the preprocessing recipe
needed to reconstruct model matrices for inference.

## Current Entry Point

Training currently loads from `training/main.R`:

```r
source(file.path("training", "main.R"))
manifest <- FAST_treatment_ML(...)
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
