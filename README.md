# FAST Track 2.2

This repository is organized around two related user surfaces:

- [training](training/README.md): fit ENET/XGB treatment-response models and
  emit reproducible artifacts plus exported model JSON packages.
- [inference](inference/README.md): score exported model JSON packages on
  labeled future datasets and summarize cross-trial validation.

The shared artifact contract is:

- `manifest.json`
- `models/reports/preprocessing.csv`
- follow-up-specific model artifacts under `models/FU*/`
- exported model packages under `models/exported_models/`

`models/` contains artifacts for model review and export.
`data/` contains subject-level training artifacts such as model matrices,
predictions, fold assignments, and row maps.

## Current Entry Point

Training currently loads from `training/main.R`:

```r
source(file.path("training", "main.R"))
manifest <- FAST_treatment_ML(...)
exported <- FAST_export_models(manifest)
```

Inference loads from `inference/main.R`:

```r
source(file.path("inference", "main.R"))
result <- FAST_evaluate(
  pheno,
  omics,
  model_path = "runs/run_a/models/exported_models/model_id.json"
)
bulk <- FAST_bulk_evaluate(
  pheno,
  omics,
  models_dir = "runs/run_a/models/exported_models",
  output_dir = "runs/run_a_on_run_b"
)
```

## Documentation

- [Training README](training/README.md)
- [Inference README](inference/README.md)
- [Training input and output schemas](training/INPUTS_OUTPUTS.md)
- [Training code walkthrough](training/CODE_WALKTHROUGH.md)
- [Inference input and output schemas](inference/INPUTS_OUTPUTS.md)
- [Inference code walkthrough](inference/CODE_WALKTHROUGH.md)

## Tests

```bash
Rscript tests/run_tests.R
```
