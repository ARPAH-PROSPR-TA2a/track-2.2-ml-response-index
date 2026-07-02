# FAST Treatment ML Inference

Inference scores exported FAST treatment models on a labeled cohort and reports validation against observed treatment status. It uses the same input validation, omics change definition, feature naming, and training-derived preprocessing recipe as the training pipeline, but only deployable features are required: omics changes and retained `FEMALE`. Training-only ENET adjustment covariates are not required on future datasets.

Run from the repository root:

``` r
source(file.path("inference", "main.R"))

result <- FAST_evaluate(
  pheno = pheno,
  omics = omics,
  model_path = "runs/trial_a_treatment_ml/models/exported_models/model_id.json",
  output_dir = "runs/trial_a_on_trial_b"
)
```

`FAST_evaluate()` evaluates one self-contained JSON model package emitted by `FAST_export_models()`. `pheno` must include observed `TREATMENT_GROUP` labels.

To evaluate every successful model package in a directory:

``` r
bulk <- FAST_bulk_evaluate(
  pheno = pheno,
  omics = omics,
  models_dir = "runs/trial_a_treatment_ml/models/exported_models",
  output_dir = "runs/trial_a_bulk_on_trial_b"
)
```

Bulk evaluation scans the directory for exported model JSON packages, evaluates only packages with `successful == true`, writes per-model predictions/validation files, and writes a shareable `validation_summary.csv`.

Validation uses the training success threshold `TRAINING_CV_AUC >= 0.8`. For successful models, inference reports validation AUC and the logistic association:

``` text
TREATMENT_GROUP ~ PREDICTED_PROB
```

## Documentation

-   [Input and output schemas](INPUTS_OUTPUTS.md)
-   [Code walkthrough](CODE_WALKTHROUGH.md)
-   [Training inputs and outputs](../training/INPUTS_OUTPUTS.md)
