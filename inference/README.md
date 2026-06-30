# FAST Treatment ML Inference

The inference pipeline scores saved ENET/XGB artifacts on a labeled cohort.
The preferred cross-trial API evaluates exported model JSON packages.

Run from the repository root:

```r
source(file.path("inference", "main.R"))

result <- FAST_evaluate(
  pheno = pheno,
  omics = omics,
  model_path = "runs/trial_a_treatment_ml/models/exported_models/model_id.json",
  output_dir = "runs/trial_a_on_trial_b"
)
```

`FAST_evaluate()` consumes one self-contained JSON model package emitted by
`FAST_export_models()`. `pheno` must include observed
`TREATMENT_GROUP` labels.

Evaluation writes:

```text
output_dir/
  predictions.csv
  validation.csv
```

`validation.csv` contains:

- `FU`
- `MODEL`
- `TRAINING_CV_AUC`
- `SUCCESS_AUC_THRESHOLD`
- `SUCCESSFUL`
- `N`, `N_CONTROL`, `N_TREATMENT`
- `AUC`
- `LOGIT_BETA`, `LOGIT_OR`, `LOGIT_P`
- `VALIDATED_P05`

`SUCCESSFUL` is hard-coded as `TRAINING_CV_AUC >= 0.8`. Logistic validation
statistics are computed only for successful models, using:

```text
TREATMENT_GROUP ~ PREDICTED_PROB
```

To evaluate all successful models in an exported model directory:

```r
bulk <- FAST_bulk_evaluate(
  pheno = pheno,
  omics = omics,
  models_dir = "runs/trial_a_treatment_ml/models/exported_models",
  output_dir = "runs/trial_a_bulk_on_trial_b"
)
```

Bulk evaluation only evaluates rows where `SUCCESSFUL == TRUE` in
`exported_models.csv` and writes:

```text
output_dir/
  validation_summary.csv
  models/
    <model_id>/
      predictions.csv
      validation.csv
```

`validation_summary.csv` contains one row per evaluated model and uses relative
paths for the per-model prediction and validation files. Columns are:

- `MODEL_ID`
- `MODEL_PATH`
- `FU`
- `MODEL`
- `TRAINING_CV_AUC`
- `SUCCESS_AUC_THRESHOLD`
- `SUCCESSFUL`
- `N`, `N_CONTROL`, `N_TREATMENT`
- `AUC`
- `LOGIT_BETA`, `LOGIT_OR`, `LOGIT_P`
- `VALIDATED_P05`
- `PREDICTIONS_PATH`
- `VALIDATION_PATH`

The older manifest-based replay API remains available for whole-run scoring and
testing:

```r
result <- FAST_treatment_predict(
  pheno = pheno,
  omics = omics,
  manifest_path = "runs/trial_a_treatment_ml/manifest.json",
  models = c("enet", "xgb"),
  output_dir = "runs/trial_a_inference"
)
```

The inference surface consumes artifacts emitted by training:

- `manifest.json`
- `models/reports/preprocessing.csv`
- `models/FU*/xgb/model.json`
- `models/FU*/enet/weights.csv`

Manifest-based predictions are returned by follow-up and model family:

```r
result$predictions$FU1$enet
result$predictions$FU1$xgb
```

Each prediction table contains:

- `SUBJECT_ID`
- `FU`
- `TREATMENT_GROUP`
- `PREDICTED_PROB`

When `output_dir` is supplied, inference writes:

```text
output_dir/
  validation.csv
  predictions/
    FU1/
      enet.csv
      xgb.csv
```

The test suite verifies that inference replay reproduces the training matrices
and saved predictions for ENET and XGB, and that validation outputs are written.
