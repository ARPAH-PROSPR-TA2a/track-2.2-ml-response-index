# FAST Treatment ML Inference

The inference pipeline replays a training run's preprocessing recipe and scores
saved ENET/XGB artifacts on a labeled cohort.

Run from the repository root:

```r
source(file.path("inference", "main.R"))

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

`pheno` must include observed `TREATMENT_GROUP` labels. Predictions are returned
by follow-up and model family:

```r
result$predictions$FU1$enet
result$predictions$FU1$xgb
```

Each prediction table contains:

- `SUBJECT_ID`
- `FU`
- `TREATMENT_GROUP`
- `PREDICTED_PROB`

Inference also returns `result$validation`, one row per scored follow-up/model:

- `FU`
- `MODEL`
- `TRAINING_CV_AUC`
- `CATALOG_AUC_THRESHOLD`
- `CATALOGED`
- `N`, `N_CONTROL`, `N_TREATMENT`
- `AUC`
- `LOGIT_BETA`, `LOGIT_OR`, `LOGIT_P`
- `VALIDATED_P05`

`CATALOGED` is hard-coded as `TRAINING_CV_AUC >= 0.8`. Logistic validation
statistics are computed only for cataloged models, using:

```text
TREATMENT_GROUP ~ PREDICTED_PROB
```

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
