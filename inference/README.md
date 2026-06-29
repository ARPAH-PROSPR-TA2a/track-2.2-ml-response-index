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
  models = c("enet", "xgb")
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

The test suite verifies that inference replay reproduces the training matrices
and saved predictions for ENET and XGB.
