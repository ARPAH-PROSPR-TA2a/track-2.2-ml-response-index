# FAST Treatment ML: Inputs and Outputs

## Inputs

### `pheno`

Required columns:

| Column | Requirement |
|:---|:---|
| `SAMPLE_ID` | Unique sample identifier |
| `SUBJECT_ID` | Subject identifier shared across visits |
| `FU` | Consecutive integers starting at `0` |
| `TREATMENT_GROUP` | Binary assigned treatment, `0/1` |
| `FEMALE` | Binary sex indicator, `0/1` |

Variables named through `additional_covariates` must be numeric, factor, or
logical. Request `FEMALE` explicitly to include sex in the models.

### `omics`

| Column | Requirement |
|:---|:---|
| `ANALYTE_NAME` | Feature/probe identifier |
| Sample columns | Numeric measurements named by `pheno$SAMPLE_ID` |

### Model Features

For every nonzero follow-up:

```text
omics_change = omics at FU k - omics at FU 0
```

- ENET: omics changes plus all requested covariates.
- XGB: omics changes plus `FEMALE` when requested.

No baseline or level omics features are included. Processing uses training-set
medians, variance filtering, centers, and scales.

DNAm is restricted to the reliable probe list in
`Data/FAST_epicv1_epicv2_sugden_TruD_probe_list.rds`.

## Output Directory

```text
output_dir/
  manifest.json
  FU1/
    enet_train.csv.gz
    enet_test.csv.gz
    xgb_train.csv.gz
    xgb_test.csv.gz
    subjects.csv
    xgb_folds.csv
    preprocessing.csv
    enet/
      metrics.csv
      predictions.csv
      weights.csv
    xgb/
      metrics.csv
      predictions.csv
      importance.csv
      tuning.csv
      model.json
```

One directory is produced per modelable follow-up. There are no modeling
strata.

### Model Matrices

`enet_train.csv.gz`, `enet_test.csv.gz`, `xgb_train.csv.gz`, and
`xgb_test.csv.gz` are the exact numeric matrices passed to each worker.

Column prefixes identify provenance:

- `omics::<analyte>`
- `covariate::<encoded covariate>`

### `subjects.csv`

| Column | Meaning |
|:---|:---|
| `SUBJECT_ID` | Subject identifier |
| `FU` | Follow-up modeled |
| `SET` | `train` or `test` |
| `TREATMENT_GROUP` | Outcome |
| `ENET_FOLD_ID` | ENET CV fold for training rows; blank for test rows |

### `xgb_folds.csv`

Contains the repeated XGB cross-validation assignments:

| Column | Meaning |
|:---|:---|
| `SUBJECT_ID` | Training subject identifier |
| `REPEAT` | XGB CV repeat number |
| `FOLD_ID` | Treatment-stratified fold within that repeat |

### `preprocessing.csv`

Contains one row per retained ENET feature:

| Column | Meaning |
|:---|:---|
| `FEATURE_NAME` | Prefixed model-matrix column name |
| `FEATURE_TYPE` | `omics` or `covariate` |
| `MEDIAN` | Training-set imputation median |
| `CENTER` | Training-set centering value |
| `SCALE` | Training-set scaling value |

Removed features are omitted.

### `manifest.json`

Contains the run settings, requested models and covariates, follow-up entries,
and paths to every emitted artifact.

## ENET Output

`enet/metrics.csv` contains:

- `CV_AUC`
- `TEST_AUC`
- `INSAMPLE_AUC`
- `LAMBDA`
- `LAMBDA_1SE`
- `ALPHA`
- `N_FEATURES`
- `N_NONZERO`
- `N_UNPENALIZED`

`enet/predictions.csv` contains:

- `SET`
- `SUBJECT_ID`
- `FU`
- `TREATMENT_GROUP`
- `PREDICTED_PROB`

`enet/weights.csv` contains the intercept, selected nonzero omics coefficients,
and all unpenalized covariate coefficients:

- `FEATURE_NAME`
- `WEIGHT`
- `FEATURE_TYPE`

Together with the model-ready matrix, these weights reproduce ENET predictions
without an R model object. The selected lambda remains in `metrics.csv`.

## XGB Output

`xgb/metrics.csv` contains mean repeated-CV AUC, CV AUC SD, repeat count, test
and in-sample AUCs, median best iteration, feature count, and selected XGBoost
parameters.

`xgb/predictions.csv` uses the same schema as ENET predictions.

`xgb/importance.csv` contains:

- `FEATURE_NAME`
- `GAIN`
- `FEATURE_TYPE`

`xgb/tuning.csv` contains mean and SD CV AUC, per-repeat AUCs and best
iterations, median best iteration, and parameters for every Optuna trial. When
tuning is disabled, it contains one row for the fixed parameter set.

`xgb/model.json` is the fitted XGBoost model.

## Reports

`FAST_treatment_ML_reports()` returns:

```r
list(
  pheno_summary = ...,
  variable_summaries = ...,
  randomization_reports = ...
)
```
