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
logical. `FEMALE` is included in the models automatically and does not need to
be listed in `additional_covariates`.

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

- ENET: omics changes plus `FEMALE` and all requested covariates.
- XGB: omics changes plus `FEMALE`.

No baseline or level omics features are included. Processing removes all-missing
training features, then uses training-set medians, variance filtering, centers,
and scales. This removes `FEMALE` when it has zero training variance.

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
    cohort.csv
    change_summary.csv
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

### `cohort.csv`

Contains one row each for the FU-specific eligible, training, and test cohorts:

| Column | Meaning |
|:---|:---|
| `FU` | Follow-up modeled |
| `SET` | `eligible`, `train`, or `test` |
| `N_SUBJECTS` | Number of subjects |
| `N_CONTROL` | Control-arm subjects |
| `N_TREATMENT` | Treatment-arm subjects |
| `N_MALE` | Male subjects |
| `N_FEMALE` | Female subjects |

Eligibility requires both baseline and the modeled follow-up after input-level
covariate filtering.

### `change_summary.csv`

Summarizes raw, pre-imputation omics changes by set and treatment arm:

```text
omics(FU k) - omics(FU 0)
```

| Column | Meaning |
|:---|:---|
| `FU` | Follow-up modeled |
| `SET` | `eligible`, `train`, or `test` |
| `TREATMENT_GROUP` | `0` or `1` |
| `ANALYTE_NAME` | Omics feature |
| `N_SUBJECTS` | Subjects in the set and treatment arm |
| `N_NONMISSING` | Observed change values |
| `MEAN`, `MEDIAN`, `SD`, `MIN`, `MAX` | Raw change-score statistics |

### `preprocessing.csv`

Contains one row per candidate omics or encoded covariate feature:

| Column | Meaning |
|:---|:---|
| `FEATURE_NAME` | Prefixed model-matrix column name |
| `FEATURE_TYPE` | `omics` or `covariate` |
| `STATUS` | `retained`, `all_missing_training`, or `zero_variance_training` |
| `MEDIAN` | Training-set imputation median |
| `CENTER` | Training-set centering value |
| `SCALE` | Training-set scaling value |
| `IN_ENET` | Whether the final ENET matrix contains the feature |
| `IN_XGB` | Whether the final XGB matrix contains the feature |

Unavailable transformation values are blank for removed features.

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

`LAMBDA` is `cv.glmnet()`'s deviance-selected `lambda.min`. `CV_AUC` is the
pooled out-of-fold AUC at that lambda. `LAMBDA_1SE` is recorded for reference
but is not used for the final model.

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
iterations, median best iteration, and parameters for every Optuna trial. XGB
requires at least 10 trials; fewer than 30 produce a limited-search warning.

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
