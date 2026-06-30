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
  models/
    reports/
      cohort.csv
      change_summary.csv
      preprocessing.csv
    FU1/
      enet/
        metrics.csv
        weights.csv
      xgb/
        metrics.csv
        importance.csv
        tuning.csv
        model.json
  data/
    FU1/
      enet_train.csv.gz
      xgb_train.csv.gz
      subjects.csv
      xgb_folds.csv
      enet/
        predictions.csv
      xgb/
        predictions.csv
```

One directory is produced per modelable follow-up under both `models/` and
`data/`. The `models/` tree contains aggregate reports and fitted model
artifacts; `data/` contains subject-level training artifacts. There are no
modeling strata.

### Model Matrices

`data/FU*/enet_train.csv.gz` and `data/FU*/xgb_train.csv.gz` are the exact
numeric matrices passed to each worker.

Column prefixes identify provenance:

- `omics::<analyte>`
- `covariate::<encoded covariate>`

### `data/FU*/subjects.csv`

| Column | Meaning |
|:---|:---|
| `SUBJECT_ID` | Subject identifier |
| `FU` | Follow-up modeled |
| `TREATMENT_GROUP` | Outcome |
| `ENET_FOLD_ID` | ENET CV fold |

### `data/FU*/xgb_folds.csv`

Contains the repeated XGB cross-validation assignments:

| Column | Meaning |
|:---|:---|
| `SUBJECT_ID` | Training subject identifier |
| `REPEAT` | XGB CV repeat number |
| `FOLD_ID` | Treatment-stratified fold within that repeat |

### `models/reports/cohort.csv`

Contains one row each for the eligible and training cohorts for every modelable
follow-up. These rows contain the same subjects.

| Column | Meaning |
|:---|:---|
| `FU` | Follow-up modeled |
| `SET` | `eligible` or `train` |
| `N_SUBJECTS` | Number of subjects |
| `N_CONTROL` | Control-arm subjects |
| `N_TREATMENT` | Treatment-arm subjects |
| `N_MALE` | Male subjects |
| `N_FEMALE` | Female subjects |

Eligibility requires both baseline and the modeled follow-up after input-level
covariate filtering.

### `models/reports/change_summary.csv`

Summarizes raw, pre-imputation omics changes by set and treatment arm:

```text
omics(FU k) - omics(FU 0)
```

| Column | Meaning |
|:---|:---|
| `FU` | Follow-up modeled |
| `SET` | `eligible` or `train` |
| `TREATMENT_GROUP` | `0` or `1` |
| `ANALYTE_NAME` | Omics feature |
| `N_SUBJECTS` | Subjects in the set and treatment arm |
| `N_NONMISSING` | Observed change values |
| `MEAN`, `MEDIAN`, `SD`, `MIN`, `MAX` | Raw change-score statistics |

### `models/reports/preprocessing.csv`

Contains one row per candidate omics or encoded covariate feature. This file is
both an audit table and the serialized preprocessing recipe for reconstructing
model matrices on future data.

| Column | Meaning |
|:---|:---|
| `FU` | Follow-up modeled |
| `FEATURE_NAME` | Prefixed model-matrix column name |
| `FEATURE_TYPE` | `omics` or `covariate` |
| `STATUS` | `retained`, `all_missing_training`, or `zero_variance_training` |
| `MEDIAN` | Training-set imputation median |
| `CENTER` | Training-set centering value |
| `SCALE` | Training-set scaling value |
| `IN_ENET` | Whether the final ENET matrix contains the feature |
| `IN_XGB` | Whether the final XGB matrix contains the feature |

Unavailable transformation values are blank for removed features. For a given
follow-up and model, filtering to retained rows with `IN_ENET` or `IN_XGB`
defines the model feature set. The filtered row order is the required matrix
column order.

### `manifest.json`

Contains the run settings, requested models and covariates, follow-up entries,
and paths to every emitted artifact.

## ENET Output

`enet/metrics.csv` contains:

- `CV_AUC`
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

`data/FU*/enet/predictions.csv` contains:

- `SUBJECT_ID`
- `FU`
- `TREATMENT_GROUP`
- `PREDICTED_PROB`

`models/FU*/enet/weights.csv` contains the intercept, selected nonzero omics
coefficients, and all unpenalized covariate coefficients:

- `FEATURE_NAME`
- `WEIGHT`
- `FEATURE_TYPE`

Together with the model-ready matrix, these weights reproduce ENET predictions
without an R model object. The selected lambda remains in `metrics.csv`.

## XGB Output

`xgb/metrics.csv` contains mean repeated-CV AUC, CV AUC SD, repeat count,
in-sample AUC, median best iteration, feature count, and selected XGBoost
parameters.

`data/FU*/xgb/predictions.csv` uses the same schema as ENET predictions.

`xgb/importance.csv` contains:

- `FEATURE_NAME`
- `GAIN`
- `FEATURE_TYPE`

`xgb/tuning.csv` contains mean and SD CV AUC, per-repeat AUCs and best
iterations, median best iteration, and parameters for every Optuna trial. XGB
requires at least 10 trials; fewer than 30 produce a limited-search warning.

`models/FU*/xgb/model.json` is the fitted XGBoost model.

## Exported Model Packages

`FAST_export_models()` writes one JSON package for each fitted model. Each
package is self-contained and can be evaluated without the training manifest.
Models with `CV_AUC >= 0.8` have `successful = true`; lower-AUC models are still
exported with `successful = false`.

Default export location:

```text
output_dir/
  models/
    exported_models/
      exported_models.csv
      <model_id>.json
```

`exported_models.csv` contains:

- `MODEL_ID`
- `FU`
- `MODEL`
- `TRAINING_CV_AUC`
- `SUCCESS_AUC_THRESHOLD`
- `SUCCESSFUL`
- `PATH`

Common fields include:

- `schema_version`
- `model_id`
- `family`
- `fu`
- `target`
- `omics_type`
- `training_cv_auc`
- `success_auc_threshold`
- `successful`
- `covariates`
- `training_metrics`
- `preprocessing`

ENET packages additionally contain `weights`. XGB packages additionally contain
`xgb_model_json`.
