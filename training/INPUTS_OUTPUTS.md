---
output:
  html_document: default
  pdf_document: default
---
# FAST Treatment ML: Inputs and Outputs

This document is a schema reference for `FAST_treatment_ML()` and
`FAST_export_models()`. It describes what the training pipeline expects, what it
writes during fitting, and what is packaged for inference.

## Inputs

### `pheno`

Required columns:

| Column | Requirement |
|:---|:---|
| `SAMPLE_ID` | Unique sample identifier; must match omics sample columns |
| `SUBJECT_ID` | Subject identifier shared across visits |
| `FU` | Consecutive integer visit index beginning at `0` |
| `TREATMENT_GROUP` | Binary assigned treatment, `0/1` |
| `FEMALE` | Binary sex indicator, `0/1` |

`additional_covariates`, when supplied, names extra phenotype columns used only
as ENET training adjustment covariates. Numeric, integer, unordered factor, and
logical columns are accepted. Character and ordered-factor columns are rejected.
`FEMALE` is added automatically and should not be listed in
`additional_covariates`.

Unordered factors use treatment contrasts. For `k` observed levels, the first
declared factor level is the reference and the other levels produce `k - 1`
numeric columns. Logical columns use `FALSE` as the reference and produce one
`<name>TRUE` column when both values are observed. Callers should construct
factors with explicit `levels` or use `relevel()` when the reference category
matters.

Encoding is performed separately for each follow-up after unused factor levels
are dropped. A factor with only one observed level in that follow-up is dropped
with a message. Encoded model-matrix names are normalized with
`make.names(..., unique = TRUE)` before the `covariate::` prefix is added.

Rows missing requested additional covariates are removed before modeling.
Subjects must have a baseline sample and at least one follow-up sample.

### `omics`

| Column | Requirement |
|:---|:---|
| `ANALYTE_NAME` | Feature/probe identifier |
| Sample columns | Numeric measurements named by `pheno$SAMPLE_ID` |

For DNAm runs, the pipeline restricts features to the reliable probe list in
`Data/FAST_epicv1_epicv2_sugden_TruD_probe_list.rds`.

### Model Features

Each nonzero follow-up is modeled independently. Omics features are baseline to
follow-up changes:

```text
omics_change = omics at FU k - omics at FU 0
```

Model-ready feature sets:

| Model | Training features | Deployment features |
|:---|:---|:---|
| ENET | Retained omics changes, retained `FEMALE`, retained encoded `additional_covariates` | Retained omics changes and retained `FEMALE` |
| XGB | Retained omics changes and retained `FEMALE` | Retained omics changes and retained `FEMALE` |

All candidate features use training-only preprocessing: all-missing removal,
median imputation, zero-variance filtering, centering, and scaling. A retained
training-only ENET adjustment covariate is omitted during inference, equivalent
to setting its standardized value to zero.

## Run Output Layout

`FAST_treatment_ML(..., output_dir = "path")` writes:

```text
output_dir/
  manifest.json
  data/
    FU1/
      subjects.csv
      enet_train.csv.gz
      xgb_train.csv.gz
      xgb_folds.csv
      enet/
        predictions.csv
      xgb/
        predictions.csv
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
```

The `FU*` directory repeats for each modelable nonzero follow-up. `data/`
contains model-ready matrices and subject-level outputs. `models/` contains
reports and fitted model artifacts. `manifest.json` records run settings,
requested models/covariates, follow-up entries, and paths to emitted artifacts.

`models/exported_models/` is created only when `FAST_export_models()` is called.

## `data/`: Model-Ready Data

### `data/FU*/subjects.csv`

Shared row map for the model matrices and prediction files.

| Column | Meaning |
|:---|:---|
| `SUBJECT_ID` | Subject represented by the matrix row |
| `FU` | Follow-up modeled |
| `TREATMENT_GROUP` | Binary prediction target |
| `ENET_FOLD_ID` | Treatment-stratified ENET CV fold |

### `data/FU*/enet_train.csv.gz`

Exact numeric matrix passed to ENET. Columns are in the same order as retained
`IN_ENET` rows in `models/reports/preprocessing.csv`.

Feature names are prefixed:

| Prefix | Meaning |
|:---|:---|
| `omics::` | Baseline-to-follow-up omics change |
| `covariate::` | Preprocessed phenotype covariate |

ENET receives retained omics features, retained `FEMALE`, and retained encoded
additional covariates. Covariate columns are unpenalized in ENET. For example,
a factor `site` with levels `A`, `B`, and `C` produces
`covariate::siteB` and `covariate::siteC` when `A` is the reference.

### `data/FU*/xgb_train.csv.gz`

Exact numeric matrix passed to XGB. Columns are in the same order as retained
`IN_XGB` rows in `models/reports/preprocessing.csv`.

XGB receives retained omics features and retained `FEMALE` only. Other requested
covariates are intentionally absent.

### `data/FU*/xgb_folds.csv`

Repeated treatment-stratified XGB CV assignments.

| Column | Meaning |
|:---|:---|
| `SUBJECT_ID` | Training subject identifier |
| `REPEAT` | XGB CV repeat number |
| `FOLD_ID` | Fold within that repeat |

### `data/FU*/{enet,xgb}/predictions.csv`

Training-cohort predictions from the final fitted model.

| Column | Meaning |
|:---|:---|
| `SUBJECT_ID` | Training subject identifier |
| `FU` | Follow-up modeled |
| `TREATMENT_GROUP` | Binary target |
| `PREDICTED_PROB` | Predicted probability for treatment group `1` |

## `models/`: Reports and Fitted Artifacts

### `models/reports/cohort.csv`

One row per modelable follow-up.

| Column | Meaning |
|:---|:---|
| `FU` | Follow-up modeled |
| `N_SUBJECTS` | Number of modeled subjects |
| `N_CONTROL` | Control-arm subjects |
| `N_TREATMENT` | Treatment-arm subjects |
| `N_MALE` | Male subjects |
| `N_FEMALE` | Female subjects |

### `models/reports/change_summary.csv`

Raw, pre-imputation omics change summaries by follow-up, treatment arm, and
analyte.

| Column | Meaning |
|:---|:---|
| `FU` | Follow-up modeled |
| `TREATMENT_GROUP` | `0` or `1` |
| `ANALYTE_NAME` | Omics feature |
| `N_SUBJECTS` | Subjects in the treatment arm |
| `N_NONMISSING` | Observed change values |
| `MEAN`, `MEDIAN`, `SD`, `MIN`, `MAX` | Raw change-score statistics |

### `models/reports/preprocessing.csv`

Stacked audit table and serialized preprocessing recipe. It has one row per
candidate omics feature or encoded covariate for each modeled follow-up.

| Column | Meaning |
|:---|:---|
| `FU` | Follow-up modeled |
| `FEATURE_NAME` | Prefixed model-matrix column name |
| `FEATURE_TYPE` | `omics` or `covariate` |
| `STATUS` | `retained`, `all_missing_training`, or `zero_variance_training` |
| `MEDIAN` | Training-set imputation median |
| `CENTER` | Training-set centering value |
| `SCALE` | Training-set scaling value |
| `IN_ENET` | Whether the ENET training matrix contains the feature |
| `IN_XGB` | Whether the XGB training matrix contains the feature |
| `DEPLOYABLE` | Whether inference uses the feature on future datasets |

Rows removed before scaling have blank transformation values as appropriate. For
deployment, inference uses retained rows with `DEPLOYABLE = TRUE` and the model
flag set. Factor covariates appear as one row per encoded treatment-contrast
column. A single-level factor produces no encoded column and therefore has no
row in this table for that follow-up.

### `models/FU*/enet/metrics.csv`

One row per fitted ENET model.

| Column | Meaning |
|:---|:---|
| `CV_AUC` | Pooled out-of-fold AUC at `lambda.min` |
| `INSAMPLE_AUC` | AUC on the training matrix |
| `LAMBDA` | `cv.glmnet()` `lambda.min`, selected using cross-validated AUC |
| `LAMBDA_1SE` | `cv.glmnet()` one-standard-error lambda, recorded only |
| `ALPHA` | Elastic-net mixing parameter |
| `N_FEATURES` | Number of ENET matrix columns |
| `N_NONZERO` | Number of non-intercept nonzero coefficients |
| `N_UNPENALIZED` | Number of covariate columns fit with `penalty.factor = 0` |

### `models/FU*/enet/weights.csv`

Scoring table for ENET.

| Column | Meaning |
|:---|:---|
| `FEATURE_NAME` | `(Intercept)` or prefixed feature name |
| `WEIGHT` | Fitted coefficient at `LAMBDA` |
| `FEATURE_TYPE` | `intercept`, `omics`, or `covariate` |

This file includes the intercept, selected nonzero omics coefficients, and all
unpenalized covariate coefficients.

### `models/FU*/xgb/metrics.csv`

One row per fitted XGB model.

| Column | Meaning |
|:---|:---|
| `CV_AUC` | Mean AUC across repeated CV runs for the selected trial |
| `CV_AUC_SD` | Standard deviation of repeated CV AUCs |
| `CV_REPEATS` | Number of CV repeats |
| `INSAMPLE_AUC` | AUC on the training matrix |
| `BEST_ITERATION` | Median best boosting round used for final fitting |
| `N_FEATURES` | Number of XGB matrix columns |
| `PARAM_*` | Selected XGBoost parameters |

### `models/FU*/xgb/importance.csv`

Gain-based XGBoost feature importance.

| Column | Meaning |
|:---|:---|
| `FEATURE_NAME` | Prefixed feature name |
| `GAIN` | XGBoost gain importance |
| `FEATURE_TYPE` | `omics` or `covariate` |

### `models/FU*/xgb/tuning.csv`

One row per Optuna trial. It records the trial parameters, mean repeated-CV AUC,
CV AUC SD, per-repeat AUCs, per-repeat best iterations, and median best
iteration. XGB requires at least 10 trials; fewer than 30 produce a
limited-search warning.

### `models/FU*/xgb/model.json`

Serialized fitted XGBoost model written by `xgboost`. This is the worker-level
model artifact used before export.

## Exported Model Packages

`FAST_export_models()` writes one self-contained JSON package per fitted
`(follow-up, model)` and an index file:

```text
output_dir/
  models/
    exported_models/
      exported_models.csv
      <model_id>.json
```

`exported_models.csv` contains:

| Column | Meaning |
|:---|:---|
| `MODEL_ID` | Stable package identifier |
| `FU` | Follow-up modeled |
| `MODEL` | `enet` or `xgb` |
| `TRAINING_CV_AUC` | Training CV AUC copied from model metrics |
| `SUCCESS_AUC_THRESHOLD` | Current success threshold, `0.8` |
| `SUCCESSFUL` | `TRUE` when `TRAINING_CV_AUC >= SUCCESS_AUC_THRESHOLD` |
| `PATH` | Path to the JSON package |

Each `<model_id>.json` has this shape. The example below is ENET; XGB packages
use `xgb_model_json` instead of `weights`.

```json
{
  "schema_version": "1.0",
  "model_id": "...",
  "family": "enet",
  "fu": 1,
  "target": "TREATMENT_GROUP",
  "omics_type": "Proteomics",
  "feature_mode": "change",
  "training_cv_auc": 0.91,
  "success_auc_threshold": 0.8,
  "successful": true,
  "source": {
    "output_dir": "...",
    "manifest_path": "..."
  },
  "covariates": {
    "additional_covariates": ["age"],
    "model_covariates": ["FEMALE", "age"]
  },
  "training_metrics": [{ "...": "..." }],
  "preprocessing": [{ "...": "..." }],
  "weights": [{ "...": "..." }]
}
```

Common fields:

| Field | Meaning |
|:---|:---|
| `schema_version` | Export package schema version |
| `model_id` | Identifier used for the JSON filename and export index |
| `family` | `enet` or `xgb` |
| `fu` | Required follow-up level |
| `target` | Prediction target, currently `TREATMENT_GROUP` |
| `omics_type` | Training omics type |
| `feature_mode` | Feature construction mode, currently `change` |
| `training_cv_auc` | Model CV AUC from training |
| `success_auc_threshold` | Threshold used to set `successful` |
| `successful` | Whether the model passed the training CV AUC threshold |
| `source` | Original training output directory and manifest path |
| `covariates` | Recorded training covariate settings |
| `training_metrics` | The model's `metrics.csv` row |
| `preprocessing` | Retained model-specific preprocessing rows in matrix order |

ENET packages additionally contain:

| Field | Meaning |
|:---|:---|
| `weights` | The ENET `weights.csv` table embedded in JSON |

XGB packages additionally contain:

| Field | Meaning |
|:---|:---|
| `xgb_model_json` | The serialized XGBoost model JSON embedded as a string |

The package `preprocessing` rows include all retained training rows for that
model. During inference, only rows with `DEPLOYABLE = TRUE` are reconstructed.
For ENET, omitted nondeployable covariates are treated as standardized zero.
