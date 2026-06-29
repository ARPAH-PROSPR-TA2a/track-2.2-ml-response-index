# FAST Omics ML / Treatment ML: Code Walkthrough

This walkthrough documents the current Track 2.2 implementation. The pipeline
predicts randomized treatment assignment from baseline-to-follow-up omics
changes. Its modeling target is `TREATMENT_GROUP`.

The public function is `FAST_treatment_ML()`: it prepares model-ready datasets,
fits ENET and/or XGB, writes all artifacts to disk, and returns a manifest.

## Table of Contents

1. [Overview](#overview)
2. [Public Function and Arguments](#public-function-and-arguments)
3. [Step 1: Validate and Prepare Inputs](#step-1-validate-and-prepare-inputs)
4. [Step 2: Select Follow-Ups and Split Subjects](#step-2-select-follow-ups-and-split-subjects)
5. [Step 3: Construct Change Features](#step-3-construct-change-features)
6. [Step 4: Preprocess Features](#step-4-preprocess-features)
7. [Step 5: Build Model-Specific Matrices](#step-5-build-model-specific-matrices)
8. [Step 6: Fit ENET](#step-6-fit-enet)
9. [Step 7: Fit XGB](#step-7-fit-xgb)
10. [Step 8: Write Outputs, Reports, and Manifest](#step-8-write-outputs-reports-and-manifest)

---

## Overview

`FAST_treatment_ML()` is the single public entry point for the Track 2.2 ML
pipeline. It takes phenotype and omics data, validates and aligns them, creates
one modeling dataset per nonzero follow-up, converts each dataset into
baseline-to-follow-up change features, preprocesses features using training-only
parameters, fits ENET and/or XGB, writes model artifacts plus reports, and
returns a manifest.

```text
FAST_treatment_ML()
├── validate arguments
├── validate pheno/omics and apply DNAm probe filtering
├── for each nonzero FU:
│   ├── require baseline + FU
│   ├── validate treatment arms and requested CV folds
│   ├── compute FU - baseline omics changes
│   ├── preprocess using train-only parameters
│   ├── build ENET and XGB matrices
│   ├── write model-ready inputs
│   ├── fit requested models
│   └── collect report rows
├── write models/reports/
├── write manifest.json
└── return manifest
```

Each section below follows that execution order. The complete external input
and output contract is in `training/INPUTS_OUTPUTS.md`.

## Public Function and Arguments

File: `training/main.R`

```r
FAST_treatment_ML <- function(
  pheno,
  omics,
  omics_type = "Proteomics",
  additional_covariates = NULL,
  models = c("enet", "xgb"),
  output_dir = NULL,
  enet_cv_folds = 10L,
  xgb_cv_folds = 10L,
  xgb_cv_repeats = 3L,
  xgb_n_trials = 50L,
  n_cores = NULL,
  python_bin = NULL,
  seed = 1L
)
```

`FAST_treatment_ML()` accepts these arguments:

| Argument | Default | Role |
|:---|:---|:---|
| `pheno` | Required | Sample-level phenotype and treatment data. |
| `omics` | Required | Analyte-by-sample numeric measurements. |
| `omics_type` | `"Proteomics"` | Selects `"Proteomics"`, `"Metabolomics"`, or `"DNAm"` handling. It does not transform input values. |
| `additional_covariates` | `NULL` | Optional phenotype columns used as ENET model covariates. `FEMALE` is already required and is added to both models automatically. Rows missing a requested covariate are removed during validation. |
| `models` | `c("enet", "xgb")` | Selects one or both model workers. |
| `output_dir` | `NULL` | Selects the artifact directory; `NULL` creates a timestamped directory under `runs/`. |
| `enet_cv_folds` | `10L` | Requests the ENET fold count. |
| `xgb_cv_folds` | `10L` | Requests the XGB fold count per repeat. |
| `xgb_cv_repeats` | `3L` | Sets the number of independent XGB fold assignments. |
| `xgb_n_trials` | `50L` | Sets the required Optuna search size. XGB requires at least `10`; values below `30` produce a warning. |
| `n_cores` | `NULL` | Sets XGBoost threads; `NULL` uses one fewer than the detected core count. |
| `python_bin` | `NULL` | Selects the XGB Python executable; `NULL` uses `python3`. |
| `seed` | `1L` | Controls folds, model fitting, and tuning. |

Argument validation happens before the pipeline touches data. `.validate_ml_args()`
checks model names, fold counts, XGB repeat count, seed, runtime settings, and
XGB trial count.

## Step 1: Validate and Prepare Inputs

Input preparation is handled by `.prepare_inputs()`.

```text
.prepare_inputs()
├── .validate_omics_type()
├── .validate_pheno()
├── .validate_omics()
└── [DNAm only]
    ├── load full and reliable probe lists
    ├── .validate_dnam_probe_coverage()
    └── .subset_omics(..., reliable_probes)
```

### Phenotype Data

`pheno` must be a data frame or matrix with these columns:

| Column | Requirement |
|:---|:---|
| `SAMPLE_ID` | Globally unique sample identifier |
| `SUBJECT_ID` | Subject identifier shared across visits |
| `FU` | Consecutive integer visit index beginning at `0` |
| `TREATMENT_GROUP` | Binary `0/1`; both arms must be present |
| `FEMALE` | Binary `0/1` |

Variables named in `additional_covariates` must exist and be numeric, factor,
or logical.

### Omics Data

`omics` must be a data frame or matrix with:

| Column | Requirement |
|:---|:---|
| `ANALYTE_NAME` | One feature identifier per row |
| Sample columns | Numeric measurements named by `pheno$SAMPLE_ID` |

Extra samples are allowed on either side. The validator retains the
intersection of phenotype sample IDs and omics column names.

### Phenotype Validation

`.validate_pheno()`:

- Converts matrix input to a data frame.
- Requires all five core columns.
- Requires `FU` values to be consecutive integers `0, 1, ..., max(FU)`.
- Requires baseline and at least one `FU == 1` sample.
- Converts `FU`, `FEMALE`, and `TREATMENT_GROUP` to factors when needed.
- Rejects duplicate `SAMPLE_ID` values.
- For duplicate `SUBJECT_ID`/`FU` pairs, keeps the first row with a warning.
- Drops rows with missing requested covariates.
- Retains only subjects with a baseline and at least one follow-up.
- Returns the validated phenotype table used by this ML pipeline.

### Omics Validation

`.validate_omics()`:

- Converts matrix input to a data frame.
- Requires `ANALYTE_NAME`.
- Requires every measurement column to be numeric.
- Intersects omics columns with validated phenotype sample IDs.
- Warns when analytes contain missing values or near-zero variance.
- Returns the validated omics table used by this ML pipeline.

Missing values and near-zero variance are reported here but handled later using
training-only preprocessing for each follow-up.

### DNAm Handling

When `omics_type == "DNAm"`, `.prepare_inputs()` loads:

```text
Data/FAST_epicv1_epicv2_probe_list.rds
Data/FAST_epicv1_epicv2_sugden_TruD_probe_list.rds
```

It validates that the input contains probes from both lists, then restricts
the ML run to the reliable Sugden/TruD probe list.

---

## Step 2: Select Follow-Ups and Validate Folds

After validation, `FAST_treatment_ML()` enters `.run_ml_disk()`, which finds
every nonzero `FU` level and processes each modelable follow-up independently.
Each follow-up uses all subjects with both baseline and that specific follow-up.

Before model preparation, `.run_ml_disk()` validates that both treatment arms
are present and that the smaller arm has at least as many subjects as each
requested fold count for the selected models. If a request is too large, the run
stops with a message naming the offending argument and suggesting the largest
valid fold count.

### Cross-Validation Folds

`.stratified_subject_folds()` assigns subjects to folds separately within each
treatment arm.

ENET receives one fold assignment, written to `subjects.csv` as
`ENET_FOLD_ID`. XGB receives `xgb_cv_repeats` independently seeded assignments,
written to `xgb_folds.csv` as `SUBJECT_ID`, `REPEAT`, and `FOLD_ID`.

The fold assignment and model fitting for a follow-up use `seed + fu_level`.

---

## Step 3: Construct Change Features

`.prepare_fu_change_dataset()` builds one dataset for each nonzero follow-up.

For each subject with both visits:

```text
omics_change(FU k) = omics(FU k) - omics(FU 0)
```

The matrix operation is:

```r
change_matrix <- t(followup_values - baseline_values)
```

Rows are subjects and columns are analytes. Subject order is explicitly:

```text
all retained training subjects
```

Treatment labels come from the follow-up phenotype rows. Preparation is skipped
for a follow-up when:

- no subjects remain after requiring both visits;
- the training subjects lack one treatment arm; or
- stratified CV folds cannot be created.

Baseline omics levels are not model features. Only baseline-to-follow-up change
values enter the models.

---

## Step 4: Preprocess Features

Omics and requested covariates use the same four-stage preprocessing sequence:

1. Remove all-missing training features.
2. Median imputation.
3. Zero-variance filtering.
4. Centering and scaling.

All learned values come from the training set.

### All-Missing Features

`.drop_all_missing_train()` removes columns with no observed training values.

### Median Imputation

`.impute_train_median()` computes one median per training column and applies it
to missing values.

### Variance Filtering

`.drop_zero_variance_train()` retains columns with training variance greater
than `1e-12`.

### Scaling

`.scale_train()` computes training means and standard deviations. Missing or
zero standard deviations are replaced by `1`.

### Preprocessing Recipe

For each candidate omics or encoded covariate feature, `preprocessing.csv`
records:

| Column | Meaning |
|:---|:---|
| `FEATURE_NAME` | Final prefixed model-matrix column |
| `FEATURE_TYPE` | `omics` or `covariate` |
| `STATUS` | `retained`, `all_missing_training`, or `zero_variance_training` |
| `MEDIAN` | Training imputation median |
| `CENTER` | Training mean after imputation/filtering |
| `SCALE` | Training standard deviation |
| `IN_ENET` | Whether ENET receives the feature |
| `IN_XGB` | Whether XGB receives the feature |

Removed features remain in the audit with their removal status and unavailable
transformation values left blank. XGB uses retained omics features plus retained
`FEMALE`; other retained covariates are ENET-only.

The file is also the preprocessing recipe for future inference. For a given
follow-up, the retained rows where `IN_ENET` or `IN_XGB` is true define that
model's feature set, and their row order is the required matrix column order.

---

## Step 5: Build Model-Specific Matrices

Feature names encode provenance:

```text
omics::<ANALYTE_NAME>
covariate::<model.matrix column>
```

### ENET

ENET receives:

```text
all retained omics changes + FEMALE + all requested additional covariates
```

Factor and logical covariates are expanded with `model.matrix(~ . - 1)`.

### XGB

XGB receives:

```text
all retained omics changes
+ FEMALE
```

Other requested covariates are intentionally excluded from XGB. Consequently,
the ENET and XGB matrices can have different column counts. `FEMALE` is subject
to the same training-only variance filtering as every other covariate.

The `IN_ENET` and `IN_XGB` columns in `preprocessing.csv` identify the retained
feature set used by each model, in model-matrix order.

---

## Step 6: Fit ENET

File: `R/write_training_artifacts.R`

`.run_enet_worker()` fits a binomial elastic net using `glmnet`.

```r
glmnet::cv.glmnet(
  x = enet_x_train,
  y = y_train,
  family = "binomial",
  alpha = 0.5,
  foldid = enet_foldid,
  type.measure = "deviance",
  standardize = FALSE,
  penalty.factor = penalty_factor,
  keep = TRUE
)
```

The prepared matrices are already standardized, so `glmnet` standardization is
disabled.

Penalty factors are:

| Feature type | Penalty factor |
|:---|---:|
| Omics | `1` |
| Requested covariate | `0` |

Covariates are therefore included as unpenalized adjustment features.

### Lambda Selection

`cv.glmnet()` selects `lambda.min`, the lambda with the lowest mean
cross-validated binomial deviance. The worker then calculates pooled
out-of-fold AUC from `fit.preval` at that selected lambda for reporting; AUC
does not influence lambda selection. `lambda.1se` is retained as a reference
value but is not used to fit the final model.

### ENET Outputs

- `metrics.csv`: pooled out-of-fold, test, and in-sample AUC; selected
  `lambda.min`; `lambda.1se`; alpha; feature counts.
- `predictions.csv`: training-subject probabilities with subject IDs and labels.
- `weights.csv`: intercept, nonzero omics coefficients, and every unpenalized
  covariate coefficient.

The saved weights and model-ready matrix reproduce the saved probabilities
without serializing the R model object.

---

## Step 7: Fit XGB

Files: `R/write_training_artifacts.R`, `training/scripts/run_xgb.py`

R launches the worker with:

```text
<python_bin> training/scripts/run_xgb.py
  --data-dir <data/FU directory>
  --out-dir <models/FU directory>/xgb
  --predictions-dir <data/FU directory>/xgb
  --seed <seed + FU>
  --nthread <n_cores>
  --n-trials <xgb_n_trials>
```

`python_bin` defaults to `python3`. Callers can provide another executable name
or path through `FAST_treatment_ML(python_bin = ...)`.

The Python worker reads:

- `xgb_train.csv.gz`
- `subjects.csv`
- `xgb_folds.csv`

It requires `numpy`, `pandas`, `scikit-learn`, `xgboost`, and `optuna`.

### Training and Search Settings

All XGB runs use:

```text
objective          binary:logistic
eval_metric        auc
```

Optuna tunes:

| Parameter | Search range |
|:---|:---|
| `max_depth` | Integer `1` to `4` |
| `eta` | Log-uniform `0.005` to `0.08` |
| `min_child_weight` | Log-uniform `2` to `25` |
| `subsample` | Uniform `0.55` to `0.95` |
| `colsample_bytree` | Uniform `0.25` to `0.85` |
| `lambda` | Log-uniform `1` to `100` |
| `alpha` | Log-uniform `0.01` to `20` |

For each trial, `xgb.cv()` runs once per predefined repeat using up to 500
boosting rounds and early stopping after 25 rounds without improvement. The
trial is scored by mean best AUC across repeats. CV AUC SD is retained as a
stability diagnostic, and the median repeat-specific best iteration is used to
fit the final model.

By default, `xgb_n_trials = 50` and `xgb_cv_repeats = 3`, producing 150 XGBoost
CV evaluations per follow-up. At least 10 trials are required; fewer than 30
produce a warning because the search is limited.

### XGB Outputs

- `metrics.csv`: mean and SD repeated-CV AUC, repeat count, test and in-sample
  AUC, median best iteration, feature count, and selected parameters.
- `predictions.csv`: training-subject probabilities.
- `importance.csv`: gain importance and feature type.
- `tuning.csv`: trial-level mean and SD AUC, per-repeat AUC and best iteration,
  median best iteration, and parameter values.
- `model.json`: serialized XGBoost model.

If the Python process exits nonzero, R stops with the worker exit status while
the Python traceback remains visible in the console.

---

## Step 8: Write Outputs, Reports, and Manifest

The pipeline writes exact model inputs before fitting either model. It also
accumulates per-follow-up report rows during preparation and writes one
`models/reports/` directory after all modelable follow-ups finish.

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

The `FU*` structure repeats under both `models/` and `data/` for every
modelable nonzero follow-up. The `models/reports/` files stack rows across all
modelable follow-ups.

### `data/FU*/subjects.csv`

This is the shared row map for model matrices:

| Column | Meaning |
|:---|:---|
| `SUBJECT_ID` | Subject represented by the matrix row |
| `FU` | Follow-up modeled |
| `TREATMENT_GROUP` | Binary prediction target |
| `ENET_FOLD_ID` | ENET CV fold |

### `data/FU*/xgb_folds.csv`

This file preserves every repeated XGB fold assignment:

| Column | Meaning |
|:---|:---|
| `SUBJECT_ID` | Training subject |
| `REPEAT` | Repeat number |
| `FOLD_ID` | Fold within the repeat |

### `models/reports/cohort.csv`

For each modeled follow-up, this file records subject, treatment-arm, and sex
counts for:

- `eligible`: subjects with baseline and the follow-up;
- `train`: subjects assigned to model training.

### `models/reports/change_summary.csv`

This file summarizes the raw, pre-imputation change matrix by `eligible` and
`train` set and by treatment arm. Each analyte row records subject
count, nonmissing count, mean, median, SD, minimum, and maximum.

### `models/reports/preprocessing.csv`

This is the stacked preprocessing audit described above. The `FU` column
identifies the follow-up that produced each feature row.

### Manifest

`manifest.json` records run settings, follow-up directories, canonical input
artifacts, and requested model outputs. The returned R manifest adds
`manifest_path` after the JSON file is written.

Follow-ups that cannot produce a valid prepared dataset are stored as `NULL` in
the in-memory manifest.

---
