# FAST Omics ML / Treatment ML: Code Walkthrough

This walkthrough documents the current Track 2.2 implementation. The pipeline
predicts randomized treatment assignment from baseline-to-follow-up omics
changes. Its modeling target is `TREATMENT_GROUP`.

The current public functions are:

- `FAST_treatment_ML()`: prepares model-ready datasets, fits ENET and/or XGB,
  writes all artifacts to disk, and returns a manifest.
- `FAST_treatment_ML_reports()`: runs the shared validation stack and returns
  descriptive and randomization reports without fitting models.

## Table of Contents

1. [Overview and Public API](#overview-and-public-api)
2. [Shared Input Preparation](#shared-input-preparation)
3. [Modeling Pipeline](#modeling-pipeline)
4. [Subject Splitting and CV Folds](#subject-splitting-and-cv-folds)
5. [Follow-Up Feature Construction](#follow-up-feature-construction)
6. [Preprocessing](#preprocessing)
7. [Model-Specific Feature Sets](#model-specific-feature-sets)
8. [ENET Worker](#enet-worker)
9. [XGB Worker](#xgb-worker)
10. [Disk Artifacts and Manifest](#disk-artifacts-and-manifest)
11. [Reporting Pipeline](#reporting-pipeline)

---

## Overview and Public API

Both public functions begin with the same input preparation:

```text
pheno + omics
      │
      ▼
.prepare_inputs()
      │
      ├── FAST_treatment_ML()
      │     └── split, construct features, preprocess, fit, write artifacts
      │
      └── FAST_treatment_ML_reports()
            └── descriptive summaries and randomization reports
```

The shared arguments and preparation path are documented once in
[Shared Input Preparation](#shared-input-preparation). The complete external
input and output contract is in `INPUTS_OUTPUTS.md`.

### `FAST_treatment_ML()`

File: `main.R`

```r
FAST_treatment_ML <- function(
  pheno,
  omics,
  omics_type = "Proteomics",
  additional_covariates = NULL,
  models = c("enet", "xgb"),
  output_dir = NULL,
  test_frac = 0.2,
  enet_cv_folds = 5L,
  xgb_cv_folds = 5L,
  xgb_cv_repeats = 3L,
  xgb_n_trials = 50L,
  n_cores = NULL,
  python_bin = NULL,
  seed = 1L
)
```

This is the modeling entry point. After shared input preparation, it calls
`.run_ml_disk()` to process each follow-up, fit the requested models, write the
artifacts, and return the run manifest.

### `FAST_treatment_ML_reports()`

File: `main.R`

```r
FAST_treatment_ML_reports <- function(
  pheno,
  omics,
  omics_type = "Proteomics",
  additional_covariates = NULL
)
```

This is the reporting entry point. After the same shared input preparation, it
calls `.generate_reports()`. It does not split subjects, create model matrices,
or fit models.

---

## Shared Input Preparation

Both public functions accept the same four input arguments:

| Argument | Default | Role |
|:---|:---|:---|
| `pheno` | Required | Sample-level phenotype and treatment data. |
| `omics` | Required | Analyte-by-sample numeric measurements. |
| `omics_type` | `"Proteomics"` | Selects `"Proteomics"`, `"Metabolomics"`, or `"DNAm"` handling. It does not transform input values. |
| `additional_covariates` | `NULL` | Optional phenotype columns used as model covariates or reporting variables. `FEMALE` is already required and is added to both models automatically. Rows missing a requested covariate are removed during validation. |

The complete user-facing contract is in `INPUTS_OUTPUTS.md`. The details below
focus on how the code prepares those inputs.

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

Both public functions call `.prepare_inputs()`.

```text
.prepare_inputs()
├── .validate_omics_type()
├── .validate_pheno()
├── .validate_omics()
└── [DNAm only]
    ├── load full and reliable probe lists
    ├── .validate_dnam_probe_coverage()
    └── .subset_omics_list(..., reliable_probes)
```

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
- Creates `all`, `male`, and `female` phenotype subsets for reporting.

The returned `requires_mixed_effects` flag is inherited from shared validation
code but is not used by this ML pipeline.

### Omics Validation

`.validate_omics()`:

- Converts matrix input to a data frame.
- Requires `ANALYTE_NAME`.
- Requires every measurement column to be numeric.
- Intersects omics columns with validated phenotype sample IDs.
- Warns when analytes contain missing values or near-zero variance.
- Creates all-subject and sex-specific omics subsets.

Missing values and near-zero variance are reported here but handled later using
training-only preprocessing for each follow-up.

### DNAm Handling

When `omics_type == "DNAm"`, `.prepare_inputs()` loads:

```text
Data/FAST_epicv1_epicv2_probe_list.rds
Data/FAST_epicv1_epicv2_sugden_TruD_probe_list.rds
```

It validates that the input contains probes from both lists, then restricts
modeling and reporting to the reliable Sugden/TruD probe list.

---

## Modeling Pipeline

The modeling-only controls are:

| Argument | Default | Role |
|:---|:---|:---|
| `models` | `c("enet", "xgb")` | Selects one or both model workers. |
| `output_dir` | `NULL` | Selects the artifact directory; `NULL` creates a timestamped directory under `runs/`. |
| `test_frac` | `0.2` | Sets the treatment-stratified held-out fraction. |
| `enet_cv_folds` | `5L` | Requests the ENET fold count. |
| `xgb_cv_folds` | `5L` | Requests the XGB fold count per repeat. |
| `xgb_cv_repeats` | `3L` | Sets the number of independent XGB fold assignments. |
| `xgb_n_trials` | `50L` | Sets the required Optuna search size. XGB requires at least `10`; values below `30` produce a warning. |
| `n_cores` | `NULL` | Sets XGBoost threads; `NULL` uses one fewer than the detected core count. |
| `python_bin` | `NULL` | Selects the XGB Python executable; `NULL` uses `python3`. |
| `seed` | `1L` | Controls splitting, folds, model fitting, and tuning. |

`FAST_treatment_ML()` validates these controls, resolves their defaults, runs
shared preparation, and enters `.run_ml_disk()`.

```text
FAST_treatment_ML()
│
├── validate arguments and resolve defaults
├── .prepare_inputs()
│
└── .run_ml_disk()
    ├── find every nonzero FU level
    │
    └── for each FU k
        ├── retain subjects with baseline and FU k
        ├── create an FU-specific stratified train/test split
        ├── .prepare_fu_change_dataset()
        ├── write exact model-ready matrices and metadata
        ├── [enet requested] .run_enet_worker()
        ├── [xgb requested] .run_xgb_worker()
        └── add artifact paths to the FU manifest
    │
    └── write manifest.json
```

Each follow-up receives its own train/test split. The split is created only from
subjects with both baseline and that specific follow-up, so treatment balance
and test-set credibility are assessed against the cohort actually available for
that model.

---

## Subject Splitting and CV Folds

### Train/Test Split

For each follow-up, `.stratified_subject_split()` works on one row per eligible
`SUBJECT_ID`. Eligibility requires both baseline and the follow-up being
modeled.

For each treatment arm:

```r
n_test <- max(1L, floor(n_subjects_in_arm * test_frac))
```

It samples test subjects separately within treatment arms, then defines all
remaining subjects as training subjects. A subject therefore cannot appear in
both sets.

Before splitting, `.run_ml_disk()` validates that both treatment arms are
present and that each arm will retain at least:

```r
max(enet_cv_folds, xgb_cv_folds)
```

training subjects after the test allocation. If this condition fails, that
follow-up is skipped while other follow-ups may continue.

### Cross-Validation Folds

`.stratified_subject_folds()` assigns training subjects to folds separately
within each treatment arm.

ENET receives one fold assignment, written to `subjects.csv` as
`ENET_FOLD_ID`. XGB receives `xgb_cv_repeats` independently seeded assignments,
written to `xgb_folds.csv` as `SUBJECT_ID`, `REPEAT`, and `FOLD_ID`.

The train/test split, fold assignment, and model fitting for a follow-up use
`seed + fu_level`.

---

## Follow-Up Feature Construction

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
all retained training subjects, then all retained test subjects
```

Treatment labels come from the follow-up phenotype rows. Preparation is skipped
for a follow-up when:

- no training or test subjects remain after requiring both visits;
- either train or test lacks one treatment arm; or
- stratified CV folds cannot be created.

Baseline omics levels are not model features. Only baseline-to-follow-up change
values enter the models.

---

## Preprocessing

Omics and requested covariates use the same four-stage preprocessing sequence:

1. Remove all-missing training features.
2. Median imputation.
3. Zero-variance filtering.
4. Centering and scaling.

All learned values come from the training set.

### All-Missing Features

`.drop_all_missing_train()` removes columns with no observed training values.
The same columns are removed from the test set, even if test values are
available.

### Median Imputation

`.impute_train_median()` computes one median per training column and applies it
to missing values in both train and test.

### Variance Filtering

`.drop_zero_variance_train()` retains columns with training variance greater
than `1e-12`. Test-set variance does not affect feature retention.

### Scaling

`.scale_train_test()` computes training means and standard deviations. The same
values transform train and test. Missing or zero standard deviations are
replaced by `1`.

### Preprocessing Audit

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

---

## Model-Specific Feature Sets

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
feature set used by each model.

---

## ENET Worker

File: `ml_helpers.R`

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
- `predictions.csv`: train and test probabilities with subject IDs and labels.
- `weights.csv`: intercept, nonzero omics coefficients, and every unpenalized
  covariate coefficient.

The saved weights and model-ready matrix reproduce the saved probabilities
without serializing the R model object.

---

## XGB Worker

Files: `ml_helpers.R`, `scripts/run_xgb.py`

R launches the worker with:

```text
<python_bin> scripts/run_xgb.py
  --data-dir <FU directory>
  --out-dir <FU directory>/xgb
  --seed <seed + FU>
  --nthread <n_cores>
  --n-trials <xgb_n_trials>
```

`python_bin` defaults to `python3`. Callers can provide another executable name
or path through `FAST_treatment_ML(python_bin = ...)`.

The Python worker reads:

- `xgb_train.csv.gz`
- `xgb_test.csv.gz`
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
- `predictions.csv`: train and test probabilities.
- `importance.csv`: gain importance and feature type.
- `tuning.csv`: trial-level mean and SD AUC, per-repeat AUC and best iteration,
  median best iteration, and parameter values.
- `model.json`: serialized XGBoost model.

If the Python process exits nonzero, R stops with the worker exit status while
the Python traceback remains visible in the console.

---

## Disk Artifacts and Manifest

The pipeline writes the exact model inputs before fitting either model.

```text
output_dir/
  manifest.json
  reports/
    cohort.csv
    change_summary.csv
    preprocessing.csv
  FU1/
    enet_train.csv.gz
    enet_test.csv.gz
    xgb_train.csv.gz
    xgb_test.csv.gz
    subjects.csv
    xgb_folds.csv
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

The `FU*` structure repeats for every modelable nonzero follow-up. The
top-level `reports/` files stack rows across all modelable follow-ups.

### `subjects.csv`

This is the shared row map for model matrices:

| Column | Meaning |
|:---|:---|
| `SUBJECT_ID` | Subject represented by the matrix row |
| `FU` | Follow-up modeled |
| `SET` | `train` or `test` |
| `TREATMENT_GROUP` | Binary prediction target |
| `ENET_FOLD_ID` | ENET CV fold for training rows; missing for test rows |

### `xgb_folds.csv`

This file preserves every repeated XGB fold assignment:

| Column | Meaning |
|:---|:---|
| `SUBJECT_ID` | Training subject |
| `REPEAT` | Repeat number |
| `FOLD_ID` | Fold within the repeat |

### `reports/cohort.csv`

For each modeled follow-up, this file records subject, treatment-arm, and sex
counts for:

- `eligible`: subjects with baseline and the follow-up;
- `train`: subjects assigned to model training;
- `test`: subjects assigned to held-out evaluation.

### `reports/change_summary.csv`

This file summarizes the raw, pre-imputation change matrix by `eligible`,
`train`, and `test` set and by treatment arm. Each analyte row records subject
count, nonmissing count, mean, median, SD, minimum, and maximum.

### `reports/preprocessing.csv`

This is the stacked preprocessing audit described above. The `FU` column
identifies the follow-up that produced each feature row.

### Manifest

`manifest.json` records run settings, follow-up directories, canonical input
artifacts, and requested model outputs. The returned R manifest adds
`manifest_path` after the JSON file is written.

Follow-ups that cannot produce a valid prepared dataset are stored as `NULL` in
the in-memory manifest.

---

## Reporting Pipeline

`FAST_treatment_ML_reports()` returns:

```r
list(
  pheno_summary = ...,
  variable_summaries = list(all = ..., male = ..., female = ...),
  randomization_reports = list(
    omics_report = ...,
    covariate_report = ...
  )
)
```

### Phenotype Summary

`.create_pheno_data_report()` reports subject, treatment-arm, and sample counts
for each `FU`/`FEMALE` cell.

### Variable Summaries

For each available reporting stratum and each `FU`/treatment cell:

- omics summaries include nonmissing count, mean, median, SD, minimum, and
  maximum per analyte;
- requested covariate summaries depend on whether the variable is numeric,
  factor, or logical.

### Randomization Reports

Randomization reports use baseline data from the all-subject dataset:

- Omics: Welch t-test per analyte, treatment-minus-control mean difference,
  Cohen's d, SE, raw p-value, and global BH correction.
- Covariates: Welch t-test for numeric variables; Fisher or chi-squared tests
  for factor/logical variables depending on cell counts.

Reporting strata do not imply modeling strata.
