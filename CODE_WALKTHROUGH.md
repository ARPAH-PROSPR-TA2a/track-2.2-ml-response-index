# FAST Omics ML / Treatment ML: Code Walkthrough

This walkthrough documents the current Track 2.2 implementation. The pipeline
predicts randomized treatment assignment from baseline-to-follow-up omics
changes. Its modeling target is `TREATMENT_GROUP`.

The current public functions are:

- `FAST_treatment_ML()`: prepares model-ready datasets, fits ENET and/or XGB,
  writes all artifacts to disk, and returns a manifest.
- `FAST_treatment_ML_reports()`: runs the shared validation stack and returns
  descriptive and randomization reports without fitting models.

The pipeline is all-subject only. Reporting may produce male and female
summaries, but modeling is not sex-stratified.

## Table of Contents

1. [Main Functions](#main-functions)
2. [Accepted Inputs](#accepted-inputs)
3. [Validation Flow](#validation-flow)
4. [High-Level Modeling Flow](#high-level-modeling-flow)
5. [Subject Splitting and CV Folds](#subject-splitting-and-cv-folds)
6. [Follow-Up Feature Construction](#follow-up-feature-construction)
7. [Preprocessing](#preprocessing)
8. [Model-Specific Feature Sets](#model-specific-feature-sets)
9. [ENET Worker](#enet-worker)
10. [XGB Worker](#xgb-worker)
11. [Disk Artifacts and Manifest](#disk-artifacts-and-manifest)
12. [DNAm Handling](#dnam-handling)
13. [Reporting Pipeline](#reporting-pipeline)
14. [Tests](#tests)

---

## Main Functions

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

Arguments:

| Argument | Default | Meaning |
|:---|:---|:---|
| `pheno` | Required | Phenotype data frame or matrix. Rows represent samples. It must contain `SAMPLE_ID`, `SUBJECT_ID`, `FU`, `TREATMENT_GROUP`, and `FEMALE`. Subjects must have a baseline and at least one follow-up. |
| `omics` | Required | Omics data frame or matrix. Rows represent analytes, `ANALYTE_NAME` identifies each row, and remaining numeric columns are named by `pheno$SAMPLE_ID`. |
| `omics_type` | `"Proteomics"` | One of `"Proteomics"`, `"Metabolomics"`, or `"DNAm"`. It controls preprocessing reminders and DNAm probe filtering; it does not transform input values. |
| `additional_covariates` | `NULL` | Character vector naming phenotype columns to include as covariates. ENET receives all requested covariates as unpenalized features. XGB receives only `FEMALE`, and only when it is explicitly requested here. Rows missing a requested covariate are removed during validation. |
| `models` | `c("enet", "xgb")` | Character vector selecting models to fit. Accepted values are `"enet"` and `"xgb"`. One or both may be requested; duplicates are removed. |
| `output_dir` | `NULL` | Directory for model matrices, metadata, model outputs, and `manifest.json`. When `NULL`, a timestamped directory under `runs/` is created. An existing directory is reused; the function does not clear it first. |
| `test_frac` | `0.2` | Fraction of subjects in each treatment arm assigned to the held-out test set. Must be greater than `0` and less than `0.5`. At least one subject per arm is assigned to test while at least one remains in training. |
| `enet_cv_folds` | `5L` | Requested number of treatment-stratified ENET folds within the training set. Must be at least `2`. The actual count is reduced when the smaller training arm contains fewer subjects. |
| `xgb_cv_folds` | `5L` | Requested number of treatment-stratified XGB folds in each repeat. Must be at least `2`. Each repeat receives an independently seeded fold assignment. |
| `xgb_cv_repeats` | `3L` | Number of independent XGB CV fold assignments used to evaluate each parameter set. Must be at least `1`. Trial selection uses mean AUC across repeats and records repeat variability. |
| `xgb_n_trials` | `50L` | Number of Optuna hyperparameter trials for XGB. Each trial is evaluated across all XGB CV repeats. `0` skips Optuna and evaluates the fixed parameters with the same repeated-CV scheme. |
| `n_cores` | `NULL` | Thread count passed to XGBoost. When `NULL`, defaults to one fewer than the detected core count, with a minimum of one. ENET fitting is not parallelized by this argument. |
| `python_bin` | `NULL` | Python executable used to launch the internal XGB worker. When `NULL`, `"python3"` is used. The selected environment must contain the required Python packages. This argument is irrelevant when XGB is not requested. |
| `seed` | `1L` | Reproducibility seed. It controls the train/test split, CV fold assignment, ENET fitting, XGBoost, and Optuna sampling. Follow-up-specific work uses `seed + fu_level`. |

What it does:

1. Validates model names and modeling arguments.
2. Defaults `n_cores` to one fewer than the detected core count.
3. Creates a timestamped `output_dir` when none is supplied.
4. Defaults `python_bin` to `python3`.
5. Runs shared phenotype and omics validation.
6. Applies DNAm probe restriction when requested.
7. Calls `.run_ml_disk()` to prepare each follow-up, fit requested models, and
   write artifacts.
8. Returns the run manifest.

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

Arguments:

| Argument | Default | Meaning |
|:---|:---|:---|
| `pheno` | Required | The same phenotype structure accepted by `FAST_treatment_ML()`. Validation, covariate filtering, and subject-completeness requirements are identical. |
| `omics` | Required | The same analyte-by-sample omics structure accepted by `FAST_treatment_ML()`. Samples are intersected with the validated phenotype data. |
| `omics_type` | `"Proteomics"` | One of `"Proteomics"`, `"Metabolomics"`, or `"DNAm"`. DNAm reports use the same reliable-probe restriction as modeling. |
| `additional_covariates` | `NULL` | Character vector naming phenotype variables to summarize and assess for baseline treatment-arm balance. Rows missing these variables are removed through the shared validation path. |

This function uses the same validation and DNAm handling as the modeling
pipeline, then calls `.generate_reports()`. It does not create model matrices,
split subjects, or fit models.

---

## Accepted Inputs

The full user-facing contract is in `INPUTS_OUTPUTS.md`.

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

### Modeling Arguments

Important constraints:

- `models` contains one or both of `enet`, `xgb`.
- `test_frac` is greater than `0` and less than `0.5`.
- `enet_cv_folds` and `xgb_cv_folds` are at least `2`.
- `xgb_cv_repeats` is at least `1`.
- `n_cores` is `NULL` or at least `1`.
- `xgb_n_trials` is non-negative.

---

## Validation Flow

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

---

## High-Level Modeling Flow

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

There are no male/female model loops.

---

## Subject Splitting and CV Folds

### Train/Test Split

For each follow-up, `.stratified_subject_split()` works on one row per eligible
`SUBJECT_ID`. Eligibility requires both baseline and the follow-up being
modeled.

For each treatment arm:

```r
n_test <- max(1L, floor(n_subjects_in_arm * test_frac))
n_test <- min(n_test, n_subjects_in_arm - 1L)
```

It samples test subjects separately within treatment arms, then defines all
remaining subjects as training subjects. A subject therefore cannot appear in
both sets.

`.run_ml_disk()` requires at least:

```r
max(2L, enet_cv_folds, xgb_cv_folds)
```

eligible subjects per treatment arm before creating the split. If this
condition fails, that follow-up is skipped while other follow-ups may continue.

### Cross-Validation Folds

`.stratified_subject_folds()` assigns training subjects to folds separately
within each treatment arm.

For each model, the requested fold count is reduced to the size of the smaller
treatment arm:

```r
actual_folds <- min(requested_folds, min(class_counts))
```

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

Omics and requested covariates use the same three-stage preprocessing sequence:

1. Median imputation.
2. Zero-variance filtering.
3. Centering and scaling.

All learned values come from the training set.

### Median Imputation

`.impute_train_median()` computes one median per training column and applies it
to missing values in both train and test. An all-missing training column is
temporarily assigned median `0`, then removed by variance filtering.

### Variance Filtering

`.drop_zero_variance_train()` retains columns with training variance greater
than `1e-12`. Test-set variance does not affect feature retention.

### Scaling

`.scale_train_test()` computes training means and standard deviations. The same
values transform train and test. Missing or zero standard deviations are
replaced by `1`.

### Preprocessing Metadata

For each retained ENET feature, `preprocessing.csv` records:

| Column | Meaning |
|:---|:---|
| `FEATURE_NAME` | Final prefixed model-matrix column |
| `FEATURE_TYPE` | `omics` or `covariate` |
| `MEDIAN` | Training imputation median |
| `CENTER` | Training mean after imputation/filtering |
| `SCALE` | Training standard deviation |

Removed features are omitted.

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
all retained omics changes + all requested additional covariates
```

Factor and logical covariates are expanded with `model.matrix(~ . - 1)`.

### XGB

XGB receives:

```text
all retained omics changes
+ FEMALE only when FEMALE was explicitly requested
```

Other requested covariates are intentionally excluded from XGB. Consequently,
the ENET and XGB matrices can have different column counts.

`preprocessing.csv` describes the complete ENET feature set. Every XGB feature
is a subset of those rows.

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

Although `cv.glmnet()` is run with deviance, the worker calculates AUC for each
lambda from `fit.preval` and selects the lambda with the highest CV AUC. If all
CV AUC values are unavailable, it falls back to `lambda.min`.

### ENET Outputs

- `metrics.csv`: CV, test, and in-sample AUC; selected lambda; `lambda.1se`;
  alpha; feature counts.
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

It requires `numpy`, `pandas`, `scikit-learn`, and `xgboost`. `optuna` is
required only when `xgb_n_trials > 0`.

### Fixed Defaults

Without Optuna tuning, XGB starts with:

```text
objective          binary:logistic
eval_metric        auc
max_depth          2
eta                0.03
min_child_weight   5
subsample          0.8
colsample_bytree   0.6
lambda             10.0
alpha              1.0
```

For each parameter set, `xgb.cv()` runs once per predefined repeat using up to
500 boosting rounds and early stopping after 25 rounds without improvement.
The parameter set is scored by mean best AUC across repeats. CV AUC SD is
retained as a stability diagnostic, and the median repeat-specific best
iteration is used to fit the final model.

### Optuna Tuning

By default, `xgb_n_trials = 50` and `xgb_cv_repeats = 3`, producing 150 XGBoost
CV evaluations per follow-up. Each trial tunes:

- `max_depth`
- `eta`
- `min_child_weight`
- `subsample`
- `colsample_bytree`
- `lambda`
- `alpha`

The objective is maximum cross-validated AUC. The best parameters and best
iteration are used for the final model.

### XGB Outputs

- `metrics.csv`: mean and SD repeated-CV AUC, repeat count, test and in-sample
  AUC, median best iteration, feature count, and selected parameters.
- `predictions.csv`: train and test probabilities.
- `importance.csv`: gain importance and feature type.
- `tuning.csv`: trial-level mean and SD AUC, per-repeat AUC and best iteration,
  median best iteration, and parameter values. With tuning disabled, it has one
  row for the fixed parameter set.
- `model.json`: serialized XGBoost model.

If the Python process exits nonzero, R stops with the worker exit status while
the Python traceback remains visible in the console.

---

## Disk Artifacts and Manifest

The pipeline writes the exact model inputs before fitting either model.

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

The structure repeats for every modelable nonzero follow-up.

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

### Manifest

`manifest.json` records run settings, follow-up directories, canonical input
artifacts, and requested model outputs. The returned R manifest adds
`manifest_path` after the JSON file is written.

Follow-ups that cannot produce a valid prepared dataset are stored as `NULL` in
the in-memory manifest.

---

## DNAm Handling

When `omics_type == "DNAm"`, `.prepare_inputs()` loads:

```text
Data/FAST_epicv1_epicv2_probe_list.rds
Data/FAST_epicv1_epicv2_sugden_TruD_probe_list.rds
```

It validates that the input contains probes from both lists, then restricts
modeling and reporting to the reliable Sugden/TruD probe list.

Unlike the WAS tracks, this ML pipeline does not fit a full-probe model and then
apply a filtered multiple-testing correction. Probe filtering happens before
feature construction.

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

---

## Tests

`run_tests.R` executes:

1. `test_feature_correctness.R`
2. `test_validation.R`
3. `test_end_to_end.R`

Coverage includes:

- subject and outcome alignment;
- train/test isolation;
- follow-up-minus-baseline calculation;
- training-only imputation, filtering, and scaling;
- ENET/XGB feature boundaries;
- validation failures;
- artifact schemas and row counts;
- ENET probability reconstruction from saved weights;
- XGB model serialization and tuning output;
- report return structure.

Run:

```bash
Rscript run_tests.R
```

The test process must use a Python environment containing the dependencies
needed by `scripts/run_xgb.py`. A larger inspectable smoke run is available in
`test_demo_run.R`.
