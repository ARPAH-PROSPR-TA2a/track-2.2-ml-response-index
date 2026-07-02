# FAST Treatment ML: Code Walkthrough

This document explains how the training code moves data through the pipeline.
It is intentionally lighter than the schema reference. For exact input columns,
artifact layouts, and exported JSON fields, see
`training/INPUTS_OUTPUTS.md`.

The training target is randomized treatment assignment:

```text
TREATMENT_GROUP
```

Models use baseline-to-follow-up omics changes, not baseline levels.

## Entrypoints

File: `training/main.R`

```r
FAST_treatment_ML(
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

`FAST_treatment_ML()` validates inputs, builds one modeling dataset per nonzero
follow-up, fits ENET and/or XGB, writes artifacts, and returns a manifest.

File: `R/write_training_artifacts.R`

```r
FAST_export_models(manifest, output_dir = NULL)
```

`FAST_export_models()` packages each fitted `(follow-up, model)` as a
self-contained JSON file for inference.

## File Map

| File | Role |
|:---|:---|
| `training/main.R` | Public training API and argument validation |
| `R/validate_inputs.R` | Phenotype, omics, omics-type, and DNAm validation |
| `R/build_training_features.R` | Follow-up datasets, change features, preprocessing, folds |
| `R/write_training_artifacts.R` | Disk outputs, ENET fitting, XGB worker launch, export packages |
| `training/scripts/run_xgb.py` | XGBoost tuning, fitting, and XGB artifacts |

## High-Level Flow

```text
FAST_treatment_ML()
├── validate function arguments
├── validate pheno and omics
├── apply DNAm probe validation/restriction when omics_type == "DNAm"
├── for each nonzero FU:
│   ├── require subjects with baseline + this FU
│   ├── validate treatment arms and CV fold counts
│   ├── compute follow-up minus baseline omics changes
│   ├── preprocess omics and covariates using this training cohort
│   ├── build ENET and XGB matrices
│   ├── write model-ready data
│   ├── fit requested models
│   └── collect report rows
├── write stacked reports
├── write manifest.json
└── return manifest
```

`FAST_export_models()` is separate. It reads the manifest and model artifacts,
then writes portable JSON packages under `models/exported_models/`.

## 1. Argument and Input Validation

`FAST_treatment_ML()` first calls `.validate_ml_args()` in `training/main.R`.
That function checks model names, fold counts, XGB repeat count, XGB trial count,
seed, and runtime settings before touching the input data.

Input preparation then flows through `.prepare_inputs()` in
`R/validate_inputs.R`:

```text
.prepare_inputs()
├── .validate_omics_type()
├── .validate_pheno()
├── .validate_omics()
└── [DNAm only] validate probe coverage and subset to reliable probes
```

Important phenotype behavior:

- `SAMPLE_ID` must be globally unique.
- `FU` must be consecutive integers starting at `0`.
- `FEMALE` and `TREATMENT_GROUP` must be binary `0/1`.
- Duplicate `SUBJECT_ID`/`FU` rows are reduced to the first row with a warning.
- `additional_covariates` must be numeric.
- Rows missing requested additional covariates are dropped.
- Subjects without both baseline and at least one follow-up are dropped.

Additional covariates are training-only ENET adjustment covariates. Raw
categorical covariates are rejected; callers should one-hot encode them upstream
if they want numeric indicators included.

Omics validation requires `ANALYTE_NAME`, numeric sample columns, and overlap
with validated phenotype sample IDs. Missing values and near-zero variance are
reported here, but the actual imputation and filtering happen per follow-up.

## 2. Follow-Up Loop

The main run body is `.run_ml_disk()` in `R/write_training_artifacts.R`.

For each nonzero `FU`, it:

1. Finds subjects with both baseline and that follow-up.
2. Checks that both treatment arms are present.
3. Checks that each requested CV fold count is no larger than the smaller arm.
4. Calls `.prepare_fu_change_dataset()`.
5. Writes prepared matrices and fold files.
6. Fits requested models.
7. Stores report rows and paths in the manifest.

Each follow-up is independent. A subject can contribute to one follow-up and be
absent from another if that subject is missing the corresponding visit.

## 3. Change Features

`.make_followup_change_data()` constructs the omics matrix for one follow-up.
For each subject:

```text
omics_change(FU k) = omics(FU k) - omics(FU 0)
```

The implementation is a matrix subtraction:

```r
change_matrix <- t(followup_values - baseline_values)
```

Rows are subjects. Columns are analytes. Treatment labels come from the
follow-up phenotype rows.

Baseline omics levels are not model features.

## 4. Preprocessing

Preprocessing is in `R/build_training_features.R`. Omics features and covariate
features go through the same numeric sequence:

1. Drop all-missing training columns.
2. Median-impute missing values.
3. Drop zero-variance training columns.
4. Center and scale retained columns.

All learned values come from the current follow-up's training cohort.

`.make_preprocessing_table()` writes the audit and replay recipe for every
candidate feature. The key flags are:

| Flag | Meaning |
|:---|:---|
| `IN_ENET` | Feature is retained in the ENET training matrix |
| `IN_XGB` | Feature is retained in the XGB training matrix |
| `DEPLOYABLE` | Feature is reconstructed during inference |

`DEPLOYABLE` is true for omics features and `covariate::FEMALE`. Other
additional covariates can adjust ENET during training, but future datasets are
not expected to contain them.

## 5. Model Matrices

`.prepare_fu_change_dataset()` creates model-specific matrices.

ENET receives:

```text
retained omics changes
+ retained FEMALE
+ retained numeric additional_covariates
```

XGB receives:

```text
retained omics changes
+ retained FEMALE
```

Feature names are prefixed:

```text
omics::<ANALYTE_NAME>
covariate::<covariate name>
```

`FEMALE` is added automatically to the model covariate list. If `FEMALE` has
zero variance within a follow-up, preprocessing removes it from both model
matrices.

During inference, only deployable retained rows are reconstructed. For ENET,
omitted non-`FEMALE` adjustment covariates are treated as standardized zero.

## 6. Cross-Validation Folds

`.stratified_subject_folds()` assigns folds separately within each treatment
arm.

ENET gets one fold assignment, stored as `ENET_FOLD_ID` in `subjects.csv`.
XGB gets `xgb_cv_repeats` independent fold assignments, stored in
`xgb_folds.csv`.

The follow-up seed is:

```text
seed + fu_level
```

XGB repeats use incremented seeds from that follow-up seed.

## 7. ENET Fitting

File: `R/write_training_artifacts.R`

`.run_enet_worker()` fits a binomial elastic net with `glmnet::cv.glmnet()`:

```r
glmnet::cv.glmnet(
  x = dataset$enet_x_train,
  y = dataset$y_train,
  family = "binomial",
  alpha = 0.5,
  foldid = dataset$enet_foldid,
  type.measure = "deviance",
  standardize = FALSE,
  penalty.factor = penalty_factor,
  keep = TRUE
)
```

The matrix is already standardized, so `standardize = FALSE`.

Penalty factors:

| Feature | Penalty |
|:---|---:|
| Omics | `1` |
| Covariates | `0` |

`lambda.min` is selected by cross-validated binomial deviance. `CV_AUC` is then
computed from out-of-fold predictions at `lambda.min`; AUC is reported but does
not select the lambda.

The worker writes metrics, training-cohort predictions, and a compact
`weights.csv` with the intercept, selected nonzero omics coefficients, and all
unpenalized covariate coefficients.

## 8. XGB Fitting

Files: `R/write_training_artifacts.R`, `training/scripts/run_xgb.py`

R launches the Python worker with the prepared follow-up directory, output
directory, prediction directory, seed, thread count, and trial count.

The Python worker reads:

```text
xgb_train.csv.gz
subjects.csv
xgb_folds.csv
```

It uses Optuna to tune a small XGBoost search space for binary logistic
classification with AUC evaluation. Each trial runs XGB CV once per predefined
repeat, with early stopping. The trial score is the mean best AUC across
repeats. The final model is fit using the selected parameters and the median
best iteration.

Search bounds:

| Parameter | Bounds |
|:---|:---|
| `max_depth` | Integer `1` to `3` |
| `eta` | Log-uniform `0.005` to `0.08` |
| `min_child_weight` | Log-uniform `5` to `50` |
| `subsample` | Uniform `0.6` to `0.9` |
| `colsample_bytree` | Uniform `0.1` to `0.6` |
| `lambda` | Log-uniform `5` to `200` |
| `alpha` | Log-uniform `0.05` to `30` |
| `gamma` | Uniform `0` to `5` |

The worker writes metrics, training-cohort predictions, gain importance, tuning
history, and the fitted XGBoost `model.json`.

## 9. Outputs and Manifest

Training writes exact model inputs before fitting models, then writes model
artifacts after each worker finishes. Stacked reports are written under
`models/reports/` after all modelable follow-ups complete.

`manifest.json` records:

- run settings;
- requested models and covariates;
- report paths;
- follow-up-specific data directories;
- follow-up-specific model artifact paths.

The returned R manifest is the same object plus `manifest_path`.

The complete artifact schema is in `training/INPUTS_OUTPUTS.md`.

## 10. Exported Model Packages

`FAST_export_models()` reads the manifest and writes one package per fitted
model. Each package embeds:

- model identity and source run metadata;
- training metrics and `successful = CV_AUC >= 0.8`;
- retained model-specific preprocessing rows;
- ENET weights or XGBoost model JSON;
- recorded covariate settings.

The package is the preferred artifact for inference. It avoids requiring the
whole training directory and keeps the preprocessing recipe beside the scoring
payload.

For exact JSON fields, see `training/INPUTS_OUTPUTS.md`.

## Design Invariants

- Each follow-up is modeled independently.
- All preprocessing parameters are learned within the follow-up training cohort.
- Model matrix column order is defined by preprocessing rows.
- XGB never receives additional covariates beyond retained `FEMALE`.
- Additional covariates are numeric, ENET-only, and training-only for inference.
- Export success is based on training `CV_AUC >= 0.8`; unsuccessful models are
  still exported and marked accordingly.
