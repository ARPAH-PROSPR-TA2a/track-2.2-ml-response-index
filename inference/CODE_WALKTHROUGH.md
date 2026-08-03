# FAST Inference Code Walkthrough

Inference reuses the training validation and feature-reconstruction rules. It
does not relearn preprocessing parameters. Instead, it rebuilds the model matrix
from the saved training recipe, scores the saved model, and validates the score
against observed treatment labels in the new cohort.

The preferred entry points are:

```r
FAST_evaluate()
FAST_bulk_evaluate()
```

## Overview

```text
pheno + omics
  │
  ├── validate inputs using training-compatible rules
  │
  ├── load exported model package
  │
  ├── rebuild one follow-up model matrix
  │     ├── baseline-to-follow-up omics changes
  │     ├── training medians
  │     ├── training centers/scales
  │     └── retained feature order
  │
  ├── score model
  │     ├── ENET: weights.csv-equivalent table
  │     └── XGB: embedded model JSON
  │
  └── validate predictions
        ├── AUC
        └── TREATMENT_GROUP ~ PREDICTED_PROB
```

## File Loading

File: `inference/main.R`

```r
source(file.path("R", "validate_inputs.R"))
source(file.path("R", "build_training_features.R"))
source(file.path("R", "replay_preprocessing.R"))
source(file.path("R", "score_models.R"))
```

Inference intentionally depends on the same input validation and preprocessing
replay helpers as training. That keeps training and validation cohorts aligned:
the same sample identifiers, follow-up encoding, omics change definition,
covariate handling, and feature naming conventions are used.

## `FAST_evaluate()`

File: `inference/main.R`

```r
FAST_evaluate <- function(
  pheno,
  omics,
  model_path,
  output_dir = NULL,
  return_matrix = FALSE
)
```

This is the single-model evaluation API. It expects one exported JSON model
package from `FAST_export_models()`.

### Step 1: Load the Package

```r
model_package <- .load_exported_model_package(model_path)
```

The loader checks:

- JSON support is available through `jsonlite`;
- `schema_version` is `"1.0"`;
- `family` is either `enet` or `xgb`.

The package carries the follow-up, model family, omics type, covariate lists,
model-specific preprocessing rows, training CV AUC, and the fitted model
payload.

### Step 2: Validate Inputs

```r
inputs <- .prepare_inputs(
  pheno = pheno,
  omics = omics,
  omics_type = model_package$omics_type,
  additional_covariates = NULL
)
```

This uses the same validation path as training. Inference additionally depends
on observed `TREATMENT_GROUP`, because validation statistics require labels.
Training-only `additional_covariates` are not required for inference. Raw
analyte names are mapped to the same internal XGB-safe form used during
training before the model matrix is rebuilt.

### Step 3: Rebuild the Model Matrix

```r
replay <- .build_inference_matrix(
  pheno_df = inputs$pheno,
  omics_df = inputs$omics,
  preprocessing = preprocessing,
  fu_level = as.integer(model_package$fu),
  model = model_package$family,
  model_covariates = model_covariates
)
```

The replay step differs from training in one key way: it does not decide which
features to retain. It applies the retained feature list and learned parameters
from training, then keeps only deployable rows. This includes:

- omics change features, defined as follow-up minus baseline;
- feature removal decisions;
- median imputation values;
- centering and scaling values;
- model-specific feature order.

Deployable rows are omics features and retained `FEMALE`. Other ENET covariates
are training-only adjustment features and are omitted from inference, equivalent
to standardized value zero.

If `return_matrix = TRUE`, this matrix is returned for audit/debugging.

### Step 4: Score the Model

```r
pred <- .score_exported_model_package(model_package, replay$matrix)
```

Scoring dispatches by model family:

- ENET uses the embedded weight table and computes the logistic response.
- XGB writes the embedded JSON model to a temporary file and scores with
  `xgboost`.

The prediction table is:

```r
data.frame(
  SUBJECT_ID = replay$subject_ids,
  FU = as.integer(model_package$fu),
  TREATMENT_GROUP = replay$y,
  PREDICTED_PROB = pred
)
```

### Step 5: Validate the Score

File: `R/score_models.R`

```r
validation <- .validation_row(
  fu_level = as.integer(model_package$fu),
  model_name = model_package$family,
  predictions = predictions,
  training_cv_auc = as.numeric(model_package$training_cv_auc)
)
```

Validation first checks the training success rule:

```text
TRAINING_CV_AUC >= 0.8
```

Only successful models get validation AUC and logistic association statistics.
The logistic model is:

```text
TREATMENT_GROUP ~ PREDICTED_PROB
```

`VALIDATED_P05` is `TRUE` when the score association p-value is below `0.05`.

### Step 6: Write Outputs

If `output_dir` is supplied, `FAST_evaluate()` writes:

```text
output_dir/
  predictions.csv
  validation.csv
```

Otherwise results are returned in memory only.

## `FAST_bulk_evaluate()`

File: `inference/main.R`

```r
FAST_bulk_evaluate <- function(
  pheno,
  omics,
  models_dir,
  output_dir,
  only_successful = TRUE
)
```

Bulk evaluation is a thin loop over single-model evaluation. It expects a
directory of exported model JSON packages.

### Step 1: Find Model Packages

```r
model_paths <- .list_exported_model_packages(models_dir)
```

Only `*.json` files in `models_dir` are considered. `exported_models.csv` may
exist as an export audit artifact, but bulk evaluation does not read it.

### Step 2: Keep Successful Models

```r
if (only_successful && !isTRUE(model_package$successful)) next
```

Failed training models are not evaluated in bulk by default. They can still
remain in the directory, with `successful == false`, for audit and transfer.

### Step 3: Evaluate Each Model

For each selected package, bulk evaluation calls the same package-level helper
used by `FAST_evaluate()`:

```r
.evaluate_model_package(
  pheno = pheno,
  omics = omics,
  model_package = model_package,
  output_dir = model_output_dir
)
```

Errors are rethrown with the model ID and path, because package failures are
expected to be actionable by the caller.

### Step 4: Write the Summary

Bulk evaluation writes:

```text
output_dir/
  validation_summary.csv
  models/
    <MODEL_ID>/
      predictions.csv
      validation.csv
```

`validation_summary.csv` is the per-model validation row plus the model ID,
model package path, and relative paths to the per-model output files.

If no models are successful, the summary is still written with the expected
columns and zero rows.

## Shared Behavior With Training

Inference shares these concepts with training:

- input validation rules;
- baseline-to-follow-up omics change calculation;
- exact feature naming conventions;
- training-derived preprocessing recipe;
- ENET and XGB scoring semantics;
- CV AUC as the training success criterion.

Inference differs in these ways:

- it requires observed labels for validation;
- it does not fit models;
- it does not learn preprocessing parameters;
- it evaluates only successful models in bulk;
- it reports validation association p-values for downstream cross-trial review.
