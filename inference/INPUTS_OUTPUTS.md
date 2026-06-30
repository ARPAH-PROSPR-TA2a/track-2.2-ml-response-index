# FAST Inference Inputs and Outputs

Inference uses the same validated phenotype and omics table conventions as
training. The important difference is that inference data must also include
observed `TREATMENT_GROUP` labels, because validation reports compare model
scores against treatment status.

The preferred API evaluates self-contained model JSON packages written by
`FAST_export_models()`. The exported package schema is owned by
`training/INPUTS_OUTPUTS.md`; this document focuses on what inference consumes
from those packages and what inference writes.

## Functions

### `FAST_evaluate()`

```r
FAST_evaluate(
  pheno,
  omics,
  model_path,
  output_dir = NULL,
  return_matrix = FALSE
)
```

Inputs:

| Argument | Meaning |
|:---|:---|
| `pheno` | Phenotype table using the same columns and validation rules as training |
| `omics` | Long omics table using the same columns and validation rules as training |
| `model_path` | Path to one exported model JSON package |
| `output_dir` | Optional directory for `predictions.csv` and `validation.csv` |
| `return_matrix` | If `TRUE`, return the reconstructed model matrix in memory |

`pheno` must include observed `TREATMENT_GROUP`. The evaluated cohort is the
set of subjects with baseline and the package follow-up.

Returns a list with:

| Element | Meaning |
|:---|:---|
| `model_id` | Exported model identifier |
| `predictions` | Per-subject prediction table |
| `validation` | One-row validation table |
| `model` | Parsed exported model package |
| `output_dir` | Normalized output directory, or `NULL` |
| `output_files` | Written file paths, or `NULL` |
| `matrix` | Reconstructed matrix, only when `return_matrix = TRUE` |

### `FAST_bulk_evaluate()`

```r
FAST_bulk_evaluate(
  pheno,
  omics,
  models_dir,
  output_dir
)
```

Inputs:

| Argument | Meaning |
|:---|:---|
| `pheno` | Labeled phenotype table |
| `omics` | Omics table |
| `models_dir` | Directory containing exported model JSONs and `exported_models.csv` |
| `output_dir` | Directory for bulk validation outputs |

Bulk evaluation reads `models_dir/exported_models.csv`, evaluates only rows
where `SUCCESSFUL == TRUE`, and stops if an expected model package cannot be
scored.

Returns a list with:

| Element | Meaning |
|:---|:---|
| `validation` | Bulk validation summary |
| `output_dir` | Normalized output directory |
| `output_files$validation_summary` | Path to `validation_summary.csv` |

### `FAST_treatment_predict()`

```r
FAST_treatment_predict(
  pheno,
  omics,
  manifest_path,
  models = NULL,
  followups = NULL,
  omics_type = NULL,
  enet_mode = "reproduce_training",
  output_dir = NULL,
  return_matrices = FALSE
)
```

This is the older manifest-based replay API. It scores a full training manifest
using model artifacts and `models/reports/preprocessing.csv`. It remains useful
for tests and whole-run replay, but cross-trial validation should generally use
exported model JSON packages.

## Single-Model Output

When `FAST_evaluate(..., output_dir = "path")` is used:

```text
output_dir/
  predictions.csv
  validation.csv
```

### `predictions.csv`

| Column | Meaning |
|:---|:---|
| `SUBJECT_ID` | Evaluated subject |
| `FU` | Follow-up required by the model package |
| `TREATMENT_GROUP` | Observed treatment label |
| `PREDICTED_PROB` | Model score/probability for treatment status |

### `validation.csv`

| Column | Meaning |
|:---|:---|
| `FU` | Follow-up evaluated |
| `MODEL` | `enet` or `xgb` |
| `TRAINING_CV_AUC` | CV AUC recorded in the exported model package |
| `SUCCESS_AUC_THRESHOLD` | Training success threshold, currently `0.8` |
| `SUCCESSFUL` | `TRUE` when `TRAINING_CV_AUC >= 0.8` |
| `N` | Evaluated subjects |
| `N_CONTROL` | Evaluated control subjects |
| `N_TREATMENT` | Evaluated treatment subjects |
| `AUC` | Validation AUC, computed only for successful models |
| `LOGIT_BETA` | Coefficient for `PREDICTED_PROB` in the validation logistic model |
| `LOGIT_OR` | `exp(LOGIT_BETA)` |
| `LOGIT_P` | Wald p-value for `PREDICTED_PROB` |
| `VALIDATED_P05` | `TRUE` when `LOGIT_P < 0.05` |

Validation uses:

```text
TREATMENT_GROUP ~ PREDICTED_PROB
```

If the training model did not pass the CV AUC threshold, `AUC`, `LOGIT_BETA`,
`LOGIT_OR`, and `LOGIT_P` are written as missing values.

## Bulk Output

When `FAST_bulk_evaluate(..., output_dir = "path")` is used:

```text
output_dir/
  validation_summary.csv
  models/
    <MODEL_ID>/
      predictions.csv
      validation.csv
```

Each per-model directory contains the same `predictions.csv` and
`validation.csv` files described above.

### `validation_summary.csv`

| Column | Meaning |
|:---|:---|
| `MODEL_ID` | Exported model identifier |
| `MODEL_PATH` | Model package path, relative to bulk `output_dir` when possible |
| `FU` | Follow-up evaluated |
| `MODEL` | `enet` or `xgb` |
| `TRAINING_CV_AUC` | CV AUC recorded in the exported model package |
| `SUCCESS_AUC_THRESHOLD` | Training success threshold, currently `0.8` |
| `SUCCESSFUL` | Always `TRUE` for evaluated rows |
| `N` | Evaluated subjects |
| `N_CONTROL` | Evaluated control subjects |
| `N_TREATMENT` | Evaluated treatment subjects |
| `AUC` | Validation AUC |
| `LOGIT_BETA` | Validation logistic coefficient |
| `LOGIT_OR` | Validation odds ratio |
| `LOGIT_P` | Validation p-value |
| `VALIDATED_P05` | `TRUE` when `LOGIT_P < 0.05` |
| `PREDICTIONS_PATH` | Per-model predictions path |
| `VALIDATION_PATH` | Per-model validation path |

If no exported models are successful, `validation_summary.csv` is written with
the schema above and zero rows.

## Manifest Replay Output

When `FAST_treatment_predict(..., output_dir = "path")` is used:

```text
output_dir/
  validation.csv
  predictions/
    FU1/
      enet.csv
      xgb.csv
```

`validation.csv` has the same validation schema as single-model evaluation.
Each prediction file has the same prediction schema as `predictions.csv`.
