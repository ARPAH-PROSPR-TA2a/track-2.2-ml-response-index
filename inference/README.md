# FAST Treatment ML Cross-Trial Validation

Cross-trial validation is a separate Track 2.2 tool. It takes an exported model,
applies it to another labeled trial, and reports how well the
treatment-associated signal transfers. It does not refit or tune the model.

The validation trial uses the same `pheno` and `omics` table layout as training.
Its phenotype table must include observed `TREATMENT_GROUP` labels. Additional
covariates used to adjust ENET during training are not needed here.

Run commands from the repository root.

## Validate Every Successful Model

Use this when evaluating a directory of exported models:

```r
source(file.path("inference", "main.R"))

validation <- FAST_bulk_evaluate(
  pheno = trial_b_pheno,
  omics = trial_b_omics,
  models_dir = "runs/trial_a_treatment_ml/models/exported_models",
  output_dir = "runs/trial_a_on_trial_b"
)

read.csv(validation$output_files$validation_summary)
```

`FAST_bulk_evaluate()` finds the exported model packages in `models_dir` and,
by default, evaluates those whose training CV AUC was at least 0.8. It writes:

- `validation_summary.csv`, with one row per evaluated model; and
- a `models/<MODEL_ID>/` directory containing that model's predictions and
  validation result.

If no exported model passed the training threshold, the summary is still
written with the expected columns but has zero rows.

## Validate One Model

Use `FAST_evaluate()` when you want to inspect one exported JSON package:

```r
source(file.path("inference", "main.R"))

result <- FAST_evaluate(
  pheno = trial_b_pheno,
  omics = trial_b_omics,
  model_path = "runs/trial_a_treatment_ml/models/exported_models/model_id.json",
  output_dir = "runs/trial_a_on_trial_b/model_id"
)

result$validation
```

The result includes one probability per subject and a one-row validation
summary. The summary reports the new trial's AUC and the association between the
model score and observed treatment status. Exact columns and statistical details
are in the [input/output reference](INPUTS_OUTPUTS.md).

## Requirements and Help

Validation runs entirely in R. ENET validation needs `jsonlite` and `pROC`;
XGBoost validation also needs the R `xgboost` package. It does not use the
Python environment that trained XGBoost.

To check those R packages before loading data:

```r
source(file.path("inference", "main.R"))
FAST_check_R("validation")
```

See the root [requirements](../README.md#requirements) and
[troubleshooting guide](../README.md#troubleshooting) for setup help. For exact
input columns, output files, and model-package behavior, use the documents
below.

## Documentation

- [Input and output reference](INPUTS_OUTPUTS.md)
- [Code walkthrough](CODE_WALKTHROUGH.md)
- [Training overview](../training/README.md)
- [Training inputs and outputs](../training/INPUTS_OUTPUTS.md)
