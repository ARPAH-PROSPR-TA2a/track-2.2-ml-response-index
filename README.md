# Treatment ML Pipeline

Track 2.2 predicts assigned treatment from omics change signatures in randomized
trial datasets. The target is `TREATMENT_GROUP`.

The pipeline exposes:

- `FAST_treatment_ML()` for model fitting.
- `FAST_treatment_ML_reports()` for QC/data summaries and randomization reports.

## Quick Example

```r
source("main.R")

manifest <- FAST_treatment_ML(
  pheno = pheno,
  omics = omics,
  omics_type = "Proteomics",
  additional_covariates = c("agebl", "mbmi"),
  models = c("enet", "xgb"),
  output_dir = "runs/trial_a_treatment_ml",
  test_frac = 0.2,
  enet_cv_folds = 5,
  xgb_cv_folds = 5,
  xgb_cv_repeats = 3,
  xgb_n_trials = 50,
  n_cores = 8,
  seed = 1
)

read.csv("runs/trial_a_treatment_ml/FU1/enet/metrics.csv")
```

## Modeling

One all-subject model is fit per nonzero follow-up. There are no sex-stratified
models. Omics features are:

```text
omics(FU k) - omics(FU 0)
```

Train/test splits and CV folds are subject-level and stratified by treatment.

ENET receives omics changes, `FEMALE`, and every requested
`additional_covariates` variable. XGB receives omics changes plus `FEMALE`;
other covariates are excluded from XGB. If `FEMALE` has zero variance in the
training set for a follow-up, preprocessing removes it from both models.

Feature columns use `omics::` and `covariate::` prefixes.

## Python

XGB uses `python3` from `PATH` by default. Activate the intended Python
environment before starting R, or pass its interpreter explicitly:

```r
FAST_treatment_ML(
  ...,
  python_bin = "/path/to/python"
)
```

The Python environment must provide `numpy`, `pandas`, `scikit-learn`,
`xgboost`, and `optuna`.

## Models

### ENET

`glmnet::cv.glmnet()` fits a binomial elastic-net model. Omics coefficients are
penalized; `FEMALE` and requested covariates are included with
`penalty.factor = 0`. The final model uses `lambda.min`, selected by minimum
cross-validated binomial deviance; out-of-fold AUC is reported at that lambda.
Outputs include metrics, subject predictions, and coefficients sufficient to
reproduce predictions without a saved R model object.

### XGB

XGB runs through `scripts/run_xgb.py`. Each parameter set is evaluated across
`xgb_cv_repeats` independent stratified fold assignments and scored by mean CV
AUC. By default, 50 Optuna trials are evaluated across three 5-fold repeats.
XGB always uses Optuna tuning: at least 10 trials are required, and fewer than
30 produce a limited-search warning. Outputs include metrics and selected
parameters, predictions, feature importance, tuning history, and the fitted
JSON model.

## DNAm

DNAm inputs are restricted before modeling to:

```text
Data/FAST_epicv1_epicv2_sugden_TruD_probe_list.rds
```

## Output

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

The model-ready matrices contain exactly the columns consumed by each model.
`subjects.csv` combines outcomes, train/test assignment, and ENET fold
assignment. `xgb_folds.csv` stores every repeated XGB fold assignment. Run
settings and artifact paths are stored once in `manifest.json`.

See `INPUTS_OUTPUTS.md` for schemas.

## Tests

```bash
Rscript run_tests.R
```

The suite requires both ENET and XGB and leaves the integration artifacts at:

```text
test_outputs/track22_integration/
```

For a larger demo-scale smoke run with inspectable, analysis-like outputs:

```bash
Rscript test_demo_run.R
```

This writes to:

```text
test_outputs/track22_demo/
```

## Reports

```r
reports <- FAST_treatment_ML_reports(
  pheno = pheno,
  omics = omics,
  omics_type = "Proteomics",
  additional_covariates = c("agebl", "mbmi")
)
```

Reports return `pheno_summary`, `variable_summaries`, and
`randomization_reports`. Reporting may summarize by sex; modeling is not
sex-stratified.

## Documentation

- [INPUTS_OUTPUTS.md](INPUTS_OUTPUTS.md): input requirements and artifact schemas.
- [CODE_WALKTHROUGH.md](CODE_WALKTHROUGH.md): implementation flow and model details.
