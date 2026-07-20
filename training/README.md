# FAST Treatment ML Training

Track 2.2 predicts assigned treatment from omics change signatures in
randomized trial datasets. The target is `TREATMENT_GROUP`.

The training pipeline exposes:

- `FAST_treatment_ML()` for model fitting and ML-specific reports.
- `FAST_export_models()` for self-contained JSON model packages.

## Quick Example

Run from the repository root:

```r
source(file.path("training", "main.R"))

manifest <- FAST_treatment_ML(
  pheno = pheno,
  omics = omics,
  omics_type = "Proteomics",
  additional_covariates = c("agebl", "mbmi"),
  models = c("enet", "xgb"),
  output_dir = "runs/trial_a_treatment_ml",
  enet_cv_folds = 10,
  xgb_cv_folds = 10,
  xgb_cv_repeats = 3,
  xgb_n_trials = 50,
  n_cores = 8,
  seed = 1
)

read.csv("runs/trial_a_treatment_ml/models/FU1/enet/metrics.csv")

exported <- FAST_export_models(manifest)
```

## Modeling

One model is fit per nonzero follow-up. There are no sex-stratified
models. Omics features are:

```text
omics(FU k) - omics(FU 0)
```

CV folds are subject-level and stratified by treatment. Each follow-up must have
at least as many subjects in each treatment arm as the requested fold count.

ENET receives omics changes, `FEMALE`, and every requested
`additional_covariates` variable as training-only adjustment features. Numeric,
integer, unordered factor, and logical covariates are accepted. Character and
ordered-factor covariates are rejected. XGB receives omics changes plus
`FEMALE`; other covariates are excluded from XGB. If `FEMALE` has zero variance
in the training set for a follow-up, preprocessing removes it from both models.

### Categorical Covariates

Unordered factors are converted to treatment contrasts independently within
each follow-up. A factor with `k` observed levels produces `k - 1` columns, with
the first declared factor level as the reference. Set factor levels explicitly
when the reference matters:

```r
pheno$site <- factor(pheno$site, levels = c("A", "B", "C"))
```

This produces `covariate::siteB` and `covariate::siteC`; level `A` is the
reference. Logical covariates use `FALSE` as the reference and produce a
`<name>TRUE` column. Unused factor levels are dropped for each follow-up, and a
factor with only one observed level is dropped with a message.

Encoded column names are normalized with `make.names(..., unique = TRUE)`.
Categorical columns are centered and scaled like numeric covariates. They are
unpenalized ENET adjustment features and are intentionally excluded from XGB
and deployment scoring.

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

The Python environment must have `numpy`, `pandas`, `scikit-learn`,
`xgboost`, and `optuna` installed.

## Models

Models with `CV_AUC >= 0.8` are marked successful for cross-trial validation.
`FAST_export_models()` writes one self-contained JSON package per fitted model
under `models/exported_models/`; each package records whether it passed the
success threshold.

### ENET

`glmnet::cv.glmnet()` fits a binomial elastic-net model. Omics coefficients are
penalized; `FEMALE` and requested covariates are included with
`penalty.factor = 0`. The final model uses `lambda.min`, selected with
cross-validated AUC (`type.measure = "auc"`); pooled out-of-fold AUC is reported
at that lambda.
Outputs include metrics, subject predictions, and coefficients sufficient to
reproduce deployment predictions without a saved R model object. Deployment
scoring omits training-only adjustment covariates, equivalent to setting their
standardized values to zero.

### XGB

XGB runs through `training/scripts/run_xgb.py`. Each parameter set is evaluated across
`xgb_cv_repeats` independent stratified fold assignments and scored by mean CV
AUC. By default, 50 Optuna trials are evaluated across three 10-fold repeats.
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

The model-ready matrices contain exactly the columns consumed by each model.
`subjects.csv` combines outcomes and ENET fold assignment. `xgb_folds.csv`
stores every repeated XGB fold assignment. These subject-level artifacts are
written under `data/`.
Reports under `models/reports/` stack information across follow-ups:
`cohort.csv` records modeled cohort counts, `change_summary.csv` describes raw
change scores by treatment arm, and `preprocessing.csv` audits feature
transformations and removals while recording the preprocessing recipe needed to
reconstruct model matrices. Run settings and artifact paths are stored once in
`manifest.json`.

Exported model JSON packages contain the model-specific preprocessing recipe,
training metrics, covariates, a `successful` flag, and either ENET weights or
embedded XGBoost model JSON. They are the preferred input for single-model
cross-trial evaluation. They are written only after `FAST_export_models()` is
called:

```text
output_dir/
  models/
    exported_models/
      exported_models.csv
      <model_id>.json
```

See [INPUTS_OUTPUTS.md](INPUTS_OUTPUTS.md) for schemas.

## Tests

```bash
Rscript tests/run_tests.R
```

The suite requires both ENET and XGB and leaves the integration artifacts at:

```text
tests/test_outputs/track22_integration/
```

## Documentation

- [INPUTS_OUTPUTS.md](INPUTS_OUTPUTS.md): input requirements and artifact schemas.
- [CODE_WALKTHROUGH.md](CODE_WALKTHROUGH.md): implementation flow and model details.
