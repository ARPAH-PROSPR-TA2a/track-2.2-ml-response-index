# FAST Treatment ML Training

Training is one Track 2.2 tool. It uses a labeled randomized trial to learn
whether treatment and control subjects can be distinguished by how their omics
measurements changed after baseline.

The training code:

1. checks the phenotype and omics inputs;
2. calculates follow-up minus baseline change for each omics feature;
3. fits one ENET and one XGBoost model for each follow-up; and
4. writes review files and exports portable model packages for another trial.

## Run Training

Run from the repository root. The
[self-contained quick start](../README.md#quick-start) uses bundled example data
and fits both models. With your own `pheno` and `omics` objects, the main call is:

```r
source(file.path("training", "main.R"))

manifest <- FAST_treatment_ML(
  pheno = pheno,
  omics = omics,
  omics_type = "Proteomics",
  models = c("enet", "xgb"),
  output_dir = "runs/trial_a_treatment_ml",
  python_bin = "/full/path/to/python"
)

exported <- FAST_export_models(manifest)
exported
```

`FAST_treatment_ML()` fits the requested models and returns a manifest of the
run. `FAST_export_models()` then writes one self-contained JSON package per
fitted model.

## Inputs at a Glance

`pheno` has one row per sample and identifies the subject, visit, assigned
treatment, and `FEMALE`. `omics` has one row per analyte and one numeric column
per sample. Omics sample-column names must match `pheno$SAMPLE_ID`.

If you already ran Track 1.1.1, use the same original `pheno` and `omics` input
tables, not its results. See the [root guidance](../README.md#if-you-already-ran-track-111)
and the [exact input requirements](INPUTS_OUTPUTS.md#inputs).

## What the Models Use

Each nonzero follow-up is modeled separately. Both models use the change from
baseline:

```text
omics at follow-up - omics at baseline
```

ENET uses the omics changes, `FEMALE`, and any requested
`additional_covariates`. XGBoost uses the omics changes and `FEMALE`. There are
no sex-stratified models.

Cross-validation folds are assigned by subject and balanced by treatment arm.
At every follow-up, each arm must have at least as many usable subjects as the
requested fold count.

Models with training CV AUC of at least 0.8 are marked successful for
cross-trial validation. All fitted models can be exported, but bulk validation
uses successful packages by default.

For exact covariate rules, preprocessing, tuning, and model behavior, use the
[input/output reference](INPUTS_OUTPUTS.md) and
[code walkthrough](CODE_WALKTHROUGH.md).

## Check R and Python

Check the R packages needed for training:

```r
FAST_check_R("training")
```

ENET runs entirely in R. When XGBoost is requested, R starts the selected
Python executable and runs `training/scripts/run_xgb.py`. That environment must
contain `numpy`, `pandas`, `scikit-learn`, `xgboost`, and `optuna`.

`FAST_treatment_ML()` checks the selected environment before starting XGBoost.
For a standalone setup check, source `training/main.R` and run:

```r
FAST_check_python("/full/path/to/python")
```

The checks report R and Python readiness separately. If either stops, follow its
next steps or see the root
[requirements](../README.md#requirements) and
[troubleshooting guide](../README.md#troubleshooting) for installation help.

## Where Results Go

Inside `output_dir`:

- `manifest.json` records the run settings and artifact paths.
- `models/reports/` contains cohort, change, preprocessing, and analyte-name
  reports.
- `models/FU*/` contains ENET and XGBoost results for each follow-up.
- `data/FU*/` contains model-ready matrices, folds, subject maps, and
  predictions.
- `models/exported_models/` is created by `FAST_export_models()` and contains
  the portable model packages.

DNAm inputs are restricted to the reliable probe list in `Data/` before models
are fit. The [input/output reference](INPUTS_OUTPUTS.md) gives the complete file
layout and schemas.

## Documentation

- [Input and output reference](INPUTS_OUTPUTS.md): exact data requirements and
  artifact schemas.
- [Code walkthrough](CODE_WALKTHROUGH.md): implementation flow and technical
  model details.
- [Cross-trial validation overview](../inference/README.md): how to use the
  exported models.

## Tests

From the repository root:

```bash
Rscript tests/run_tests.R
```

The full suite requires the ENET, XGBoost, export, and validation dependencies
listed in the [root README](../README.md#requirements).
