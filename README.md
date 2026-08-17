# FAST Track 2.2

FAST Track 2.2 implements machine learning to predict treatment status from
omics measurements.

You can use it to:

- [Train models](training/README.md) on one trial and export them for later use.
- [Test existing models](inference/README.md) in another trial to see how well
  they transfer.

Run all commands from the repository root: the folder containing `training/`,
`inference/`, and `R/`.

## Try the Example

The included [example](Examples/run_quickstart.R) checks your R and Python setup,
fits ENET and XGBoost models using simulated data, and exports the fitted models:

```bash
Rscript Examples/run_quickstart.R
```

FAST uses `python3` by default. To use a different Python installation, provide
its path:

```bash
Rscript Examples/run_quickstart.R "/full/path/to/python"
```

The example saves its results in a new `runs/quickstart_<timestamp>/` folder. It
uses 30 XGBoost tuning trials and one cross-validation repeat so that it finishes
more quickly. These reduced settings are for checking your setup, not for a
scientific analysis. See [Examples/README.md](Examples/README.md) for the files
it creates.

## If You Already Ran Track 1.1.1

Reuse the same `pheno` and `omics` tables that you supplied to
[Track 1.1.1](https://github.com/ARPAH-PROSPR-TA2a/track-1.1.1-single-analyte).
Do not use its results or reports as inputs. Proteomics, Metabolomics, and DNAm
tables that worked in Track 1.1.1 do not need to be reshaped.

Before you begin, note these Track 2.2 requirements:

- `omics_type` must be `"Proteomics"`, `"Metabolomics"`, or `"DNAm"`.
- For each follow-up, both treatment groups need at least as many usable subjects
  as the requested number of cross-validation folds. The default is 10 folds.
- The standard training run includes XGBoost and therefore requires Python.

See the [training input reference](training/INPUTS_OUTPUTS.md) for the full list
of required columns and accepted data types.

## Requirements

The standard Track 2.2 run fits both ENET and XGBoost, so it uses R and Python.
The packages needed for each task are:

| Task | R packages | Python packages |
|:---|:---|:---|
| Train ENET | `glmnet`, `jsonlite`, `pROC` | None |
| Train XGBoost | `jsonlite` | `numpy`, `pandas`, `scikit-learn`, `xgboost`, `optuna` |
| Validate an ENET model | `jsonlite`, `pROC` | None |
| Validate an XGBoost model | `jsonlite`, `pROC`, `xgboost` | None |

To use both training and cross-trial validation, install these four packages in
the R installation you will use to run FAST:

```r
install.packages(c("glmnet", "jsonlite", "pROC", "xgboost"))
```

Install these five packages using the same Python installation you will give to
FAST:

```bash
"/full/path/to/python" -m pip install numpy pandas scikit-learn xgboost optuna
```

XGBoost training uses Python. Cross-trial validation runs in R, even for an
XGBoost model. The R and Python `xgboost` packages must be installed separately.

### Check Your Setup

Check R and Python without starting a model run:

```r
source(file.path("training", "main.R"))
FAST_check_R()
FAST_check_python("/full/path/to/python")
```

The two checks report R and Python readiness separately. `FAST_check_R()` checks
everything needed for training and validation. To check only one task, use
`FAST_check_R("training")` or `FAST_check_R("validation")`.

If either check finds a problem, follow the instructions it prints or see
[Troubleshooting](#troubleshooting). FAST also checks Python automatically before
XGBoost training.

## Run Your Own Data

### Train and export models

```r
source(file.path("training", "main.R"))

# Use the Python installation checked above.
python_bin <- "/full/path/to/python"

manifest <- FAST_treatment_ML(
  pheno = trial_a_pheno,
  omics = trial_a_omics,
  omics_type = "Proteomics",
  models = c("enet", "xgb"),
  output_dir = "runs/trial_a_treatment_ml",
  python_bin = python_bin
)

exported <- FAST_export_models(manifest)
exported
```

FAST prints where it saved the training results and exported model files.

### Validate models in another trial

You can validate models created earlier or supplied by a collaborator. You do
not need to run training in the same R session.

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

By default, FAST validates only models whose training CV AUC was at least 0.8.

## Find Your Results

The main result files are:

- `runs/trial_a_treatment_ml/manifest.json`: the settings and file locations for
  a training run.
- `runs/trial_a_treatment_ml/models/reports/`: summary tables describing the
  included subjects, changes from baseline, preprocessing, and analyte names.
- `runs/trial_a_treatment_ml/models/FU*/`: ENET and XGBoost results for each
  follow-up visit.
- `runs/trial_a_treatment_ml/models/exported_models/`: the model files used for
  validation, plus `exported_models.csv`.
- `runs/trial_a_on_trial_b/validation_summary.csv`: the cross-trial validation
  results.

The [training](training/INPUTS_OUTPUTS.md) and
[validation](inference/INPUTS_OUTPUTS.md) references describe every file and
column in detail.

## Troubleshooting

- **FAST cannot find a file it needs.** Run `getwd()` in R and confirm that you
  are in the repository root, not in `training/`, `inference/`, or `Examples/`.

- **FAST says that an R package is missing.** Run `FAST_check_R()` to check both
  training and validation. Install missing packages in the same R installation
  you use to run FAST. The command is under [Requirements](#requirements).

- **FAST cannot find Python.** On macOS or Linux, open a terminal and run
  `which python3` and `which python`; on Windows, run `where python`. Pass the
  Python path you want to use to `FAST_check_python()`. You can also check
  from R with `Sys.which(c("python3", "python"))`. FAST does not silently switch
  environments. If no Python is found, use the official
  [Python downloads](https://www.python.org/downloads/).

- **FAST says that Python packages are missing.** Install them using the same
  Python path supplied as `python_bin`. The command is under
  [Requirements](#requirements). If needed, see the official guides to
  [virtual environments](https://docs.python.org/3/tutorial/venv.html) and
  [installing packages](https://packaging.python.org/en/latest/tutorials/installing-packages/).

- **FAST says there are not enough subjects.** Count usable subjects separately
  in each treatment group at each follow-up. The smaller group determines the
  largest possible number of cross-validation folds. ENET requires at least four
  folds, so both groups need at least four usable subjects.

- **Phenotype and omics samples do not match.** Each `SAMPLE_ID` in `pheno`
  should exactly match a sample-column name in `omics`, including capitalization
  and punctuation.

- **`validation_summary.csv` is empty.** Check
  `models/exported_models/exported_models.csv`. By default, FAST skips models
  whose training CV AUC was below 0.8.

## More Detail

- [Training overview](training/README.md)
- [Cross-trial validation overview](inference/README.md)
- [Training inputs and outputs](training/INPUTS_OUTPUTS.md)
- [Training code walkthrough](training/CODE_WALKTHROUGH.md)
- [Validation inputs and outputs](inference/INPUTS_OUTPUTS.md)
- [Validation code walkthrough](inference/CODE_WALKTHROUGH.md)

## Tests

The full test suite runs ENET, XGBoost, model export, and validation. It requires
all packages listed above:

```bash
Rscript tests/run_tests.R
```
