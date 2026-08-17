# Quick Start

This example is a short environment and training smoke test. It checks the R
and Python environments separately, fits both ENET and XGBoost to bundled
synthetic Proteomics data, exports both models, and checks the resulting file
structure. Its reduced tuning settings are not scientific analysis defaults.

Run from the repository root. To use `python3` from `PATH`:

```bash
Rscript Examples/run_quickstart.R
```

To use a specific Python environment, pass its interpreter exactly as you would
to `FAST_treatment_ML()`:

```bash
Rscript Examples/run_quickstart.R /path/to/python
```

The script first runs `FAST_check_R("training")`. The XGBoost preflight then
checks that the selected Python can import `numpy`, `pandas`, `sklearn`,
`xgboost`, and `optuna`. You can run either check separately after sourcing
`training/main.R`:

```r
FAST_check_R("training")
FAST_check_python("/path/to/python")
```

The smoke test uses one follow-up, five CV folds, one XGB CV repeat, 30 Optuna
trials, one core, and no additional covariates. The single XGB CV repeat keeps
this example short; scientific analyses should use the documented default of
three repeats.

Each run writes to a fresh `runs/quickstart_YYYYMMDD_HHMMSS/` directory, adding a
numeric suffix if needed. Key outputs are:

```text
manifest.json
models/FU1/enet/metrics.csv
models/FU1/enet/weights.csv
models/FU1/xgb/metrics.csv
models/FU1/xgb/tuning.csv
models/FU1/xgb/model.json
models/exported_models/exported_models.csv
models/exported_models/<model_id>.json
```

Successful completion means that both model families ran and two portable model
packages were exported. The printed AUC values come from synthetic data and are
not scientific results.
