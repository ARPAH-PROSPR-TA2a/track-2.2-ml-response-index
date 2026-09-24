# CALERIE DNAm beta training on the cloud VM

Run as the existing `rstudio` user. This uses the same `~/FAST` layout as Track 3.3.
If SSH signs you in as `ubuntu`, first run `sudo -iu rstudio`; the RStudio terminal
already uses the RStudio account.
The runner fits treatment-prediction models from **follow-up minus baseline beta
values**, separately at 12 months (`FU1`) and 24 months (`FU2`). It does not apply
INF/MetS outcome restrictions. Each visit uses its own baseline-paired cohort.

## Prepare the checkout and environment

Use the checkout at `~/FAST/GitHub/track-2.2`. Get these changes into that server
checkout before launching. For a first checkout, once the runner is available in
the remote repository:

```bash
mkdir -p "$HOME/FAST/GitHub"
git clone https://github.com/ARPAH-PROSPR-TA2a/track-2.2-ml-response-index.git "$HOME/FAST/GitHub/track-2.2"
```

In the VM's RStudio console, install any missing R dependencies:

```r
install.packages(c("glmnet", "jsonlite", "pROC", "haven"),
                 repos = "https://cloud.r-project.org")
```

In the VM terminal, create the Python environment and check both runtimes with
the bundled synthetic quickstart:

```bash
python3 -m venv "$HOME/FAST/Envs/track22"
"$HOME/FAST/Envs/track22/bin/python" -m pip install numpy pandas scikit-learn xgboost optuna
cd "$HOME/FAST/GitHub/track-2.2"
Rscript Examples/run_quickstart.R "$HOME/FAST/Envs/track22/bin/python"
```

The quickstart trains and exports small synthetic models under `runs/`. Continue
after it reports `Quick start complete.` If `venv` or `tmux` is unavailable,
install those VM prerequisites before proceeding.

## Review the runner settings

Edit the paths and settings at the top of `run_DNAm_betas.R` before launch:

| Setting | Default |
|:---|:---|
| Repository | `~/FAST/GitHub/track-2.2` |
| Beta matrix | `~/FAST/Data/CALERIE/Raw/DNAm/GRSet_fully_filtered_bmiq_chunk.rds` |
| Phenotypes | `~/FAST/Data/CALERIE/Raw/DNAm/CALERIE_CPR_processed_pheno.rds` |
| Python | `~/FAST/Envs/track22/bin/python` |
| Results | `~/FAST/Outputs/2.2/DNAm_betas_2.2_Output` |
| XGBoost threads | `7` on the 8-vCPU machine |
| ENET CV | `10` folds |
| XGBoost tuning | `100` trials, `10` folds, `3` repeats |
| Additional covariates | `NULL`; both models include retained `FEMALE` |

This runner uses 100 trials per follow-up and retains the pipeline's 10-fold,
3-repeat XGBoost CV. Seven threads leave headroom on the 8-vCPU machine; adjust
this setting if the allocation changes. Adding covariates also requires loading
those fields into the runner's `pheno` table.

Inputs must already be available on the VM. The runner requires a numeric
CpG-by-sample matrix with complete, finite beta values in `[0, 1]`, and phenotype
columns `Barcode`, `Participant_ID`, `fu`, `CR`, and `female`. Visit codes are
`0/1/2`; treatment and sex codes are `0/1`. It checks sample matching, duplicates,
and both follow-up cohorts before fitting. The core pipeline selects the reliable
DNAm probe set from the supplied matrix.

## Launch and monitor

Start a persistent Bash session:

```bash
tmux new -s track22 bash
```

Inside that session:

```bash
cd "$HOME/FAST/GitHub/track-2.2"
mkdir -p "$HOME/FAST/Outputs/2.2"
set -o pipefail
console_log="$HOME/FAST/Outputs/2.2/console_DNAm_betas_$(date +%Y%m%d_%H%M%S).log"
Rscript training/CALERIE/DNAm/run_DNAm_betas.R 2>&1 | tee "$console_log"
run_status=$?
printf 'Run exit status: %s\nConsole log: %s\n' "$run_status" "$console_log"
```

`pipefail` ensures an R error produces a failed pipeline status even when `tee`
succeeds. A zero status and the runner's `COMPLETE` message indicate completion.
Detach with **Ctrl-b, then d**; reconnect with `tmux attach -t track22`.

The console log captures stdout/stderr, including errors. The result directory
also contains a concise `run.log` and `provenance.rds`, alongside the
pipeline's `manifest.json`, `data/`, and `models/`. Exported model packages and their
index are under `models/exported_models/`; a complete run exports four models.

The current XGBoost worker is quiet during tuning; a lack of new log lines alone
does not mean it has stalled. This runner has **no checkpoint/resume** and stops
if the result directory already exists. For a rerun, choose a new `out_dir` at the
top of the script; training starts from the beginning. The runner creates its
result directory, so do not create `DNAm_betas_2.2_Output` before launching.
