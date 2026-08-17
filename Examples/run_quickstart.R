args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L || (length(args) == 1L && !nzchar(args[[1L]]))) {
  stop(
    "Usage: Rscript Examples/run_quickstart.R [PYTHON_EXECUTABLE]",
    call. = FALSE
  )
}

required_project_files <- c(
  file.path("training", "main.R"),
  file.path("Examples", "ExampleData", "pheno_quickstart.rds"),
  file.path("Examples", "ExampleData", "proteomics_quickstart.rds")
)
if (!all(file.exists(required_project_files))) {
  stop("Run this script from the Track 2.2 repository root.", call. = FALSE)
}

python_bin <- if (length(args) == 1L) args[[1L]] else NULL

source(file.path("training", "main.R"))

FAST_check_R("training")

pheno <- readRDS(file.path("Examples", "ExampleData", "pheno_quickstart.rds"))
omics <- readRDS(file.path("Examples", "ExampleData", "proteomics_quickstart.rds"))

assert <- function(condition, message) {
  if (!isTRUE(condition)) {
    stop("Quick-start check failed: ", message, call. = FALSE)
  }
}

assert(
  identical(
    names(pheno),
    c("SAMPLE_ID", "SUBJECT_ID", "FU", "TREATMENT_GROUP", "FEMALE")
  ),
  "unexpected phenotype columns"
)
assert(nrow(pheno) == 240L, "unexpected phenotype row count")
assert(nrow(omics) == 24L, "unexpected omics row count")
assert(
  identical(names(omics)[-1L], pheno$SAMPLE_ID),
  "omics sample columns do not match phenotype sample IDs"
)

run_stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_stem <- file.path("runs", paste0("quickstart_", run_stamp))
output_dir <- output_stem
suffix <- 1L
while (file.exists(output_dir)) {
  output_dir <- paste0(output_stem, "_", suffix)
  suffix <- suffix + 1L
}

message(
  "Running the bundled smoke test with 30 XGB trials and one CV repeat. ",
  "Use the documented defaults for scientific analyses."
)

manifest <- FAST_treatment_ML(
  pheno = pheno,
  omics = omics,
  omics_type = "Proteomics",
  additional_covariates = NULL,
  models = c("enet", "xgb"),
  output_dir = output_dir,
  enet_cv_folds = 5L,
  xgb_cv_folds = 5L,
  xgb_cv_repeats = 1L,
  xgb_n_trials = 30L,
  n_cores = 1L,
  python_bin = python_bin,
  seed = 2202L
)

exported <- FAST_export_models(manifest)

expected_files <- c(
  file.path(output_dir, "manifest.json"),
  file.path(output_dir, "models", "reports", "cohort.csv"),
  file.path(output_dir, "models", "FU1", "enet", "metrics.csv"),
  file.path(output_dir, "models", "FU1", "enet", "weights.csv"),
  file.path(output_dir, "models", "FU1", "xgb", "metrics.csv"),
  file.path(output_dir, "models", "FU1", "xgb", "importance.csv"),
  file.path(output_dir, "models", "FU1", "xgb", "tuning.csv"),
  file.path(output_dir, "models", "FU1", "xgb", "model.json"),
  file.path(output_dir, "models", "exported_models", "exported_models.csv")
)
assert(all(file.exists(expected_files)), "one or more expected output files are missing")
assert(identical(names(manifest$followups), "FU1"), "expected exactly one follow-up")
assert(nrow(exported$models) == 2L, "expected two exported model packages")
assert(
  setequal(as.character(exported$models$MODEL), c("enet", "xgb")),
  "exported model families are incomplete"
)
assert(all(file.exists(exported$models$PATH)), "an exported model package is missing")

enet_metrics <- read.csv(
  file.path(output_dir, "models", "FU1", "enet", "metrics.csv"),
  stringsAsFactors = FALSE
)
xgb_metrics <- read.csv(
  file.path(output_dir, "models", "FU1", "xgb", "metrics.csv"),
  stringsAsFactors = FALSE
)
xgb_tuning <- read.csv(
  file.path(output_dir, "models", "FU1", "xgb", "tuning.csv"),
  stringsAsFactors = FALSE
)

assert(
  nrow(enet_metrics) == 1L && "CV_AUC" %in% names(enet_metrics),
  "unexpected ENET metrics schema"
)
assert(
  nrow(xgb_metrics) == 1L &&
    all(c("CV_AUC", "CV_REPEATS") %in% names(xgb_metrics)),
  "unexpected XGB metrics schema"
)
assert(nrow(xgb_tuning) == 30L, "expected one XGB tuning row per smoke-test trial")

reported_auc <- c(enet_metrics$CV_AUC, xgb_metrics$CV_AUC)
assert(
  all(is.finite(reported_auc)) && all(reported_auc >= 0 & reported_auc <= 1),
  "reported CV AUC values must be finite and between zero and one"
)

cat("\nQuick start complete.\n")
cat("Output directory: ", manifest$output_dir, "\n", sep = "")
print(data.frame(
  MODEL = c("enet", "xgb"),
  CV_AUC = reported_auc,
  row.names = NULL
))
