# -----------------------------
# Paths and settings: edit these for your cloud machine
# -----------------------------
pipeline_repo <- path.expand("~/FAST/GitHub/track-2.2")
omics_raw_path <- path.expand("~/FAST/Data/CALERIE/Raw/DNAm/GRSet_fully_filtered_bmiq_chunk.rds")
pheno_raw_path <- path.expand("~/FAST/Data/CALERIE/Raw/DNAm/CALERIE_CPR_processed_pheno.rds")
python_bin <- path.expand("~/FAST/Envs/track22/bin/python")
out_dir <- path.expand("~/FAST/Outputs/2.2/DNAm_betas_2.2_Output")

models <- c("enet", "xgb")
n_cores <- 7L # XGBoost threads; leave headroom on the 8-vCPU cloud machine.
enet_cv_folds <- 10L
xgb_cv_folds <- 10L
xgb_cv_repeats <- 3L
xgb_n_trials <- 100L
seed <- 2202L
# Both models use methylation changes and FEMALE. No extra adjustment covariates.
additional_covariates <- NULL

# Keep BLAS serial, but allow XGBoost to use its explicit thread allocation.
Sys.setenv(
  OMP_NUM_THREADS = as.character(n_cores), OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1"
)
setwd(pipeline_repo)
source(file.path("training", "main.R"))
FAST_check_R("training")
if (!requireNamespace("haven", quietly = TRUE)) stop("Install the R package 'haven'.")
if ("xgb" %in% models) python_bin <- FAST_check_python(python_bin)
models <- .validate_ml_args(
  models, enet_cv_folds, xgb_cv_folds, xgb_cv_repeats, seed, n_cores, xgb_n_trials
)
stopifnot(file.exists(omics_raw_path), file.exists(pheno_raw_path))
if (file.exists(out_dir)) stop("Output directory already exists; choose a new out_dir.")
dir.create(out_dir, recursive = TRUE)

log_file <- file.path(out_dir, "run.log")
log_status <- function(text) {
  line <- paste0("[2.2] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", text)
  message(line)
  cat(line, "\n", file = log_file, append = TRUE, sep = "")
}
started <- Sys.time()
log_status(paste("START CALERIE DNAm beta training; output:", out_dir))
log_status(sprintf(
  "XGB: %d threads, %d trials, %d folds, %d repeats; ENET: %d folds",
  n_cores, xgb_n_trials, xgb_cv_folds, xgb_cv_repeats, enet_cv_folds
))

# -----------------------------
# Read raw inputs and check beta scale
# -----------------------------
log_status("Loading raw beta matrix and phenotype table")
omics_raw <- readRDS(omics_raw_path)
pheno_raw <- readRDS(pheno_raw_path)
if (is.data.frame(omics_raw)) {
  if (!all(vapply(omics_raw, is.numeric, logical(1)))) {
    stop("Raw DNAm data-frame columns must all be numeric beta measurements.")
  }
  omics_raw <- as.matrix(omics_raw)
}
if (!is.matrix(omics_raw) || !is.numeric(omics_raw)) {
  stop("Expected the raw DNAm file to contain a numeric CpG-by-sample matrix or data frame.")
}
valid_names <- function(x) {
  length(x) > 0L && !anyNA(x) && all(nzchar(trimws(x))) && !anyDuplicated(x)
}
if (!valid_names(rownames(omics_raw)) || !valid_names(colnames(omics_raw))) {
  stop("Raw beta matrix needs unique, nonblank probe and sample names.")
}
beta_range <- range(omics_raw)
if (!all(is.finite(beta_range)) || beta_range[1] < 0 || beta_range[2] > 1) {
  stop("Expected complete, finite DNAm beta values in [0, 1].")
}
log_status(sprintf(
  "Loaded %d probes and %d samples; beta range %.6f to %.6f",
  nrow(omics_raw), ncol(omics_raw), beta_range[1], beta_range[2]
))

# -----------------------------
# Build the treatment-prediction cohort (no INF/MetS outcome filtering)
# -----------------------------
required_raw <- c("Barcode", "Participant_ID", "fu", "CR", "female")
if (!is.data.frame(pheno_raw) || !all(required_raw %in% names(pheno_raw))) {
  stop("Phenotype input must contain: ", paste(required_raw, collapse = ", "))
}
# Preserve underlying 0/1 and 0/1/2 codes from haven-labelled columns.
# For factors, parse the displayed codes, never the internal level numbers.
read_code <- function(x, allowed, label) {
  raw <- haven::zap_labels(x)
  code <- suppressWarnings(as.numeric(as.character(raw)))
  if (any(!is.na(raw) & (is.na(code) | !code %in% allowed))) {
    stop(label, " must use codes ", paste(allowed, collapse = "/"), ".")
  }
  code
}
pheno <- data.frame(
  SAMPLE_ID = as.character(pheno_raw$Barcode),
  SUBJECT_ID = as.character(pheno_raw$Participant_ID),
  FU = read_code(pheno_raw$fu, 0:2, "fu"),
  TREATMENT_GROUP = read_code(pheno_raw$CR, 0:1, "CR"),
  FEMALE = read_code(pheno_raw$female, 0:1, "female"),
  stringsAsFactors = FALSE
)
complete <- complete.cases(pheno) & nzchar(trimws(pheno$SAMPLE_ID)) &
  nzchar(trimws(pheno$SUBJECT_ID))
log_status(sprintf("Excluding %d phenotype rows missing required values", sum(!complete)))
pheno <- pheno[complete, , drop = FALSE]
matched <- pheno$SAMPLE_ID %in% colnames(omics_raw)
log_status(sprintf(
  "Sample matching: excluding %d phenotype rows without DNAm and %d DNAm-only columns",
  sum(!matched), sum(!colnames(omics_raw) %in% pheno$SAMPLE_ID)
))
pheno <- pheno[matched, , drop = FALSE]
if (anyDuplicated(pheno$SAMPLE_ID) || anyDuplicated(pheno[c("SUBJECT_ID", "FU")])) {
  stop("Duplicate usable sample or subject/visit rows; resolve replicates before training.")
}
for (column in c("TREATMENT_GROUP", "FEMALE")) {
  consistent <- vapply(split(pheno[[column]], pheno$SUBJECT_ID), function(x) {
    length(unique(x)) == 1L
  }, logical(1))
  if (!all(consistent)) stop(column, " is inconsistent across visits for a subject.")
}
baseline_subjects <- pheno$SUBJECT_ID[pheno$FU == 0L]
paired_subjects <- intersect(baseline_subjects, pheno$SUBJECT_ID[pheno$FU > 0L])
log_status(sprintf(
  "Excluding %d subjects without matched baseline and at least one follow-up",
  length(setdiff(unique(pheno$SUBJECT_ID), paired_subjects))
))
pheno <- pheno[pheno$SUBJECT_ID %in% paired_subjects, , drop = FALSE]
pheno$FU <- factor(pheno$FU, levels = 0:2)
pheno$TREATMENT_GROUP <- factor(pheno$TREATMENT_GROUP, levels = 0:1)
pheno$FEMALE <- factor(pheno$FEMALE, levels = 0:1)

# Preflight both visits before fitting either; each visit keeps its own cohort.
fold_counts <- c()
if ("enet" %in% models) fold_counts["enet_cv_folds"] <- enet_cv_folds
if ("xgb" %in% models) fold_counts["xgb_cv_folds"] <- xgb_cv_folds
for (fu in 1:2) {
  visit_pheno <- pheno[pheno$FU == fu, , drop = FALSE]
  issue <- .validate_followup_cohort(visit_pheno, fold_counts)
  if (!is.null(issue)) stop("FU", fu, ": ", issue)
  counts <- table(visit_pheno$TREATMENT_GROUP)
  log_status(sprintf(
    "FU%d (%d months): %d baseline-paired subjects; control=%d, treatment=%d",
    fu, fu * 12L, nrow(visit_pheno), counts[1], counts[2]
  ))
}

# -----------------------------
# Build aligned omics input; the pipeline selects reliable probes
# -----------------------------
reliable_probes <- readRDS("Data/FAST_epicv1_epicv2_sugden_TruD_probe_list.rds")
n_reliable <- sum(reliable_probes %in% rownames(omics_raw))
if (n_reliable == 0L) stop("No reliable DNAm probes found in the input.")
log_status(sprintf("Reliable probe coverage: %d of %d", n_reliable, length(reliable_probes)))
omics <- data.frame(
  ANALYTE_NAME = rownames(omics_raw),
  omics_raw[, pheno$SAMPLE_ID, drop = FALSE], check.names = FALSE
)
stopifnot(identical(names(omics)[-1L], pheno$SAMPLE_ID))
rm(omics_raw, pheno_raw)
invisible(gc())

provenance <- list(
  started = started, measurement_scale = "beta", feature_mode = "follow-up minus baseline",
  inputs = c(omics = omics_raw_path, pheno = pheno_raw_path),
  pipeline_repo = pipeline_repo, python_bin = python_bin, models = models,
  n_cores = n_cores, enet_cv_folds = enet_cv_folds, xgb_cv_folds = xgb_cv_folds,
  xgb_cv_repeats = xgb_cv_repeats, xgb_n_trials = xgb_n_trials, seed = seed,
  additional_covariates = additional_covariates, n_reliable_probes = n_reliable,
  visits = c(baseline = 0L, month12 = 1L, month24 = 2L), session_info = sessionInfo()
)
saveRDS(provenance, file.path(out_dir, "provenance.rds"))

# -----------------------------
# Train and export models
# -----------------------------
log_status(paste("Training", paste(models, collapse = "/"), "on beta changes at FU1 and FU2"))
manifest <- FAST_treatment_ML(
  pheno = pheno, omics = omics, omics_type = "DNAm",
  additional_covariates = additional_covariates, models = models, output_dir = out_dir,
  enet_cv_folds = enet_cv_folds, xgb_cv_folds = xgb_cv_folds,
  xgb_cv_repeats = xgb_cv_repeats, xgb_n_trials = xgb_n_trials,
  n_cores = n_cores, python_bin = python_bin, seed = seed
)
log_status("Training complete; exporting fitted models")
exported <- FAST_export_models(manifest)
stopifnot(setequal(names(manifest$followups), c("FU1", "FU2")))
stopifnot(nrow(exported$models) == 2L * length(models), all(file.exists(exported$models$PATH)))
provenance$finished <- Sys.time()
provenance$elapsed_seconds <- as.numeric(difftime(provenance$finished, started, units = "secs"))
saveRDS(provenance, file.path(out_dir, "provenance.rds"))
log_status(sprintf("COMPLETE: %d exported models in %.1f minutes", nrow(exported$models),
                   provenance$elapsed_seconds / 60))
print(exported$models[c("FU", "MODEL", "TRAINING_CV_AUC", "SUCCESSFUL")])
