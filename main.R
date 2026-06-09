source("validation_helpers.R")
source("reporting_helpers.R")
source("feature_helpers.R")
source("ml_helpers.R")


.validate_ml_args <- function(models, test_frac, cv_folds, seed, n_cores, xgb_n_trials) {
  allowed_models <- c("enet", "xgb")
  if (!is.character(models) || length(models) == 0L ||
      any(!models %in% allowed_models)) {
    stop("models must contain one or more of: ", paste(allowed_models, collapse = ", "))
  }

  if (!is.numeric(test_frac) || length(test_frac) != 1L ||
      test_frac <= 0 || test_frac >= 0.5) {
    stop("test_frac must be a single numeric value greater than 0 and less than 0.5.")
  }

  if (!is.numeric(cv_folds) || length(cv_folds) != 1L || cv_folds < 2L) {
    stop("cv_folds must be a single integer >= 2.")
  }

  if (!is.numeric(seed) || length(seed) != 1L) {
    stop("seed must be a single numeric value.")
  }

  if (!is.null(n_cores) && (!is.numeric(n_cores) || length(n_cores) != 1L || n_cores < 1L)) {
    stop("n_cores must be NULL or a single integer >= 1.")
  }

  if (!is.numeric(xgb_n_trials) || length(xgb_n_trials) != 1L || xgb_n_trials < 0L) {
    stop("xgb_n_trials must be a single integer >= 0.")
  }

  unique(models)
}


.prepare_inputs <- function(pheno, omics, omics_type, additional_covariates = NULL) {
  .validate_omics_type(omics_type)

  pheno_list <- .validate_pheno(pheno, additional_covariates)
  omics_list <- .validate_omics(omics, pheno_list)

  if (omics_type == "DNAm") {
    full_probes <- readRDS("Data/FAST_epicv1_epicv2_probe_list.rds")
    reliable_probes <- readRDS("Data/FAST_epicv1_epicv2_sugden_TruD_probe_list.rds")
    .validate_dnam_probe_coverage(full_probes, reliable_probes, omics_list$all$ANALYTE_NAME)
    omics_list <- .subset_omics_list(omics_list, reliable_probes)
    message("DNAm: restricted analysis to ", nrow(omics_list$all), " reliable probes.")
  }

  list(pheno_list = pheno_list, omics_list = omics_list)
}


FAST_treatment_ML <- function(pheno,
                              omics,
                              omics_type = "Proteomics",
                              additional_covariates = NULL,
                              models = c("enet", "xgb"),
                              output_dir = NULL,
                              test_frac = 0.2,
                              cv_folds = 5L,
                              seed = 1L,
                              n_cores = NULL,
                              python_bin = NULL,
                              xgb_n_trials = 0L) {

  models <- .validate_ml_args(models, test_frac, cv_folds, seed, n_cores, xgb_n_trials)
  if (is.null(n_cores)) {
    n_cores <- max(1L, parallel::detectCores() - 1L)
  }
  if (is.null(output_dir)) {
    run_stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    output_dir <- file.path("runs", paste0("treatment_ml_", run_stamp))
  }
  if (is.null(python_bin)) {
    python_bin <- "python3"
  }

  inputs <- .prepare_inputs(pheno, omics, omics_type, additional_covariates)

  manifest <- .run_ml_disk(
    pheno_df = inputs$pheno_list$all,
    omics_df = inputs$omics_list$all,
    additional_covariates = additional_covariates,
    models = models,
    output_dir = output_dir,
    test_frac = test_frac,
    cv_folds = as.integer(cv_folds),
    seed = as.integer(seed),
    n_cores = as.integer(n_cores),
    python_bin = python_bin,
    xgb_n_trials = as.integer(xgb_n_trials)
  )

  manifest
}


FAST_treatment_ML_reports <- function(pheno,
                                      omics,
                                      omics_type = "Proteomics",
                                      additional_covariates = NULL) {
  inputs <- .prepare_inputs(pheno, omics, omics_type, additional_covariates)
  .generate_reports(inputs$pheno_list, inputs$omics_list, additional_covariates)
}
