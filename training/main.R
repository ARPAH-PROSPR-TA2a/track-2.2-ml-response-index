source(file.path("R", "validate_inputs.R"))
source(file.path("R", "build_training_features.R"))
source(file.path("R", "write_training_artifacts.R"))


.validate_ml_args <- function(models, enet_cv_folds, xgb_cv_folds,
                              xgb_cv_repeats, seed, n_cores, xgb_n_trials) {
  allowed_models <- c("enet", "xgb")
  if (!is.character(models) || length(models) == 0L ||
      any(!models %in% allowed_models)) {
    stop("models must contain one or more of: ", paste(allowed_models, collapse = ", "))
  }

  if (!is.numeric(enet_cv_folds) || length(enet_cv_folds) != 1L || enet_cv_folds < 2L) {
    stop("enet_cv_folds must be a single integer >= 2.")
  }

  if (!is.numeric(xgb_cv_folds) || length(xgb_cv_folds) != 1L || xgb_cv_folds < 2L) {
    stop("xgb_cv_folds must be a single integer >= 2.")
  }

  if (!is.numeric(xgb_cv_repeats) || length(xgb_cv_repeats) != 1L || xgb_cv_repeats < 1L) {
    stop("xgb_cv_repeats must be a single integer >= 1.")
  }

  if (!is.numeric(seed) || length(seed) != 1L) {
    stop("seed must be a single numeric value.")
  }

  if (!is.null(n_cores) && (!is.numeric(n_cores) || length(n_cores) != 1L || n_cores < 1L)) {
    stop("n_cores must be NULL or a single integer >= 1.")
  }

  if ("xgb" %in% models &&
      (!is.numeric(xgb_n_trials) || length(xgb_n_trials) != 1L ||
       !is.finite(xgb_n_trials) || xgb_n_trials != as.integer(xgb_n_trials) ||
       xgb_n_trials < 10L)) {
    stop("xgb_n_trials must be a single integer >= 10 when XGB is requested.")
  }
  if ("xgb" %in% models && xgb_n_trials < 30L) {
    warning(
      "xgb_n_trials below 30 provides a limited hyperparameter search; ",
      "30 or more trials are recommended."
    )
  }

  unique(models)
}


FAST_treatment_ML <- function(pheno,
                              omics,
                              omics_type = "Proteomics",
                              additional_covariates = NULL,
                              models = c("enet", "xgb"),
                              output_dir = NULL,
                              enet_cv_folds = 10L,
                              xgb_cv_folds = 10L,
                              xgb_cv_repeats = 3L,
                              xgb_n_trials = 50L,
                              n_cores = NULL,
                              python_bin = NULL,
                              seed = 1L) {

  models <- .validate_ml_args(
    models, enet_cv_folds, xgb_cv_folds, xgb_cv_repeats,
    seed, n_cores, xgb_n_trials
  )
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

  model_covariates <- unique(c("FEMALE", additional_covariates))
  inputs <- .prepare_inputs(pheno, omics, omics_type, additional_covariates)

  manifest <- .run_ml_disk(
    pheno_df = inputs$pheno,
    omics_df = inputs$omics,
    additional_covariates = additional_covariates,
    model_covariates = model_covariates,
    models = models,
    output_dir = output_dir,
    omics_type = omics_type,
    enet_cv_folds = as.integer(enet_cv_folds),
    xgb_cv_folds = as.integer(xgb_cv_folds),
    xgb_cv_repeats = as.integer(xgb_cv_repeats),
    seed = as.integer(seed),
    n_cores = as.integer(n_cores),
    python_bin = python_bin,
    xgb_n_trials = as.integer(xgb_n_trials)
  )

  manifest
}
