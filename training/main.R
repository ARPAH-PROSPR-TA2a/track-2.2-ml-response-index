source(file.path("R", "check_environment.R"))
source(file.path("R", "validate_inputs.R"))
source(file.path("R", "build_training_features.R"))
source(file.path("R", "write_training_artifacts.R"))


.python_setup_help <- function() {
  paste(
    "What to do next:",
    "  1. In a terminal, run `which python3` and `which python`.",
    "     On Windows, run `where python`.",
    "  2. Copy the path for the Python you intend to use, then run:",
    "     FAST_check_python(\"/full/path/to/python\")",
    "  3. Open the Troubleshooting section in README.md for setup help.",
    sep = "\n"
  )
}


FAST_check_python <- function(python_bin = NULL) {
  if (is.null(python_bin)) {
    python_bin <- "python3"
  }
  if (!is.character(python_bin) || length(python_bin) != 1L ||
      is.na(python_bin) || !nzchar(python_bin)) {
    stop("python_bin must be one non-empty character string.", call. = FALSE)
  }

  supplied_path <- grepl("[/\\\\]", python_bin)
  # Expand `~`, but do not resolve symlinks: a venv's Python may be a symlink.
  python_command <- if (supplied_path) path.expand(python_bin) else python_bin
  if (supplied_path) {
    if (!file.exists(python_command)) {
      stop(
        "Python executable was not found: '", python_bin, "'.\n\n",
        .python_setup_help(),
        call. = FALSE
      )
    }
    if (dir.exists(python_command)) {
      stop(
        "python_bin points to a directory, not a Python executable: '",
        python_bin, "'.\n\n", .python_setup_help(),
        call. = FALSE
      )
    }
    if (file.access(python_command, mode = 1L) != 0L) {
      stop(
        "Python executable is not executable: '", python_bin, "'.\n\n",
        .python_setup_help(),
        call. = FALSE
      )
    }
  } else if (!nzchar(Sys.which(python_bin))) {
    stop(
      "Python executable '", python_bin, "' was not found on PATH.\n\n",
      .python_setup_help(),
      call. = FALSE
    )
  }

  checker <- file.path("training", "scripts", "check_python.py")
  if (!file.exists(checker)) {
    stop(
      "Could not find ", checker,
      ". Run FAST from the Track 2.2 repository root.",
      call. = FALSE
    )
  }

  output <- tryCatch(
    suppressWarnings(system2(
      python_command,
      args = shQuote(checker),
      stdout = TRUE,
      stderr = TRUE,
      timeout = 60L
    )),
    error = identity
  )
  if (inherits(output, "error")) {
    stop(
      "Python executable could not be run: '", python_bin, "'. ",
      conditionMessage(output), "\n\n", .python_setup_help(),
      call. = FALSE
    )
  }

  status <- attr(output, "status")
  if (is.null(status)) {
    status <- 0L
  }
  output <- output[nzchar(output)]
  if (status %in% c(10L, 11L, 12L)) {
    stop(
      paste(output, collapse = "\n"),
      "\n\n", .python_setup_help(),
      call. = FALSE
    )
  }
  if (!identical(as.integer(status), 0L)) {
    stop(
      "Python executable failed the FAST environment check with exit status ",
      status, ": '", python_bin, "'.",
      if (length(output) > 0L) paste0("\n", paste(output, collapse = "\n")) else "",
      "\n\n", .python_setup_help(),
      call. = FALSE
    )
  }

  if (length(output) == 0L) {
    stop(
      "Python executable returned no FAST environment-check output: '",
      python_bin, "'.\n\n", .python_setup_help(),
      call. = FALSE
    )
  }
  message(paste(output, collapse = "\n"))

  invisible(python_command)
}


.validate_ml_args <- function(models, enet_cv_folds, xgb_cv_folds,
                              xgb_cv_repeats, seed, n_cores, xgb_n_trials) {
  allowed_models <- c("enet", "xgb")
  if (!is.character(models) || length(models) == 0L ||
      any(!models %in% allowed_models)) {
    stop("models must contain one or more of: ", paste(allowed_models, collapse = ", "))
  }

  if ("enet" %in% models &&
      (!is.numeric(enet_cv_folds) || length(enet_cv_folds) != 1L ||
       enet_cv_folds < 4L)) {
    stop("enet_cv_folds must be a single integer >= 4.")
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
  if ("xgb" %in% models) {
    python_bin <- FAST_check_python(python_bin)
  }
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
    analyte_name_map = inputs$analyte_name_map,
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

  message("Training complete. Results: ", manifest$output_dir)
  manifest
}
