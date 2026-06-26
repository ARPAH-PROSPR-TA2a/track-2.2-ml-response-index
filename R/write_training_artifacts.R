.write_matrix_csv_gz <- function(x, path) {
  con <- gzfile(path, open = "wt")
  on.exit(close(con), add = TRUE)
  write.csv(as.data.frame(x, check.names = FALSE), con, row.names = FALSE)
}


.write_csv <- function(x, path) {
  utils::write.csv(x, path, row.names = FALSE)
}


.write_json_simple <- function(x, path) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required to write ML run config files.")
  }
  jsonlite::write_json(x, path, auto_unbox = TRUE, pretty = TRUE, null = "null")
}


.write_prepared_dataset <- function(dataset, fu_dir) {
  dir.create(fu_dir, recursive = TRUE, showWarnings = FALSE)

  .write_matrix_csv_gz(dataset$enet_x_train, file.path(fu_dir, "enet_train.csv.gz"))
  .write_matrix_csv_gz(dataset$xgb_x_train, file.path(fu_dir, "xgb_train.csv.gz"))

  .write_csv(
    data.frame(
      SUBJECT_ID = dataset$subject_ids_train,
      FU = dataset$fu_level,
      SET = "train",
      TREATMENT_GROUP = dataset$y_train,
      ENET_FOLD_ID = dataset$enet_foldid,
      stringsAsFactors = FALSE
    ),
    file.path(fu_dir, "subjects.csv")
  )
  .write_csv(dataset$xgb_folds, file.path(fu_dir, "xgb_folds.csv"))
}


.run_enet_worker <- function(dataset, out_dir, seed = 1L, alpha = 0.5) {
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("Package 'glmnet' is required for ENET models.")
  }

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  penalty_factor <- ifelse(
    startsWith(colnames(dataset$enet_x_train), "covariate::"),
    0,
    1
  )

  set.seed(seed)
  cv_fit <- glmnet::cv.glmnet(
    x = dataset$enet_x_train,
    y = dataset$y_train,
    family = "binomial",
    alpha = alpha,
    foldid = dataset$enet_foldid,
    type.measure = "deviance",
    standardize = FALSE,
    penalty.factor = penalty_factor,
    keep = TRUE
  )

  lambda <- cv_fit$lambda.min
  lambda_idx <- which.min(abs(cv_fit$lambda - lambda))
  cv_auc <- .safe_auc(dataset$y_train, cv_fit$fit.preval[, lambda_idx])

  fit <- cv_fit$glmnet.fit
  train_pred <- as.numeric(predict(fit, newx = dataset$enet_x_train, s = lambda, type = "response"))

  coefs <- as.matrix(stats::coef(fit, s = lambda))
  weights <- data.frame(
    FEATURE_NAME = rownames(coefs),
    WEIGHT = as.numeric(coefs[, 1]),
    stringsAsFactors = FALSE
  )
  weights$FEATURE_TYPE <- ifelse(
    weights$FEATURE_NAME == "(Intercept)",
    "intercept",
    ifelse(startsWith(weights$FEATURE_NAME, "omics::"), "omics", "covariate")
  )
  weights <- weights[
    weights$WEIGHT != 0 |
      weights$FEATURE_NAME == "(Intercept)" |
      weights$FEATURE_TYPE == "covariate",
    ,
    drop = FALSE
  ]
  row.names(weights) <- NULL

  metrics <- data.frame(
    CV_AUC = cv_auc,
    INSAMPLE_AUC = .safe_auc(dataset$y_train, train_pred),
    LAMBDA = lambda,
    LAMBDA_1SE = cv_fit$lambda.1se,
    ALPHA = alpha,
    N_FEATURES = ncol(dataset$enet_x_train),
    N_NONZERO = sum(
      weights$FEATURE_NAME != "(Intercept)" & weights$WEIGHT != 0
    ),
    N_UNPENALIZED = sum(penalty_factor == 0),
    stringsAsFactors = FALSE
  )

  predictions <- data.frame(
    SET = "train",
    SUBJECT_ID = dataset$subject_ids_train,
    FU = dataset$fu_level,
    TREATMENT_GROUP = dataset$y_train,
    PREDICTED_PROB = train_pred,
    stringsAsFactors = FALSE
  )

  .write_csv(metrics, file.path(out_dir, "metrics.csv"))
  .write_csv(predictions, file.path(out_dir, "predictions.csv"))
  .write_csv(weights, file.path(out_dir, "weights.csv"))

  list(
    metrics = file.path(out_dir, "metrics.csv"),
    predictions = file.path(out_dir, "predictions.csv"),
    weights = file.path(out_dir, "weights.csv")
  )
}


.run_xgb_worker <- function(fu_dir, out_dir, python_bin = "python3",
                            n_trials = 50L, seed = 1L, n_cores = 1L) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  args <- c(
    file.path("training", "scripts", "run_xgb.py"),
    "--data-dir", fu_dir,
    "--out-dir", out_dir,
    "--seed", as.character(seed),
    "--nthread", as.character(max(1L, as.integer(n_cores))),
    "--n-trials", as.character(as.integer(n_trials))
  )

  status <- system2(python_bin, args = args)
  if (!identical(status, 0L)) {
    stop("XGB worker failed with exit status ", status, ".")
  }

  list(
    metrics = file.path(out_dir, "metrics.csv"),
    predictions = file.path(out_dir, "predictions.csv"),
    importance = file.path(out_dir, "importance.csv"),
    tuning = file.path(out_dir, "tuning.csv"),
    model = file.path(out_dir, "model.json")
  )
}


.run_ml_disk <- function(pheno_df, omics_df, additional_covariates = NULL,
                         model_covariates = "FEMALE",
                         models = c("enet", "xgb"), output_dir,
                         omics_type = NULL,
                         enet_cv_folds = 10L,
                         xgb_cv_folds = 10L, xgb_cv_repeats = 3L, seed = 1L,
                         n_cores = 1L, python_bin = "python3",
                         xgb_n_trials = 50L) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  manifest <- list(
    output_dir = normalizePath(output_dir, mustWork = FALSE),
    target = "TREATMENT_GROUP",
    omics_type = omics_type,
    feature_mode = "change",
    requested_models = models,
    enet_cv_folds = enet_cv_folds,
    xgb_cv_folds = xgb_cv_folds,
    xgb_cv_repeats = xgb_cv_repeats,
    xgb_n_trials = xgb_n_trials,
    seed = seed,
    additional_covariates = additional_covariates,
    model_covariates = model_covariates,
    reports = list(
      cohort = file.path(output_dir, "reports", "cohort.csv"),
      change_summary = file.path(output_dir, "reports", "change_summary.csv"),
      preprocessing = file.path(output_dir, "reports", "preprocessing.csv")
    ),
    followups = list()
  )
  report_tables <- list(
    cohort = list(),
    change_summary = list(),
    preprocessing = list()
  )

  fu_levels <- sort(unique(as.integer(as.character(pheno_df$FU))))
  fu_levels <- fu_levels[fu_levels != 0L]

  for (fu_level in fu_levels) {
    fu_key <- paste0("FU", fu_level)
    fu_dir <- file.path(output_dir, fu_key)
    fu_num <- as.integer(as.character(pheno_df$FU))
    baseline_subjects <- pheno_df$SUBJECT_ID[fu_num == 0L]
    followup_subjects <- pheno_df$SUBJECT_ID[fu_num == fu_level]
    complete_subjects <- intersect(baseline_subjects, followup_subjects)
    split_pheno <- pheno_df[
      fu_num == fu_level & pheno_df$SUBJECT_ID %in% complete_subjects,
      ,
      drop = FALSE
    ]

    requested_fold_counts <- c()
    if ("enet" %in% models) requested_fold_counts["enet_cv_folds"] <- enet_cv_folds
    if ("xgb" %in% models) requested_fold_counts["xgb_cv_folds"] <- xgb_cv_folds
    cohort_issue <- .validate_followup_cohort(split_pheno, requested_fold_counts)
    if (!is.null(cohort_issue)) {
      stop(fu_key, ": ", cohort_issue, ".")
    }

    message(fu_key, ": preparing model-ready datasets.")
    prepared <- .prepare_fu_change_dataset(
      pheno_df = pheno_df,
      omics_df = omics_df,
      fu_level = fu_level,
      model_covariates = model_covariates,
      enet_cv_folds = enet_cv_folds,
      xgb_cv_folds = xgb_cv_folds,
      xgb_cv_repeats = xgb_cv_repeats,
      seed = seed + fu_level
    )
    if (is.null(prepared)) {
      manifest$followups[[fu_key]] <- NULL
      next
    }

    .write_prepared_dataset(prepared, fu_dir)
    report_tables$cohort[[fu_key]] <- prepared$cohort
    report_tables$change_summary[[fu_key]] <- prepared$change_summary
    report_tables$preprocessing[[fu_key]] <- prepared$preprocessing
    message(
      fu_key, ": prepared ", length(prepared$subject_ids_train), " training subjects."
    )

    fu_manifest <- list(
      data_dir = normalizePath(fu_dir, mustWork = FALSE),
      artifacts = list(
        enet_train = file.path(fu_dir, "enet_train.csv.gz"),
        xgb_train = file.path(fu_dir, "xgb_train.csv.gz"),
        subjects = file.path(fu_dir, "subjects.csv"),
        xgb_folds = file.path(fu_dir, "xgb_folds.csv")
      ),
      models = list()
    )

    if ("enet" %in% models) {
      message(fu_key, ": fitting ENET.")
      fu_manifest$models$enet <- .run_enet_worker(
        prepared,
        out_dir = file.path(fu_dir, "enet"),
        seed = seed + fu_level
      )
      message(fu_key, ": ENET complete.")
    }

    if ("xgb" %in% models) {
      message(
        fu_key, ": fitting XGB with ", xgb_n_trials, " Optuna trials."
      )
      fu_manifest$models$xgb <- .run_xgb_worker(
        fu_dir = fu_dir,
        out_dir = file.path(fu_dir, "xgb"),
        python_bin = python_bin,
        n_trials = xgb_n_trials,
        seed = seed + fu_level,
        n_cores = n_cores
      )
      message(fu_key, ": XGB complete.")
    }

    manifest$followups[[fu_key]] <- fu_manifest
  }

  reports_dir <- file.path(output_dir, "reports")
  dir.create(reports_dir, recursive = TRUE, showWarnings = FALSE)
  for (report_name in names(report_tables)) {
    if (length(report_tables[[report_name]]) > 0L) {
      .write_csv(
        do.call(rbind, report_tables[[report_name]]),
        manifest$reports[[report_name]]
      )
    }
  }

  manifest_path <- file.path(output_dir, "manifest.json")
  .write_json_simple(manifest, manifest_path)
  manifest$manifest_path <- manifest_path
  manifest
}
