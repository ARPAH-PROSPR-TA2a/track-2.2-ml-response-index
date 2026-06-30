source(file.path("R", "validate_inputs.R"))
source(file.path("R", "build_training_features.R"))
source(file.path("R", "replay_preprocessing.R"))
source(file.path("R", "score_models.R"))


.write_inference_outputs <- function(predictions, validation, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  validation_path <- file.path(output_dir, "validation.csv")
  utils::write.csv(validation, validation_path, row.names = FALSE)

  prediction_paths <- list()
  for (fu_key in names(predictions)) {
    prediction_paths[[fu_key]] <- list()
    fu_dir <- file.path(output_dir, "predictions", fu_key)
    dir.create(fu_dir, recursive = TRUE, showWarnings = FALSE)

    for (model_name in names(predictions[[fu_key]])) {
      path <- file.path(fu_dir, paste0(model_name, ".csv"))
      utils::write.csv(predictions[[fu_key]][[model_name]], path, row.names = FALSE)
      prediction_paths[[fu_key]][[model_name]] <- path
    }
  }

  list(
    validation = validation_path,
    predictions = prediction_paths
  )
}


.write_single_model_evaluation_outputs <- function(predictions, validation, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  prediction_path <- file.path(output_dir, "predictions.csv")
  validation_path <- file.path(output_dir, "validation.csv")
  utils::write.csv(predictions, prediction_path, row.names = FALSE)
  utils::write.csv(validation, validation_path, row.names = FALSE)

  list(
    predictions = prediction_path,
    validation = validation_path
  )
}


.load_exported_model_package <- function(model_path) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required to read model JSON files.")
  }
  package <- jsonlite::fromJSON(model_path, simplifyVector = TRUE)
  if (!identical(package$schema_version, "1.0")) {
    stop("Unsupported model package schema version: ", package$schema_version)
  }
  if (!package$family %in% c("enet", "xgb")) {
    stop("Unsupported model family: ", package$family)
  }
  package
}


.package_vector <- function(x) {
  if (is.null(x) || length(x) == 0L) {
    return(character(0))
  }
  as.character(unlist(x, use.names = FALSE))
}


.score_exported_model_package <- function(model_package, x) {
  if (model_package$family == "enet") {
    .score_enet_weight_table(x, as.data.frame(model_package$weights, stringsAsFactors = FALSE))
  } else {
    .score_xgb_model_json(x, model_package$xgb_model_json)
  }
}


.relative_path <- function(path, base_dir) {
  path <- normalizePath(path, mustWork = FALSE)
  base_dir <- normalizePath(base_dir, mustWork = FALSE)
  prefix <- paste0(base_dir, .Platform$file.sep)
  if (startsWith(path, prefix)) {
    return(sub(paste0("^", gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", prefix)), "", path))
  }
  path
}


FAST_treatment_predict <- function(pheno,
                                   omics,
                                   manifest_path,
                                   models = NULL,
                                   followups = NULL,
                                   omics_type = NULL,
                                   enet_mode = "reproduce_training",
                                   output_dir = NULL,
                                   return_matrices = FALSE) {
  enet_mode <- match.arg(enet_mode, "reproduce_training")
  manifest <- .load_manifest(manifest_path)
  if (is.null(omics_type)) {
    omics_type <- manifest$omics_type
  }
  if (is.null(omics_type)) {
    stop("omics_type must be supplied for manifests that do not record it.")
  }

  requested_models <- if (is.null(models)) .manifest_vector(manifest$requested_models) else models
  requested_models <- match.arg(requested_models, c("enet", "xgb"), several.ok = TRUE)
  additional_covariates <- .manifest_vector(manifest$additional_covariates)
  model_covariates <- .manifest_vector(manifest$model_covariates)
  if (length(model_covariates) == 0L) {
    model_covariates <- unique(c("FEMALE", additional_covariates))
  }

  inputs <- .prepare_inputs(
    pheno = pheno,
    omics = omics,
    omics_type = omics_type,
    additional_covariates = additional_covariates
  )
  preprocessing <- read.csv(manifest$reports$preprocessing, stringsAsFactors = FALSE)

  followup_names <- names(manifest$followups)
  if (!is.null(followups)) {
    followup_names <- ifelse(
      grepl("^FU", as.character(followups)),
      as.character(followups),
      paste0("FU", as.integer(followups))
    )
  }

  predictions <- list()
  matrices <- list()
  validation_rows <- list()
  validation_index <- 1L

  for (fu_key in followup_names) {
    fu_manifest <- manifest$followups[[fu_key]]
    if (is.null(fu_manifest)) next

    fu_level <- as.integer(sub("^FU", "", fu_key))
    predictions[[fu_key]] <- list()
    matrices[[fu_key]] <- list()

    for (model_name in requested_models) {
      model_manifest <- fu_manifest$models[[model_name]]
      if (is.null(model_manifest)) next
      training_cv_auc <- .model_cv_auc(model_manifest)

      replay <- .build_inference_matrix(
        pheno_df = inputs$pheno,
        omics_df = inputs$omics,
        preprocessing = preprocessing,
        fu_level = fu_level,
        model = model_name,
        model_covariates = model_covariates
      )
      matrices[[fu_key]][[model_name]] <- replay$matrix

      pred <- if (model_name == "enet") {
        .score_enet_weights(replay$matrix, model_manifest$weights)
      } else {
        .score_xgb_model(replay$matrix, model_manifest$model)
      }

      predictions[[fu_key]][[model_name]] <- data.frame(
        SUBJECT_ID = replay$subject_ids,
        FU = fu_level,
        TREATMENT_GROUP = replay$y,
        PREDICTED_PROB = pred,
        stringsAsFactors = FALSE
      )
      validation_rows[[validation_index]] <- .validation_row(
        fu_level = fu_level,
        model_name = model_name,
        predictions = predictions[[fu_key]][[model_name]],
        training_cv_auc = training_cv_auc
      )
      validation_index <- validation_index + 1L
    }
  }

  validation <- if (length(validation_rows) > 0L) {
    do.call(rbind, validation_rows)
  } else {
    data.frame()
  }

  output_files <- NULL
  normalized_output_dir <- NULL
  if (!is.null(output_dir)) {
    output_files <- .write_inference_outputs(predictions, validation, output_dir)
    normalized_output_dir <- normalizePath(output_dir, mustWork = FALSE)
  }

  result <- list(
    predictions = predictions,
    validation = validation,
    manifest = manifest,
    output_dir = normalized_output_dir,
    output_files = output_files
  )
  if (return_matrices) {
    result$matrices <- matrices
  }
  result
}


FAST_evaluate <- function(pheno,
                          omics,
                          model_path,
                          output_dir = NULL,
                          return_matrix = FALSE) {
  model_package <- .load_exported_model_package(model_path)
  additional_covariates <- .package_vector(model_package$covariates$additional_covariates)
  model_covariates <- .package_vector(model_package$covariates$model_covariates)
  if (length(model_covariates) == 0L) {
    model_covariates <- unique(c("FEMALE", additional_covariates))
  }

  inputs <- .prepare_inputs(
    pheno = pheno,
    omics = omics,
    omics_type = model_package$omics_type,
    additional_covariates = additional_covariates
  )

  preprocessing <- as.data.frame(model_package$preprocessing, stringsAsFactors = FALSE)
  replay <- .build_inference_matrix(
    pheno_df = inputs$pheno,
    omics_df = inputs$omics,
    preprocessing = preprocessing,
    fu_level = as.integer(model_package$fu),
    model = model_package$family,
    model_covariates = model_covariates
  )
  pred <- .score_exported_model_package(model_package, replay$matrix)

  predictions <- data.frame(
    SUBJECT_ID = replay$subject_ids,
    FU = as.integer(model_package$fu),
    TREATMENT_GROUP = replay$y,
    PREDICTED_PROB = pred,
    stringsAsFactors = FALSE
  )
  validation <- .validation_row(
    fu_level = as.integer(model_package$fu),
    model_name = model_package$family,
    predictions = predictions,
    training_cv_auc = as.numeric(model_package$training_cv_auc)
  )

  output_files <- NULL
  normalized_output_dir <- NULL
  if (!is.null(output_dir)) {
    output_files <- .write_single_model_evaluation_outputs(
      predictions,
      validation,
      output_dir
    )
    normalized_output_dir <- normalizePath(output_dir, mustWork = FALSE)
  }

  result <- list(
    model_id = model_package$model_id,
    predictions = predictions,
    validation = validation,
    model = model_package,
    output_dir = normalized_output_dir,
    output_files = output_files
  )
  if (return_matrix) {
    result$matrix <- replay$matrix
  }
  result
}


FAST_bulk_evaluate <- function(pheno,
                               omics,
                               models_dir,
                               output_dir) {
  index_path <- file.path(models_dir, "exported_models.csv")
  if (!file.exists(index_path)) {
    stop("models_dir must contain exported_models.csv.")
  }

  index <- read.csv(index_path, stringsAsFactors = FALSE)
  required_cols <- c("MODEL_ID", "PATH", "SUCCESSFUL")
  missing_cols <- setdiff(required_cols, names(index))
  if (length(missing_cols) > 0L) {
    stop("exported_models.csv is missing column(s): ", paste(missing_cols, collapse = ", "))
  }

  successful <- index[index$SUCCESSFUL, , drop = FALSE]
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  rows <- list()
  for (i in seq_len(nrow(successful))) {
    model_row <- successful[i, , drop = FALSE]
    model_path <- model_row$PATH
    if (!file.exists(model_path)) {
      model_path <- file.path(models_dir, basename(model_path))
    }
    if (!file.exists(model_path)) {
      stop("Model package does not exist for ", model_row$MODEL_ID, ": ", model_row$PATH)
    }

    model_output_dir <- file.path(output_dir, "models", model_row$MODEL_ID)
    evaluated <- tryCatch(
      FAST_evaluate(
        pheno = pheno,
        omics = omics,
        model_path = model_path,
        output_dir = model_output_dir
      ),
      error = function(e) {
        stop(
          "Failed evaluating model ", model_row$MODEL_ID,
          " (", model_path, "): ", conditionMessage(e),
          call. = FALSE
        )
      }
    )

    validation <- evaluated$validation
    validation$MODEL_ID <- model_row$MODEL_ID
    validation$MODEL_PATH <- .relative_path(model_path, output_dir)
    validation$PREDICTIONS_PATH <- .relative_path(evaluated$output_files$predictions, output_dir)
    validation$VALIDATION_PATH <- .relative_path(evaluated$output_files$validation, output_dir)
    rows[[i]] <- validation[
      ,
      c(
        "MODEL_ID", "MODEL_PATH", "FU", "MODEL", "TRAINING_CV_AUC",
        "SUCCESS_AUC_THRESHOLD", "SUCCESSFUL", "N", "N_CONTROL",
        "N_TREATMENT", "AUC", "LOGIT_BETA", "LOGIT_OR", "LOGIT_P",
        "VALIDATED_P05", "PREDICTIONS_PATH", "VALIDATION_PATH"
      ),
      drop = FALSE
    ]
  }

  summary <- if (length(rows) > 0L) {
    do.call(rbind, rows)
  } else {
    data.frame(
      MODEL_ID = character(0),
      MODEL_PATH = character(0),
      FU = integer(0),
      MODEL = character(0),
      TRAINING_CV_AUC = numeric(0),
      SUCCESS_AUC_THRESHOLD = numeric(0),
      SUCCESSFUL = logical(0),
      N = integer(0),
      N_CONTROL = integer(0),
      N_TREATMENT = integer(0),
      AUC = numeric(0),
      LOGIT_BETA = numeric(0),
      LOGIT_OR = numeric(0),
      LOGIT_P = numeric(0),
      VALIDATED_P05 = logical(0),
      PREDICTIONS_PATH = character(0),
      VALIDATION_PATH = character(0),
      stringsAsFactors = FALSE
    )
  }

  summary_path <- file.path(output_dir, "validation_summary.csv")
  utils::write.csv(summary, summary_path, row.names = FALSE)

  list(
    validation = summary,
    output_dir = normalizePath(output_dir, mustWork = FALSE),
    output_files = list(validation_summary = summary_path)
  )
}
