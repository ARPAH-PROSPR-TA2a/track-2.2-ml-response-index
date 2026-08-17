source(file.path("R", "check_environment.R"))
source(file.path("R", "validate_inputs.R"))
source(file.path("R", "build_training_features.R"))
source(file.path("R", "replay_preprocessing.R"))
source(file.path("R", "score_models.R"))


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


.list_exported_model_packages <- function(models_dir) {
  if (!dir.exists(models_dir)) {
    stop("models_dir does not exist: ", models_dir)
  }

  model_paths <- list.files(
    models_dir,
    pattern = "\\.json$",
    full.names = TRUE,
    recursive = FALSE
  )
  if (length(model_paths) == 0L) {
    stop("models_dir must contain exported model JSON packages.")
  }

  model_paths
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


.evaluate_model_package <- function(pheno,
                                    omics,
                                    model_package,
                                    output_dir = NULL,
                                    return_matrix = FALSE) {
  additional_covariates <- .package_vector(model_package$covariates$additional_covariates)
  model_covariates <- .package_vector(model_package$covariates$model_covariates)
  if (length(model_covariates) == 0L) {
    model_covariates <- unique(c("FEMALE", additional_covariates))
  }

  inputs <- .prepare_inputs(
    pheno = pheno,
    omics = omics,
    omics_type = model_package$omics_type,
    additional_covariates = NULL
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


FAST_evaluate <- function(pheno,
                          omics,
                          model_path,
                          output_dir = NULL,
                          return_matrix = FALSE) {
  model_package <- .load_exported_model_package(model_path)
  result <- .evaluate_model_package(
    pheno = pheno,
    omics = omics,
    model_package = model_package,
    output_dir = output_dir,
    return_matrix = return_matrix
  )
  message(
    "Validation complete for ", result$model_id,
    if (is.null(result$output_dir)) "." else paste0(". Results: ", result$output_dir)
  )
  result
}


FAST_bulk_evaluate <- function(pheno,
                               omics,
                               models_dir,
                               output_dir,
                               only_successful = TRUE) {
  model_paths <- .list_exported_model_packages(models_dir)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  rows <- list()
  row_index <- 1L
  for (model_path in model_paths) {
    model_package <- tryCatch(
      .load_exported_model_package(model_path),
      error = function(e) {
        stop(
          "Failed loading model package (", model_path, "): ",
          conditionMessage(e),
          call. = FALSE
        )
      }
    )

    if (only_successful && !isTRUE(model_package$successful)) {
      next
    }

    model_output_dir <- file.path(output_dir, "models", model_package$model_id)
    evaluate_one <- function() {
      .evaluate_model_package(
        pheno = pheno,
        omics = omics,
        model_package = model_package,
        output_dir = model_output_dir
      )
    }
    evaluated <- tryCatch(
      if (length(rows) == 0L) {
        evaluate_one()
      } else {
        suppressMessages(evaluate_one())
      },
      error = function(e) {
        stop(
          "Failed evaluating model ", model_package$model_id,
          " (", model_path, "): ", conditionMessage(e),
          call. = FALSE
        )
      }
    )

    validation <- evaluated$validation
    validation$MODEL_ID <- model_package$model_id
    validation$MODEL_PATH <- .relative_path(model_path, output_dir)
    validation$PREDICTIONS_PATH <- .relative_path(evaluated$output_files$predictions, output_dir)
    validation$VALIDATION_PATH <- .relative_path(evaluated$output_files$validation, output_dir)
    rows[[row_index]] <- validation[
      ,
      c(
        "MODEL_ID", "MODEL_PATH", "FU", "MODEL", "TRAINING_CV_AUC",
        "SUCCESS_AUC_THRESHOLD", "SUCCESSFUL", "N", "N_CONTROL",
        "N_TREATMENT", "AUC", "LOGIT_BETA", "LOGIT_OR", "LOGIT_P",
        "VALIDATED_P05", "PREDICTIONS_PATH", "VALIDATION_PATH"
      ),
      drop = FALSE
    ]
    row_index <- row_index + 1L
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

  message(
    "Cross-trial validation complete: ",
    if (nrow(summary) == 0L) {
      "no eligible models were evaluated. "
    } else if (nrow(summary) == 1L) {
      "1 model evaluated. "
    } else {
      paste0(nrow(summary), " models evaluated. ")
    },
    "Summary: ", normalizePath(summary_path, mustWork = FALSE)
  )

  list(
    validation = summary,
    output_dir = normalizePath(output_dir, mustWork = FALSE),
    output_files = list(validation_summary = summary_path)
  )
}
