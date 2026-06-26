.manifest_vector <- function(x) {
  if (is.null(x)) {
    return(character(0))
  }
  as.character(unlist(x, use.names = FALSE))
}


.load_manifest <- function(manifest_path) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required to read manifest files.")
  }
  jsonlite::fromJSON(manifest_path, simplifyVector = FALSE)
}


.model_preprocessing <- function(preprocessing, fu_level, model) {
  flag <- paste0("IN_", toupper(model))
  if (!flag %in% names(preprocessing)) {
    stop("preprocessing.csv is missing column ", flag, ".")
  }

  preprocessing[
    preprocessing$FU == fu_level &
      preprocessing$STATUS == "retained" &
      preprocessing[[flag]],
    ,
    drop = FALSE
  ]
}


.empty_matrix <- function(n) {
  matrix(numeric(0), nrow = n, ncol = 0)
}


.apply_preprocessing_recipe <- function(raw_matrix, feature_rows, prefix) {
  if (nrow(feature_rows) == 0L) {
    return(.empty_matrix(nrow(raw_matrix)))
  }

  raw_names <- sub(paste0("^", prefix), "", feature_rows$FEATURE_NAME)
  missing_features <- setdiff(raw_names, colnames(raw_matrix))
  if (length(missing_features) > 0L) {
    stop(
      "Missing raw feature(s) needed for preprocessing replay: ",
      paste(missing_features, collapse = ", ")
    )
  }

  out <- raw_matrix[, raw_names, drop = FALSE]
  for (j in seq_len(ncol(out))) {
    missing <- is.na(out[, j])
    if (any(missing)) {
      out[missing, j] <- feature_rows$MEDIAN[j]
    }
    out[, j] <- (out[, j] - feature_rows$CENTER[j]) / feature_rows$SCALE[j]
  }
  colnames(out) <- feature_rows$FEATURE_NAME
  out
}


.make_covariate_model_matrix <- function(pheno_followup, model_covariates) {
  if (length(model_covariates) == 0L) {
    return(.empty_matrix(nrow(pheno_followup)))
  }

  missing_covariates <- setdiff(model_covariates, names(pheno_followup))
  if (length(missing_covariates) > 0L) {
    stop("Missing covariate(s): ", paste(missing_covariates, collapse = ", "))
  }

  combined <- pheno_followup[, model_covariates, drop = FALSE]
  for (col in names(combined)) {
    if (col == "FEMALE") {
      combined[[col]] <- .as_binary_numeric(combined[[col]])
    } else if (is.logical(combined[[col]])) {
      combined[[col]] <- factor(combined[[col]], levels = c(FALSE, TRUE))
    }
  }

  model.matrix(~ . - 1, data = combined)
}


.build_inference_matrix <- function(pheno_df, omics_df, preprocessing,
                                    fu_level, model, model_covariates) {
  change_data <- .make_followup_change_data(pheno_df, omics_df, fu_level)
  if (is.null(change_data)) {
    stop("FU", fu_level, ": no subjects after requiring baseline and follow-up.")
  }

  feature_rows <- .model_preprocessing(preprocessing, fu_level, model)
  omics_rows <- feature_rows[feature_rows$FEATURE_TYPE == "omics", , drop = FALSE]
  covariate_rows <- feature_rows[feature_rows$FEATURE_TYPE == "covariate", , drop = FALSE]

  omics_matrix <- .apply_preprocessing_recipe(
    change_data$change_matrix,
    omics_rows,
    "omics::"
  )
  covariate_raw <- .make_covariate_model_matrix(
    change_data$pheno_followup,
    model_covariates
  )
  covariate_matrix <- .apply_preprocessing_recipe(
    covariate_raw,
    covariate_rows,
    "covariate::"
  )

  matrix_parts <- cbind(omics_matrix, covariate_matrix)
  matrix_parts <- matrix_parts[, feature_rows$FEATURE_NAME, drop = FALSE]

  list(
    matrix = matrix_parts,
    subject_ids = change_data$subject_ids,
    fu = fu_level,
    y = change_data$y
  )
}
