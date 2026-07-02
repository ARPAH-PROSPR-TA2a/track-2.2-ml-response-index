run_inference_equivalence_tests <- function(fixture, manifest) {
  cat("\nInference Equivalence\n")
  cat("=====================\n")

  exported <- FAST_export_models(
    manifest,
    output_dir = file.path(dirname(manifest$manifest_path), "exported_models")
  )
  .expect_true(
    nrow(exported$models) == sum(vapply(manifest$followups, function(fu_manifest) {
      if (is.null(fu_manifest)) return(0L)
      length(fu_manifest$models)
    }, integer(1))) &&
      all(exported$models$SUCCESSFUL ==
            (exported$models$TRAINING_CV_AUC >= 0.8)) &&
      all(file.exists(exported$models$PATH)) &&
      file.exists(exported$manifest),
    "model JSON packages are exported with success flags"
  )

  evaluated_models <- list()
  unlink(exported$manifest)
  bulk <- suppressWarnings(FAST_bulk_evaluate(
    pheno = fixture$pheno,
    omics = fixture$omics,
    models_dir = exported$output_dir,
    output_dir = file.path(dirname(manifest$manifest_path), "bulk_evaluation")
  ))
  .expect_true(
    nrow(bulk$validation) == sum(exported$models$SUCCESSFUL) &&
      all(bulk$validation$SUCCESSFUL) &&
      !file.exists(exported$manifest) &&
      file.exists(bulk$output_files$validation_summary) &&
      all(file.exists(file.path(dirname(bulk$output_files$validation_summary), bulk$validation$PREDICTIONS_PATH))) &&
      all(file.exists(file.path(dirname(bulk$output_files$validation_summary), bulk$validation$VALIDATION_PATH))),
    "bulk evaluation scores successful model packages without an export index"
  )

  for (export_index in seq_len(nrow(exported$models))) {
    exported_model <- exported$models[export_index, , drop = FALSE]
    evaluated <- suppressWarnings(FAST_evaluate(
      pheno = fixture$pheno,
      omics = fixture$omics,
      model_path = exported_model$PATH,
      output_dir = file.path(
        dirname(manifest$manifest_path),
        "single_model_evaluation",
        exported_model$MODEL_ID
      ),
      return_matrix = TRUE
    ))
    evaluated_models[[exported_model$MODEL_ID]] <- evaluated
    fu_key <- paste0("FU", exported_model$FU)
    fu_manifest <- manifest$followups[[fu_key]]
    train_matrix_path <- fu_manifest$artifacts[[paste0(exported_model$MODEL, "_train")]]
    train_matrix <- read.csv(train_matrix_path, check.names = FALSE)
    deployable_train_matrix <- train_matrix[
      ,
      names(train_matrix)[
        !startsWith(names(train_matrix), "covariate::") |
          names(train_matrix) == "covariate::FEMALE"
      ],
      drop = FALSE
    ]
    saved_predictions <- read.csv(
      fu_manifest$models[[exported_model$MODEL]]$predictions,
      stringsAsFactors = FALSE
    )

    .expect_equal(
      names(evaluated$validation),
      c(
        "FU", "MODEL", "TRAINING_CV_AUC", "SUCCESS_AUC_THRESHOLD",
        "SUCCESSFUL", "N", "N_CONTROL", "N_TREATMENT", "AUC",
        "LOGIT_BETA", "LOGIT_OR", "LOGIT_P", "VALIDATED_P05"
      ),
      paste(exported_model$MODEL_ID, "single-model validation schema")
    )
    .expect_equal(
      colnames(evaluated$matrix),
      names(deployable_train_matrix),
      paste(exported_model$MODEL_ID, "inference matrix columns match deployable training columns")
    )
    .expect_equal(
      as.matrix(evaluated$matrix),
      as.matrix(deployable_train_matrix),
      paste(exported_model$MODEL_ID, "inference matrix values match deployable training columns"),
      tolerance = 1e-10
    )
    .expect_equal(
      names(evaluated$predictions),
      c("SUBJECT_ID", "FU", "TREATMENT_GROUP", "PREDICTED_PROB"),
      paste(exported_model$MODEL_ID, "single-model prediction schema")
    )
    .expect_equal(
      evaluated$predictions$SUBJECT_ID,
      saved_predictions$SUBJECT_ID,
      paste(exported_model$MODEL_ID, "single-model prediction subjects match training")
    )
    .expect_true(
      evaluated$validation$SUCCESSFUL == exported_model$SUCCESSFUL &&
        evaluated$validation$SUCCESS_AUC_THRESHOLD == 0.8 &&
        file.exists(evaluated$output_files$predictions) &&
        file.exists(evaluated$output_files$validation),
      paste(exported_model$MODEL_ID, "single-model package validation outputs are written")
    )
    .expect_equal(
      evaluated$validation$TRAINING_CV_AUC,
      exported_model$TRAINING_CV_AUC,
      paste(exported_model$MODEL_ID, "single-model package training CV AUC"),
      tolerance = 1e-12
    )
    if (exported_model$MODEL == "xgb") {
      .expect_equal(
        evaluated$predictions$PREDICTED_PROB,
        saved_predictions$PREDICTED_PROB,
        paste(exported_model$MODEL_ID, "single-model XGB predictions match training artifact"),
        tolerance = 1e-7
      )
    } else {
      .expect_true(
        !any(startsWith(colnames(evaluated$matrix), "covariate::") &
               colnames(evaluated$matrix) != "covariate::FEMALE"),
        paste(exported_model$MODEL_ID, "single-model ENET omits training-only adjustment covariates")
      )
    }
    .expect_true(
      evaluated$validation$N == nrow(evaluated$predictions) &&
        evaluated$validation$N_CONTROL == sum(evaluated$predictions$TREATMENT_GROUP == 0L) &&
        evaluated$validation$N_TREATMENT == sum(evaluated$predictions$TREATMENT_GROUP == 1L),
      paste(exported_model$MODEL_ID, "single-model validation counts match predictions")
    )
    if (isTRUE(evaluated$validation$SUCCESSFUL)) {
      direct_fit <- stats::glm(
        TREATMENT_GROUP ~ PREDICTED_PROB,
        data = evaluated$predictions,
        family = stats::binomial()
      )
      direct_p <- summary(direct_fit)$coefficients["PREDICTED_PROB", "Pr(>|z|)"]
      .expect_equal(
        evaluated$validation$LOGIT_P,
        direct_p,
        paste(exported_model$MODEL_ID, "single-model logit p-value"),
        tolerance = 1e-10
      )
    }
  }

  for (bulk_index in seq_len(nrow(bulk$validation))) {
    bulk_row <- bulk$validation[bulk_index, , drop = FALSE]
    evaluated <- evaluated_models[[bulk_row$MODEL_ID]]
    .expect_equal(
      bulk_row[
        ,
        c(
          "FU", "MODEL", "TRAINING_CV_AUC", "SUCCESS_AUC_THRESHOLD",
          "SUCCESSFUL", "N", "N_CONTROL", "N_TREATMENT", "AUC",
          "LOGIT_BETA", "LOGIT_OR", "LOGIT_P", "VALIDATED_P05"
        ),
        drop = FALSE
      ],
      evaluated$validation,
      paste(bulk_row$MODEL_ID, "bulk validation matches single-model evaluation"),
      tolerance = 1e-10
    )
  }
}
