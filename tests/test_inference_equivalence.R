run_inference_equivalence_tests <- function(fixture, manifest) {
  cat("\nInference Equivalence\n")
  cat("=====================\n")

  inference <- suppressWarnings(FAST_treatment_predict(
    pheno = fixture$pheno,
    omics = fixture$omics,
    manifest_path = manifest$manifest_path,
    models = c("enet", "xgb"),
    omics_type = "Proteomics",
    enet_mode = "reproduce_training",
    output_dir = file.path(dirname(manifest$manifest_path), "inference"),
    return_matrices = TRUE
  ))

  .expect_equal(
    names(inference$validation),
    c(
      "FU", "MODEL", "TRAINING_CV_AUC", "SUCCESS_AUC_THRESHOLD",
      "SUCCESSFUL", "N", "N_CONTROL", "N_TREATMENT", "AUC",
      "LOGIT_BETA", "LOGIT_OR", "LOGIT_P", "VALIDATED_P05"
    ),
    "inference validation schema"
  )
  .expect_true(
    file.exists(inference$output_files$validation) &&
      all(inference$validation$SUCCESS_AUC_THRESHOLD == 0.8) &&
      all(inference$validation$SUCCESSFUL ==
            (inference$validation$TRAINING_CV_AUC >= 0.8)),
    "inference validation outputs are written and gated by training CV AUC"
  )

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

  bulk <- suppressWarnings(FAST_bulk_evaluate(
    pheno = fixture$pheno,
    omics = fixture$omics,
    models_dir = exported$output_dir,
    output_dir = file.path(dirname(manifest$manifest_path), "bulk_evaluation")
  ))
  .expect_true(
    nrow(bulk$validation) == sum(exported$models$SUCCESSFUL) &&
      all(bulk$validation$SUCCESSFUL) &&
      file.exists(bulk$output_files$validation_summary) &&
      all(file.exists(file.path(dirname(bulk$output_files$validation_summary), bulk$validation$PREDICTIONS_PATH))) &&
      all(file.exists(file.path(dirname(bulk$output_files$validation_summary), bulk$validation$VALIDATION_PATH))),
    "bulk evaluation scores successful exported models and writes a shareable summary"
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
    fu_key <- paste0("FU", exported_model$FU)
    manifest_predictions <- inference$predictions[[fu_key]][[exported_model$MODEL]]

    .expect_equal(
      evaluated$predictions$PREDICTED_PROB,
      manifest_predictions$PREDICTED_PROB,
      paste(exported_model$MODEL_ID, "single-model package predictions match manifest inference"),
      tolerance = if (exported_model$MODEL == "xgb") 1e-7 else 1e-10
    )
    .expect_true(
      evaluated$validation$SUCCESSFUL == exported_model$SUCCESSFUL &&
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
  }

  for (fu_key in names(manifest$followups)) {
    fu_manifest <- manifest$followups[[fu_key]]
    if (is.null(fu_manifest)) next

    fu_level <- as.integer(sub("^FU", "", fu_key))

    for (model_name in c("enet", "xgb")) {
      if (is.null(fu_manifest$models[[model_name]])) next

      train_matrix_path <- fu_manifest$artifacts[[paste0(model_name, "_train")]]
      train_matrix <- read.csv(train_matrix_path, check.names = FALSE)
      replayed_matrix <- inference$matrices[[fu_key]][[model_name]]

      .expect_equal(
        colnames(replayed_matrix),
        names(train_matrix),
        paste(fu_key, model_name, "inference matrix columns match training")
      )
      .expect_equal(
        as.matrix(replayed_matrix),
        as.matrix(train_matrix),
        paste(fu_key, model_name, "inference matrix values match training"),
        tolerance = 1e-10
      )

      saved_predictions <- read.csv(
        fu_manifest$models[[model_name]]$predictions,
        stringsAsFactors = FALSE
      )
      inferred_predictions <- inference$predictions[[fu_key]][[model_name]]
      validation_row <- inference$validation[
        inference$validation$FU == fu_level &
          inference$validation$MODEL == model_name,
        ,
        drop = FALSE
      ]

      .expect_equal(
        names(inferred_predictions),
        c("SUBJECT_ID", "FU", "TREATMENT_GROUP", "PREDICTED_PROB"),
        paste(fu_key, model_name, "inference prediction schema")
      )

      .expect_equal(
        inferred_predictions$SUBJECT_ID,
        saved_predictions$SUBJECT_ID,
        paste(fu_key, model_name, "inference prediction subjects match training")
      )
      .expect_equal(
        inferred_predictions$PREDICTED_PROB,
        saved_predictions$PREDICTED_PROB,
        paste(fu_key, model_name, "inference predictions match training"),
        tolerance = if (model_name == "xgb") 1e-7 else 1e-10
      )
      .expect_true(
        nrow(validation_row) == 1L &&
          validation_row$N == nrow(inferred_predictions) &&
          validation_row$N_CONTROL == sum(inferred_predictions$TREATMENT_GROUP == 0L) &&
          validation_row$N_TREATMENT == sum(inferred_predictions$TREATMENT_GROUP == 1L) &&
          file.exists(inference$output_files$predictions[[fu_key]][[model_name]]),
        paste(fu_key, model_name, "inference validation row and prediction file")
      )
      if (isTRUE(validation_row$SUCCESSFUL)) {
        direct_fit <- stats::glm(
          TREATMENT_GROUP ~ PREDICTED_PROB,
          data = inferred_predictions,
          family = stats::binomial()
        )
        direct_p <- summary(direct_fit)$coefficients["PREDICTED_PROB", "Pr(>|z|)"]
        .expect_equal(
          validation_row$LOGIT_P,
          direct_p,
          paste(fu_key, model_name, "inference logit p-value"),
          tolerance = 1e-10
        )
      }
    }
  }
}
