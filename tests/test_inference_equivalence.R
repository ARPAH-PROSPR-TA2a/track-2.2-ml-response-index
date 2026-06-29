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
      "FU", "MODEL", "TRAINING_CV_AUC", "CATALOG_AUC_THRESHOLD",
      "CATALOGED", "N", "N_CONTROL", "N_TREATMENT", "AUC",
      "LOGIT_BETA", "LOGIT_OR", "LOGIT_P", "VALIDATED_P05"
    ),
    "inference validation schema"
  )
  .expect_true(
    file.exists(inference$output_files$validation) &&
      all(inference$validation$CATALOG_AUC_THRESHOLD == 0.8) &&
      all(inference$validation$CATALOGED ==
            (inference$validation$TRAINING_CV_AUC >= 0.8)),
    "inference validation outputs are written and gated by training CV AUC"
  )

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
      if (isTRUE(validation_row$CATALOGED)) {
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
