source(file.path("R", "validate_inputs.R"))
source(file.path("R", "build_training_features.R"))
source(file.path("R", "replay_preprocessing.R"))
source(file.path("R", "score_models.R"))


FAST_treatment_predict <- function(pheno,
                                   omics,
                                   manifest_path,
                                   models = NULL,
                                   followups = NULL,
                                   omics_type = NULL,
                                   enet_mode = "reproduce_training") {
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
  pred_index <- 1L

  for (fu_key in followup_names) {
    fu_manifest <- manifest$followups[[fu_key]]
    if (is.null(fu_manifest)) next

    fu_level <- as.integer(sub("^FU", "", fu_key))
    matrices[[fu_key]] <- list()

    for (model_name in requested_models) {
      model_manifest <- fu_manifest$models[[model_name]]
      if (is.null(model_manifest)) next

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

      predictions[[pred_index]] <- data.frame(
        SET = "train",
        SUBJECT_ID = replay$subject_ids,
        FU = fu_level,
        MODEL = model_name,
        TREATMENT_GROUP = replay$y,
        PREDICTED_PROB = pred,
        stringsAsFactors = FALSE
      )
      pred_index <- pred_index + 1L
    }
  }

  list(
    predictions = do.call(rbind, predictions),
    matrices = matrices,
    manifest = manifest
  )
}
