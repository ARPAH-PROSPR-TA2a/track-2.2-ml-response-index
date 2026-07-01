.score_enet_weight_table <- function(x, weights) {
  intercept <- weights$WEIGHT[weights$FEATURE_NAME == "(Intercept)"]
  if (length(intercept) != 1L) {
    stop("ENET weights must contain exactly one intercept row.")
  }

  feature_weights <- weights[weights$FEATURE_NAME != "(Intercept)", , drop = FALSE]
  missing_features <- setdiff(feature_weights$FEATURE_NAME, colnames(x))
  if (length(missing_features) > 0L) {
    omitted_adjustment <- startsWith(missing_features, "covariate::") &
      missing_features != "covariate::FEMALE"
    if (!all(omitted_adjustment)) {
      stop(
        "ENET matrix missing feature(s): ",
        paste(missing_features[!omitted_adjustment], collapse = ", ")
      )
    }
    feature_weights <- feature_weights[
      feature_weights$FEATURE_NAME %in% colnames(x),
      ,
      drop = FALSE
    ]
  }

  linear_predictor <- rep(intercept, nrow(x))
  if (nrow(feature_weights) > 0L) {
    linear_predictor <- linear_predictor +
      as.numeric(as.matrix(x[, feature_weights$FEATURE_NAME, drop = FALSE]) %*% feature_weights$WEIGHT)
  }

  stats::plogis(linear_predictor)
}


.score_enet_weights <- function(x, weights_path) {
  weights <- read.csv(weights_path, stringsAsFactors = FALSE)
  .score_enet_weight_table(x, weights)
}


.score_xgb_model <- function(x, model_path) {
  if (!requireNamespace("xgboost", quietly = TRUE)) {
    stop("Package 'xgboost' is required for XGB inference.")
  }

  model <- xgboost::xgb.load(model_path)
  dmat <- xgboost::xgb.DMatrix(as.matrix(x))
  as.numeric(stats::predict(model, dmat))
}


.score_xgb_model_json <- function(x, model_json) {
  model_path <- tempfile(fileext = ".json")
  on.exit(unlink(model_path), add = TRUE)
  writeLines(model_json, model_path, useBytes = TRUE)
  .score_xgb_model(x, model_path)
}


.success_auc_threshold <- function() {
  0.8
}


.model_cv_auc <- function(model_manifest) {
  metrics <- read.csv(model_manifest$metrics, stringsAsFactors = FALSE)
  if (!"CV_AUC" %in% names(metrics) || nrow(metrics) != 1L) {
    stop("Model metrics must contain exactly one CV_AUC value.")
  }
  as.numeric(metrics$CV_AUC)
}


.validation_row <- function(fu_level, model_name, predictions, training_cv_auc) {
  successful <- is.finite(training_cv_auc) &&
    training_cv_auc >= .success_auc_threshold()

  y <- as.integer(predictions$TREATMENT_GROUP)
  n_control <- sum(y == 0L)
  n_treatment <- sum(y == 1L)

  auc <- NA_real_
  logit_beta <- NA_real_
  logit_or <- NA_real_
  logit_p <- NA_real_

  if (successful && n_control > 0L && n_treatment > 0L) {
    auc <- .safe_auc(y, predictions$PREDICTED_PROB)
    fit <- stats::glm(
      TREATMENT_GROUP ~ PREDICTED_PROB,
      data = predictions,
      family = stats::binomial()
    )
    coef_table <- summary(fit)$coefficients
    if ("PREDICTED_PROB" %in% rownames(coef_table)) {
      logit_beta <- unname(coef_table["PREDICTED_PROB", "Estimate"])
      logit_or <- exp(logit_beta)
      logit_p <- unname(coef_table["PREDICTED_PROB", "Pr(>|z|)"])
    }
  }

  data.frame(
    FU = fu_level,
    MODEL = model_name,
    TRAINING_CV_AUC = training_cv_auc,
    SUCCESS_AUC_THRESHOLD = .success_auc_threshold(),
    SUCCESSFUL = successful,
    N = nrow(predictions),
    N_CONTROL = n_control,
    N_TREATMENT = n_treatment,
    AUC = auc,
    LOGIT_BETA = logit_beta,
    LOGIT_OR = logit_or,
    LOGIT_P = logit_p,
    VALIDATED_P05 = is.finite(logit_p) && logit_p < 0.05,
    stringsAsFactors = FALSE
  )
}
