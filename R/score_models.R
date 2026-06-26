.score_enet_weights <- function(x, weights_path) {
  weights <- read.csv(weights_path, stringsAsFactors = FALSE)
  intercept <- weights$WEIGHT[weights$FEATURE_NAME == "(Intercept)"]
  if (length(intercept) != 1L) {
    stop("ENET weights must contain exactly one intercept row.")
  }

  feature_weights <- weights[weights$FEATURE_NAME != "(Intercept)", , drop = FALSE]
  missing_features <- setdiff(feature_weights$FEATURE_NAME, colnames(x))
  if (length(missing_features) > 0L) {
    stop("ENET matrix missing feature(s): ", paste(missing_features, collapse = ", "))
  }

  linear_predictor <- rep(intercept, nrow(x))
  if (nrow(feature_weights) > 0L) {
    linear_predictor <- linear_predictor +
      as.numeric(as.matrix(x[, feature_weights$FEATURE_NAME, drop = FALSE]) %*% feature_weights$WEIGHT)
  }

  stats::plogis(linear_predictor)
}


.score_xgb_model <- function(x, model_path) {
  if (!requireNamespace("xgboost", quietly = TRUE)) {
    stop("Package 'xgboost' is required for XGB inference.")
  }

  model <- xgboost::xgb.load(model_path)
  dmat <- xgboost::xgb.DMatrix(as.matrix(x))
  as.numeric(stats::predict(model, dmat))
}
