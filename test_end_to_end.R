source("main.R")
source("test_helpers.R")

cat("\nEnd-to-end pipeline\n")
cat("===================\n")

fixture <- .make_simulated_fixture()
output_dir <- file.path(getwd(), "test_outputs", "track22_integration")
if (dir.exists(output_dir)) {
  unlink(output_dir, recursive = TRUE)
}

manifest <- suppressWarnings(FAST_treatment_ML(
  pheno = fixture$pheno,
  omics = fixture$omics,
  omics_type = "Proteomics",
  additional_covariates = c("FEMALE", "age", "bmi", "site", "smoker"),
  models = c("enet", "xgb"),
  output_dir = output_dir,
  test_frac = 0.2,
  cv_folds = 5L,
  seed = 2202L,
  n_cores = 1L,
  xgb_n_trials = 10L,
  python_bin = "python"
))

fu1 <- manifest$followups$FU1
.expect_true(!is.null(fu1), "all-subject FU1 run completes")
.expect_true(
  identical(normalizePath(output_dir), manifest$output_dir),
  "test outputs persist in the working directory"
)
.expect_true(file.exists(manifest$manifest_path), "manifest is written")
.expect_true(!is.null(jsonlite::fromJSON(manifest$manifest_path)), "manifest JSON is parseable")

artifact_paths <- unlist(fu1$artifacts, use.names = FALSE)
.expect_true(all(file.exists(artifact_paths)), "canonical artifacts are written")

subjects <- read.csv(fu1$artifacts$subjects, stringsAsFactors = FALSE)
enet_train <- read.csv(fu1$artifacts$enet_train, check.names = FALSE)
enet_test <- read.csv(fu1$artifacts$enet_test, check.names = FALSE)
xgb_train <- read.csv(fu1$artifacts$xgb_train, check.names = FALSE)
xgb_test <- read.csv(fu1$artifacts$xgb_test, check.names = FALSE)
preprocessing <- read.csv(fu1$artifacts$preprocessing, stringsAsFactors = FALSE)

.expect_true(
  length(intersect(
    subjects$SUBJECT_ID[subjects$SET == "train"],
    subjects$SUBJECT_ID[subjects$SET == "test"]
  )) == 0L,
  "integration split has no subject overlap"
)
.expect_true(
  nrow(enet_train) == sum(subjects$SET == "train") &&
    nrow(enet_test) == sum(subjects$SET == "test") &&
    nrow(xgb_train) == sum(subjects$SET == "train") &&
    nrow(xgb_test) == sum(subjects$SET == "test") &&
    all(!is.na(subjects$FOLD_ID[subjects$SET == "train"])) &&
    all(is.na(subjects$FOLD_ID[subjects$SET == "test"])),
  "canonical artifact row counts agree"
)
.expect_true(
  any(startsWith(names(enet_train), "covariate::age")) &&
    "covariate::FEMALE" %in% names(xgb_train) &&
    sum(startsWith(names(xgb_train), "covariate::FEMALE")) == 1L &&
    !any(startsWith(names(xgb_train), "covariate::age")),
  "model-ready matrices preserve feature boundaries"
)
.expect_true(
  identical(
    names(preprocessing),
    c("FEATURE_NAME", "FEATURE_TYPE", "MEDIAN", "CENTER", "SCALE")
  ) &&
    identical(preprocessing$FEATURE_NAME, names(enet_train)) &&
    all(is.finite(as.matrix(preprocessing[, c("MEDIAN", "CENTER", "SCALE")]))) &&
    all(names(xgb_train) %in% preprocessing$FEATURE_NAME),
  "preprocessing CSV describes retained model features"
)

for (model_name in c("enet", "xgb")) {
  model <- fu1$models[[model_name]]
  .expect_true(
    all(file.exists(unlist(model, use.names = FALSE))),
    paste(model_name, "declared outputs are written")
  )
  metrics <- read.csv(model$metrics)
  .expect_true(
    nrow(metrics) == 1L &&
      all(c("CV_AUC", "TEST_AUC", "INSAMPLE_AUC") %in% names(metrics)),
    paste(model_name, "metrics schema")
  )
  .assert_probability_file(model$predictions, nrow(subjects))
}

enet_weights <- read.csv(fu1$models$enet$weights, stringsAsFactors = FALSE)
enet_metrics <- read.csv(fu1$models$enet$metrics, stringsAsFactors = FALSE)
.expect_true(
  identical(
    names(enet_weights),
    c("FEATURE_NAME", "WEIGHT", "FEATURE_TYPE")
  ) &&
    all(
      names(enet_train)[startsWith(names(enet_train), "covariate::")] %in%
        enet_weights$FEATURE_NAME
    ),
  "ENET weights are consolidated"
)
.expect_true(
  enet_metrics$N_UNPENALIZED ==
    sum(startsWith(names(enet_train), "covariate::")),
  "ENET reports all covariates as unpenalized"
)

enet_predictions <- read.csv(fu1$models$enet$predictions, stringsAsFactors = FALSE)
intercept <- enet_weights$WEIGHT[enet_weights$FEATURE_NAME == "(Intercept)"]
feature_weights <- enet_weights[enet_weights$FEATURE_NAME != "(Intercept)", ]
score_matrix <- rbind(enet_train, enet_test)
linear_predictor <- rep(intercept, nrow(score_matrix))
if (nrow(feature_weights) > 0L) {
  linear_predictor <- linear_predictor +
    as.numeric(
      as.matrix(score_matrix[, feature_weights$FEATURE_NAME, drop = FALSE]) %*%
        feature_weights$WEIGHT
    )
}
reconstructed_prob <- stats::plogis(linear_predictor)
.expect_equal(
  reconstructed_prob,
  enet_predictions$PREDICTED_PROB,
  "ENET weights reproduce saved probabilities",
  tolerance = 1e-10
)
.expect_true(
  !is.null(jsonlite::fromJSON(fu1$models$xgb$model, simplifyVector = FALSE)),
  "XGB model JSON reloads"
)
.expect_true(
  file.exists(fu1$models$xgb$tuning),
  "XGB tuning output is written"
)
.expect_true(
  any(startsWith(names(read.csv(fu1$models$xgb$metrics)), "PARAM_")),
  "XGB selected parameters are recorded in metrics"
)

reports <- suppressWarnings(FAST_treatment_ML_reports(
  pheno = fixture$pheno,
  omics = fixture$omics,
  omics_type = "Proteomics",
  additional_covariates = c("FEMALE", "age", "bmi", "site", "smoker")
))
.expect_true(
  identical(
    names(reports),
    c("pheno_summary", "variable_summaries", "randomization_reports")
  ),
  "reporting returns the documented top-level structure"
)
.expect_true(
  !is.null(reports$randomization_reports$omics_report) &&
    !is.null(reports$randomization_reports$covariate_report),
  "randomization reports are populated"
)
