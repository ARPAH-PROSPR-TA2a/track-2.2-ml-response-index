source("main.R")
source("test_helpers.R")

cat("\nEnd-to-end pipeline\n")
cat("===================\n")

fixture <- .make_simulated_fixture(followups = 1:2)
output_dir <- file.path(getwd(), "test_outputs", "track22_integration")
if (dir.exists(output_dir)) {
  unlink(output_dir, recursive = TRUE)
}

manifest <- suppressWarnings(FAST_treatment_ML(
  pheno = fixture$pheno,
  omics = fixture$omics,
  omics_type = "Proteomics",
  additional_covariates = c("age", "bmi", "site", "smoker"),
  models = c("enet", "xgb"),
  output_dir = output_dir,
  test_frac = 0.2,
  enet_cv_folds = 5L,
  xgb_cv_folds = 5L,
  xgb_cv_repeats = 3L,
  seed = 2202L,
  n_cores = 1L,
  xgb_n_trials = 10L,
  python_bin = "python"
))

fu1 <- manifest$followups$FU1
fu2 <- manifest$followups$FU2
.expect_true(!is.null(fu1) && !is.null(fu2), "all-subject FU1 and FU2 runs complete")
.expect_true(
  identical(normalizePath(output_dir), manifest$output_dir),
  "test outputs persist in the working directory"
)
.expect_true(file.exists(manifest$manifest_path), "manifest is written")
.expect_true(!is.null(jsonlite::fromJSON(manifest$manifest_path)), "manifest JSON is parseable")

artifact_paths <- unlist(fu1$artifacts, use.names = FALSE)
.expect_true(all(file.exists(artifact_paths)), "canonical artifacts are written")

subjects <- read.csv(fu1$artifacts$subjects, stringsAsFactors = FALSE)
xgb_folds <- read.csv(fu1$artifacts$xgb_folds, stringsAsFactors = FALSE)
enet_train <- read.csv(fu1$artifacts$enet_train, check.names = FALSE)
enet_test <- read.csv(fu1$artifacts$enet_test, check.names = FALSE)
xgb_train <- read.csv(fu1$artifacts$xgb_train, check.names = FALSE)
xgb_test <- read.csv(fu1$artifacts$xgb_test, check.names = FALSE)
cohort <- read.csv(manifest$reports$cohort, stringsAsFactors = FALSE)
change_summary <- read.csv(manifest$reports$change_summary, stringsAsFactors = FALSE)
preprocessing <- read.csv(manifest$reports$preprocessing, stringsAsFactors = FALSE)

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
    all(!is.na(subjects$ENET_FOLD_ID[subjects$SET == "train"])) &&
    all(is.na(subjects$ENET_FOLD_ID[subjects$SET == "test"])),
  "canonical artifact row counts agree"
)
.expect_true(
  identical(names(xgb_folds), c("SUBJECT_ID", "REPEAT", "FOLD_ID")) &&
    nrow(xgb_folds) == sum(subjects$SET == "train") * 3L &&
    identical(sort(unique(xgb_folds$REPEAT)), 1:3) &&
    all(table(xgb_folds$SUBJECT_ID) == 3L),
  "XGB fold artifact contains every training subject in every repeat"
)
.expect_true(
  any(startsWith(names(enet_train), "covariate::age")) &&
    "covariate::FEMALE" %in% names(xgb_train) &&
    sum(startsWith(names(xgb_train), "covariate::FEMALE")) == 1L &&
    !any(startsWith(names(xgb_train), "covariate::age")),
  "model-ready matrices preserve feature boundaries"
)
.expect_true(
  all(file.exists(unlist(manifest$reports, use.names = FALSE))),
  "consolidated report artifacts are written"
)
.expect_true(
  all(c("FU1", "FU2") %in% names(manifest$followups)),
  "manifest records both modeled follow-ups"
)
.expect_true(
  identical(
    names(cohort),
    c(
      "FU", "SET", "N_SUBJECTS", "N_CONTROL", "N_TREATMENT",
      "N_MALE", "N_FEMALE"
    )
  ) &&
    identical(sort(unique(cohort$FU)), 1:2) &&
    all(table(cohort$FU) == 3L) &&
    all(c("eligible", "train", "test") %in% cohort$SET) &&
    cohort$N_SUBJECTS[cohort$FU == 1L & cohort$SET == "eligible"] == nrow(subjects) &&
    cohort$N_SUBJECTS[cohort$FU == 1L & cohort$SET == "train"] == sum(subjects$SET == "train") &&
    cohort$N_SUBJECTS[cohort$FU == 1L & cohort$SET == "test"] == sum(subjects$SET == "test"),
  "cohort artifact stacks eligible, train, and test subjects across follow-ups"
)
.expect_true(
  identical(
    names(change_summary),
    c(
      "FU", "SET", "TREATMENT_GROUP", "ANALYTE_NAME", "N_SUBJECTS",
      "N_NONMISSING", "MEAN", "MEDIAN", "SD", "MIN", "MAX"
    )
  ) &&
    identical(sort(unique(change_summary$FU)), 1:2) &&
    all(c("eligible", "train", "test") %in% change_summary$SET) &&
    all(0:1 %in% change_summary$TREATMENT_GROUP) &&
    all(change_summary$ANALYTE_NAME %in% fixture$omics$ANALYTE_NAME),
  "change summary artifact stacks raw changes by set and treatment"
)
.expect_true(
  identical(
    names(preprocessing),
    c(
      "FU", "FEATURE_NAME", "FEATURE_TYPE", "STATUS", "MEDIAN", "CENTER",
      "SCALE", "IN_ENET", "IN_XGB"
    )
  ) &&
    identical(sort(unique(preprocessing$FU)), 1:2) &&
    identical(
      preprocessing$FEATURE_NAME[preprocessing$FU == 1L & preprocessing$IN_ENET],
      names(enet_train)
    ) &&
    identical(
      preprocessing$FEATURE_NAME[preprocessing$FU == 1L & preprocessing$IN_XGB],
      names(xgb_train)
    ) &&
    all(
      is.finite(as.matrix(
        preprocessing[
          preprocessing$STATUS == "retained",
          c("MEDIAN", "CENTER", "SCALE")
        ]
      ))
    ),
  "preprocessing CSV audits retained and removed model features"
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
  file.exists(fu1$models$xgb$tuning) &&
    all(
      c("CV_AUC_MEAN", "CV_AUC_SD", "BEST_ITERATION_MEDIAN") %in%
        names(read.csv(fu1$models$xgb$tuning))
    ),
  "XGB repeated-CV tuning output is written"
)
.expect_true(
  {
    xgb_metrics <- read.csv(fu1$models$xgb$metrics)
    any(startsWith(names(xgb_metrics), "PARAM_")) &&
      xgb_metrics$CV_REPEATS == 3L &&
      is.finite(xgb_metrics$CV_AUC_SD)
  },
  "XGB selected parameters and repeated-CV metrics are recorded"
)
