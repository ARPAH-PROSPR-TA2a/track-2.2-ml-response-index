source("main.R")
source("test_helpers.R")

cat("\nValidation\n")
cat("==========\n")

fixture <- .make_simulated_fixture(n_subjects = 12L, n_features = 6L)

missing_column <- fixture$pheno
missing_column$SAMPLE_ID <- NULL
.expect_error(
  suppressWarnings(.validate_pheno(missing_column)),
  "Missing required columns: SAMPLE_ID",
  "missing required phenotype columns are rejected"
)

invalid_treatment <- fixture$pheno
invalid_treatment$TREATMENT_GROUP[1] <- 2L
.expect_error(
  suppressWarnings(.validate_pheno(invalid_treatment)),
  "TREATMENT_GROUP must contain only values 0 or 1",
  "invalid treatment values are rejected"
)

duplicate_sample <- fixture$pheno
duplicate_sample$SAMPLE_ID[2] <- duplicate_sample$SAMPLE_ID[1]
.expect_error(
  suppressWarnings(.validate_pheno(duplicate_sample)),
  "SAMPLE_ID contains duplicate values",
  "duplicate sample IDs are rejected"
)

validated_pheno <- suppressWarnings(.validate_pheno(fixture$pheno))
no_overlap <- fixture$omics
names(no_overlap)[-1] <- paste0("OTHER_", seq_len(ncol(no_overlap) - 1L))
.expect_error(
  suppressWarnings(.validate_omics(no_overlap, validated_pheno)),
  "No overlap between omics column names and pheno SAMPLE_IDs",
  "omics without shared samples are rejected"
)

.expect_error(
  .validate_ml_args(
    models = "xgb",
    test_frac = 0.2,
    enet_cv_folds = 5L,
    xgb_cv_folds = 5L,
    xgb_cv_repeats = 3L,
    seed = 1L,
    n_cores = 1L,
    xgb_n_trials = 9L
  ),
  "xgb_n_trials must be a single integer >= 10 when XGB is requested",
  "XGB rejects fewer than 10 tuning trials"
)

trial_warning <- NULL
invisible(withCallingHandlers(
  .validate_ml_args(
    models = "xgb",
    test_frac = 0.2,
    enet_cv_folds = 5L,
    xgb_cv_folds = 5L,
    xgb_cv_repeats = 3L,
    seed = 1L,
    n_cores = 1L,
    xgb_n_trials = 20L
  ),
  warning = function(w) {
    trial_warning <<- conditionMessage(w)
    invokeRestart("muffleWarning")
  }
))
.expect_true(
  identical(
    trial_warning,
    "xgb_n_trials below 30 provides a limited hyperparameter search; 30 or more trials are recommended."
  ),
  "XGB warns when fewer than 30 tuning trials are requested"
)

.expect_equal(
  suppressWarnings(.validate_ml_args(
    models = "enet",
    test_frac = 0.2,
    enet_cv_folds = 5L,
    xgb_cv_folds = 5L,
    xgb_cv_repeats = 3L,
    seed = 1L,
    n_cores = 1L,
    xgb_n_trials = 0L
  )),
  "enet",
  "ENET-only runs ignore the XGB trial setting"
)

too_small <- data.frame(
  SUBJECT_ID = paste0("S", 1:6),
  TREATMENT_GROUP = factor(c(0, 0, 0, 1, 1, 1), levels = 0:1)
)
.expect_equal(
  .validate_followup_cohort(
    too_small,
    test_frac = 0.2,
    min_train_per_class = 3L
  ),
  "each treatment arm must retain at least 3 training subjects after the test split",
  "undersized follow-up cohorts are rejected before splitting"
)

large_enough <- data.frame(
  SUBJECT_ID = paste0("S", 1:8),
  TREATMENT_GROUP = factor(c(0, 0, 0, 0, 1, 1, 1, 1), levels = 0:1)
)
.expect_true(
  is.null(.validate_followup_cohort(
    large_enough,
    test_frac = 0.2,
    min_train_per_class = 3L
  )),
  "follow-up cohorts with enough post-split training subjects pass validation"
)

validated_split <- .stratified_subject_split(
  large_enough,
  test_frac = 0.2,
  seed = 1L
)
.expect_equal(
  length(validated_split$test_subjects),
  2L,
  "validated stratified split assigns one test subject per arm"
)
