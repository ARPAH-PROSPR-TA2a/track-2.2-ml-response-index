source("main.R")
source("test_helpers.R")

cat("\nFeature correctness\n")
cat("===================\n")

fixture <- .make_exact_fixture()
prepared <- .prepare_fu_change_dataset(
  pheno_df = fixture$pheno,
  omics_df = fixture$omics,
  fu_level = 1L,
  split = fixture$split,
  additional_covariates = c("FEMALE", "age"),
  cv_folds = 2L,
  seed = 1L
)

.expect_equal(prepared$subject_ids_train, paste0("S", 1:4), "training subjects stay aligned")
.expect_equal(prepared$subject_ids_test, paste0("S", 5:6), "test subjects stay aligned")
.expect_equal(prepared$y_train, c(0L, 1L, 0L, 1L), "training outcomes stay aligned")
.expect_true(
  length(intersect(prepared$subject_ids_train, prepared$subject_ids_test)) == 0L,
  "train and test subjects do not overlap"
)

feature_idx <- match("omics::feature_change", colnames(prepared$xgb_x_train))
expected_train <- as.numeric(scale(c(1, 3, 5, 7)))
expected_test <- (c(9, 11) - mean(c(1, 3, 5, 7))) / sd(c(1, 3, 5, 7))
.expect_equal(
  prepared$xgb_x_train[, feature_idx],
  expected_train,
  "follow-up minus baseline is computed and scaled correctly"
)
.expect_equal(
  prepared$xgb_x_test[, feature_idx],
  expected_test,
  "test changes use training scaling parameters"
)

missing_idx <- match("omics::missing_change", colnames(prepared$xgb_x_train))
.expect_equal(
  prepared$preprocessing$MEDIAN[
    prepared$preprocessing$FEATURE_NAME == "omics::missing_change"
  ],
  4,
  "imputation median comes from training subjects"
)
.expect_equal(
  prepared$xgb_x_train[1, missing_idx],
  0,
  "missing training change is median-imputed"
)
.expect_true(
  !"omics::constant_change" %in% colnames(prepared$xgb_x_train),
  "training-zero-variance omics features are removed"
)
.expect_true(
  any(startsWith(colnames(prepared$enet_x_train), "covariate::")),
  "ENET includes requested covariates"
)
.expect_true(
  "covariate::FEMALE" %in% colnames(prepared$xgb_x_train) &&
    sum(startsWith(colnames(prepared$xgb_x_train), "covariate::FEMALE")) == 1L &&
    !any(startsWith(colnames(prepared$xgb_x_train), "covariate::age")),
  "XGB includes one FEMALE column and excludes other covariates"
)
.expect_equal(
  prepared$preprocessing$FEATURE_NAME,
  colnames(prepared$enet_x_train),
  "preprocessing rows match retained ENET features"
)

missing_visit_pheno <- fixture$pheno[
  !(fixture$pheno$SUBJECT_ID == "S4" & fixture$pheno$FU == 1),
]
prepared_missing <- .prepare_fu_change_dataset(
  pheno_df = missing_visit_pheno,
  omics_df = fixture$omics,
  fu_level = 1L,
  split = list(
    train_subjects = paste0("S", 1:6),
    test_subjects = paste0("S", 7:8)
  ),
  cv_folds = 2L,
  seed = 1L
)
.expect_true(
  !"S4" %in% prepared_missing$subject_ids_train &&
    all(c("S1", "S2", "S3", "S5", "S6") %in% prepared_missing$subject_ids_train),
  "a missing visit excludes only that subject from the follow-up dataset"
)
