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
  model_covariates = "age",
  enet_cv_folds = 2L,
  xgb_cv_folds = 2L,
  xgb_cv_repeats = 3L,
  seed = 1L
)

.expect_equal(prepared$subject_ids_train, paste0("S", 1:4), "training subjects stay aligned")
.expect_equal(prepared$subject_ids_test, paste0("S", 5:6), "test subjects stay aligned")
.expect_equal(prepared$y_train, c(0L, 1L, 0L, 1L), "training outcomes stay aligned")
.expect_true(
  length(intersect(prepared$subject_ids_train, prepared$subject_ids_test)) == 0L,
  "train and test subjects do not overlap"
)

missing_only_in_train <- .drop_all_missing_train(
  x_train = matrix(
    c(NA, NA, 1, 2),
    nrow = 2,
    dimnames = list(NULL, c("all_missing", "retained"))
  ),
  x_test = matrix(
    c(10, 11, 3, 4),
    nrow = 2,
    dimnames = list(NULL, c("all_missing", "retained"))
  )
)
.expect_equal(
  colnames(missing_only_in_train$train),
  "retained",
  "all-missing training features are removed before imputation"
)
.expect_equal(
  colnames(missing_only_in_train$test),
  "retained",
  "test observations do not retain an all-missing training feature"
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
  "covariate::FEMALE" %in% colnames(prepared$enet_x_train) &&
    any(startsWith(colnames(prepared$enet_x_train), "covariate::age")),
  "ENET includes FEMALE and requested additional covariates"
)
.expect_true(
  "covariate::FEMALE" %in% colnames(prepared$xgb_x_train) &&
    sum(startsWith(colnames(prepared$xgb_x_train), "covariate::FEMALE")) == 1L &&
    !any(startsWith(colnames(prepared$xgb_x_train), "covariate::age")),
  "XGB includes FEMALE by default and excludes other covariates"
)
.expect_true(
  {
    duplicate_female <- .prepare_fu_change_dataset(
      pheno_df = fixture$pheno,
      omics_df = fixture$omics,
      fu_level = 1L,
      split = fixture$split,
      model_covariates = c("FEMALE", "FEMALE", "age"),
      enet_cv_folds = 2L,
      xgb_cv_folds = 2L,
      xgb_cv_repeats = 1L,
      seed = 1L
    )
    sum(startsWith(colnames(duplicate_female$enet_x_train), "covariate::FEMALE")) == 1L &&
      sum(startsWith(colnames(duplicate_female$xgb_x_train), "covariate::FEMALE")) == 1L
  },
  "explicit FEMALE entries do not create duplicate model features"
)

constant_female_pheno <- fixture$pheno
constant_female_pheno$FEMALE[
  constant_female_pheno$SUBJECT_ID %in% fixture$split$train_subjects
] <- factor(0L, levels = 0:1)
constant_female <- .prepare_fu_change_dataset(
  pheno_df = constant_female_pheno,
  omics_df = fixture$omics,
  fu_level = 1L,
  split = fixture$split,
  model_covariates = "age",
  enet_cv_folds = 2L,
  xgb_cv_folds = 2L,
  xgb_cv_repeats = 1L,
  seed = 1L
)
.expect_true(
  !any(startsWith(colnames(constant_female$enet_x_train), "covariate::FEMALE")) &&
    !any(startsWith(colnames(constant_female$xgb_x_train), "covariate::FEMALE")),
  "training-zero-variance FEMALE is removed from both models"
)
.expect_equal(
  prepared$preprocessing$FEATURE_NAME,
  colnames(prepared$enet_x_train),
  "preprocessing rows match retained ENET features"
)
.expect_true(
  length(prepared$enet_foldid) == length(prepared$subject_ids_train) &&
    nrow(prepared$xgb_folds) == length(prepared$subject_ids_train) * 3L &&
    identical(sort(unique(prepared$xgb_folds$REPEAT)), 1:3),
  "ENET and repeated XGB folds are prepared separately"
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
  enet_cv_folds = 2L,
  xgb_cv_folds = 2L,
  xgb_cv_repeats = 3L,
  seed = 1L
)
.expect_true(
  !"S4" %in% prepared_missing$subject_ids_train &&
    all(c("S1", "S2", "S3", "S5", "S6") %in% prepared_missing$subject_ids_train),
  "a missing visit excludes only that subject from the follow-up dataset"
)
