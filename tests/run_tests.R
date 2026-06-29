args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) == 1L) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg)))
  repo_dir <- if (basename(script_dir) == "tests") dirname(script_dir) else script_dir
  setwd(repo_dir)
}

cat("Track 2.2 test suite\n")
cat("====================\n")

options(track22.tests.passed = 0L)
source(file.path("training", "main.R"))
source(file.path("inference", "main.R"))
source(file.path("tests", "helpers.R"))
source(file.path("tests", "test_inference_equivalence.R"))

required_r <- c("glmnet", "jsonlite", "pROC", "xgboost")
missing_r <- required_r[!vapply(required_r, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_r) > 0L) {
  stop("Missing required R packages: ", paste(missing_r, collapse = ", "), call. = FALSE)
}
cat("PASS R dependencies are ready\n")

.test_python_bin <- function() {
  repo_python <- file.path(getwd(), ".track22-python", "bin", "python")
  if (file.exists(repo_python)) {
    return(repo_python)
  }

  python3 <- Sys.which("python3")
  if (nzchar(python3)) {
    return(unname(python3))
  }

  python <- Sys.which("python")
  if (nzchar(python)) {
    return(unname(python))
  }

  "python3"
}

run_validation_tests <- function() {
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
  .expect_true(
    is.data.frame(validated_pheno) &&
      !any(c("all", "male", "female", "requires_mixed_effects") %in%
             names(validated_pheno)),
    "phenotype validation returns the ML phenotype table"
  )

  single_sex_pheno <- fixture$pheno
  single_sex_pheno$FU <- factor(single_sex_pheno$FU)
  single_sex_pheno$TREATMENT_GROUP <- factor(single_sex_pheno$TREATMENT_GROUP, levels = 0:1)
  single_sex_pheno$FEMALE <- factor(1L, levels = 0:1)
  single_sex_warnings <- character(0)
  single_sex_validated <- withCallingHandlers(
    .validate_pheno(single_sex_pheno),
    warning = function(w) {
      single_sex_warnings <<- c(single_sex_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  .expect_true(
    is.data.frame(single_sex_validated) &&
      !any(grepl("Male and female subsets", single_sex_warnings, fixed = TRUE)),
    "single-sex phenotype input is accepted without legacy subset warnings"
  )

  no_overlap <- fixture$omics
  names(no_overlap)[-1] <- paste0("OTHER_", seq_len(ncol(no_overlap) - 1L))
  .expect_error(
    suppressWarnings(.validate_omics(no_overlap, validated_pheno)),
    "No overlap between omics column names and pheno SAMPLE_IDs",
    "omics without shared samples are rejected"
  )

  duplicate_analytes <- fixture$omics
  duplicate_analytes$ANALYTE_NAME[2] <- duplicate_analytes$ANALYTE_NAME[1]
  .expect_error(
    suppressWarnings(.validate_omics(duplicate_analytes, validated_pheno)),
    "ANALYTE_NAME contains duplicate values",
    "duplicate analyte names are rejected"
  )

  validated_omics <- suppressWarnings(.validate_omics(fixture$omics, validated_pheno))
  .expect_true(
    is.data.frame(validated_omics) &&
      identical(names(validated_omics)[1], "ANALYTE_NAME") &&
      !any(c("all", "male", "female") %in% names(validated_omics)),
    "omics validation returns the ML omics table"
  )

  .expect_error(
    .validate_ml_args(
      models = "xgb",
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
      c(enet_cv_folds = 4L, xgb_cv_folds = 2L)
    ),
    "requested CV folds exceed the smaller treatment arm (3 subjects): enet_cv_folds=4. Use 3 folds or fewer",
    "undersized follow-up cohorts reject too many CV folds"
  )

  large_enough <- data.frame(
    SUBJECT_ID = paste0("S", 1:8),
    TREATMENT_GROUP = factor(c(0, 0, 0, 0, 1, 1, 1, 1), levels = 0:1)
  )
  .expect_true(
    is.null(.validate_followup_cohort(
      large_enough,
      c(enet_cv_folds = 4L, xgb_cv_folds = 3L)
    )),
    "follow-up cohorts with enough subjects pass fold validation"
  )
}

run_feature_tests <- function() {
  cat("\nFeature Preparation\n")
  cat("===================\n")

  fixture <- .make_exact_fixture()
  prepared <- .prepare_fu_change_dataset(
    pheno_df = fixture$pheno,
    omics_df = fixture$omics,
    fu_level = 1L,
    model_covariates = "age",
    enet_cv_folds = 2L,
    xgb_cv_folds = 2L,
    xgb_cv_repeats = 3L,
    seed = 1L
  )

  .expect_equal(prepared$subject_ids_train, paste0("S", 1:8), "training subjects stay aligned")
  .expect_equal(prepared$y_train, rep(0:1, 4L), "training outcomes stay aligned")

  missing_only_in_train <- .drop_all_missing_train(
    x_train = matrix(
      c(NA, NA, 1, 2),
      nrow = 2,
      dimnames = list(NULL, c("all_missing", "retained"))
    )
  )
  .expect_equal(
    colnames(missing_only_in_train$train),
    "retained",
    "all-missing training features are removed before imputation"
  )
  feature_idx <- match("omics::feature_change", colnames(prepared$xgb_x_train))
  expected_train <- as.numeric(scale(seq(1, 15, by = 2)))
  .expect_equal(
    prepared$xgb_x_train[, feature_idx],
    expected_train,
    "follow-up minus baseline is computed and scaled correctly"
  )

  missing_idx <- match("omics::missing_change", colnames(prepared$xgb_x_train))
  .expect_equal(
    prepared$preprocessing$MEDIAN[
      prepared$preprocessing$FEATURE_NAME == "omics::missing_change"
    ],
    8,
    "imputation median comes from training subjects"
  )
  .expect_equal(
    is.finite(prepared$xgb_x_train[1, missing_idx]),
    TRUE,
    "missing training change is median-imputed"
  )
  .expect_true(
    !"omics::constant_change" %in% colnames(prepared$xgb_x_train),
    "training-zero-variance omics features are removed"
  )
  .expect_true(
    !"omics::all_missing_training" %in% colnames(prepared$xgb_x_train),
    "all-missing training omics features are removed"
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

  female_prefix_pheno <- fixture$pheno
  female_prefix_pheno$FEMALE_SCORE <- seq_len(nrow(female_prefix_pheno))
  female_prefix <- .prepare_fu_change_dataset(
    pheno_df = female_prefix_pheno,
    omics_df = fixture$omics,
    fu_level = 1L,
    model_covariates = c("age", "FEMALE_SCORE"),
    enet_cv_folds = 2L,
    xgb_cv_folds = 2L,
    xgb_cv_repeats = 1L,
    seed = 1L
  )
  .expect_true(
    any(startsWith(colnames(female_prefix$enet_x_train), "covariate::FEMALE_SCORE")) &&
      !any(startsWith(colnames(female_prefix$xgb_x_train), "covariate::FEMALE_SCORE")) &&
      identical(
        female_prefix$preprocessing$FEATURE_NAME[
          female_prefix$preprocessing$IN_XGB &
            female_prefix$preprocessing$FEATURE_TYPE == "covariate"
        ],
        "covariate::FEMALE"
      ),
    "only exact FEMALE covariate enters XGB"
  )

  constant_female_pheno <- fixture$pheno
  constant_female_pheno$FEMALE <- factor(0L, levels = 0:1)
  constant_female <- .prepare_fu_change_dataset(
    pheno_df = constant_female_pheno,
    omics_df = fixture$omics,
    fu_level = 1L,
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
    prepared$preprocessing$FEATURE_NAME[
      prepared$preprocessing$STATUS == "retained" & prepared$preprocessing$IN_ENET
    ],
    colnames(prepared$enet_x_train),
    "preprocessing recipe identifies retained ENET features in matrix order"
  )
  .expect_equal(
    prepared$preprocessing$FEATURE_NAME[
      prepared$preprocessing$STATUS == "retained" & prepared$preprocessing$IN_XGB
    ],
    colnames(prepared$xgb_x_train),
    "preprocessing recipe identifies retained XGB features in matrix order"
  )
  .expect_equal(
    prepared$preprocessing$STATUS[
      prepared$preprocessing$FEATURE_NAME == "omics::constant_change"
    ],
    "zero_variance_training",
    "preprocessing audit records zero-variance removals"
  )
  .expect_equal(
    prepared$preprocessing$STATUS[
      prepared$preprocessing$FEATURE_NAME == "omics::all_missing_training"
    ],
    "all_missing_training",
    "preprocessing audit records all-missing removals"
  )
  .expect_true(
    all(is.na(prepared$preprocessing[
      prepared$preprocessing$STATUS == "all_missing_training",
      c("MEDIAN", "CENTER", "SCALE")
    ])),
    "all-missing preprocessing rows have no learned parameters"
  )
  .expect_equal(
    prepared$cohort$SET,
    c("eligible", "train"),
    "cohort report contains eligible and train rows"
  )
  .expect_equal(
    prepared$cohort$N_SUBJECTS,
    c(8L, 8L),
    "cohort report counts aligned subjects by set"
  )
  .expect_equal(
    prepared$cohort[, c("N_CONTROL", "N_TREATMENT", "N_MALE", "N_FEMALE")],
    data.frame(
      N_CONTROL = c(4L, 4L),
      N_TREATMENT = c(4L, 4L),
      N_MALE = c(4L, 4L),
      N_FEMALE = c(4L, 4L)
    ),
    "cohort report counts treatment and sex by set"
  )

  train_control_feature <- prepared$change_summary[
    prepared$change_summary$SET == "train" &
      prepared$change_summary$TREATMENT_GROUP == 0L &
      prepared$change_summary$ANALYTE_NAME == "feature_change",
  ]
  .expect_equal(
    train_control_feature[, c("N_SUBJECTS", "N_NONMISSING", "MEAN", "MEDIAN", "SD", "MIN", "MAX")],
    data.frame(
      N_SUBJECTS = 4L,
      N_NONMISSING = 4L,
      MEAN = 7,
      MEDIAN = 7,
      SD = sqrt(80 / 3),
      MIN = 1,
      MAX = 13
    ),
    "change summary uses raw pre-imputation changes"
  )
  .expect_true(
    {
      all_missing_train_rows <- prepared$change_summary[
        prepared$change_summary$SET == "train" &
          prepared$change_summary$ANALYTE_NAME == "all_missing_training",
      ]
      all(all_missing_train_rows$N_NONMISSING == 0L) &&
        all(is.na(all_missing_train_rows$MEAN))
    },
    "change summary preserves raw missingness"
  )
  .expect_true(
    length(prepared$enet_foldid) == length(prepared$subject_ids_train) &&
      nrow(prepared$xgb_folds) == length(prepared$subject_ids_train) * 3L &&
      identical(sort(unique(prepared$xgb_folds$REPEAT)), 1:3),
    "ENET and repeated XGB folds are prepared separately"
  )

  artifact_dir <- tempfile("track22_prepared_")
  .write_prepared_dataset(prepared, artifact_dir)
  .expect_true(
    all(!file.exists(file.path(
      artifact_dir,
      c("cohort.csv", "change_summary.csv", "preprocessing.csv")
    ))),
    "prepared dataset writer leaves consolidated reports to the run writer"
  )

  missing_visit_pheno <- fixture$pheno[
    !(fixture$pheno$SUBJECT_ID == "S4" & fixture$pheno$FU == 1),
  ]
  prepared_missing <- .prepare_fu_change_dataset(
    pheno_df = missing_visit_pheno,
    omics_df = fixture$omics,
    fu_level = 1L,
    enet_cv_folds = 2L,
    xgb_cv_folds = 2L,
    xgb_cv_repeats = 3L,
    seed = 1L
  )
  .expect_true(
    !"S4" %in% prepared_missing$subject_ids_train &&
      all(c("S1", "S2", "S3", "S5", "S6", "S7", "S8") %in% prepared_missing$subject_ids_train),
    "a missing visit excludes only that subject from the follow-up dataset"
  )
}

run_end_to_end_tests <- function() {
  cat("\nEnd-to-End Pipeline\n")
  cat("===================\n")

  fixture <- .make_simulated_fixture(followups = 1:2)
  output_dir <- file.path(getwd(), "tests", "test_outputs", "track22_integration")
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
    enet_cv_folds = 5L,
    xgb_cv_folds = 5L,
    xgb_cv_repeats = 3L,
    seed = 2202L,
    n_cores = 1L,
    xgb_n_trials = 10L,
    python_bin = .test_python_bin()
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
  .expect_true(
    all(startsWith(artifact_paths, file.path(manifest$data_dir, "FU1"))) &&
      all(startsWith(
        unlist(fu1$models, use.names = FALSE)[
          !grepl("predictions\\.csv$", unlist(fu1$models, use.names = FALSE))
        ],
        file.path(manifest$models_dir, "FU1")
      )) &&
      all(startsWith(
        unlist(fu1$models, use.names = FALSE)[
          grepl("predictions\\.csv$", unlist(fu1$models, use.names = FALSE))
        ],
        file.path(manifest$data_dir, "FU1")
      )),
    "manifest separates model artifacts from subject-level data artifacts"
  )

  subjects <- read.csv(fu1$artifacts$subjects, stringsAsFactors = FALSE)
  xgb_folds <- read.csv(fu1$artifacts$xgb_folds, stringsAsFactors = FALSE)
  enet_train <- read.csv(fu1$artifacts$enet_train, check.names = FALSE)
  xgb_train <- read.csv(fu1$artifacts$xgb_train, check.names = FALSE)
  cohort <- read.csv(manifest$reports$cohort, stringsAsFactors = FALSE)
  change_summary <- read.csv(manifest$reports$change_summary, stringsAsFactors = FALSE)
  preprocessing <- read.csv(manifest$reports$preprocessing, stringsAsFactors = FALSE)

  .expect_true(
    identical(names(subjects), c("SUBJECT_ID", "FU", "TREATMENT_GROUP", "ENET_FOLD_ID")) &&
      nrow(enet_train) == nrow(subjects) &&
      nrow(xgb_train) == nrow(subjects) &&
      all(!is.na(subjects$ENET_FOLD_ID)),
    "canonical artifact row counts agree"
  )
  .expect_true(
    identical(names(xgb_folds), c("SUBJECT_ID", "REPEAT", "FOLD_ID")) &&
      nrow(xgb_folds) == nrow(subjects) * 3L &&
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
      all(table(cohort$FU) == 2L) &&
      all(c("eligible", "train") %in% cohort$SET) &&
      cohort$N_SUBJECTS[cohort$FU == 1L & cohort$SET == "eligible"] == nrow(subjects) &&
      cohort$N_SUBJECTS[cohort$FU == 1L & cohort$SET == "train"] == nrow(subjects),
    "cohort artifact stacks eligible and train subjects across follow-ups"
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
      all(c("eligible", "train") %in% change_summary$SET) &&
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

  for (fu_key in names(manifest$followups)) {
    fu_manifest <- manifest$followups[[fu_key]]
    if (is.null(fu_manifest)) next

    fu_level <- as.integer(sub("^FU", "", fu_key))
    fu_preprocessing <- preprocessing[preprocessing$FU == fu_level, , drop = FALSE]
    enet_columns <- names(read.csv(fu_manifest$artifacts$enet_train, check.names = FALSE))
    xgb_columns <- names(read.csv(fu_manifest$artifacts$xgb_train, check.names = FALSE))

    .expect_equal(
      fu_preprocessing$FEATURE_NAME[
        fu_preprocessing$STATUS == "retained" & fu_preprocessing$IN_ENET
      ],
      enet_columns,
      paste(fu_key, "preprocessing order defines ENET matrix columns")
    )
    .expect_equal(
      fu_preprocessing$FEATURE_NAME[
        fu_preprocessing$STATUS == "retained" & fu_preprocessing$IN_XGB
      ],
      xgb_columns,
      paste(fu_key, "preprocessing order defines XGB matrix columns")
    )
  }

  for (model_name in c("enet", "xgb")) {
    model <- fu1$models[[model_name]]
    .expect_true(
      all(file.exists(unlist(model, use.names = FALSE))),
      paste(model_name, "declared outputs are written")
    )
    metrics <- read.csv(model$metrics)
    .expect_true(
      nrow(metrics) == 1L &&
        all(c("CV_AUC", "INSAMPLE_AUC") %in% names(metrics)) &&
        !"TEST_AUC" %in% names(metrics),
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
  score_matrix <- enet_train
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

  list(fixture = fixture, manifest = manifest)
}

completed_sections <- 0L
result <- tryCatch({
  run_validation_tests()
  completed_sections <- completed_sections + 1L

  run_feature_tests()
  completed_sections <- completed_sections + 1L

  integration <- run_end_to_end_tests()
  completed_sections <- completed_sections + 1L

  run_inference_equivalence_tests(integration$fixture, integration$manifest)
  completed_sections <- completed_sections + 1L

  NULL
}, error = identity)

cat("\nTest summary\n")
cat("============\n")
cat("Sections:         ", completed_sections, "/4\n", sep = "")
cat("Assertions passed:", getOption("track22.tests.passed", 0L), "\n")
cat("Models exercised: ENET, XGB, Optuna, inference replay\n")
cat("Result:           ", if (is.null(result)) "PASS" else "FAIL", "\n", sep = "")

if (!is.null(result)) {
  stop(conditionMessage(result), call. = FALSE)
}
