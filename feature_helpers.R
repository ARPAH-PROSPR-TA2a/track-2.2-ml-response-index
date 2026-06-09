.as_binary_numeric <- function(x) {
  out <- suppressWarnings(as.integer(as.character(x)))
  if (any(is.na(out)) || !all(out %in% 0:1)) {
    stop("TREATMENT_GROUP must be coercible to 0/1.")
  }
  out
}


.stratified_subject_split <- function(pheno_df, test_frac = 0.2, seed = 1L,
                                      min_per_class = 2L) {
  subject_rows <- pheno_df[!duplicated(pheno_df$SUBJECT_ID), c("SUBJECT_ID", "TREATMENT_GROUP")]
  subject_rows$Y <- .as_binary_numeric(subject_rows$TREATMENT_GROUP)

  class_counts <- table(subject_rows$Y)
  if (!all(c("0", "1") %in% names(class_counts)) ||
      any(class_counts[c("0", "1")] < min_per_class)) {
    warning("Skipping split: not enough subjects in both treatment arms.")
    return(NULL)
  }

  set.seed(seed)
  test_subjects <- character(0)
  for (class_value in c(0L, 1L)) {
    ids <- subject_rows$SUBJECT_ID[subject_rows$Y == class_value]
    n_test <- max(1L, floor(length(ids) * test_frac))
    n_test <- min(n_test, length(ids) - 1L)
    test_subjects <- c(test_subjects, sample(ids, n_test))
  }

  train_subjects <- setdiff(subject_rows$SUBJECT_ID, test_subjects)

  list(
    train_subjects = train_subjects,
    test_subjects = test_subjects
  )
}


.stratified_subject_folds <- function(subject_ids, y, cv_folds = 5L, seed = 1L) {
  if (length(subject_ids) != length(y)) {
    stop("subject_ids and y must have the same length.")
  }

  class_counts <- table(y)
  if (!all(c("0", "1") %in% names(class_counts))) {
    stop("CV requires both treatment arms.")
  }

  max_folds <- min(as.integer(class_counts[c("0", "1")]))
  cv_folds <- min(as.integer(cv_folds), max_folds)
  if (cv_folds < 2L) {
    warning("Skipping CV: fewer than two subjects in at least one treatment arm.")
    return(NULL)
  }

  set.seed(seed)
  foldid <- integer(length(subject_ids))
  for (class_value in c(0L, 1L)) {
    idx <- which(y == class_value)
    idx <- sample(idx)
    foldid[idx] <- rep(seq_len(cv_folds), length.out = length(idx))
  }

  foldid
}


.safe_auc <- function(y, pred) {
  if (length(unique(y[!is.na(y)])) < 2L || length(unique(pred[!is.na(pred)])) < 2L) {
    return(NA_real_)
  }

  as.numeric(pROC::auc(pROC::roc(y, pred, quiet = TRUE, direction = "<")))
}


.scale_train_test <- function(x_train, x_test) {
  centers <- colMeans(x_train, na.rm = TRUE)
  scales <- apply(x_train, 2, sd, na.rm = TRUE)
  scales[is.na(scales) | scales == 0] <- 1

  x_train_scaled <- sweep(x_train, 2, centers, "-")
  x_train_scaled <- sweep(x_train_scaled, 2, scales, "/")
  x_test_scaled <- sweep(x_test, 2, centers, "-")
  x_test_scaled <- sweep(x_test_scaled, 2, scales, "/")

  list(
    train = x_train_scaled,
    test = x_test_scaled,
    center = centers,
    scale = scales
  )
}


.impute_train_median <- function(x_train, x_test) {
  medians <- apply(x_train, 2, median, na.rm = TRUE)
  all_na <- is.na(medians)
  medians[all_na] <- 0

  for (j in seq_len(ncol(x_train))) {
    x_train[is.na(x_train[, j]), j] <- medians[j]
    x_test[is.na(x_test[, j]), j] <- medians[j]
  }

  list(train = x_train, test = x_test, medians = medians, all_na = all_na)
}


.drop_zero_variance_train <- function(x_train, x_test) {
  vars <- apply(x_train, 2, var, na.rm = TRUE)
  keep <- !is.na(vars) & vars > 1e-12

  list(
    train = x_train[, keep, drop = FALSE],
    test = x_test[, keep, drop = FALSE],
    keep = keep
  )
}


.make_preprocessing_table <- function(feature_type, prefix, medians, keep, center, scale) {
  retained_names <- names(keep)[keep]
  data.frame(
    FEATURE_NAME = paste0(prefix, retained_names),
    FEATURE_TYPE = feature_type,
    MEDIAN = as.numeric(medians[retained_names]),
    CENTER = as.numeric(center[retained_names]),
    SCALE = as.numeric(scale[retained_names]),
    stringsAsFactors = FALSE
  )
}


.prepare_covariate_matrices <- function(train_pheno, test_pheno, additional_covariates = NULL) {
  if (is.null(additional_covariates) || length(additional_covariates) == 0L) {
    return(list(
      train = NULL,
      test = NULL,
      feature_names = character(0),
      metadata = NULL
    ))
  }

  train_cov <- train_pheno[, additional_covariates, drop = FALSE]
  test_cov <- test_pheno[, additional_covariates, drop = FALSE]
  combined <- rbind(train_cov, test_cov)

  for (col in names(combined)) {
    if (col == "FEMALE") {
      combined[[col]] <- .as_binary_numeric(combined[[col]])
    } else if (is.logical(combined[[col]])) {
      combined[[col]] <- factor(combined[[col]], levels = c(FALSE, TRUE))
    }
  }

  mm <- model.matrix(~ . - 1, data = combined)
  train_mm <- mm[seq_len(nrow(train_cov)), , drop = FALSE]
  test_mm <- mm[(nrow(train_cov) + 1L):nrow(mm), , drop = FALSE]

  imputed <- .impute_train_median(train_mm, test_mm)
  dropped <- .drop_zero_variance_train(imputed$train, imputed$test)
  scaled <- .scale_train_test(dropped$train, dropped$test)

  list(
    train = scaled$train,
    test = scaled$test,
    feature_names = colnames(scaled$train),
    preprocessing = .make_preprocessing_table(
      feature_type = "covariate",
      prefix = "covariate::",
      medians = imputed$medians,
      keep = dropped$keep,
      center = scaled$center,
      scale = scaled$scale
    )
  )
}


.prepare_fu_change_dataset <- function(pheno_df, omics_df, fu_level, split,
                                       additional_covariates = NULL,
                                       cv_folds = 5L,
                                       seed = 1L) {
  fu_num <- as.integer(as.character(pheno_df$FU))
  pheno_baseline <- pheno_df[fu_num == 0L, ]
  pheno_followup <- pheno_df[fu_num == fu_level, ]

  complete_subjects <- intersect(pheno_baseline$SUBJECT_ID, pheno_followup$SUBJECT_ID)
  complete_subjects <- intersect(complete_subjects, c(split$train_subjects, split$test_subjects))

  train_subjects <- intersect(split$train_subjects, complete_subjects)
  test_subjects <- intersect(split$test_subjects, complete_subjects)

  if (length(train_subjects) == 0L || length(test_subjects) == 0L) {
    warning("FU", fu_level, ": no train/test subjects after requiring baseline and follow-up.")
    return(NULL)
  }

  ordered_subjects <- c(train_subjects, test_subjects)
  pheno_baseline <- pheno_baseline[match(ordered_subjects, pheno_baseline$SUBJECT_ID), ]
  pheno_followup <- pheno_followup[match(ordered_subjects, pheno_followup$SUBJECT_ID), ]

  y <- .as_binary_numeric(pheno_followup$TREATMENT_GROUP)
  train_idx <- seq_along(train_subjects)
  test_idx <- seq.int(length(train_subjects) + 1L, length(ordered_subjects))

  if (length(unique(y[train_idx])) < 2L || length(unique(y[test_idx])) < 2L) {
    warning("FU", fu_level, ": train and test sets must each contain both treatment arms.")
    return(NULL)
  }

  baseline_values <- as.matrix(omics_df[, pheno_baseline$SAMPLE_ID, drop = FALSE])
  followup_values <- as.matrix(omics_df[, pheno_followup$SAMPLE_ID, drop = FALSE])
  change_matrix <- t(followup_values - baseline_values)
  colnames(change_matrix) <- omics_df$ANALYTE_NAME

  omics_train <- change_matrix[train_idx, , drop = FALSE]
  omics_test <- change_matrix[test_idx, , drop = FALSE]
  imputed <- .impute_train_median(omics_train, omics_test)
  dropped <- .drop_zero_variance_train(imputed$train, imputed$test)
  scaled <- .scale_train_test(dropped$train, dropped$test)
  omics_preprocessing <- .make_preprocessing_table(
    feature_type = "omics",
    prefix = "omics::",
    medians = imputed$medians,
    keep = dropped$keep,
    center = scaled$center,
    scale = scaled$scale
  )
  colnames(scaled$train) <- paste0("omics::", colnames(scaled$train))
  colnames(scaled$test) <- colnames(scaled$train)

  train_pheno <- pheno_followup[train_idx, , drop = FALSE]
  test_pheno <- pheno_followup[test_idx, , drop = FALSE]
  covariates <- .prepare_covariate_matrices(train_pheno, test_pheno, additional_covariates)
  if (!is.null(covariates$train)) {
    colnames(covariates$train) <- paste0("covariate::", colnames(covariates$train))
    colnames(covariates$test) <- colnames(covariates$train)
  }

  enet_train <- scaled$train
  enet_test <- scaled$test

  if (!is.null(covariates$train)) {
    enet_train <- cbind(enet_train, covariates$train)
    enet_test <- cbind(enet_test, covariates$test)
  }

  xgb_covariates <- NULL
  if ("FEMALE" %in% additional_covariates) {
    female_columns <- startsWith(colnames(covariates$train), "covariate::FEMALE")
    xgb_covariates <- list(
      train = covariates$train[, female_columns, drop = FALSE],
      test = covariates$test[, female_columns, drop = FALSE]
    )
  }
  xgb_train <- scaled$train
  xgb_test <- scaled$test
  if (!is.null(xgb_covariates)) {
    xgb_train <- cbind(xgb_train, xgb_covariates$train)
    xgb_test <- cbind(xgb_test, xgb_covariates$test)
  }

  foldid <- .stratified_subject_folds(train_subjects, y[train_idx], cv_folds = cv_folds, seed = seed)
  if (is.null(foldid)) {
    warning("FU", fu_level, ": unable to create stratified CV folds.")
    return(NULL)
  }

  list(
    fu_level = fu_level,
    subject_ids_train = train_subjects,
    subject_ids_test = test_subjects,
    y_train = y[train_idx],
    y_test = y[test_idx],
    foldid = foldid,
    enet_x_train = enet_train,
    enet_x_test = enet_test,
    xgb_x_train = xgb_train,
    xgb_x_test = xgb_test,
    preprocessing = rbind(omics_preprocessing, covariates$preprocessing)
  )
}
