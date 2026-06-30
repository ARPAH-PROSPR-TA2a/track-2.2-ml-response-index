.as_binary_numeric <- function(x) {
  out <- suppressWarnings(as.integer(as.character(x)))
  if (any(is.na(out)) || !all(out %in% 0:1)) {
    stop("TREATMENT_GROUP must be coercible to 0/1.")
  }
  out
}


.validate_followup_cohort <- function(pheno_df, fold_counts) {
  subject_rows <- pheno_df[!duplicated(pheno_df$SUBJECT_ID), c("SUBJECT_ID", "TREATMENT_GROUP")]
  subject_rows$Y <- .as_binary_numeric(subject_rows$TREATMENT_GROUP)

  class_counts <- table(subject_rows$Y)
  if (!all(c("0", "1") %in% names(class_counts))) {
    return("eligible subjects must include both treatment arms")
  }

  class_counts <- as.integer(class_counts[c("0", "1")])
  min_class_count <- min(class_counts)
  too_many <- fold_counts[fold_counts > min_class_count]
  if (length(too_many) > 0L) {
    return(paste0(
      "requested CV folds exceed the smaller treatment arm (",
      min_class_count, " subjects): ",
      paste(names(too_many), too_many, sep = "=", collapse = ", "),
      ". Use ", min_class_count, " folds or fewer"
    ))
  }

  NULL
}


.stratified_subject_folds <- function(subject_ids, y, cv_folds = 10L, seed = 1L) {
  if (length(subject_ids) != length(y)) {
    stop("subject_ids and y must have the same length.")
  }

  class_counts <- table(y)
  if (!all(c("0", "1") %in% names(class_counts))) {
    stop("CV requires both treatment arms.")
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


.scale_train <- function(x_train) {
  centers <- colMeans(x_train, na.rm = TRUE)
  scales <- apply(x_train, 2, sd, na.rm = TRUE)
  scales[is.na(scales) | scales == 0] <- 1

  x_train_scaled <- sweep(x_train, 2, centers, "-")
  x_train_scaled <- sweep(x_train_scaled, 2, scales, "/")

  list(
    train = x_train_scaled,
    center = centers,
    scale = scales
  )
}


.drop_all_missing_train <- function(x_train) {
  keep <- colSums(!is.na(x_train)) > 0L

  list(
    train = x_train[, keep, drop = FALSE],
    keep = keep
  )
}


.impute_train_median <- function(x_train) {
  medians <- apply(x_train, 2, median, na.rm = TRUE)

  for (j in seq_len(ncol(x_train))) {
    x_train[is.na(x_train[, j]), j] <- medians[j]
  }

  list(train = x_train, medians = medians)
}


.drop_zero_variance_train <- function(x_train) {
  vars <- apply(x_train, 2, var, na.rm = TRUE)
  keep <- !is.na(vars) & vars > 1e-12

  list(
    train = x_train[, keep, drop = FALSE],
    keep = keep
  )
}


.make_preprocessing_table <- function(feature_type, prefix, original_names,
                                      nonmissing_keep, medians, variance_keep,
                                      center, scale, xgb_feature_names,
                                      fu_level) {
  status <- rep("all_missing_training", length(original_names))
  names(status) <- original_names

  nonmissing_names <- names(nonmissing_keep)[nonmissing_keep]
  retained_names <- names(variance_keep)[variance_keep]
  status[nonmissing_names] <- "zero_variance_training"
  status[retained_names] <- "retained"

  median_values <- rep(NA_real_, length(original_names))
  center_values <- rep(NA_real_, length(original_names))
  scale_values <- rep(NA_real_, length(original_names))
  names(median_values) <- names(center_values) <- names(scale_values) <- original_names
  median_values[names(medians)] <- medians
  center_values[names(center)] <- center
  scale_values[names(scale)] <- scale

  feature_names <- paste0(prefix, original_names)
  data.frame(
    FU = fu_level,
    FEATURE_NAME = feature_names,
    FEATURE_TYPE = feature_type,
    STATUS = unname(status),
    MEDIAN = unname(median_values),
    CENTER = unname(center_values),
    SCALE = unname(scale_values),
    IN_ENET = status == "retained",
    IN_XGB = feature_names %in% xgb_feature_names,
    stringsAsFactors = FALSE
  )
}


.make_cohort_report <- function(pheno_followup) {
  treatment <- .as_binary_numeric(pheno_followup$TREATMENT_GROUP)
  female <- .as_binary_numeric(pheno_followup$FEMALE)
  data.frame(
    FU = as.integer(as.character(pheno_followup$FU[1])),
    N_SUBJECTS = nrow(pheno_followup),
    N_CONTROL = sum(treatment == 0L),
    N_TREATMENT = sum(treatment == 1L),
    N_MALE = sum(female == 0L),
    N_FEMALE = sum(female == 1L),
    stringsAsFactors = FALSE
  )
}


.summarize_change_values <- function(values) {
  observed <- values[!is.na(values)]
  if (length(observed) == 0L) {
    return(c(
      N_NONMISSING = 0,
      MEAN = NA_real_,
      MEDIAN = NA_real_,
      SD = NA_real_,
      MIN = NA_real_,
      MAX = NA_real_
    ))
  }

  c(
    N_NONMISSING = length(observed),
    MEAN = mean(observed),
    MEDIAN = median(observed),
    SD = if (length(observed) > 1L) stats::sd(observed) else NA_real_,
    MIN = min(observed),
    MAX = max(observed)
  )
}


.make_change_summary <- function(change_matrix, y, fu_level) {
  rows <- list()
  row_index <- 1L
  for (treatment_group in 0:1) {
    group_idx <- which(y == treatment_group)
    if (length(group_idx) == 0L) next

    summaries <- t(apply(
      change_matrix[group_idx, , drop = FALSE],
      2,
      .summarize_change_values
    ))
    rows[[row_index]] <- data.frame(
      FU = fu_level,
      TREATMENT_GROUP = treatment_group,
      ANALYTE_NAME = colnames(change_matrix),
      N_SUBJECTS = length(group_idx),
      N_NONMISSING = as.integer(summaries[, "N_NONMISSING"]),
      MEAN = summaries[, "MEAN"],
      MEDIAN = summaries[, "MEDIAN"],
      SD = summaries[, "SD"],
      MIN = summaries[, "MIN"],
      MAX = summaries[, "MAX"],
      stringsAsFactors = FALSE
    )
    row_index <- row_index + 1L
  }

  do.call(rbind, rows)
}


.make_followup_change_data <- function(pheno_df, omics_df, fu_level) {
  fu_num <- as.integer(as.character(pheno_df$FU))
  pheno_baseline <- pheno_df[fu_num == 0L, ]
  pheno_followup <- pheno_df[fu_num == fu_level, ]

  complete_subjects <- intersect(pheno_baseline$SUBJECT_ID, pheno_followup$SUBJECT_ID)
  subject_ids <- intersect(pheno_followup$SUBJECT_ID, complete_subjects)

  if (length(subject_ids) == 0L) {
    return(NULL)
  }

  pheno_baseline <- pheno_baseline[match(subject_ids, pheno_baseline$SUBJECT_ID), ]
  pheno_followup <- pheno_followup[match(subject_ids, pheno_followup$SUBJECT_ID), ]

  baseline_values <- as.matrix(omics_df[, pheno_baseline$SAMPLE_ID, drop = FALSE])
  followup_values <- as.matrix(omics_df[, pheno_followup$SAMPLE_ID, drop = FALSE])
  change_matrix <- t(followup_values - baseline_values)
  colnames(change_matrix) <- omics_df$ANALYTE_NAME

  list(
    fu_level = fu_level,
    subject_ids = subject_ids,
    pheno_followup = pheno_followup,
    y = .as_binary_numeric(pheno_followup$TREATMENT_GROUP),
    change_matrix = change_matrix
  )
}


.prepare_covariate_matrices <- function(train_pheno, additional_covariates = NULL) {
  if (is.null(additional_covariates) || length(additional_covariates) == 0L) {
    return(list(
      train = NULL,
      test = NULL,
      feature_names = character(0),
      metadata = NULL
    ))
  }

  train_cov <- train_pheno[, additional_covariates, drop = FALSE]
  combined <- train_cov

  for (col in names(combined)) {
    if (col == "FEMALE") {
      combined[[col]] <- .as_binary_numeric(combined[[col]])
    } else if (is.logical(combined[[col]])) {
      combined[[col]] <- factor(combined[[col]], levels = c(FALSE, TRUE))
    }
  }

  mm <- model.matrix(~ . - 1, data = combined)
  train_mm <- mm[seq_len(nrow(train_cov)), , drop = FALSE]

  nonmissing <- .drop_all_missing_train(train_mm)
  imputed <- .impute_train_median(nonmissing$train)
  dropped <- .drop_zero_variance_train(imputed$train)
  scaled <- .scale_train(dropped$train)
  retained_names <- colnames(scaled$train)
  xgb_feature_names <- paste0(
    "covariate::",
    retained_names[retained_names == "FEMALE"]
  )

  list(
    train = scaled$train,
    feature_names = colnames(scaled$train),
    preprocessing = .make_preprocessing_table(
      feature_type = "covariate",
      prefix = "covariate::",
      original_names = colnames(train_mm),
      nonmissing_keep = nonmissing$keep,
      medians = imputed$medians,
      variance_keep = dropped$keep,
      center = scaled$center,
      scale = scaled$scale,
      xgb_feature_names = xgb_feature_names,
      fu_level = unique(as.integer(as.character(train_pheno$FU)))
    )
  )
}


.prepare_fu_change_dataset <- function(pheno_df, omics_df, fu_level,
                                       model_covariates = "FEMALE",
                                       enet_cv_folds = 10L,
                                       xgb_cv_folds = 10L,
                                       xgb_cv_repeats = 3L,
                                       seed = 1L) {
  change_data <- .make_followup_change_data(pheno_df, omics_df, fu_level)
  if (is.null(change_data)) {
    warning("FU", fu_level, ": no subjects after requiring baseline and follow-up.")
    return(NULL)
  }

  train_subjects <- change_data$subject_ids
  y <- change_data$y
  train_idx <- seq_along(train_subjects)

  if (length(unique(y[train_idx])) < 2L) {
    warning("FU", fu_level, ": training subjects must contain both treatment arms.")
    return(NULL)
  }

  omics_train <- change_data$change_matrix[train_idx, , drop = FALSE]
  nonmissing <- .drop_all_missing_train(omics_train)
  imputed <- .impute_train_median(nonmissing$train)
  dropped <- .drop_zero_variance_train(imputed$train)
  scaled <- .scale_train(dropped$train)
  retained_omics <- paste0("omics::", colnames(scaled$train))
  omics_preprocessing <- .make_preprocessing_table(
    feature_type = "omics",
    prefix = "omics::",
    original_names = colnames(omics_train),
    nonmissing_keep = nonmissing$keep,
    medians = imputed$medians,
    variance_keep = dropped$keep,
    center = scaled$center,
    scale = scaled$scale,
    xgb_feature_names = retained_omics,
    fu_level = fu_level
  )
  colnames(scaled$train) <- paste0("omics::", colnames(scaled$train))

  train_pheno <- change_data$pheno_followup[train_idx, , drop = FALSE]
  model_covariates <- unique(c("FEMALE", model_covariates))
  covariates <- .prepare_covariate_matrices(train_pheno, model_covariates)
  if (!is.null(covariates$train)) {
    colnames(covariates$train) <- paste0("covariate::", colnames(covariates$train))
  }

  enet_train <- scaled$train

  if (!is.null(covariates$train)) {
    enet_train <- cbind(enet_train, covariates$train)
  }

  female_columns <- colnames(covariates$train) == "covariate::FEMALE"
  xgb_covariates <- list(
    train = covariates$train[, female_columns, drop = FALSE]
  )
  xgb_train <- scaled$train
  if (ncol(xgb_covariates$train) > 0L) {
    xgb_train <- cbind(xgb_train, xgb_covariates$train)
  }

  enet_foldid <- .stratified_subject_folds(
    train_subjects,
    y[train_idx],
    cv_folds = enet_cv_folds,
    seed = seed
  )
  if (is.null(enet_foldid)) {
    warning("FU", fu_level, ": unable to create stratified ENET CV folds.")
    return(NULL)
  }

  xgb_fold_rows <- lapply(seq_len(xgb_cv_repeats), function(repeat_id) {
    foldid <- .stratified_subject_folds(
      train_subjects,
      y[train_idx],
      cv_folds = xgb_cv_folds,
      seed = seed + repeat_id
    )
    if (is.null(foldid)) return(NULL)
    data.frame(
      SUBJECT_ID = train_subjects,
      REPEAT = repeat_id,
      FOLD_ID = foldid,
      stringsAsFactors = FALSE
    )
  })
  if (any(vapply(xgb_fold_rows, is.null, logical(1)))) {
    warning("FU", fu_level, ": unable to create repeated stratified XGB CV folds.")
    return(NULL)
  }

  list(
    fu_level = fu_level,
    subject_ids_train = train_subjects,
    y_train = y[train_idx],
    enet_foldid = enet_foldid,
    xgb_folds = do.call(rbind, xgb_fold_rows),
    enet_x_train = enet_train,
    xgb_x_train = xgb_train,
    cohort = .make_cohort_report(change_data$pheno_followup),
    change_summary = .make_change_summary(
      change_data$change_matrix,
      y,
      fu_level
    ),
    preprocessing = rbind(omics_preprocessing, covariates$preprocessing)
  )
}
