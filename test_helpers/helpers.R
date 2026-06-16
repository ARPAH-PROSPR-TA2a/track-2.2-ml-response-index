.record_pass <- function() {
  current <- getOption("track22.tests.passed", 0L)
  options(track22.tests.passed = current + 1L)
}


.expect_true <- function(value, label) {
  if (!isTRUE(value)) {
    stop("FAIL: ", label, call. = FALSE)
  }
  .record_pass()
  cat("PASS ", label, "\n", sep = "")
}


.expect_equal <- function(actual, expected, label, tolerance = 1e-8) {
  equal <- isTRUE(all.equal(actual, expected, tolerance = tolerance, check.attributes = FALSE))
  if (!equal) {
    stop(
      "FAIL: ", label, "\nExpected: ", paste(expected, collapse = ", "),
      "\nActual: ", paste(actual, collapse = ", "),
      call. = FALSE
    )
  }
  .record_pass()
  cat("PASS ", label, "\n", sep = "")
}


.expect_error <- function(expr, pattern, label) {
  message <- tryCatch(
    {
      force(expr)
      NULL
    },
    error = function(e) conditionMessage(e)
  )
  if (is.null(message) || !grepl(pattern, message, fixed = TRUE)) {
    stop(
      "FAIL: ", label, "\nExpected error containing: ", pattern,
      "\nActual: ", if (is.null(message)) "<no error>" else message,
      call. = FALSE
    )
  }
  .record_pass()
  cat("PASS ", label, "\n", sep = "")
}


.make_exact_fixture <- function() {
  subjects <- paste0("S", 1:8)
  treatment <- rep(0:1, 4L)
  pheno <- do.call(rbind, lapply(seq_along(subjects), function(i) {
    data.frame(
      SAMPLE_ID = paste0(subjects[i], c("_B", "_F")),
      SUBJECT_ID = subjects[i],
      FU = factor(c(0L, 1L), levels = 0:1),
      TREATMENT_GROUP = factor(rep(treatment[i], 2L), levels = 0:1),
      FEMALE = factor(rep(i %% 2L, 2L), levels = 0:1),
      age = rep(seq(20, 90, by = 10)[i], 2L),
      stringsAsFactors = FALSE
    )
  }))

  feature_change <- c(1, 3, 5, 7, 9, 11, 13, 15)
  missing_change <- c(NA, 2, 4, 6, 100, 8, 10, 12)
  constant_change <- c(5, 5, 5, 5, 20, 20, 20, 20)
  all_missing_training <- c(NA, NA, NA, NA, 30, 31, 32, 33)
  sample_ids <- pheno$SAMPLE_ID

  make_values <- function(changes) {
    values <- numeric(length(sample_ids))
    values[seq(1, length(values), by = 2)] <- 10
    values[seq(2, length(values), by = 2)] <- 10 + changes
    values
  }

  omics <- data.frame(
    ANALYTE_NAME = c(
      "feature_change",
      "missing_change",
      "constant_change",
      "all_missing_training"
    ),
    stringsAsFactors = FALSE
  )
  omics[sample_ids] <- rbind(
    make_values(feature_change),
    make_values(missing_change),
    make_values(constant_change),
    make_values(all_missing_training)
  )

  list(
    pheno = pheno,
    omics = omics,
    split = list(
      train_subjects = subjects[1:4],
      test_subjects = subjects[5:6]
    )
  )
}


.make_simulated_fixture <- function(n_subjects = 36L, n_features = 24L,
                                    followups = 1L, seed = 2202L) {
  set.seed(seed)
  subjects <- sprintf("SUBJ%03d", seq_len(n_subjects))
  treatment <- rep(0:1, length.out = n_subjects)
  female <- rep(c(0L, 0L, 1L, 1L), length.out = n_subjects)
  baseline_ids <- paste0(subjects, "_FU0")

  make_pheno <- function(fu, sample_ids) {
    data.frame(
      SAMPLE_ID = sample_ids,
      SUBJECT_ID = subjects,
      FU = fu,
      TREATMENT_GROUP = treatment,
      FEMALE = female,
      age = 35 + seq_len(n_subjects) / 3,
      bmi = 22 + (seq_len(n_subjects) %% 7) / 2,
      site = factor(rep(c("A", "B", "C"), length.out = n_subjects)),
      smoker = rep(c(FALSE, TRUE, FALSE), length.out = n_subjects)
    )
  }
  followup_ids <- lapply(followups, function(fu) paste0(subjects, "_FU", fu))
  names(followup_ids) <- paste0("FU", followups)

  pheno <- do.call(rbind, c(
    list(make_pheno(0L, baseline_ids)),
    Map(make_pheno, followups, followup_ids)
  ))

  baseline <- matrix(rnorm(n_features * n_subjects), nrow = n_features)
  baseline[4, 2] <- NA_real_
  followup_values <- lapply(seq_along(followups), function(i) {
    change <- matrix(rnorm(n_features * n_subjects, sd = 0.7), nrow = n_features)
    change[1:3, treatment == 1L] <- change[1:3, treatment == 1L] + 2.5 + 0.25 * i
    change[n_features, ] <- 0
    out <- baseline + change
    out[4, 2] <- NA_real_
    out
  })

  omics_values <- do.call(cbind, c(list(baseline), followup_values))
  colnames(omics_values) <- c(baseline_ids, unlist(followup_ids, use.names = FALSE))
  omics <- data.frame(
    ANALYTE_NAME = paste0("protein_", seq_len(n_features)),
    omics_values,
    check.names = FALSE
  )

  list(pheno = pheno, omics = omics)
}


.assert_probability_file <- function(path, expected_subjects) {
  predictions <- read.csv(path, stringsAsFactors = FALSE)
  .expect_true(
    identical(
      names(predictions),
      c("SET", "SUBJECT_ID", "FU", "TREATMENT_GROUP", "PREDICTED_PROB")
    ),
    paste(basename(dirname(path)), "prediction schema")
  )
  .expect_true(
    nrow(predictions) == expected_subjects &&
      all(c("train", "test") %in% predictions$SET),
    paste(basename(dirname(path)), "prediction rows")
  )
  .expect_true(
    all(is.finite(predictions$PREDICTED_PROB)) &&
      all(predictions$PREDICTED_PROB >= 0 & predictions$PREDICTED_PROB <= 1),
    paste(basename(dirname(path)), "probabilities are valid")
  )
}
