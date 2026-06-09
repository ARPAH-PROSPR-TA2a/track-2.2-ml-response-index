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

too_small <- data.frame(
  SUBJECT_ID = paste0("S", 1:6),
  TREATMENT_GROUP = factor(c(0, 0, 0, 1, 1, 1), levels = 0:1)
)
.expect_true(
  is.null(suppressWarnings(.stratified_subject_split(
    too_small,
    test_frac = 0.2,
    seed = 1L,
    min_per_class = 5L
  ))),
  "undersized treatment arms are skipped"
)
