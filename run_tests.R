args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) == 1L) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg)))
  setwd(script_dir)
}

cat("Track 2.2 test suite\n")
cat("====================\n")

options(track22.tests.passed = 0L)
source("main.R")

required_r <- c("glmnet", "jsonlite", "pROC")
missing_r <- required_r[!vapply(required_r, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_r) > 0L) {
  stop("Missing required R packages: ", paste(missing_r, collapse = ", "), call. = FALSE)
}

test_files <- c(
  "test_feature_correctness.R",
  "test_validation.R",
  "test_end_to_end.R"
)
completed_files <- 0L

result <- tryCatch({
  cat("PASS R dependencies are ready\n")

  for (test_file in test_files) {
    source(test_file, local = new.env(parent = globalenv()))
    completed_files <- completed_files + 1L
  }

  NULL
}, error = identity)

cat("\nTest summary\n")
cat("============\n")
cat("Test files:       ", completed_files, "/", length(test_files), "\n", sep = "")
cat("Assertions passed:", getOption("track22.tests.passed", 0L), "\n")
cat("Models exercised: ENET, XGB, Optuna\n")
cat("Result:           ", if (is.null(result)) "PASS" else "FAIL", "\n", sep = "")

if (!is.null(result)) {
  stop(conditionMessage(result), call. = FALSE)
}
