source("main.R")


make_demo_trial <- function(n_subjects = 250L, n_features = 1000L,
                            followups = 1:2, seed = 220200L) {
  set.seed(seed)

  subjects <- sprintf("DEMO%04d", seq_len(n_subjects))
  treatment <- sample(rep(0:1, length.out = n_subjects))
  female <- rbinom(n_subjects, size = 1L, prob = 0.55)
  age <- round(rnorm(n_subjects, mean = 48, sd = 11), 1)
  bmi <- round(rnorm(n_subjects, mean = 27, sd = 4.5), 1)
  site <- factor(sample(c("site_a", "site_b", "site_c", "site_d"), n_subjects, replace = TRUE))
  smoker <- sample(c(FALSE, TRUE), n_subjects, replace = TRUE, prob = c(0.78, 0.22))

  fu_values <- c(0L, followups)
  pheno <- do.call(rbind, lapply(fu_values, function(fu) {
    data.frame(
      SAMPLE_ID = sprintf("%s_FU%d", subjects, fu),
      SUBJECT_ID = subjects,
      FU = factor(fu, levels = fu_values),
      TREATMENT_GROUP = factor(treatment, levels = 0:1),
      FEMALE = factor(female, levels = 0:1),
      age = age,
      bmi = bmi,
      site = site,
      smoker = smoker,
      stringsAsFactors = FALSE
    )
  }))

  n_samples <- nrow(pheno)
  sample_ids <- pheno$SAMPLE_ID
  feature_names <- sprintf("protein_%04d", seq_len(n_features))

  subject_baseline <- matrix(rnorm(n_features * n_subjects, sd = 1.0), nrow = n_features)
  values <- matrix(NA_real_, nrow = n_features, ncol = n_samples)

  for (fu in fu_values) {
    sample_idx <- which(pheno$FU == fu)
    if (fu == 0L) {
      values[, sample_idx] <- subject_baseline
      next
    }

    change <- matrix(rnorm(n_features * n_subjects, sd = 0.85), nrow = n_features)
    signal_block <- seq_len(8L)
    female_block <- 9:12
    change[signal_block, treatment == 1L] <-
      change[signal_block, treatment == 1L] + if (fu == 1L) 0.45 else 0.65
    change[female_block, female == 1L] <-
      change[female_block, female == 1L] + if (fu == 1L) 0.2 else 0.3
    values[, sample_idx] <- subject_baseline + change
  }

  missing_mask <- matrix(runif(length(values)) < 0.0002, nrow = nrow(values))
  values[missing_mask] <- NA_real_

  colnames(values) <- sample_ids
  omics <- data.frame(
    ANALYTE_NAME = feature_names,
    values,
    check.names = FALSE
  )

  list(pheno = pheno, omics = omics)
}


print_metrics <- function(output_dir) {
  manifest <- jsonlite::fromJSON(file.path(output_dir, "manifest.json"), simplifyVector = FALSE)
  for (fu in names(manifest$followups)) {
    cat("\n", fu, " metric preview\n", sep = "")
    cat(strrep("-", nchar(fu) + 15L), "\n", sep = "")
    for (model in c("enet", "xgb")) {
      metrics_path <- file.path(output_dir, fu, model, "metrics.csv")
      if (!file.exists(metrics_path)) next
      metrics <- read.csv(metrics_path, stringsAsFactors = FALSE)
      preview <- metrics[, intersect(c("CV_AUC", "TEST_AUC", "INSAMPLE_AUC", "N_FEATURES", "BEST_ITERATION"), names(metrics)), drop = FALSE]
      cat("\n", toupper(model), "\n", sep = "")
      print(preview, row.names = FALSE)
    }
  }
}


cat("Track 2.2 demo-scale smoke run\n")
cat("==============================\n\n")

output_dir <- file.path(getwd(), "test_outputs", "track22_demo")
if (dir.exists(output_dir)) {
  cat("Removing prior demo output: ", output_dir, "\n", sep = "")
  unlink(output_dir, recursive = TRUE)
}

cat("Generating synthetic trial: 250 subjects, 1000 features, FU1/FU2\n")
demo <- make_demo_trial()
cat("Phenotype rows: ", nrow(demo$pheno), "\n", sep = "")
cat("Omics features: ", nrow(demo$omics), "\n", sep = "")
cat("Output directory: ", output_dir, "\n\n", sep = "")

cat("Starting FAST_treatment_ML with ENET + XGB and 10 Optuna trials...\n")
manifest <- FAST_treatment_ML(
  pheno = demo$pheno,
  omics = demo$omics,
  omics_type = "Proteomics",
  additional_covariates = c("age", "bmi", "site", "smoker"),
  models = c("enet", "xgb"),
  output_dir = output_dir,
  test_frac = 0.2,
  enet_cv_folds = 5L,
  xgb_cv_folds = 5L,
  xgb_cv_repeats = 3L,
  seed = 220200L,
  n_cores = 2L,
  xgb_n_trials = 10L
)

cat("\nModel run complete.\n")
cat("Manifest: ", manifest$manifest_path, "\n", sep = "")
print_metrics(output_dir)

cat("\nRunning reports...\n")
reports <- FAST_treatment_ML_reports(
  pheno = demo$pheno,
  omics = demo$omics,
  omics_type = "Proteomics",
  additional_covariates = c("age", "bmi", "site", "smoker")
)

cat("Report objects:\n")
cat("  pheno_summary rows: ", nrow(reports$pheno_summary), "\n", sep = "")
cat("  variable_summaries strata: ", paste(names(reports$variable_summaries), collapse = ", "), "\n", sep = "")
cat("  randomization reports: ", paste(names(reports$randomization_reports), collapse = ", "), "\n", sep = "")

cat("\nDemo output is available for inspection at:\n")
cat(output_dir, "\n")
