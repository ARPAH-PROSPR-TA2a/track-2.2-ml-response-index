.create_pheno_data_report <- function(pheno_df) {

  # One row per (FU, FEMALE) cell. Subject-level counts (N_SUBJECTS,
  # N_CONTROL, N_TREATMENT) dedupe by SUBJECT_ID so multiple samples for
  # the same person at the same FU don't inflate the totals. N_SAMPLES
  # is row-level so the difference between people present and samples
  # collected stays visible.
  groups <- unique(pheno_df[, c("FU", "FEMALE")])
  groups <- groups[order(groups$FU, groups$FEMALE), , drop = FALSE]

  report <- data.frame(
    FU          = groups$FU,
    FEMALE      = groups$FEMALE,
    N_SUBJECTS  = NA_integer_,
    N_CONTROL   = NA_integer_,
    N_TREATMENT = NA_integer_,
    N_SAMPLES   = NA_integer_,
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(report))) {
    cell <- pheno_df[pheno_df$FU == report$FU[i] &
                       pheno_df$FEMALE == report$FEMALE[i], ]
    report$N_SUBJECTS[i]  <- length(unique(cell$SUBJECT_ID))
    report$N_CONTROL[i]   <- length(unique(cell$SUBJECT_ID[cell$TREATMENT_GROUP == 0]))
    report$N_TREATMENT[i] <- length(unique(cell$SUBJECT_ID[cell$TREATMENT_GROUP == 1]))
    report$N_SAMPLES[i]   <- nrow(cell)
  }

  rownames(report) <- NULL
  report
}


.create_omics_data_report <- function(sample_ids, omics_df) {

  analyte_names <- omics_df$ANALYTE_NAME

  omics_numeric <- omics_df[, setdiff(names(omics_df), "ANALYTE_NAME"), drop = FALSE]
  omics_numeric <- omics_numeric[, colnames(omics_numeric) %in% sample_ids, drop = FALSE]
  
  # Initialize results data.frame
  report <- data.frame(
    ANALYTE_NAME = analyte_names,
    N_NONMISSING = NA_integer_,
    MEAN = NA_real_,
    MEDIAN = NA_real_,
    SD = NA_real_,
    MIN = NA_real_,
    MAX = NA_real_,
    stringsAsFactors = FALSE
  )
  
  # Calculate per-analyte statistics
  for (i in seq_along(analyte_names)) {
    analyte_values <- as.numeric(omics_numeric[i, ])
    
    report$N_NONMISSING[i] <- sum(!is.na(analyte_values))
    report$MEAN[i] <- mean(analyte_values, na.rm = TRUE)
    report$MEDIAN[i] <- median(analyte_values, na.rm = TRUE)
    report$SD[i] <- sd(analyte_values, na.rm = TRUE)
    report$MIN[i] <- min(analyte_values, na.rm = TRUE)
    report$MAX[i] <- max(analyte_values, na.rm = TRUE)
  }
  
  return(report)
}


.create_addx_covariate_report <- function(pheno_df, covariate_names) {

   # Handle NULL or empty covariate_names
   if (is.null(covariate_names) || length(covariate_names) == 0) {
     return(NULL)
   }

   # Initialize results list to build then convert to data.frame
   results_list <- list(
     COVARIATE_NAME = character(),
     TYPE = character(),
     N_NA = integer(),
     SUMMARY = list()
   )
   
   # Process each covariate
   for (i in seq_along(covariate_names)) {
     covar_name <- covariate_names[i]
     covar_data <- pheno_df[[covar_name]]
     
     # Count NAs
     n_na <- sum(is.na(covar_data))
     
     # Determine type and generate summary
     if (is.numeric(covar_data)) {
       covar_type <- "numeric"
       summary_stats <- list(
         mean = mean(covar_data, na.rm = TRUE),
         median = median(covar_data, na.rm = TRUE),
         sd = sd(covar_data, na.rm = TRUE),
         min = min(covar_data, na.rm = TRUE),
         max = max(covar_data, na.rm = TRUE)
       )
     } else if (is.factor(covar_data)) {
       covar_type <- "factor"
       level_counts <- table(covar_data, useNA = "no")
       summary_stats <- list(
         n_levels = nlevels(covar_data),
         level_names = levels(covar_data),
         counts = as.numeric(level_counts)
       )
     } else if (is.logical(covar_data)) {
       covar_type <- "logical"
       summary_stats <- list(
         n_true = sum(covar_data == TRUE, na.rm = TRUE),
         n_false = sum(covar_data == FALSE, na.rm = TRUE)
       )
     } else {
       # Fallback for unexpected types
       covar_type <- "unknown"
       summary_stats <- list()
     }
     
     # Add to results list
     results_list$COVARIATE_NAME <- c(results_list$COVARIATE_NAME, covar_name)
     results_list$TYPE <- c(results_list$TYPE, covar_type)
     results_list$N_NA <- c(results_list$N_NA, n_na)
     results_list$SUMMARY[[i]] <- summary_stats
   }
   
   # Convert to data.frame
   report <- data.frame(
     COVARIATE_NAME = results_list$COVARIATE_NAME,
     TYPE = results_list$TYPE,
     N_NA = results_list$N_NA,
     stringsAsFactors = FALSE
   )
   
   # Add SUMMARY as a list column
   report$SUMMARY <- results_list$SUMMARY

   return(report)
}


.create_pheno_randomization_report <- function(pheno_df, additional_covariates = NULL) {

  # Restrict to baseline
  pheno_baseline <- pheno_df[pheno_df$FU == 0, ]

  control_status <- pheno_baseline$TREATMENT_GROUP
  ctrl_mask <- as.integer(as.character(control_status)) == 0
  trt_mask  <- as.integer(as.character(control_status)) == 1

  # FEMALE only if both sexes present; then additional covariates
  variables <- character(0)
  if (length(unique(pheno_baseline$FEMALE)) == 2) {
    variables <- c(variables, "FEMALE")
  }
  if (!is.null(additional_covariates)) {
    variables <- c(variables, additional_covariates)
  }

  if (length(variables) == 0) return(NULL)

  results_list <- list(
    VARIABLE          = character(),
    TYPE              = character(),
    TEST              = character(),
    STATISTIC         = numeric(),
    MIN_CELL_COUNT    = integer(),
    P_VALUE           = numeric(),
    SUMMARY_CONTROL   = list(),
    SUMMARY_TREATMENT = list()
  )

  for (var in variables) {
    col_data  <- pheno_baseline[[var]]
    ctrl_data <- col_data[ctrl_mask]
    trt_data  <- col_data[trt_mask]

    if (is.numeric(col_data)) {
      var_type  <- "numeric"
      test_name <- "t-test"
      tt <- tryCatch(
        t.test(trt_data, ctrl_data, var.equal = FALSE),
        error = function(e) NULL
      )
      statistic      <- if (!is.null(tt)) as.numeric(tt$statistic) else NA_real_
      p_value        <- if (!is.null(tt)) tt$p.value else NA_real_
      min_cell_count <- NA_integer_
      summary_ctrl   <- list(mean = mean(ctrl_data, na.rm = TRUE), sd = sd(ctrl_data, na.rm = TRUE))
      summary_trt    <- list(mean = mean(trt_data,  na.rm = TRUE), sd = sd(trt_data,  na.rm = TRUE))

    } else if (is.logical(col_data) || is.factor(col_data)) {
      var_type       <- if (is.logical(col_data)) "logical" else "factor"
      tab            <- table(col_data, control_status)
      min_cell_count <- as.integer(min(tab))

      # For tables larger than 2x2, simulate p-value in Fisher's exact
      if (min_cell_count < 5) {
        ft <- tryCatch(
          fisher.test(tab, simulate.p.value = nrow(tab) > 2),
          error = function(e) NULL
        )
        test_name <- "fisher"
        statistic <- NA_real_
        p_value   <- if (!is.null(ft)) ft$p.value else NA_real_
      } else {
        ct <- tryCatch(chisq.test(tab), error = function(e) NULL)
        test_name <- "chi-squared"
        statistic <- if (!is.null(ct)) as.numeric(ct$statistic) else NA_real_
        p_value   <- if (!is.null(ct)) ct$p.value else NA_real_
      }

      if (is.logical(col_data)) {
        summary_ctrl <- list(n_true  = sum(ctrl_data == TRUE,  na.rm = TRUE),
                             n_false = sum(ctrl_data == FALSE, na.rm = TRUE))
        summary_trt  <- list(n_true  = sum(trt_data  == TRUE,  na.rm = TRUE),
                             n_false = sum(trt_data  == FALSE, na.rm = TRUE))
      } else {
        ctrl_counts        <- as.numeric(tab[, "0"])
        trt_counts         <- as.numeric(tab[, "1"])
        names(ctrl_counts) <- rownames(tab)
        names(trt_counts)  <- rownames(tab)
        summary_ctrl <- list(counts = ctrl_counts)
        summary_trt  <- list(counts = trt_counts)
      }

    } else {
      next
    }

    results_list$VARIABLE          <- c(results_list$VARIABLE,       var)
    results_list$TYPE              <- c(results_list$TYPE,            var_type)
    results_list$TEST              <- c(results_list$TEST,            test_name)
    results_list$STATISTIC         <- c(results_list$STATISTIC,       statistic)
    results_list$MIN_CELL_COUNT    <- c(results_list$MIN_CELL_COUNT,  min_cell_count)
    results_list$P_VALUE           <- c(results_list$P_VALUE,         p_value)
    n <- length(results_list$SUMMARY_CONTROL)
    results_list$SUMMARY_CONTROL[[n + 1]]   <- summary_ctrl
    results_list$SUMMARY_TREATMENT[[n + 1]] <- summary_trt
  }

  if (length(results_list$VARIABLE) == 0) return(NULL)

  report <- data.frame(
    VARIABLE       = results_list$VARIABLE,
    TYPE           = results_list$TYPE,
    TEST           = results_list$TEST,
    STATISTIC      = results_list$STATISTIC,
    MIN_CELL_COUNT = results_list$MIN_CELL_COUNT,
    P_VALUE        = results_list$P_VALUE,
    stringsAsFactors = FALSE
  )
  report$SUMMARY_CONTROL   <- results_list$SUMMARY_CONTROL
  report$SUMMARY_TREATMENT <- results_list$SUMMARY_TREATMENT

  return(report)
}


.create_randomization_report <- function(pheno_df, omics_df) {

  # Filter to baseline (FU=0)
  pheno_baseline <- pheno_df[pheno_df$FU == 0, ]

  # Filter omics to baseline samples
  baseline_sample_ids <- pheno_baseline$SAMPLE_ID
  omics_baseline <- omics_df[, colnames(omics_df) %in% baseline_sample_ids, drop = FALSE]

  # Initialize results data.frame
  results <- data.frame(
    ANALYTE_NAME = character(),
    MEAN_DIFFERENCE = numeric(),
    COHENS_D = numeric(),
    SE = numeric(),
    P_VALUE = numeric(),
    stringsAsFactors = FALSE
  )

  # T-test per analyte
  analyte_names <- omics_df$ANALYTE_NAME

  for (i in seq_along(analyte_names)) {
    tryCatch({
      # Extract values for this analyte
      analyte_values <- as.numeric(omics_baseline[i, ])
      control_status <- pheno_baseline$TREATMENT_GROUP

      # Split by group
      values_control <- analyte_values[control_status == 0]
      values_treatment <- analyte_values[control_status == 1]

      # Remove NAs
      values_control <- na.omit(values_control)
      values_treatment <- na.omit(values_treatment)

      # Skip if insufficient samples in either group
      if (length(values_control) < 2 || length(values_treatment) < 2) {
        return()
      }

      # Calculate statistics
      n_control <- length(values_control)
      n_treatment <- length(values_treatment)
      mean_control <- mean(values_control)
      mean_treatment <- mean(values_treatment)
      var_control <- var(values_control)
      var_treatment <- var(values_treatment)

      # Mean difference (Treatment - Control)
      mean_diff <- mean_treatment - mean_control

      # Standard Error (Welch's approach - doesn't assume equal variances)
      se <- sqrt((var_control / n_control) + (var_treatment / n_treatment))

      # Cohen's d (using pooled SD)
      pooled_sd <- sqrt(((n_control - 1) * var_control + (n_treatment - 1) * var_treatment) /
                          (n_control + n_treatment - 2))
      cohens_d <- mean_diff / pooled_sd

      # Welch's t-test
      tt <- t.test(values_treatment, values_control, var.equal = FALSE)

      # Add to results
      results <- rbind(results, data.frame(
        ANALYTE_NAME = analyte_names[i],
        MEAN_DIFFERENCE = mean_diff,
        COHENS_D = cohens_d,
        SE = se,
        P_VALUE = tt$p.value,
        stringsAsFactors = FALSE
      ))

    }, error = function(e) {
      warning("Error processing analyte '", analyte_names[i], "': ", e$message)
    })
  }

  # Apply multiple testing correction (global across all analytes)
  if (nrow(results) > 0) {
    results$BH_P_VALUE <- p.adjust(results$P_VALUE, method = "BH")
    return(results)
  } else {
    return(NULL)
  }
}


.generate_reports <- function(pheno_list, omics_list, additional_covariates = NULL) {

  # --- Variable summaries: per stratum, per FU x TREATMENT_GROUP cell ---
  variable_summaries <- list(all = NULL, male = NULL, female = NULL)

  for (dataset in c("all", "male", "female")) {
    if (is.null(pheno_list[[dataset]])) next

    pheno_df <- pheno_list[[dataset]]
    omics_df <- omics_list[[dataset]]

    fu_levels <- as.character(sort(unique(as.integer(as.character(pheno_df$FU)))))
    tx_levels <- c("0", "1")

    cell_reports <- list()

    for (fu in fu_levels) {
      for (tx in tx_levels) {
        cell_key   <- paste0("_FU", fu, "_Tx", tx)
        cell_mask  <- as.character(pheno_df$FU) == fu &
                      as.character(pheno_df$TREATMENT_GROUP) == tx
        cell_pheno <- pheno_df[cell_mask, ]

        if (nrow(cell_pheno) == 0) next

        cell_reports[[paste0("omics", cell_key)]] <-
          .create_omics_data_report(cell_pheno$SAMPLE_ID, omics_df)

        if (!is.null(additional_covariates)) {
          cell_reports[[paste0("covariates", cell_key)]] <-
            .create_addx_covariate_report(cell_pheno, additional_covariates)
        }
      }
    }

    variable_summaries[[dataset]] <- cell_reports
  }

  # --- Randomization reports: study-level, not sex-stratified ---
  randomization_reports <- list(
    omics_report     = .create_randomization_report(pheno_list$all, omics_list$all),
    covariate_report = .create_pheno_randomization_report(pheno_list$all, additional_covariates)
  )

  # --- Pheno summary: study-level ---
  pheno_summary <- .create_pheno_data_report(pheno_list$all)

  return(list(
    pheno_summary         = pheno_summary,
    variable_summaries    = variable_summaries,
    randomization_reports = randomization_reports
  ))
}
