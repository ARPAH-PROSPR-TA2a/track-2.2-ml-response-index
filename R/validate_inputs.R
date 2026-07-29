# Helper function for NZV detection
.is_near_zero_variance <- function(x) {
  if (!is.numeric(x)) return(FALSE)
  if (all(is.na(x))) return(FALSE)
  var(x, na.rm = TRUE) < 1e-8
}


.validate_omics_type <- function(omics_type){

  acceptable_types <- c("DNAm", "Proteomics", "Metabolomics")

  if (!omics_type %in% acceptable_types) {
    stop(
      "Invalid omics_type '", omics_type, "'. ",
      "Must be one of: ", paste(acceptable_types, collapse = ", ")
    )
  }

  # Reminders about expected input format. These are not data checks --
  # the pipeline cannot tell whether values have been pre-processed
  # correctly -- they exist so the caller is reminded of the convention.
  if (omics_type == "DNAm") {
    message("DNAm: input should be M-values.")
  } else if (omics_type == "Metabolomics") {
    message("Metabolomics: inputs should be log2-transformed prior to analysis.")
  } else if (omics_type == "Proteomics") {
    message("Proteomics: inputs should be log2-transformed prior to analysis.")
  }
}


.validate_pheno <- function(pheno, additional_covariates = NULL) {
  
  # Step 1: Input validation and conversion
  if (is.matrix(pheno)) {
    pheno <- as.data.frame(pheno)
  } else if (!is.data.frame(pheno)) {
    stop("pheno must be a data.frame or matrix")
  }
  
  if (!is.null(additional_covariates) && !is.character(additional_covariates)) {
    stop("additional_covariates must be NULL or a character vector")
  }
  
  # Step 2: Required columns check
  required_cols <- c("SAMPLE_ID", "FU", "SUBJECT_ID", "FEMALE", "TREATMENT_GROUP")
  missing_cols <- setdiff(required_cols, names(pheno))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Step 3: Validate required columns (values, then types)
  # FU validation
  fu_num <- suppressWarnings(as.integer(as.character(pheno$FU)))
  if (any(is.na(fu_num))) {
    stop("FU must be integer-valued (e.g., 0, 1, 2, ...).")
  }
  if (any(fu_num < 0)) {
    stop("FU must be non-negative (baseline is FU == 0).")
  }
  if (!any(fu_num == 0)) {
    stop("pheno must contain at least one baseline sample (FU == 0)")
  }
  if (!any(fu_num == 1)) {
    stop("pheno must contain at least one follow-up sample (FU >= 1)")
  }
  # Enforce FU encoding as consecutive integers: 0, 1, 2, 3, ...
  # This prevents common trial encodings like months (0, 3, 6, 12).
  unique_fu <- sort(unique(fu_num))
  expected_fu <- seq.int(0, max(unique_fu))
  if (!identical(unique_fu, expected_fu)) {
    stop(
      "FU must be encoded as consecutive integers starting at 0: ",
      paste(expected_fu, collapse = ", "),
      ". Found: ", paste(unique_fu, collapse = ", "),
      ". If your raw follow-up is in months (e.g., 0,3,6,12), recode to 0,1,2,3 before running."
    )
  }
  if (!is.factor(pheno$FU)) {
    warning("FU column is not a factor. Converting to factor.")
    pheno$FU <- factor(fu_num)
  } else {
    # Normalize factor representation (ensures levels are "0","1",... in order).
    pheno$FU <- factor(fu_num)
  }
  
  # FEMALE validation
  if (!all(pheno$FEMALE %in% 0:1)) {
    stop("FEMALE must contain only values 0 or 1")
  }
  if (!is.factor(pheno$FEMALE)) {
    warning("FEMALE column is not a factor. Converting to factor.")
    pheno$FEMALE <- factor(pheno$FEMALE)
  }
  
  # TREATMENT_GROUP validation
  if (!all(pheno$TREATMENT_GROUP %in% 0:1)) {
    stop("TREATMENT_GROUP must contain only values 0 or 1")
  }
  if (!any(pheno$TREATMENT_GROUP == 0) || !any(pheno$TREATMENT_GROUP == 1)) {
    stop("pheno must contain both control (TREATMENT_GROUP == 0) and treatment (TREATMENT_GROUP == 1)")
  }
  if (!is.factor(pheno$TREATMENT_GROUP)) {
    warning("TREATMENT_GROUP column is not a factor. Converting to factor.")
    pheno$TREATMENT_GROUP <- factor(pheno$TREATMENT_GROUP)
  }
  
  # Step 4: Column-specific validation
  if (any(duplicated(pheno$SAMPLE_ID))) {
    stop("SAMPLE_ID contains duplicate values")
  }
  
  # Check SUBJECT_ID/FU pair uniqueness
  subject_fu_pairs <- paste(pheno$SUBJECT_ID, pheno$FU, sep = "_")
  if (any(duplicated(subject_fu_pairs))) {
    # Keep only first occurrence of each SUBJECT_ID/FU pair
    warning("Found duplicate SUBJECT_ID/FU pairs. Keeping first occurrence, discarding replicates.")
    pheno <- pheno[!duplicated(subject_fu_pairs), ]
  }
  
  # Step 5: Additional covariates validation and NA filtering
  if (!is.null(additional_covariates)) {
    # Check all covariates exist in pheno
    missing_addl <- setdiff(additional_covariates, names(pheno))
    if (length(missing_addl) > 0) {
      stop("Additional covariates not found in pheno: ", paste(missing_addl, collapse = ", "))
    }
  }

  # Validate each covariate type and drop rows with NA values
  for (covar in if (is.null(additional_covariates)) character(0) else additional_covariates) {
    col_data <- pheno[[covar]]

    if (is.ordered(col_data)) {
      stop(
        "Additional covariate '", covar, "' must be an unordered factor. ",
        "Ordered factors are not accepted because their default contrasts are polynomial."
      )
    }

    is_allowed_covar <- is.numeric(col_data) ||
      is.integer(col_data) ||
      is.factor(col_data) ||
      is.logical(col_data)

    if (!is_allowed_covar) {
      stop(
        "Additional covariate '", covar, "' must be numeric, integer, factor, or logical. ",
        "Character covariates are not accepted."
      )
    }

    na_rows <- is.na(col_data)
    if (any(na_rows)) {
      message("Dropping ", sum(na_rows), " sample(s) with NA values in covariate '", covar, "'.")
      pheno <- pheno[!na_rows, ]
    }
  }

  # Step 5b: Subject-level completeness — keep only subjects with both a
  # baseline (FU == 0) and at least one follow-up (FU > 0).  This ensures
  # the analysis and reports operate on exactly the same set of subjects.
  fu_num_current <- as.integer(as.character(pheno$FU))
  subjects_with_baseline <- unique(pheno$SUBJECT_ID[fu_num_current == 0])
  subjects_with_followup <- unique(pheno$SUBJECT_ID[fu_num_current > 0])
  complete_subjects <- intersect(subjects_with_baseline, subjects_with_followup)

  n_incomplete <- length(unique(pheno$SUBJECT_ID)) - length(complete_subjects)
  if (n_incomplete > 0) {
    message("Dropping ", n_incomplete, " subject(s) missing either a baseline or a follow-up sample.")
  }
  pheno <- pheno[pheno$SUBJECT_ID %in% complete_subjects, ]

  if (nrow(pheno) == 0) {
    stop("No subjects remain after requiring both a baseline and a follow-up sample with complete covariates.")
  }

  # Step 6: Column cleanup
  cols_to_keep <- c(required_cols, additional_covariates)
  cols_to_keep <- intersect(cols_to_keep, names(pheno))
  
  pheno <- pheno[, cols_to_keep, drop = FALSE]
  
  pheno
}


.validate_omics <- function(omics, pheno_df) {
  
  # Step 1: Input validation and conversion
  if (is.matrix(omics)) {
    omics <- as.data.frame(omics)
  } else if (!is.data.frame(omics)) {
    stop("omics must be a data.frame or matrix")
  }
  
  # Step 2: Extract analyte names
  if (!"ANALYTE_NAME" %in% names(omics)) {
    stop("omics must contain ANALYTE_NAME column")
  }
  
  analyte_names <- omics$ANALYTE_NAME
  if (any(duplicated(analyte_names))) {
    stop("ANALYTE_NAME contains duplicate values")
  }
  
  # Remove ANALYTE_NAME column to get numeric data
  omics_numeric <- omics[, setdiff(names(omics), "ANALYTE_NAME"), drop = FALSE]
  
  # Step 3: Validate all data is numeric
  if (!all(sapply(omics_numeric, is.numeric))) {
    stop("All columns in omics (except ANALYTE_NAME) must be numeric")
  }
  
   # Step 4: Filter to shared SAMPLE_IDs with pheno_df
   pheno_sample_ids <- pheno_df$SAMPLE_ID
   omics_sample_ids <- names(omics_numeric)
    
   shared_samples <- intersect(omics_sample_ids, pheno_sample_ids)
   
   if (length(shared_samples) == 0) {
     stop("No overlap between omics column names and pheno SAMPLE_IDs")
   }
   
   message("Found ", length(shared_samples), " samples shared between omics and pheno")
   message("  Omics only: ", length(setdiff(omics_sample_ids, pheno_sample_ids)))
   message("  Pheno only: ", length(setdiff(pheno_sample_ids, omics_sample_ids)))
  
  # Filter omics to shared samples (keep order from pheno for consistency)
  omics_numeric <- omics_numeric[, shared_samples, drop = FALSE]
  
  # Step 5: Quality checks
  n_with_na <- 0
  n_with_nzv <- 0
  
  for (i in seq_along(analyte_names)) {
    analyte_data <- as.numeric(omics_numeric[i, ])
    
    if (any(is.na(analyte_data))) {
      n_with_na <- n_with_na + 1
    }
    
    if (.is_near_zero_variance(analyte_data)) {
      n_with_nzv <- n_with_nzv + 1
    }
  }
  
  # Summmarize analyte NA/NZVs
  if (n_with_na > 0) {
    warning(n_with_na, " analytes contain NA values")
  }
  
  if (n_with_nzv > 0) {
    warning(n_with_nzv, " analytes have near-zero variance")
  }
  
  # Add ANALYTE_NAME back to the omics data.frame
  omics_all <- cbind(ANALYTE_NAME = analyte_names, omics_numeric)

  omics_all
}


# Helper function to validate DNAm probe coverage against available data
.validate_dnam_probe_coverage <- function(full_probes, filtered_probes, available_probes) {
  full_present <- sum(full_probes %in% available_probes)
  full_missing <- length(full_probes) - full_present
  if (full_present == 0) {
    stop("No probes from full probe list found in data")
  }
  if (full_missing > 0) {
    warning(sprintf("%d of %d probes from full probe list not found in data",
                    full_missing, length(full_probes)))
  }

  filtered_present <- sum(filtered_probes %in% available_probes)
  filtered_missing <- length(filtered_probes) - filtered_present
  if (filtered_present == 0) {
    stop("No probes from filtered probe list found in data")
  }
  if (filtered_missing > 0) {
    warning(sprintf("%d of %d probes from filtered probe list not found in data",
                    filtered_missing, length(filtered_probes)))
  }
}


# Helper function to subset omics data to a specific set of analytes
.subset_omics <- function(omics_df, analyte_subset) {
  if (is.null(analyte_subset)) {
    return(omics_df)
  }

  matching_analytes <- omics_df$ANALYTE_NAME %in% analyte_subset
  omics_df[matching_analytes, , drop = FALSE]
}


.prepare_inputs <- function(pheno, omics, omics_type, additional_covariates = NULL) {
  .validate_omics_type(omics_type)

  pheno_df <- .validate_pheno(pheno, additional_covariates)
  omics_df <- .validate_omics(omics, pheno_df)

  if (omics_type == "DNAm") {
    full_probes <- readRDS("Data/FAST_epicv1_epicv2_probe_list.rds")
    reliable_probes <- readRDS("Data/FAST_epicv1_epicv2_sugden_TruD_probe_list.rds")
    .validate_dnam_probe_coverage(full_probes, reliable_probes, omics_df$ANALYTE_NAME)
    omics_df <- .subset_omics(omics_df, reliable_probes)
    message("DNAm: restricted analysis to ", nrow(omics_df), " reliable probes.")
  }

  analyte_map <- .make_xgb_safe_analyte_map(omics_df$ANALYTE_NAME)
  omics_df$ANALYTE_NAME <- analyte_map$INTERNAL_ANALYTE_NAME

  list(
    pheno = pheno_df,
    omics = omics_df,
    analyte_name_map = analyte_map
  )
}


.make_xgb_safe_analyte_map <- function(analyte_names) {
  analyte_names <- as.character(analyte_names)
  
  internal <- janitor::make_clean_names(
    analyte_names,
    case = "snake",
    ascii = TRUE,
    allow_dupes = FALSE
  )

  if (anyDuplicated(internal)) {
    stop(
      "XGB-safe analyte-name mapping is not one-to-one. ",
      "Please resolve duplicate analyte names after sanitization."
    )
  }

  data.frame(
    ORIGINAL_ANALYTE_NAME = analyte_names,
    INTERNAL_ANALYTE_NAME = internal,
    WAS_MODIFIED = internal != analyte_names,
    stringsAsFactors = FALSE
  )
}

.apply_analyte_map_to_omics <- function(omics_df, analyte_map) {
  idx <- match(omics_df$ANALYTE_NAME, analyte_map$ORIGINAL_ANALYTE_NAME)
  omics_df$ANALYTE_NAME[!is.na(idx)] <- analyte_map$INTERNAL_ANALYTE_NAME[idx[!is.na(idx)]]
  omics_df
}

.restore_original_analyte_name <- function(df, analyte_map,
                                           analyte_col = "ANALYTE_NAME") {
  if (is.null(analyte_map) || nrow(analyte_map) == 0L ||
      !analyte_col %in% names(df)) {
    return(df)
  }

  idx <- match(df[[analyte_col]], analyte_map$INTERNAL_ANALYTE_NAME)
  mapped <- !is.na(idx)
  df[[analyte_col]][mapped] <- analyte_map$ORIGINAL_ANALYTE_NAME[idx[mapped]]
  df
}

.add_original_feature_name <- function(df, analyte_map,
                                       feature_col = "FEATURE_NAME") {
  if (is.null(analyte_map) || nrow(analyte_map) == 0L ||
      !feature_col %in% names(df)) {
    return(df)
  }

  feature_names <- df[[feature_col]]
  original_feature_names <- feature_names

  is_omics <- startsWith(feature_names, "omics::")
  internal_analytes <- sub("^omics::", "", feature_names[is_omics])
  idx <- match(internal_analytes, analyte_map$INTERNAL_ANALYTE_NAME)
  mapped <- !is.na(idx)

  omics_positions <- which(is_omics)
  original_feature_names[omics_positions[mapped]] <- paste0(
    "omics::",
    analyte_map$ORIGINAL_ANALYTE_NAME[idx[mapped]]
  )

  df$ORIGINAL_FEATURE_NAME <- original_feature_names
  df
}
