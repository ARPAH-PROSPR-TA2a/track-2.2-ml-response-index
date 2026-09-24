# Standalone synthetic runner check: Rscript tests/test_calerie_dnam_runner.R
# Uses temporary input/output files and the real ENET, XGBoost and export paths.
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) == 1L) {
  setwd(dirname(dirname(normalizePath(sub("^--file=", "", file_arg)))))
}
repo <- getwd()
source("training/main.R")
source("tests/helpers.R")
options(track22.tests.passed = 0L)
if (!requireNamespace("haven", quietly = TRUE)) stop("Tests require haven.")
python <- Sys.getenv("FAST_TEST_PYTHON", file.path(repo, ".track22-python/bin/python"))
if (!file.exists(python)) python <- unname(Sys.which("python3"))
runner <- parse("training/CALERIE/DNAm/run_DNAm_betas.R")
assignment <- function(expr, name) {
  is.call(expr) && identical(expr[[1L]], as.name("<-")) &&
    identical(expr[[2L]], as.name(name))
}

run_tests <- function() {
  work <- tempfile("track22-calerie-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  raw_paths <- file.path(work, c("betas.rds", "pheno.rds"))
  set.seed(2202L)
  probes <- head(readRDS("Data/FAST_epicv1_epicv2_sugden_TruD_probe_list.rds"), 8L)
  visits <- rbind(
    data.frame(subject = 1:50, fu = 0L),
    data.frame(subject = c(1:40, 49L), fu = 1L),
    data.frame(subject = 9:48, fu = 2L)
  )
  subjects <- sprintf("SIM-%02d", visits$subject)
  barcodes <- paste0("2014 chip-", subjects, "/visit", visits$fu)
  pheno <- data.frame(Barcode = barcodes, Participant_ID = subjects,
                      Time_Point = paste0("Month", visits$fu * 12L))
  pheno$fu <- haven::labelled(visits$fu, c(Baseline = 0L, Month12 = 1L, Month24 = 2L))
  pheno$CR <- haven::labelled(visits$subject %% 2L, c(Control = 0L, Treatment = 1L))
  pheno$female <- haven::labelled((visits$subject %/% 2L) %% 2L, c(Male = 0L, Female = 1L))
  # Missing required data, a pheno-only row and a missing matched baseline.
  pheno$female[visits$subject == 50L] <- NA
  pheno <- rbind(pheno, pheno[1L, ])
  pheno$Barcode[nrow(pheno)] <- "phenotype-only"
  pheno$Participant_ID[nrow(pheno)] <- "SIM-unmatched"
  beta <- matrix(runif(9L * length(barcodes), 0.2, 0.8), nrow = 9L,
                 dimnames = list(c(probes, "synthetic_nonreliable"), barcodes))
  beta <- beta[, colnames(beta) != barcodes[visits$subject == 49L & visits$fu == 0L], drop = FALSE]
  beta <- cbind(beta, "omics-only" = rep(0.5, nrow(beta)))
  beta <- beta[, sample(ncol(beta)), drop = FALSE]
  pheno <- pheno[sample(nrow(pheno)), , drop = FALSE]
  # A later assay of the same visit has distinctive values; the first raw row
  # must supply the beta deltas checked in the real training run below.
  winning_barcode <- pheno$Barcode[pheno$Participant_ID == "SIM-01" & as.numeric(pheno$fu) == 1L]
  later_replicate <- pheno[match(winning_barcode, pheno$Barcode), ]
  later_barcode <- later_replicate$Barcode <- "later-technical-replicate"
  middle_replicate <- later_replicate
  middle_replicate$Barcode <- "middle-technical-replicate"
  replicate_barcodes <- c(middle_replicate$Barcode, later_barcode)
  pheno <- rbind(pheno, middle_replicate, later_replicate)
  beta <- cbind(beta, "middle-technical-replicate" = rep(0.05, nrow(beta)),
                "later-technical-replicate" = rep(0.95, nrow(beta)))
  save_inputs <- function(b = beta, p = pheno) {
    saveRDS(b, raw_paths[1L])
    saveRDS(p, raw_paths[2L])
  }
  save_inputs()

  # Change only the script's editable settings, leaving its execution path intact.
  overrides <- list(
    pipeline_repo = repo, omics_raw_path = raw_paths[1L], pheno_raw_path = raw_paths[2L],
    python_bin = python, out_dir = file.path(work, "success"), n_cores = 1L,
    enet_cv_folds = 4L, xgb_cv_folds = 2L, xgb_cv_repeats = 1L, xgb_n_trials = 10L
  )
  expressions <- runner
  for (name in names(overrides)) {
    idx <- which(vapply(expressions, assignment, logical(1), name = name))
    stopifnot(length(idx) == 1L)
    expressions[[idx]][[3L]] <- overrides[[name]]
  }
  execution <- new.env(parent = globalenv())
  suppressWarnings(for (expr in expressions) eval(expr, execution))
  manifest <- execution$manifest
  .expect_true(setequal(names(manifest$followups), c("FU1", "FU2")), "both follow-ups trained")
  .expect_true(is.null(manifest$additional_covariates) &&
                 identical(manifest$model_covariates, "FEMALE"),
               "runner defaults to FEMALE without extra adjustment covariates")
  .expect_true(identical(names(execution$omics)[-1L], execution$pheno$SAMPLE_ID),
               "shuffled nonsyntactic sample names align without renaming")
  .expect_true(winning_barcode %in% execution$pheno$SAMPLE_ID &&
                 !any(replicate_barcodes %in% execution$pheno$SAMPLE_ID),
               "training selects the first raw assay of a replicated visit")
  preprocessing <- read.csv(manifest$reports$preprocessing, stringsAsFactors = FALSE)
  for (fu in 1:2) {
    artifacts <- manifest$followups[[paste0("FU", fu)]]$artifacts
    retained <- read.csv(artifacts$subjects, stringsAsFactors = FALSE)
    intended <- sprintf("SIM-%02d", if (fu == 1L) 1:40 else 9:48)
    expected_order <- pheno$Participant_ID[!pheno$Barcode %in% replicate_barcodes & as.numeric(pheno$fu) == fu &
                                           pheno$Participant_ID %in% intended]
    .expect_true(identical(retained$SUBJECT_ID, expected_order) && nrow(retained) == 40L,
                 paste0("FU", fu, " preserves its exact baseline-paired population and order"))
    baseline_ids <- pheno$Barcode[match(paste0(retained$SUBJECT_ID, "/0"),
                                      paste0(pheno$Participant_ID, "/", as.numeric(pheno$fu)))]
    followup_ids <- pheno$Barcode[match(paste0(retained$SUBJECT_ID, "/", fu),
                                      paste0(pheno$Participant_ID, "/", as.numeric(pheno$fu)))]
    delta <- t(beta[probes, followup_ids] - beta[probes, baseline_ids])
    for (model in c("enet", "xgb")) {
      x <- read.csv(artifacts[[paste0(model, "_train")]], check.names = FALSE)
      .expect_true(setequal(names(x), c(paste0("omics::", probes), "covariate::FEMALE")),
                   paste0("FU", fu, " ", model, " selects reliable probes and FEMALE only"))
      prep <- preprocessing[preprocessing$FU == fu, ]
      prep <- prep[match(paste0("omics::", probes), prep$FEATURE_NAME), ]
      recovered <- sweep(sweep(as.matrix(x[paste0("omics::", probes)]), 2L,
                               prep$SCALE, "*"), 2L, prep$CENTER, "+")
      .expect_equal(recovered, unname(delta),
                    paste0("FU", fu, " ", model, " trains on exact beta follow-up minus baseline"))
    }
  }
  exported <- execution$exported$models
  .expect_true(nrow(exported) == 4L && all(file.exists(exported$PATH)) &&
                 setequal(paste(exported$FU, exported$MODEL), c("1 enet", "1 xgb", "2 enet", "2 xgb")),
               "four fitted model packages exported")
  provenance <- readRDS(file.path(execution$out_dir, "provenance.rds"))
  .expect_true(identical(provenance$measurement_scale, "beta") &&
                 identical(provenance$feature_mode, "follow-up minus baseline") &&
                 identical(provenance$inputs, setNames(raw_paths, c("omics", "pheno"))) &&
                 provenance$n_reliable_probes == length(probes) &&
                 !is.null(provenance$finished) && provenance$elapsed_seconds > 0,
               "completed provenance records beta scale, input paths, coverage and timing")
  selection <- provenance$replicate_selection
  .expect_true(identical(selection$rule, "first raw row per Participant_ID/Time_Point") &&
                 identical(selection$stage, "before completeness filtering and sample matching") &&
                 identical(selection$retained_order, "raw input order") &&
                 selection$raw_rows == nrow(pheno) && selection$duplicate_groups == 1L &&
                 selection$removed_rows == 2L && selection$retained_rows == nrow(pheno) - 2L,
               "provenance records the established replicate rule and aggregate counts")
  log <- readLines(file.path(execution$out_dir, "run.log"))
  .expect_true(!any(grepl("SIM-|chip-", log)), "runner status log contains no sample or subject identifiers")

  # Exercise alternate containers and invalid inputs through real preflight,
  # without fitting more models or replacing any validation function.
  start <- which(vapply(runner, assignment, logical(1), name = "omics_raw"))
  end <- which(vapply(runner, assignment, logical(1), name = "reliable_probes")) - 1L
  preflight <- runner[start:end]
  run_preflight <- function(b = beta, p = pheno) {
    save_inputs(b, p)
    context <- list2env(list(omics_raw_path = raw_paths[1L], pheno_raw_path = raw_paths[2L],
                            models = c("enet", "xgb"), enet_cv_folds = 4L,
                            xgb_cv_folds = 2L, log_status = function(text) invisible(NULL)),
                       parent = globalenv())
    for (expr in preflight) eval(expr, context)
    context
  }
  reject <- function(label, pattern, b = beta, p = pheno) {
    .expect_error(run_preflight(b, p), pattern, label)
  }
  numeric_frame <- as.data.frame(beta, check.names = FALSE)
  rownames(numeric_frame)[nrow(numeric_frame)] <- "synthetic [non reliable]"
  numeric_matrix <- beta
  rownames(numeric_matrix) <- rownames(numeric_frame)
  frame_preflight <- run_preflight(numeric_frame)
  .expect_true(identical(frame_preflight$omics_raw, numeric_matrix),
               "numeric data.frame preserves matrix values and nonsyntactic probe/sample names")
  .expect_true(identical(frame_preflight$pheno, execution$pheno),
               "numeric data.frame retains the matrix input's exact paired cohort and ordering")
  mixed_frame <- numeric_frame
  mixed_frame[[1L]] <- as.character(mixed_frame[[1L]])
  reject("mixed numeric/character data.frame rejected without coercion", "numeric", b = mixed_frame)
  character_frame <- numeric_frame
  character_frame[] <- lapply(character_frame, as.character)
  reject("character data.frame rejected without coercion", "numeric", b = character_frame)
  invalid_code <- pheno
  invalid_code$fu <- as.numeric(invalid_code$fu)
  invalid_code$fu[1L] <- 0.5
  reject("fractional visit code rejected", "fu must use codes", p = invalid_code)
  for (value in c(NA_real_, Inf, -0.01, 1.01)) {
    invalid_beta <- beta
    invalid_beta[1L, 1L] <- value
    reject(paste("invalid beta rejected:", value), "complete, finite DNAm beta", b = invalid_beta)
  }
  reordered <- pheno[c(nrow(pheno), seq_len(nrow(pheno) - 1L)), ]
  selected <- run_preflight(p = reordered)$pheno$SAMPLE_ID
  .expect_true(later_barcode %in% selected && !winning_barcode %in% selected,
               "replicate selection follows raw row order rather than barcode sorting")
  invalid_discarded <- pheno
  invalid_discarded$fu[nrow(pheno)] <- 99L
  .expect_true(identical(run_preflight(p = invalid_discarded)$pheno, execution$pheno),
               "replicate selection precedes visit-code validation")
  first_incomplete <- pheno
  first_incomplete$female[match(winning_barcode, first_incomplete$Barcode)] <- NA
  .expect_true(!"SIM-01" %in% run_preflight(p = first_incomplete)$pheno$SUBJECT_ID,
               "missing required value in first assay does not promote a later assay")
  .expect_true(!"SIM-01" %in% run_preflight(b = beta[, colnames(beta) != winning_barcode])$pheno$SUBJECT_ID,
               "unmatched first assay does not promote a later matched assay")
  duplicate_barcode <- pheno[match(winning_barcode, pheno$Barcode), ]
  duplicate_barcode$female <- NA
  reject("duplicate raw Barcode rejected even when later row is incomplete", "Duplicate Barcode",
         p = rbind(pheno, duplicate_barcode))
  ambiguous_visit <- pheno
  ambiguous_visit$Time_Point[nrow(pheno)] <- "OtherTimePoint"
  reject("distinct Time_Point groups mapping to the same subject/FU rejected", "Duplicate usable",
         p = ambiguous_visit)
  inconsistent <- pheno
  idx <- which(inconsistent$Participant_ID == "SIM-01" & as.numeric(inconsistent$fu) == 1L)
  inconsistent$CR[idx] <- 1L - as.numeric(inconsistent$CR[idx])
  reject("contradictory within-subject treatment rejected", "TREATMENT_GROUP is inconsistent", p = inconsistent)
  small_fu2 <- pheno[as.numeric(pheno$fu) != 2L | pheno$Participant_ID %in% sprintf("SIM-%02d", 9:14), ]
  reject("insufficient FU2 fails before either model starts", "FU2: ENET needs at least 4", p = small_fu2)
  cat("CALERIE runner: ", getOption("track22.tests.passed"), " assertions passed\n", sep = "")
}

run_tests()
