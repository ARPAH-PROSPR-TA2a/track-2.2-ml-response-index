.fast_r_requirements <- list(
  training = c("glmnet", "jsonlite", "pROC"),
  validation = c("jsonlite", "pROC", "xgboost")
)


.fast_r_check_label <- function(tool) {
  switch(
    tool,
    training = "training",
    validation = "cross-trial validation",
    all = "training and cross-trial validation"
  )
}


FAST_check_R <- function(tool = c("all", "training", "validation")) {
  tool <- match.arg(tool)
  packages <- if (tool == "all") {
    unique(unlist(.fast_r_requirements, use.names = FALSE))
  } else {
    .fast_r_requirements[[tool]]
  }
  missing <- packages[
    !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  label <- .fast_r_check_label(tool)

  if (length(missing) > 0L) {
    install_command <- paste0(
      "install.packages(c(",
      paste(shQuote(missing), collapse = ", "),
      "))"
    )
    stop(
      "R is not ready for FAST ", label, ".\n",
      "Missing required R packages: ", paste(missing, collapse = ", "),
      ".\n\n",
      "What to do next:\n",
      "  1. In R, run:\n",
      "     ", install_command, "\n",
      "  2. Run FAST_check_R(\"", tool, "\") again.\n",
      "  3. Open the Troubleshooting section in README.md for more help.",
      call. = FALSE
    )
  }

  message(
    "R is ready for FAST ", label, ".\n",
    "Using: ", R.version.string
  )
  invisible(packages)
}
