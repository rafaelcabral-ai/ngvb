## The rgeneric engine reads these from its definition environment (populated by
## inla.rgeneric.define); they are not true globals. Declared to quiet R CMD check.
## `.data` is the ggplot2/rlang pronoun used in aes(); declare it so R CMD check
## does not flag "no visible binding for global variable '.data'".
utils::globalVariables(c("Dfunc", "Vinv", "rankdef", "graph.pattern",
                         "theta.initial", "logprior", "lognc", ".data"))

## INLA is a Suggests dependency (it is not on CRAN). Every entry point that
## needs it calls this first, so a missing INLA fails with an actionable message
## instead of an opaque "could not find function" error.
#' @keywords internal
.need_inla <- function() {
  if (!requireNamespace("INLA", quietly = TRUE))
    stop("ngvb2 requires the 'INLA' package, which is not on CRAN. Install it with:\n",
         "  install.packages('INLA', repos = c(getOption('repos'),\n",
         "                    INLA = 'https://inla.r-inla-download.org/R/stable'), dep = TRUE)",
         call. = FALSE)
  invisible(TRUE)
}
