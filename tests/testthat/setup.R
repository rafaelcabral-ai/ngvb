# INLA is a Suggests dependency and is not on CRAN, where installations are
# frequently broken. Nothing in the suite may touch it there: every test that
# needs INLA calls skip_if_no_inla(), which skips on CRAN first and only then
# looks for the package. INLA is never attached -- all calls are INLA::-prefixed
# -- so a broken install cannot fail the suite at load time either.
skip_if_no_inla <- function() {
  testthat::skip_on_cran()
  if (!requireNamespace("INLA", quietly = TRUE)) testthat::skip("INLA not available")
}
