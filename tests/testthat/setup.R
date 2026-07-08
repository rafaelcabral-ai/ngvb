# INLA is a Suggests dependency (not on CRAN) and is no longer attached by the
# package. Attach it for the tests that call inla() directly; individual tests
# still skip_if_not_installed("INLA") so the suite is a no-op without it.
if (requireNamespace("INLA", quietly = TRUE)) suppressMessages(library(INLA))
