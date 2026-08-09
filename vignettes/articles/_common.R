# Shared setup for the pkgdown articles.
#
# Everything under vignettes/articles/ is excluded from the package tarball via
# .Rbuildignore, so these documents are never built or run by CRAN. They are
# rendered only by the pkgdown workflow, which installs INLA explicitly. The
# guards below still let the site build when a suggested package is missing.

knitr::opts_chunk$set(collapse = TRUE, comment = "#>",
                      fig.width = 8, fig.height = 3.6, fig.align = "center",
                      dpi = 96, message = FALSE, warning = FALSE,
                      eval = requireNamespace("INLA", quietly = TRUE))

have_inla    <- requireNamespace("INLA", quietly = TRUE)
have_areal   <- have_inla && all(vapply(c("spdep", "sf", "spData"), requireNamespace,
                                        logical(1), quietly = TRUE))
have_fmesher <- have_inla && requireNamespace("fmesher", quietly = TRUE)
have_ngme2   <- have_inla && requireNamespace("ngme2", quietly = TRUE)

set.seed(1)
