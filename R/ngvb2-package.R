#' ngvb2: Latent Non-Gaussian Models via INLA and Variational Bayes
#'
#' Checks the latent Gaussian assumption ([ng.check()]) and fits latent
#' non-Gaussian models ([ngvb()]) on top of R-INLA. A single unified rgeneric
#' engine implements the conditional precision
#' \eqn{Q(\theta, V) = D(\theta)^T \mathrm{diag}(1/V) D(\theta)} and owns the
#' V-dependent normalizing constant, driven by a model registry
#' ([ngvb_operator()]) covering iid, rw1, rw2, ar1, SAR, CAR and SPDE/Matern.
#'
#' @keywords internal
"_PACKAGE"
