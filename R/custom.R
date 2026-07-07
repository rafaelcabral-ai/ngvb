## ---------------------------------------------------------------------------
## Custom operators. The fundamental object in ngvb2 is an "operator" = a
## dependency matrix D(theta) plus a constant vector h. The latent field has
## conditional precision Q(V) = D(theta)^T diag(1/V) D(theta), with V = h giving
## the Gaussian base model Q = D^T diag(1/h) D. Built-in models are auto-detected;
## anything else you define with ngvb_custom() by giving D(theta) and h.
## ---------------------------------------------------------------------------

## The rgeneric callbacks run in a separate R process to which the model function
## is serialized. Closures written at the top level reference the GLOBAL
## environment, whose contents are NOT serialized (so `object 'W' not found` in
## the subprocess). Capture the data each function actually uses into a small,
## serializable environment so custom operators work wherever they are defined.
#' @keywords internal
.ngvb_capture <- function(f) {
  if (!is.function(f)) return(f)
  vars <- tryCatch(codetools::findGlobals(f, merge = FALSE)$variables,
                   error = function(e) character(0))
  src  <- environment(f)
  e    <- new.env(parent = globalenv())
  for (v in vars)
    if (exists(v, envir = src, inherits = TRUE))
      assign(v, get(v, envir = src, inherits = TRUE), envir = e)
  environment(f) <- e
  f
}

#' Define a custom latent non-Gaussian model from a precision "square root".
#'
#' Supply the dependency matrix `D(theta)` (a function of the hyperparameters
#' INLA fits, on the internal/unconstrained scale) and the constant vector `h`.
#' The field then has conditional precision
#' `Q(V) = D(theta)^T diag(1/V) D(theta)`, and `V = h` recovers the Gaussian model.
#' This turns *any* Gaussian model whose precision factors as `Q = D^T diag(1/h) D`
#' into a latent non-Gaussian one, without hand-writing an rgeneric.
#'
#' Typical workflow: build `op`, fit the Gaussian LGM with
#' `f(idx, model = ngvb_rgeneric(op))`, then call `ng.check(fit, components = ...)`
#' and `ngvb(fit, components = ...)`.
#'
#' @param D `function(theta)` returning the sparse dependency matrix `D(theta)`,
#'   of dimension `m x n` (`m` driving-noise terms, `n` latent nodes).
#' @param h constant vector of length `m` (`E[V] = h`; `1` for discrete models).
#'   A scalar is recycled.
#' @param ntheta number of hyperparameters.
#' @param theta.initial initial hyperparameters (internal scale).
#' @param logprior `function(theta)` log prior on the internal scale (including
#'   Jacobians). Default: independent `N(0, 3)`.
#' @param rankdef rank deficiency of `D^T D` (0 = proper; positive for intrinsic
#'   models such as `m < n` difference operators).
#' @param graph optional sparsity pattern of `D^T D` (the union over `theta`);
#'   computed from `D(theta.initial)` if omitted.
#' @param lognc optional analytic normalizing constant `function(theta, Vinv)`
#'   returning `0.5 log|Q| - 0.5 n log(2 pi)`. Supply it when the generic
#'   determinant is numerically unstable (e.g. `D` near-singular, as in SAR/CAR).
#' @return An operator descriptor for use in [ngvb()] / [ng.check()] via
#'   `components`, or to build a Gaussian fit via [ngvb_rgeneric()].
#' @note `D`, `logprior` and `lognc` are evaluated in a separate R process where
#'   only base packages are attached. Data they reference (matrices, vectors) is
#'   captured automatically, but any *non-base function* must be namespace-qualified
#'   (e.g. `Matrix::Diagonal`, `stats::dnorm`) or written with base R. For instance
#'   use `exp(t)/(1 + exp(t))` rather than `plogis(t)`.
#' @seealso [ngvb()], [ng.check()], [ngvb_rgeneric()]
#' @examples
#' # First differences by hand (an RW1), showing the custom interface:
#' n  <- 8
#' D0 <- Matrix::bandSparse(n - 1, n, c(0, 1), list(rep(-1, n - 1), rep(1, n - 1)))
#' op <- ngvb_custom(D = function(theta) sqrt(exp(theta[1])) * D0,
#'                   h = 1, ntheta = 1, rankdef = 1)
#' dim(op$Dfunc(0))
#' @export
ngvb_custom <- function(D, h = 1, ntheta, theta.initial = rep(0, ntheta),
                        logprior = function(theta) sum(stats::dnorm(theta, 0, 3, log = TRUE)),
                        rankdef = 0L, graph = NULL, lognc = NULL) {
  stopifnot(is.function(D), length(ntheta) == 1L, ntheta >= 1)
  D0 <- D(theta.initial)
  m  <- nrow(D0); n <- ncol(D0)
  if (length(h) == 1L) h <- rep(h, m)
  if (length(h) != m)
    stop("ngvb2: h must have length nrow(D(theta)) = ", m, " (the number of noise terms).")
  if (is.null(graph)) graph <- Matrix::crossprod(D0)
  ## make the user functions self-contained for the rgeneric subprocess
  op <- list(type = "custom", n = n, ntheta = as.integer(ntheta),
             rankdef = as.integer(rankdef), h = h, theta.initial = theta.initial,
             Dfunc = .ngvb_capture(D), graph = graph, logprior = .ngvb_capture(logprior))
  if (!is.null(lognc)) op$lognc <- .ngvb_capture(lognc)
  op
}
