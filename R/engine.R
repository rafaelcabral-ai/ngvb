## ---------------------------------------------------------------------------
## ngvb2 unified latent engine (rgeneric).
##
## A single custom INLA latent model whose precision conditioned on the mixing
## vector V is
##        Q(theta, V) = D(theta)^T diag(1 / V) D(theta).
##
## The engine OWNS the V-dependent normalizing constant, fixing the bias of
## INLA's `generic0` (which silently drops 0.5*log|Cmatrix|, so the term
## -0.5 * sum(log V) -- which changes every variational iteration -- is lost):
##   * full-rank D  : return numeric(0); INLA computes 0.5*log|Q| from its
##                    Cholesky, correctly including the V-dependence.
##   * intrinsic D  : supply the generalized log-determinant explicitly.
##
## rgeneric callbacks run in a separate R process where Matrix is loaded but NOT
## attached, so every Matrix symbol is qualified (Matrix::, methods::as).
## ---------------------------------------------------------------------------

#' rgeneric model function implementing the ngvb conditional precision.
#'
#' Not called directly; passed to `INLA::inla.rgeneric.define()` by
#' [ngvb_rgeneric()]. Reads `Dfunc`, `Vinv`, `rankdef`, `ntheta`,
#' `theta.initial`, `logprior`, `graph.pattern` from its definition environment.
#' @keywords internal
ngvb.rgeneric.engine <- function(
    cmd = c("graph", "Q", "mu", "initial", "log.norm.const", "log.prior", "quit"),
    theta = NULL) {

  envir <- parent.env(environment())

  if (is.null(theta) || length(theta) == 0L) theta <- theta.initial

  graph <- function() {
    methods::as(graph.pattern, "TsparseMatrix")
  }

  Q <- function() {
    D <- Dfunc(theta)
    methods::as(Matrix::t(D) %*% Matrix::Diagonal(x = Vinv) %*% D, "TsparseMatrix")
  }

  mu <- function() numeric(0)

  log.norm.const <- function() {
    ## Supply the V-dependent constant EXPLICITLY. (Returning numeric(0) makes
    ## INLA drop the theta-dependent 0.5*log|Q| term -- same failure as generic0 --
    ## which biases the hyperparameter posterior.)
    ##
    ## If the operator provides an analytic `lognc(theta, Vinv)` (e.g. SAR/CAR use
    ## sum(log|1 - rho*eigen(W)|), exact and stable near the singularity), use it;
    ## otherwise fall back to a generic determinant.
    if (exists("lognc", envir = envir, inherits = FALSE)) {
      return(lognc(theta, Vinv))
    }
    D  <- Dfunc(theta)
    Qm <- Matrix::t(D) %*% Matrix::Diagonal(x = Vinv) %*% D
    nn <- nrow(Qm)
    if (rankdef == 0L) {
      ld <- tryCatch(as.numeric(Matrix::determinant(Qm, logarithm = TRUE)$modulus),
                     error = function(e) -Inf)
      if (!is.finite(ld)) return(-1e12)               # singular -> repel the optimizer
      return(0.5 * ld - 0.5 * nn * log(2 * pi))
    }
    ## intrinsic model: generalized determinant (product of non-zero eigenvalues)
    r  <- nn - rankdef
    ev <- eigen(as.matrix(Qm), symmetric = TRUE, only.values = TRUE)$values
    ev <- sort(ev, decreasing = TRUE)[seq_len(r)]
    0.5 * sum(log(ev)) - 0.5 * r * log(2 * pi)
  }

  log.prior <- function() logprior(theta)

  initial <- function() theta.initial

  quit <- function() invisible()

  do.call(match.arg(cmd), args = list())
}

#' Build an INLA rgeneric model for an ngvb operator at a fixed V.
#'
#' @param op An operator descriptor from [ngvb_operator()].
#' @param V Numeric vector of mixing variables (length = nrow of the latent
#'   field). The precision uses `diag(1/V)`; pass `op$h` for the Gaussian model.
#' @return An object usable as `f(idx, model = <this>)` in an INLA formula.
#' @export
ngvb_rgeneric <- function(op, V = op$h) {
  .need_inla()
  Vinv <- 1 / V
  args <- list(
    ngvb.rgeneric.engine,
    Dfunc         = op$Dfunc,
    Vinv          = Vinv,
    rankdef       = as.integer(op$rankdef),
    ntheta        = as.integer(op$ntheta),
    theta.initial = op$theta.initial,
    logprior      = op$logprior,
    graph.pattern = op$graph
  )
  if (!is.null(op$lognc)) args$lognc <- op$lognc      # operator-supplied analytic norm-const
  do.call(INLA::inla.rgeneric.define, args)
}

#' Conditional precision matrix Q(theta, V) for an operator (main-process helper).
#'
#' Used for validation and for computing the discrepancy statistic d_i.
#' @param op Operator descriptor.
#' @param theta Hyperparameters (internal scale). Defaults to `op$theta.initial`.
#' @param V Mixing vector. Defaults to `op$h` (Gaussian model).
#' @return A sparse precision matrix `Q(theta, V) = D(theta)^T diag(1/V) D(theta)`.
#' @export
ngvb_precision <- function(op, theta = op$theta.initial, V = op$h) {
  D <- op$Dfunc(theta)
  Matrix::t(D) %*% Matrix::Diagonal(x = 1 / V) %*% D
}
