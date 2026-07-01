## ---------------------------------------------------------------------------
## Operator registry.
##
## An "operator descriptor" is a plain list describing one latent model:
##   type          model name
##   n             latent dimension
##   ntheta        number of hyperparameters (internal scale)
##   rankdef       rank deficiency of D^T D (0 = proper)
##   h             constant weight vector (E[V] = h); Gaussian model uses V = h
##   theta.initial starting hyperparameters (internal scale)
##   Dfunc(theta)  builds D(theta) as a sparse Matrix (Matrix:: qualified, so it
##                 also works inside the rgeneric subprocess)
##   graph         structural sparsity pattern of D^T D (theta-independent)
##   logprior(theta) log prior on the internal scale (incl. Jacobians)
## ---------------------------------------------------------------------------

#' Construct an ngvb operator descriptor.
#'
#' @param type Model type: one of `"iid"`, `"rw1"`, `"rw2"`, `"ar1"`, `"sar"`,
#'   `"car"`, `"spde"`, `"ou"`.
#' @param ... Model-specific arguments (e.g. `n` for rw/ar/iid, `W` for sar/car,
#'   `spde` for spde, `loc` for ou).
#' @return An operator descriptor: a list with `Dfunc(theta)`, the constant
#'   vector `h`, `rankdef`, and prior/graph metadata.
#' @examples
#' op <- ngvb_operator("rw1", n = 10)
#' dim(op$Dfunc(0))       # D(theta): the 9 x 10 first-difference operator
#' ngvb_precision(op)     # Q = D^T diag(1/h) D  (the RW1 structure matrix)
#' @export
ngvb_operator <- function(type, ...) {
  switch(type,
    ar1 = op_ar1(...),
    rw1 = op_rw1(...),
    rw2 = op_rw2(...),
    iid = op_iid(...),
    sar = op_sar(...),
    car = op_car(...),
    spde = op_spde(...),
    ou  = op_ou(...),
    stop("ngvb2: operator type not implemented yet: ", type)
  )
}

## ---- Ornstein-Uhlenbeck: continuous-time AR1 for irregularly spaced times ---
## dX = -kappa X dt + sigma dB; marginal precision tau, rho_i = exp(-kappa dt_i).
## D is lower-bidiagonal with the exact transition-variance scaling, so D^T D is
## the OU precision (marginal covariance (1/tau) exp(-kappa |t_i - t_j|)).
## theta = (log tau, log kappa).

op_ou <- function(loc, pc.prec = c(U = 1, alpha = 0.01)) {
  loc <- sort(loc); n <- length(loc); dt <- diff(loc)
  stopifnot(n >= 2, all(dt > 0))
  lp <- .pc_prec_logprior(pc.prec[["U"]], pc.prec[["alpha"]])
  Dfunc <- function(theta) {
    tau <- exp(theta[1L]); kappa <- exp(theta[2L])
    rho <- exp(-kappa * dt)                       # length n-1
    s   <- sqrt(tau / (1 - rho^2))                # transition scaling
    Matrix::sparseMatrix(
      i = c(1:n, 2:n), j = c(1:n, 1:(n - 1L)),
      x = c(sqrt(tau), s, -rho * s), dims = c(n, n))
  }
  graph <- Matrix::bandSparse(n, n, k = c(-1L, 0L, 1L),
                              diagonals = list(rep(1, n - 1L), rep(1, n), rep(1, n - 1L)))
  logprior <- function(theta) lp(theta[1L]) + stats::dnorm(theta[2L], 0, 3, log = TRUE)
  list(type = "ou", n = n, ntheta = 2L, rankdef = 0L, h = rep(1, n),
       theta.initial = c(1, 0), Dfunc = Dfunc, graph = graph, logprior = logprior)
}

## ---- SPDE / Matern (alpha = 2):  D = kappa^2 C + G  (C, G = FEM mass, stiffness)
## h = diag(C) (mass-lumped weights), so Q(V=h) = D^T C^{-1} D = the Matern precision
## tau (kappa^4 M0 + 2 kappa^2 M1 + M2). theta = (log tau, log kappa^2). Proper (rankdef 0),
## well-conditioned for kappa>0, so the generic Cholesky determinant is fine (no lognc).

op_spde <- function(spde, mesh = NULL) {
  if (is.null(spde) && !is.null(mesh)) spde <- INLA::inla.spde2.matern(mesh)
  pin <- spde$param.inla
  M0 <- methods::as(pin$M0, "CsparseMatrix")
  M1 <- methods::as(pin$M1, "CsparseMatrix")
  M2 <- methods::as(pin$M2, "CsparseMatrix")
  n  <- nrow(M0)
  Dfunc <- function(theta) {
    tau <- exp(theta[1L]); k2 <- exp(theta[2L])
    sqrt(tau) * (k2 * M0 + M1)
  }
  list(type = "spde", n = n, ntheta = 2L, rankdef = 0L,
       h = Matrix::diag(M0),                                   # the only model with h != 1
       theta.initial = c(0, 0),
       Dfunc = Dfunc,
       graph = methods::as(abs(M0) + abs(M1) + abs(M2), "CsparseMatrix"),
       logprior = function(theta)
         stats::dnorm(theta[1L], 0, 3, log = TRUE) + stats::dnorm(theta[2L], 0, 3, log = TRUE))
}

## ---- SAR: simultaneous autoregression  D = I - rho W  (proper, full rank) ---
## Q = tau * (I - rho W)^T diag(1/V) (I - rho W). W is a (row-standardized)
## adjacency. theta = (log tau, rho-internal); rho = 2*plogis(theta2) - 1.

op_sar <- function(W, pc.prec = c(U = 1, alpha = 0.01)) {
  W <- methods::as(W, "CsparseMatrix")
  n <- nrow(W)
  eigenv <- Re(eigen(W, only.values = TRUE)$values)            # real for row-standardized W
  lp <- .pc_prec_logprior(pc.prec[["U"]], pc.prec[["alpha"]])
  ## rho in (0, 1) via the logistic map (positive spatial correlation; stable
  ## away from the singular boundary rho = 1/lambda_max). Uniform(0,1) prior.
  Dfunc <- function(theta) {
    prec <- exp(theta[1L])
    rho  <- exp(theta[2L]) / (1 + exp(theta[2L]))
    sqrt(prec) * (Matrix::Diagonal(n) - rho * W)
  }
  ## structural sparsity of D^T D (union over rho; rho = 0.5 gives the full pattern)
  graph <- Matrix::crossprod(Matrix::Diagonal(n) - 0.5 * W)
  ## analytic, stable log normalizing constant:
  ## 0.5 log|Q| - 0.5 n log(2pi), log|Q| = n log tau + 2 sum log(1-rho*lambda) - sum log V
  lognc <- function(theta, Vinv) {
    rho <- exp(theta[2L]) / (1 + exp(theta[2L]))
    0.5 * n * theta[1L] + sum(log(1 - rho * eigenv)) +
      0.5 * sum(log(Vinv)) - 0.5 * n * log(2 * pi)
  }
  logprior <- function(theta) {
    rho <- exp(theta[2L]) / (1 + exp(theta[2L]))
    lp(theta[1L]) + log(rho) + log(1 - rho)                     # Uniform(0,1) prior + Jacobian
  }
  list(type = "sar", n = n, ntheta = 2L, rankdef = 0L, h = rep(1, n),
       theta.initial = c(0, 0), Dfunc = Dfunc, graph = graph,
       lognc = lognc, logprior = logprior)
}

## ---- CAR (proper, besagproper-style):  D = I - rho B,  B = row-standardized
## adjacency from a graph. Full rank for |rho| < 1. For the intrinsic ICAR use
## op_car(..., intrinsic = TRUE): D = the edge-incidence operator (rankdef 1).

op_car <- function(W, pc.prec = c(U = 1, alpha = 0.01), intrinsic = FALSE) {
  W <- methods::as(W, "CsparseMatrix")
  n <- nrow(W)
  lp <- .pc_prec_logprior(pc.prec[["U"]], pc.prec[["alpha"]])
  if (intrinsic) {
    ## edge-incidence D0: one row per undirected edge (i<j), +1 at i, -1 at j
    ij <- Matrix::which(Matrix::tril(W != 0, -1L), arr.ind = TRUE)
    ne <- nrow(ij)
    D0 <- Matrix::sparseMatrix(i = rep(seq_len(ne), 2),
                               j = c(ij[, 1], ij[, 2]),
                               x = rep(c(1, -1), each = ne), dims = c(ne, n))
    list(type = "car", n = n, ntheta = 1L, rankdef = 1L, h = rep(1, ne),
         theta.initial = 4,
         Dfunc = function(theta) sqrt(exp(theta[1L])) * D0,
         graph = Matrix::crossprod(D0),
         logprior = function(theta) lp(theta[1L]))
  } else {
    eigenv <- Re(eigen(W, only.values = TRUE)$values)
    Dfunc <- function(theta) {
      prec <- exp(theta[1L])
      rho  <- exp(theta[2L]) / (1 + exp(theta[2L]))
      sqrt(prec) * (Matrix::Diagonal(n) - rho * W)
    }
    graph <- Matrix::crossprod(Matrix::Diagonal(n) - 0.5 * W)
    lognc <- function(theta, Vinv) {
      rho <- exp(theta[2L]) / (1 + exp(theta[2L]))
      0.5 * n * theta[1L] + sum(log(1 - rho * eigenv)) +
        0.5 * sum(log(Vinv)) - 0.5 * n * log(2 * pi)
    }
    logprior <- function(theta) {
      rho <- exp(theta[2L]) / (1 + exp(theta[2L]))
      lp(theta[1L]) + log(rho) + log(1 - rho)
    }
    list(type = "car", n = n, ntheta = 2L, rankdef = 0L, h = rep(1, n),
         theta.initial = c(0, 0), Dfunc = Dfunc, graph = graph,
         lognc = lognc, logprior = logprior)
  }
}

## ---- intrinsic random walks (RW1 / RW2) and iid -----------------------------
## D is the order-k difference operator ((n-k) x n), so D^T D is the RW(k)
## structure matrix (rank n-k). A single precision hyperparameter; D = sqrt(tau) D0.

.pc_prec_logprior <- function(U, alpha) {
  lambda <- -log(alpha) / U
  function(lprec) {                       # log prior of theta = log(tau), incl. Jacobian
    tau <- exp(lprec)
    log(lambda / 2) - 1.5 * lprec - lambda / sqrt(tau) + lprec
  }
}

op_rw1 <- function(n, pc.prec = c(U = 1, alpha = 0.01)) {
  stopifnot(n >= 3)
  D0 <- Matrix::bandSparse(n - 1L, n, k = c(0L, 1L),
                           diagonals = list(rep(-1, n - 1L), rep(1, n - 1L)))
  lp <- .pc_prec_logprior(pc.prec[["U"]], pc.prec[["alpha"]])
  list(type = "rw1", n = n, ntheta = 1L, rankdef = 1L,
       h = rep(1, n - 1L), theta.initial = 4,
       Dfunc    = function(theta) sqrt(exp(theta[1L])) * D0,
       graph    = Matrix::crossprod(D0),
       logprior = function(theta) lp(theta[1L]))
}

op_rw2 <- function(n, pc.prec = c(U = 1, alpha = 0.01)) {
  stopifnot(n >= 4)
  D0 <- Matrix::bandSparse(n - 2L, n, k = c(0L, 1L, 2L),
                           diagonals = list(rep(1, n - 2L), rep(-2, n - 2L), rep(1, n - 2L)))
  lp <- .pc_prec_logprior(pc.prec[["U"]], pc.prec[["alpha"]])
  list(type = "rw2", n = n, ntheta = 1L, rankdef = 2L,
       h = rep(1, n - 2L), theta.initial = 4,
       Dfunc    = function(theta) sqrt(exp(theta[1L])) * D0,
       graph    = Matrix::crossprod(D0),
       logprior = function(theta) lp(theta[1L]))
}

op_iid <- function(n, pc.prec = c(U = 1, alpha = 0.01)) {
  stopifnot(n >= 1)
  D0 <- Matrix::Diagonal(n)
  lp <- .pc_prec_logprior(pc.prec[["U"]], pc.prec[["alpha"]])
  list(type = "iid", n = n, ntheta = 1L, rankdef = 0L,
       h = rep(1, n), theta.initial = 4,
       Dfunc    = function(theta) sqrt(exp(theta[1L])) * D0,
       graph    = D0,
       logprior = function(theta) lp(theta[1L]))
}

## ---- AR1 -------------------------------------------------------------------
## x_i - rho x_{i-1} = innovation; D is lower-bidiagonal with the stationary
## first-row scaling, so that D^T D is the AR1 precision (marginal precision tau).

op_ar1 <- function(n, pc.prec = c(U = 1, alpha = 0.01),
                      pc.cor0 = c(U = 0.5, alpha = 0.5)) {
  stopifnot(n >= 2)

  Dfunc <- function(theta) {
    prec  <- exp(theta[1L])
    rho   <- 2 * exp(theta[2L]) / (1 + exp(theta[2L])) - 1
    preci <- prec / (1 - rho^2)
    D <- Matrix::bandSparse(n, n, k = c(-1L, 0L),
                            diagonals = list(rep(-rho, n - 1L), rep(1, n)))
    D[1L, 1L] <- sqrt(1 - rho^2)
    sqrt(preci) * D
  }

  graph <- Matrix::bandSparse(n, n, k = c(-1L, 0L, 1L),
                              diagonals = list(rep(1, n - 1L), rep(1, n), rep(1, n - 1L)))

  U1 <- pc.prec[["U"]]; a1 <- pc.prec[["alpha"]]
  U2 <- pc.cor0[["U"]]; a2 <- pc.cor0[["alpha"]]
  logprior <- function(theta) {
    lprec <- theta[1L]; tau <- exp(lprec)
    rho.i <- theta[2L]
    phi   <- 2 * exp(rho.i) / (1 + exp(rho.i)) - 1
    ## PC prior on the marginal precision (pc.prec), with d tau/d theta Jacobian
    lambda  <- -log(a1) / U1
    lp.tau  <- log(lambda / 2) - 1.5 * lprec - lambda / sqrt(tau) + lprec
    ## PC prior on the lag-1 correlation toward rho = 0 (pc.cor0)
    th      <- -log(a2) / sqrt(-log(1 - U2^2))
    s       <- sqrt(-log(1 - phi^2))
    lp.phi  <- -log(2) + log(th) - th * s + log(abs(phi)) - log(1 - phi^2) - 0.5 * log(s^2)
    ## d phi/d theta Jacobian: phi = 2*plogis(theta)-1 => log|dphi/dtheta| = log2 + rho.i - 2 log(1+e^rho.i)
    jac.phi <- log(2) + rho.i - 2 * log(1 + exp(rho.i))
    lp.tau + lp.phi + jac.phi
  }

  list(type = "ar1", n = n, ntheta = 2L, rankdef = 0L,
       h = rep(1, n), theta.initial = c(1, 1),
       Dfunc = Dfunc, graph = graph, logprior = logprior)
}
