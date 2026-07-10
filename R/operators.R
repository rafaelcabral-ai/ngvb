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
#'   `"car"`, `"spde"`, `"ou"`, `"seasonal"`, `"generic0"`.
#' @param ... Model-specific arguments (e.g. `n` for rw/ar/iid, `W` for sar/car,
#'   `loc` for ou). For `"spde"`, pass `spde = INLA::inla.spde2.pcmatern(mesh,
#'   prior.range = , prior.sigma = )` (or `mesh` + `prior.range` + `prior.sigma`
#'   and ngvb2 builds it); the non-Gaussian SPDE uses that PC prior on the
#'   practical range and marginal SD. Plain `inla.spde2.matern()` fields are
#'   rejected — only the PC-prior parameterization is supported.
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
    seasonal = op_seasonal(...),
    generic0 = op_generic0(...),
    stop("ngvb2: operator type not implemented yet: ", type)
  )
}

## ---- seasonal: sum of `season` consecutive terms is white noise ------------
## D0 is the (n-season+1) x n sliding-sum operator, so D0^T D0 is INLA's seasonal
## structure matrix (rank n-(season-1)). One precision hyperparameter.

op_seasonal <- function(n, season, pc.prec = c(U = 1, alpha = 0.01)) {
  stopifnot(season >= 2, n > season)
  m  <- n - season + 1L
  D0 <- Matrix::sparseMatrix(
    i = rep(seq_len(m), each = season),
    j = as.integer(unlist(lapply(seq_len(m), function(r) r:(r + season - 1L)))),
    x = 1, dims = c(m, n))
  lp <- .pc_prec_logprior(pc.prec[["U"]], pc.prec[["alpha"]])
  list(type = "seasonal", n = n, ntheta = 1L, rankdef = season - 1L,
       h = rep(1, m), theta.initial = 4,
       Dfunc    = function(theta) sqrt(exp(theta[1L])) * D0,
       graph    = Matrix::crossprod(D0),
       logprior = function(theta) lp(theta[1L]), prec.logprior = lp)
}

## ---- generic0: user-supplied structure matrix C, precision Q = tau C --------
## Factor C = U diag(lambda) U^T and set D0 = diag(sqrt(lambda_+)) U_+^T over the
## positive eigenvalues, so D0^T D0 = C exactly and rank(D0) = rank(C). A single
## precision scales D = sqrt(tau) D0. Handles proper and intrinsic C alike.

op_generic0 <- function(C, pc.prec = c(U = 1, alpha = 0.01)) {
  C  <- as.matrix(C); n <- nrow(C)
  e  <- eigen((C + t(C)) / 2, symmetric = TRUE)
  tol <- max(abs(e$values)) * 1e-9
  if (min(e$values) < -tol)
    warning("ngvb2: generic0 Cmatrix has negative eigenvalues (min = ",
            signif(min(e$values), 3), "); a structure matrix must be positive ",
            "semi-definite. The negative part is dropped, so D^T D != C.", call. = FALSE)
  pos <- e$values > tol
  D0 <- methods::as(Matrix::Matrix(diag(sqrt(e$values[pos]), sum(pos)) %*%
                                   t(e$vectors[, pos, drop = FALSE])), "CsparseMatrix")
  ## Graph = sparsity of Q(theta, V) = D0^T diag(1/V) D0 over ALL V. D0 is the
  ## (generally dense) eigenvector factor, so this is NOT the sparsity of C:
  ## for V != h the cancellations that make C sparse no longer hold and Q fills
  ## in. Declare the structural union pattern of D0^T D0 (crossprod of D0's
  ## incidence), which the rgeneric engine needs to be a superset of every Q(V).
  inc   <- methods::as(abs(D0) > 0, "dsparseMatrix")
  graph <- methods::as(Matrix::crossprod(inc) > 0, "CsparseMatrix")
  lp <- .pc_prec_logprior(pc.prec[["U"]], pc.prec[["alpha"]])
  list(type = "generic0", n = n, ntheta = 1L, rankdef = n - sum(pos),
       h = rep(1, sum(pos)), theta.initial = 4,
       Dfunc    = function(theta) sqrt(exp(theta[1L])) * D0,
       graph    = graph,
       logprior = function(theta) lp(theta[1L]), prec.logprior = lp)
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
       theta.initial = c(1, 0), Dfunc = Dfunc, graph = graph, logprior = logprior,
       prec.logprior = lp)
}

## ---- SPDE / Matern (alpha = 2), PC-prior parameterization ------------------
## Built ONLY from an INLA::inla.spde2.pcmatern() object, so the prior is the
## penalized-complexity prior of Fuglstad et al. (2019) on the interpretable
## (practical range, marginal SD), exactly as INLA parameterizes it. We take
## theta = (log range, log sigma) -- the same internal coordinates INLA's pcmatern
## uses -- and map them to (kappa, tau) analytically:
##   kappa = sqrt(8 nu) / range,     nu = alpha - d/2
##   tau   : sigma^2 = Gamma(nu) / (Gamma(nu+d/2) (4 pi)^{d/2} kappa^{2nu} tau^2)
## The FEM factor D = tau (kappa^2 M0 + M1) (M0 mass-lumped, M2 = M1 M0^{-1} M1)
## then satisfies Q(V=h) = D^T diag(1/h) D = tau^2 (kappa^4 M0 + 2 kappa^2 M1 + M2),
## INLA's Matern precision (verified against inla.spde2.precision() to ~1e-16).
##
## The joint PC prior factorizes (range PC prior in dimension d, exponential PC
## prior on sigma):
##   pi(range) = (d/2) l_r range^{-d/2-1} exp(-l_r range^{-d/2}),  l_r = lambda_range
##   pi(sigma) = l_s exp(-l_s sigma),                              l_s = lambda_sigma
## with l_r, l_s read straight off the pcmatern object (hyper$theta1$param =
## c(l_r, l_s, d)); this reproduces P(range<range0)=p_r, P(sigma>sigma0)=p_s.

op_spde <- function(spde = NULL, mesh = NULL, prior.range = NULL, prior.sigma = NULL,
                    alpha = 2) {
  if (is.null(spde)) {
    if (is.null(mesh) || is.null(prior.range) || is.null(prior.sigma))
      stop("ngvb2: supply a PC-prior SPDE via spde = INLA::inla.spde2.pcmatern(...), or ",
           "give mesh + prior.range + prior.sigma so ngvb2 can build one.", call. = FALSE)
    .need_inla()
    spde <- INLA::inla.spde2.pcmatern(mesh, alpha = alpha,
                                      prior.range = prior.range, prior.sigma = prior.sigma)
  }
  if (!inherits(spde, "inla.spde2"))
    stop("ngvb2: `spde` must be an inla.spde2 object from INLA::inla.spde2.pcmatern().",
         call. = FALSE)
  hy <- spde$f$hyper
  if (is.null(hy$theta1$prior) || !identical(hy$theta1$prior, "pcmatern"))
    stop("ngvb2: the SPDE operator requires a PC-prior field built with ",
         "INLA::inla.spde2.pcmatern(). The supplied object uses a '",
         if (is.null(hy$theta1$prior)) "?" else hy$theta1$prior,
         "' prior; rebuild it with inla.spde2.pcmatern(mesh, prior.range = , prior.sigma = ).",
         call. = FALSE)

  pin <- spde$param.inla
  M0 <- methods::as(pin$M0, "CsparseMatrix")
  M1 <- methods::as(pin$M1, "CsparseMatrix")
  M2 <- methods::as(pin$M2, "CsparseMatrix")
  n  <- nrow(M0)

  ## PC-prior parameters, straight from the pcmatern object: param = c(l_r, l_s, d)
  pr <- hy$theta1$param
  l_r <- pr[1L]; l_s <- pr[2L]; d <- pr[3L]
  nu  <- alpha - d / 2
  if (nu <= 0)
    stop("ngvb2: alpha - d/2 must be positive (got alpha = ", alpha, ", d = ", d, ").",
         call. = FALSE)
  ## constant in log tau (see header): 0.5[logGamma(nu) - logGamma(nu+d/2) - (d/2)log 4pi]
  ctau <- 0.5 * (lgamma(nu) - lgamma(nu + d / 2) - (d / 2) * log(4 * pi))

  ## theta = (log range, log sigma); reuse INLA's own initial values
  theta0 <- c(if (is.null(hy$theta1$initial)) 0 else hy$theta1$initial,
              if (is.null(hy$theta2$initial)) 0 else hy$theta2$initial)

  Dfunc <- function(theta) {
    logkappa <- 0.5 * log(8 * nu) - theta[1L]                       # kappa = sqrt(8 nu)/range
    logtau   <- ctau - nu * logkappa - theta[2L]                    # from the sigma relation
    exp(logtau) * (exp(2 * logkappa) * M0 + M1)                     # tau (kappa^2 M0 + M1)
  }
  logprior <- function(theta) {
    r <- theta[1L]; s <- theta[2L]                                  # log range, log sigma
    ## range PC prior (dim d) and exponential sigma PC prior, each with its Jacobian
    lr <- log(d / 2) + log(l_r) - (d / 2) * r - l_r * exp(-(d / 2) * r)
    ls <- log(l_s) + s - l_s * exp(s)
    lr + ls
  }

  op <- list(type = "spde", n = n, ntheta = 2L, rankdef = 0L,
             h = Matrix::diag(M0),                                  # mass-lumped weights (h != 1)
             theta.initial = theta0, Dfunc = Dfunc,
             graph = methods::as(abs(M0) + abs(M1) + abs(M2), "CsparseMatrix"),
             logprior = logprior,
             pc = list(lambda.range = l_r, lambda.sigma = l_s, d = d, nu = nu, alpha = alpha))

  ## Self-check: our analytic precision must equal INLA's at the initial theta.
  ## (Guards against a non-standard alpha/manifold where the D factorization differs.)
  if (requireNamespace("INLA", quietly = TRUE)) {
    Qi <- tryCatch(INLA::inla.spde2.precision(spde, theta = theta0), error = function(e) NULL)
    if (!is.null(Qi)) {
      D  <- op$Dfunc(theta0)
      Qm <- Matrix::t(D) %*% Matrix::Diagonal(x = 1 / op$h) %*% D
      rel <- max(abs(as.matrix(Qi - Qm))) / max(abs(as.matrix(Qi)))
      if (is.finite(rel) && rel > 1e-6)
        stop("ngvb2: SPDE precision reconstruction disagrees with INLA (rel. diff ",
             signif(rel, 3), "). Only the standard alpha = 2 Matern on a flat mesh is ",
             "supported for the non-Gaussian SPDE extension.", call. = FALSE)
    }
  }
  op
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
       lognc = lognc, logprior = logprior, prec.logprior = lp)
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
         logprior = function(theta) lp(theta[1L]), prec.logprior = lp)
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
         lognc = lognc, logprior = logprior, prec.logprior = lp)
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
       logprior = function(theta) lp(theta[1L]), prec.logprior = lp)
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
       logprior = function(theta) lp(theta[1L]), prec.logprior = lp)
}

op_iid <- function(n, pc.prec = c(U = 1, alpha = 0.01)) {
  stopifnot(n >= 1)
  D0 <- Matrix::Diagonal(n)
  lp <- .pc_prec_logprior(pc.prec[["U"]], pc.prec[["alpha"]])
  list(type = "iid", n = n, ntheta = 1L, rankdef = 0L,
       h = rep(1, n), theta.initial = 4,
       Dfunc    = function(theta) sqrt(exp(theta[1L])) * D0,
       graph    = D0,
       logprior = function(theta) lp(theta[1L]), prec.logprior = lp)
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

  U2 <- pc.cor0[["U"]]; a2 <- pc.cor0[["alpha"]]
  ## PC prior on the marginal precision (theta[1] = log tau); shared closure so a
  ## user prior from f() can replace exactly this piece (see ngvb_apply_user_prec_prior).
  lp <- .pc_prec_logprior(pc.prec[["U"]], pc.prec[["alpha"]])
  logprior <- function(theta) {
    rho.i <- theta[2L]
    phi   <- 2 * exp(rho.i) / (1 + exp(rho.i)) - 1
    ## PC prior on the lag-1 correlation toward rho = 0 (pc.cor0)
    th      <- -log(a2) / sqrt(-log(1 - U2^2))
    s       <- sqrt(-log(1 - phi^2))
    lp.phi  <- -log(2) + log(th) - th * s + log(abs(phi)) - log(1 - phi^2) - 0.5 * log(s^2)
    ## d phi/d theta Jacobian: phi = 2*plogis(theta)-1 => log|dphi/dtheta| = log2 + rho.i - 2 log(1+e^rho.i)
    jac.phi <- log(2) + rho.i - 2 * log(1 + exp(rho.i))
    lp(theta[1L]) + lp.phi + jac.phi
  }

  list(type = "ar1", n = n, ntheta = 2L, rankdef = 0L,
       h = rep(1, n), theta.initial = c(1, 1),
       Dfunc = Dfunc, graph = graph, logprior = logprior, prec.logprior = lp)
}
