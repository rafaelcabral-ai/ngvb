## Operator correctness: Q(theta, V=h) = D^T diag(1/h) D must equal the analytic
## precision for each model. Fast (no INLA fit).

test_that("AR1 operator equals the analytic AR1 precision", {
  n <- 12; prec <- 2; rho <- 0.6
  op <- ngvb_operator("ar1", n = n)
  theta <- c(log(prec), qlogis((rho + 1) / 2))
  Q  <- as.matrix(ngvb_precision(op, theta = theta, V = op$h))
  Sig <- outer(1:n, 1:n, function(i, j) rho^abs(i - j)) / prec
  expect_equal(Q, solve(Sig), tolerance = 1e-7, ignore_attr = TRUE)
})

test_that("RW1 / RW2 operators equal the difference-operator structure", {
  n <- 15
  op1 <- ngvb_operator("rw1", n = n)
  expect_equal(as.matrix(ngvb_precision(op1, theta = 0, V = op1$h)),
               as.matrix(crossprod(diff(diag(n), differences = 1))),
               tolerance = 1e-10, ignore_attr = TRUE)
  expect_equal(op1$rankdef, 1L)
  op2 <- ngvb_operator("rw2", n = n)
  expect_equal(as.matrix(ngvb_precision(op2, theta = 0, V = op2$h)),
               as.matrix(crossprod(diff(diag(n), differences = 2))),
               tolerance = 1e-10, ignore_attr = TRUE)
  expect_equal(op2$rankdef, 2L)
})

test_that("SAR operator equals tau (I - rho W)^T (I - rho W)", {
  set.seed(1); n <- 10
  A <- matrix(rbinom(n * n, 1, 0.3), n, n); A[lower.tri(A)] <- t(A)[lower.tri(A)]; diag(A) <- 0
  W <- diag(1 / pmax(rowSums(A), 1)) %*% A
  op <- ngvb_operator("sar", W = W)
  tau <- 1.5; rho <- plogis(0.4)
  theta <- c(log(tau), 0.4)
  Q   <- as.matrix(ngvb_precision(op, theta = theta, V = op$h))
  D   <- diag(n) - rho * W
  expect_equal(Q, tau * t(D) %*% D, tolerance = 1e-8, ignore_attr = TRUE)
})

test_that("OU operator equals the Ornstein-Uhlenbeck precision (irregular times)", {
  set.seed(1); loc <- sort(runif(15, 0, 10)); tau <- 1.5; kappa <- 0.4
  op <- ngvb_operator("ou", loc = loc)
  Q  <- as.matrix(ngvb_precision(op, theta = c(log(tau), log(kappa)), V = op$h))
  Sig <- (1 / tau) * exp(-kappa * abs(outer(loc, loc, "-")))   # OU covariance
  expect_equal(Q, solve(Sig), tolerance = 1e-6, ignore_attr = TRUE)
})

test_that("h is 1 for discrete models and diag(C) for SPDE", {
  expect_true(all(ngvb_operator("ar1", n = 8)$h == 1))
  expect_true(all(ngvb_operator("rw1", n = 8)$h == 1))
  skip_if_not_installed("fmesher")
  skip_if_no_inla()
  set.seed(1); loc <- matrix(runif(40), 20, 2)
  mesh <- fmesher::fm_mesh_2d(loc, max.edge = c(0.3, 0.6), cutoff = 0.1)
  spde <- INLA::inla.spde2.pcmatern(mesh, prior.range = c(0.3, 0.5),
                                    prior.sigma = c(1, 0.01))
  op <- ngvb_operator("spde", spde = spde)
  expect_equal(op$h, Matrix::diag(spde$param.inla$M0), ignore_attr = TRUE)
  expect_false(all(op$h == 1))
  ## PC-prior parameters read straight off the pcmatern object
  expect_equal(op$pc$lambda.range, -log(0.5) * 0.3, tolerance = 1e-6)
  expect_equal(op$pc$lambda.sigma, -log(0.01) / 1,  tolerance = 1e-6)
})

test_that("SPDE operator matches inla.spde2.precision and reproduces its PC prior", {
  skip_if_not_installed("fmesher")
  skip_if_no_inla()
  set.seed(1); loc <- matrix(runif(60), 30, 2)
  mesh <- fmesher::fm_mesh_2d(loc, max.edge = 0.3, cutoff = 0.05)
  spde <- INLA::inla.spde2.pcmatern(mesh, alpha = 2,
                                    prior.range = c(0.3, 0.5), prior.sigma = c(1, 0.01))
  op <- ngvb_operator("spde", spde = spde)
  ## the D-factored precision equals INLA's, at several (log range, log sigma)
  for (rs in list(c(0.5, 0.8), c(0.2, 1.5))) {
    th <- c(log(rs[1]), log(rs[2]))
    Qi <- INLA::inla.spde2.precision(spde, theta = th)
    D  <- op$Dfunc(th)
    Qm <- Matrix::t(D) %*% Matrix::Diagonal(x = 1 / op$h) %*% D
    expect_lt(max(abs(as.matrix(Qi - Qm))) / max(abs(as.matrix(Qi))), 1e-8)
  }
  ## the PC prior integrates to 1 and reproduces the range/sigma tail statements
  g <- seq(log(1e-3), log(60), length.out = 600); dg <- diff(g)[1]
  lp <- outer(g, g, Vectorize(function(a, b) op$logprior(c(a, b))))
  expect_equal(sum(exp(lp)) * dg * dg, 1, tolerance = 0.02)
  Prange <- sum((rowSums(exp(lp)) * dg)[exp(g) < 0.3]) * dg
  Psigma <- sum((colSums(exp(lp)) * dg)[exp(g) > 1.0]) * dg
  expect_equal(Prange, 0.5,  tolerance = 0.02)
  expect_equal(Psigma, 0.01, tolerance = 0.01)
})

test_that("SPDE operator rejects a non-PC (plain matern) object", {
  skip_if_not_installed("fmesher")
  skip_if_no_inla()
  set.seed(1); loc <- matrix(runif(40), 20, 2)
  mesh <- fmesher::fm_mesh_2d(loc, max.edge = 0.4, cutoff = 0.1)
  expect_error(ngvb_operator("spde", spde = INLA::inla.spde2.matern(mesh)),
               "pcmatern")
})
