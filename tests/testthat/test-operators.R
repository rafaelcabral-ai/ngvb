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
  set.seed(1); loc <- matrix(runif(40), 20, 2)
  mesh <- fmesher::fm_mesh_2d(loc, max.edge = c(0.3, 0.6), cutoff = 0.1)
  spde <- INLA::inla.spde2.matern(mesh)
  op <- ngvb_operator("spde", spde = spde)
  expect_equal(op$h, Matrix::diag(spde$param.inla$M0), ignore_attr = TRUE)
  expect_false(all(op$h == 1))
})
