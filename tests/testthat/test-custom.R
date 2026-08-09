## Custom operators (ngvb_custom) and auto-detection helpers. Fast (no INLA fit).

test_that("ngvb_custom builds a valid operator and captures its data", {
  Wt <- matrix(c(0, 1, 1, 0), 2, 2)
  m  <- 2
  op <- ngvb_custom(
    D = function(theta) sqrt(exp(theta[1])) * (Matrix::Diagonal(m) - 0.3 * Wt),
    h = rep(1, m), ntheta = 1)
  expect_equal(op$type, "custom")
  expect_equal(op$ntheta, 1L)
  expect_equal(op$n, 2L)
  expect_length(op$h, 2L)
  ## the data D() uses must be captured into a non-global env so it serialises
  ## into the rgeneric subprocess
  expect_false(identical(environment(op$Dfunc), globalenv()))
  expect_true(all(c("Wt", "m") %in% ls(environment(op$Dfunc))))
})

test_that("ngvb_custom recycles a scalar h and validates its length", {
  D3 <- Matrix::Diagonal(3)
  op <- ngvb_custom(D = function(theta) D3, h = 1, ntheta = 1)
  expect_equal(op$h, rep(1, 3))
  expect_error(ngvb_custom(D = function(theta) D3, h = c(1, 2), ntheta = 1),
               "length")
})

test_that("ngvb_custom SAR operator equals tau (I - rho W)^T (I - rho W)", {
  set.seed(1); n <- 6
  A <- matrix(stats::rbinom(n * n, 1, 0.4), n, n); A <- (A + t(A)) > 0; diag(A) <- FALSE
  W <- diag(1 / pmax(rowSums(A), 1)) %*% A
  op <- ngvb_custom(
    D = function(theta) sqrt(exp(theta[1])) *
      (Matrix::Diagonal(n) - (exp(theta[2]) / (1 + exp(theta[2]))) * W),
    h = rep(1, n), ntheta = 2)
  th  <- c(log(1.5), 0.3); rho <- exp(0.3) / (1 + exp(0.3))
  Q   <- as.matrix(ngvb_precision(op, theta = th, V = op$h))
  D   <- diag(n) - rho * W
  expect_equal(Q, 1.5 * t(D) %*% D, tolerance = 1e-8, ignore_attr = TRUE)
})

test_that("ngvb_find_f_model extracts the model expression from a formula", {
  f <- y ~ 1 + z + f(s, model = spde) + f(t, model = "ar1")
  expect_equal(deparse(ngvb:::ngvb_find_f_model(f, "s")), "spde")
  expect_equal(deparse(ngvb:::ngvb_find_f_model(f, "t")), "\"ar1\"")
  expect_null(ngvb:::ngvb_find_f_model(f, "nope"))
})

test_that("SPDE operator uses lumped mass for h and factors the Matern precision", {
  skip_if_not_installed("fmesher")
  skip_if_no_inla()
  set.seed(1); loc <- matrix(runif(40), 20, 2)
  mesh <- fmesher::fm_mesh_2d(loc, max.edge = c(0.3, 0.6), cutoff = 0.1)
  spde <- INLA::inla.spde2.pcmatern(mesh, prior.range = c(0.5, 0.5),
                                    prior.sigma = c(1, 0.01))
  op   <- ngvb_operator("spde", spde = spde)
  expect_true(Matrix::isDiagonal(spde$param.inla$M0))          # lumped mass -> h = diag(M0)
  expect_equal(op$h, Matrix::diag(spde$param.inla$M0), ignore_attr = TRUE)
  ## theta = (log range, log sigma); Q = tau^2 (kappa^4 M0 + 2 kappa^2 M1 + M2)
  range <- 0.5; sigma <- 1.2
  kappa <- sqrt(8) / range; tau <- 1 / (sqrt(4 * pi) * kappa * sigma)
  Q  <- as.matrix(ngvb_precision(op, theta = c(log(range), log(sigma)), V = op$h))
  M0 <- spde$param.inla$M0; M1 <- spde$param.inla$M1; M2 <- spde$param.inla$M2
  expect_equal(Q, as.matrix(tau^2 * (kappa^4 * M0 + 2 * kappa^2 * M1 + M2)),
               tolerance = 1e-8, ignore_attr = TRUE)
})
