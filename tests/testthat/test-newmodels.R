test_that("generic0 operator reconstructs the structure matrix exactly", {
  set.seed(1); m <- 8; B <- matrix(rnorm(m * m), m); C <- crossprod(B) + diag(m)
  op <- ngvb_operator("generic0", C = C)
  expect_equal(op$rankdef, 0L)
  D0 <- op$Dfunc(0)
  expect_equal(as.matrix(Matrix::crossprod(D0)), C, tolerance = 1e-8)
})

test_that("generic0 handles a rank-deficient (intrinsic) structure matrix", {
  D1 <- diff(diag(6)); C <- crossprod(D1)          # RW1 structure, rank 5
  op <- ngvb_operator("generic0", C = C)
  expect_equal(op$rankdef, 1L)
  expect_equal(as.matrix(Matrix::crossprod(op$Dfunc(0))), as.matrix(C), tolerance = 1e-8)
})

test_that("seasonal operator matches INLA's seasonal structure (short period)", {
  skip_if_not_installed("INLA")
  n <- 12; s <- 4
  op <- ngvb_operator("seasonal", n = n, season = s)
  Qours <- as.matrix(Matrix::crossprod(op$Dfunc(0)))
  set.seed(1); y <- rnorm(n)
  fit <- INLA::inla(y ~ -1 + f(t, model = "seasonal", season.length = s, constr = FALSE,
                               hyper = list(prec = list(initial = 0, fixed = TRUE))),
                    data = data.frame(y = y, t = 1:n),
                    control.compute = list(config = TRUE))
  Qinla <- as.matrix(fit$misc$configs$config[[1]]$Qprior)[1:n, 1:n]
  Qinla <- Qinla + t(Qinla) - diag(diag(Qinla))   # INLA stores Qprior upper-triangular
  expect_equal(Qours, Qinla, tolerance = 1e-6)
})

test_that("seasonal component is auto-detected end to end", {
  skip_if_not_installed("INLA")
  n <- 36; set.seed(1)
  y <- as.numeric(arima.sim(list(ar = 0.3), n)) + rep(sin(2 * pi * (1:12) / 12), length.out = n)
  LGM <- INLA::inla(y ~ -1 + f(t, model = "seasonal", season.length = 12),
                    data = data.frame(y = y, t = 1:n), control.compute = list(config = TRUE))
  op <- ngvb_detect_operator(LGM, "t")
  expect_identical(op$type, "seasonal")
  expect_identical(op$rankdef, 11L)
})

test_that("bayes.factor returns a well-formed estimate and favours the right model", {
  skip_if_not_installed("INLA")
  ## clear non-Gaussian iid signal: most effects ~ N(0,1), one strong outlier
  set.seed(1); ng <- 40; g <- rnorm(ng); g[10] <- 8
  y <- g + rnorm(ng, sd = 0.3)
  LGM <- INLA::inla(y ~ -1 + f(i, model = "iid"), data = data.frame(y = y, i = 1:ng),
                    control.compute = list(config = TRUE))
  LnGM <- ngvb(LGM, iter = 12, verbose = FALSE)
  bf <- bayes.factor(LnGM, n.samples = 25, seed = 1)
  expect_true(is.finite(bf$log10BF))
  expect_true(bf$BF > 0)
  expect_true(bf$ess >= 1 && bf$ess <= 25)
  expect_true(bf$log10BF > 0.3)        # substantial evidence for the outlier model

  ## Gaussian null: no outlier -> should not favour the non-Gaussian model
  set.seed(2); y0 <- rnorm(ng)
  LGM0 <- INLA::inla(y ~ -1 + f(i, model = "iid"), data = data.frame(y = y0, i = 1:ng),
                     control.compute = list(config = TRUE))
  bf0 <- bayes.factor(ngvb(LGM0, iter = 12, verbose = FALSE), n.samples = 25, seed = 1)
  expect_true(bf0$log10BF < bf$log10BF)
})
