test_that("loggamma precision log-prior matches Gamma(shape,rate) with the log-Jacobian", {
  a <- 0.7; b <- 0.03; lp <- ngvb2:::.loggamma_prec_logprior(a, b)
  for (theta in c(-2, 0, 1.3)) {
    ref <- stats::dgamma(exp(theta), shape = a, rate = b, log = TRUE) + theta
    expect_equal(lp(theta), ref)
  }
})

test_that("normal precision log-prior is Gaussian on theta with (mean, precision) params", {
  lp <- ngvb2:::.normal_prec_logprior(0.5, 4)          # sd = 1/2 on the log-precision scale
  for (theta in c(-1, 0.5, 2)) {
    expect_equal(lp(theta), stats::dnorm(theta, 0.5, 0.5, log = TRUE))
  }
})

test_that("prior-family dispatch maps pc.prec/loggamma/normal and rejects the rest", {
  expect_true(ngvb2:::.prec_logprior_from_hyper(list(prior = "pc.prec",  param = c(1, 0.01)))$mappable)
  expect_true(ngvb2:::.prec_logprior_from_hyper(list(prior = "loggamma", param = c(1, 5e-5)))$mappable)
  expect_true(ngvb2:::.prec_logprior_from_hyper(list(prior = "normal",   param = c(0, 1)))$mappable)
  expect_false(ngvb2:::.prec_logprior_from_hyper(list(prior = "flat",    param = numeric(0)))$mappable)
  expect_false(ngvb2:::.prec_logprior_from_hyper(list(prior = "loggamma", param = c(1, 1),
                                                      fixed = TRUE))$mappable)
})

test_that("hyper specs are read from fit$all.hyper, including INLA defaults", {
  skip_if_not_installed("INLA")
  set.seed(2); n <- 20
  d <- data.frame(y = rnorm(n), i = 1:n, j = 1:n)
  fit <- INLA::inla(y ~ -1 +
                      f(i, model = "iid",
                        hyper = list(prec = list(prior = "loggamma", param = c(0.1, 0.01)))) +
                      f(j, model = "iid"),
                    data = d, control.compute = list(config = TRUE))
  hy <- ngvb2:::ngvb_hyper_from_fit(fit, c("i", "j"))
  pi <- ngvb2:::.find_prec_entry(hy$i)$spec
  pj <- ngvb2:::.find_prec_entry(hy$j)$spec
  expect_equal(pi$prior, "loggamma", ignore_attr = TRUE)
  expect_equal(pi$param, c(0.1, 0.01), ignore_attr = TRUE)
  expect_equal(pj$prior, "loggamma", ignore_attr = TRUE)   # INLA's iid default is recorded too
  expect_equal(unname(pj$param[1]), 1, ignore_attr = TRUE)
})

test_that("applying a mapped prior replaces exactly the precision piece of logprior", {
  ## single-hyperparameter operator: logprior becomes the mapped prior outright
  op   <- ngvb_operator("rw1", n = 8)
  hy   <- list(prec = list(prior = "loggamma", param = c(0.1, 0.01)))
  op2  <- suppressMessages(ngvb2:::ngvb_apply_user_prec_prior(op, hy, "s", verbose = FALSE))
  want <- ngvb2:::.loggamma_prec_logprior(0.1, 0.01)
  for (th in c(-1, 0.7, 2)) expect_equal(op2$logprior(th), want(th))

  ## two-hyperparameter operator (ar1): the correlation prior must survive untouched
  op.ar  <- ngvb_operator("ar1", n = 12)
  op.ar2 <- suppressMessages(ngvb2:::ngvb_apply_user_prec_prior(
    op.ar, list(prec = list(prior = "pc.prec", param = c(2, 0.01))), "t", verbose = FALSE))
  th <- c(0.3, 0.4); newp <- ngvb2:::.pc_prec_logprior(2, 0.01)
  expect_equal(op.ar2$logprior(th),
               op.ar$logprior(th) - ngvb2:::.pc_prec_logprior(1, 0.01)(th[1]) + newp(th[1]))
})

test_that("an unmappable family and a secondary-hyperparameter prior both warn and fall back", {
  op0 <- ngvb_operator("iid", n = 5)
  expect_warning(
    op1 <- ngvb2:::ngvb_apply_user_prec_prior(
      op0, list(prec = list(prior = "flat", param = numeric(0))), "z", verbose = FALSE),
    "was dropped")
  expect_identical(op1$logprior(0.5), op0$logprior(0.5))   # default PC prior retained

  expect_warning(
    ngvb2:::ngvb_apply_user_prec_prior(
      ngvb_operator("ar1", n = 6),
      list(rho = list(prior = "pc.cor0", param = c(0.5, 0.5))), "t", verbose = FALSE),
    "not carried over")
})

test_that("all.hyper-style secondary hyperparameters inform rather than warn", {
  op.ar <- ngvb_operator("ar1", n = 6)
  hyper <- list(
    theta1 = list(name = "log precision", short.name = "prec",
                  prior = "loggamma", param = c(1, 5e-5), fixed = FALSE),
    theta2 = list(name = "logit lag one correlation", short.name = "rho",
                  prior = "normal", param = c(0, 0.15), fixed = FALSE),
    theta3 = list(name = "mean", short.name = "mean",
                  prior = "normal", param = c(0, 1), fixed = TRUE))
  expect_no_warning(
    expect_message(
      ngvb2:::ngvb_apply_user_prec_prior(op.ar, hyper, "t", verbose = TRUE),
      "engine's default prior"))
})

test_that("a NULL hyper record warns loudly instead of silently keeping the default", {
  op0 <- ngvb_operator("iid", n = 5)
  expect_warning(
    ngvb2:::ngvb_apply_user_prec_prior(op0, NULL, "z", verbose = FALSE),
    "no hyperparameter record")
})

test_that("ngvb() carries a loggamma f() prior into the fit and changes it", {
  skip_if_not_installed("INLA")
  set.seed(1); n <- 25; g <- rnorm(n); g[7] <- 6; y <- g + rnorm(n, sd = 0.3)
  d <- data.frame(y = y, i = 1:n)
  L0 <- INLA::inla(y ~ -1 + f(i, model = "iid"), data = d,
                   control.compute = list(config = TRUE))
  Lg <- INLA::inla(y ~ -1 + f(i, model = "iid",
                              hyper = list(prec = list(prior = "loggamma", param = c(20, 40)))),
                   data = d, control.compute = list(config = TRUE))
  expect_message(ngvb(Lg, iter = 3, verbose = TRUE), "carried the loggamma")
  n0 <- suppressMessages(ngvb(L0, iter = 4, verbose = FALSE))
  ng <- suppressMessages(ngvb(Lg, iter = 4, verbose = FALSE))
  expect_false(isTRUE(all.equal(n0$eta[["i"]], ng$eta[["i"]])))
})

test_that("REGRESSION: the prior is carried even when ngvb() runs inside a function", {
  ## The old formula-reparsing extraction silently fell back to the default
  ## prior in exactly this scenario; all.hyper has no such failure mode.
  skip_if_not_installed("INLA")
  set.seed(3); n <- 25
  d <- data.frame(y = rnorm(n), i = 1:n)
  run <- function(shape, rate) {
    LGM <- INLA::inla(y ~ -1 + f(i, model = "iid",
                                 hyper = list(prec = list(prior = "loggamma",
                                                          param = c(shape, rate)))),
                      data = d, control.compute = list(config = TRUE))
    suppressMessages(ngvb(LGM, iter = 1, verbose = FALSE))
  }
  LnGM <- run(1000, 0.01)
  p <- 1e5                                             # the Gamma(1000, 0.01) mean
  expect_equal(LnGM$ops$i$logprior(log(p)),
               stats::dgamma(p, shape = 1000, rate = 0.01, log = TRUE) + log(p))
})
