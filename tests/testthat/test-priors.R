test_that("loggamma precision log-prior matches Gamma(shape,rate) with the log-Jacobian", {
  a <- 0.7; b <- 0.03; lp <- ngvb2:::.loggamma_prec_logprior(a, b)
  for (theta in c(-2, 0, 1.3)) {
    ref <- stats::dgamma(exp(theta), shape = a, rate = b, log = TRUE) + theta
    expect_equal(lp(theta), ref)
  }
})

test_that("prior-family dispatch maps pc.prec/loggamma and rejects the rest", {
  expect_true(ngvb2:::.prec_logprior_from_hyper(list(prior = "pc.prec",  param = c(1, 0.01)))$mappable)
  expect_true(ngvb2:::.prec_logprior_from_hyper(list(prior = "loggamma", param = c(1, 5e-5)))$mappable)
  expect_false(ngvb2:::.prec_logprior_from_hyper(list(prior = "gaussian", param = c(0, 1)))$mappable)
  expect_false(ngvb2:::.prec_logprior_from_hyper(list(prior = "loggamma", param = c(1, 1),
                                                      fixed = TRUE))$mappable)
})

test_that("f() hyper is extracted from a formula by component name", {
  f <- y ~ 1 + x + f(s, model = "besag", graph = W,
                     hyper = list(prec = list(prior = "loggamma", param = c(0.1, 0.01)))) +
    f(t, model = "rw1")
  hy <- ngvb2:::ngvb_extract_f_hyper(f, c("s", "t"))
  expect_equal(hy$s$prec$prior, "loggamma")
  expect_equal(hy$s$prec$param, c(0.1, 0.01))
  expect_null(hy$t)                                  # no hyper on t
})

test_that("applying a mapped prior replaces exactly the precision piece of logprior", {
  ## single-hyperparameter operator: logprior becomes the mapped prior outright
  op   <- ngvb_operator("rw1", n = 8)
  hy   <- list(prec = list(prior = "loggamma", param = c(0.1, 0.01)))
  op2  <- ngvb2:::ngvb_apply_user_prec_prior(op, hy, "s", verbose = FALSE)
  want <- ngvb2:::.loggamma_prec_logprior(0.1, 0.01)
  for (th in c(-1, 0.7, 2)) expect_equal(op2$logprior(th), want(th))

  ## two-hyperparameter operator (ar1): the correlation prior must survive untouched
  op.ar  <- ngvb_operator("ar1", n = 12)
  op.ar2 <- ngvb2:::ngvb_apply_user_prec_prior(
    op.ar, list(prec = list(prior = "pc.prec", param = c(2, 0.01))), "t", verbose = FALSE)
  th <- c(0.3, 0.4); newp <- ngvb2:::.pc_prec_logprior(2, 0.01)
  expect_equal(op.ar2$logprior(th),
               op.ar$logprior(th) - ngvb2:::.pc_prec_logprior(1, 0.01)(th[1]) + newp(th[1]))
})

test_that("an unmappable family and a secondary-hyperparameter prior both warn and fall back", {
  op0 <- ngvb_operator("iid", n = 5)
  expect_warning(
    op1 <- ngvb2:::ngvb_apply_user_prec_prior(
      op0, list(prec = list(prior = "gaussian", param = c(0, 1))), "z", verbose = FALSE),
    "was dropped")
  expect_identical(op1$logprior(0.5), op0$logprior(0.5))   # default PC prior retained

  expect_warning(
    ngvb2:::ngvb_apply_user_prec_prior(
      ngvb_operator("ar1", n = 6),
      list(rho = list(prior = "pc.cor0", param = c(0.5, 0.5))), "t", verbose = FALSE),
    "not carried over")
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
  n0 <- ngvb(L0, iter = 4, verbose = FALSE)
  ng <- suppressMessages(ngvb(Lg, iter = 4, verbose = FALSE))
  expect_false(isTRUE(all.equal(n0$eta[["i"]], ng$eta[["i"]])))
})
