## End-to-end checks that require INLA. Slower; skipped on CRAN.

test_that("engine reproduces native ar1 at V = h (Gaussian limit)", {
  skip_if_no_inla()
  set.seed(7); N <- 200
  x <- as.numeric(arima.sim(n = N, model = list(ar = 0.8)))
  y <- x + rnorm(N, sd = 0.2)
  dat <- data.frame(y = y, i = 1:N)
  hyper <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)),
                rho  = list(prior = "pc.cor0", param = c(0.5, 0.5)))
  fit_nat <- INLA::inla(y ~ -1 + f(i, model = "ar1", hyper = hyper), data = dat,
                        control.compute = list(config = TRUE))
  op <- ngvb_operator("ar1", n = N)
  fit_eng <- INLA::inla(y ~ -1 + f(i, model = ngvb_rgeneric(op, V = op$h)), data = dat,
                        control.compute = list(config = TRUE))
  mn <- fit_nat$summary.random$i$mean; me <- fit_eng$summary.random$i$mean
  expect_gt(cor(mn, me), 0.999)         # latent posterior essentially identical
})

test_that("ngvb(fit) auto-detects ar1 and recovers planted outliers", {
  skip_if_no_inla()
  set.seed(42); N <- 100
  innov <- rnorm(N, sd = 0.7); innov[c(30, 70)] <- c(7, -7)
  x <- numeric(N); x[1] <- innov[1]; for (i in 2:N) x[i] <- 0.5 * x[i-1] + innov[i]
  y <- x + rnorm(N, sd = 0.3)
  LGM <- INLA::inla(y ~ -1 + f(i, model = "ar1"), data = data.frame(y = y, i = 1:N),
                    control.compute = list(config = TRUE))
  res <- ngvb(LGM, iter = 6, verbose = FALSE)
  expect_equal(res$ops$i$type, "ar1")            # auto-detected
  expect_gt(res$eta[["i"]], 0.3)                 # non-Gaussianity detected
  top <- order(res$V$i, decreasing = TRUE)[1:2]  # the two jumps are the top-2 V
  expect_setequal(top, c(30, 70))
})

test_that("ngvb(fit) handles multiple components (block-diagonal, one fit)", {
  skip_if_no_inla()
  set.seed(1); ng <- 30; nrep <- 5; gi <- rep(1:ng, each = nrep)
  b0 <- rnorm(ng); b0[c(5, 20)] <- c(5, -5); b1 <- rnorm(ng, sd = 0.5); tt <- rnorm(ng * nrep)
  yy <- b0[gi] + b1[gi] * tt + rnorm(ng * nrep, sd = 0.3)
  LGM <- INLA::inla(y ~ 1 + tt + f(g, model = "iid") + f(g2, tt, model = "iid"),
                    data = data.frame(y = yy, g = gi, g2 = gi, tt = tt),
                    control.compute = list(config = TRUE))
  res <- ngvb(LGM, iter = 4, verbose = FALSE)
  expect_named(res$eta, c("g", "g2"))
  expect_gt(res$V$g[5], median(res$V$g))         # outlying intercept flagged
  expect_gt(res$V$g[20], median(res$V$g))
})

test_that("ng.check reproduces the reference diagnostic on RW1", {
  skip_if_no_inla()
  set.seed(3); N <- 100
  x <- cumsum(rnorm(N, sd = 0.3)); x[50:N] <- x[50:N] + 6; x[80:N] <- x[80:N] - 9
  y <- x + rnorm(N, sd = 0.4)
  LGM <- INLA::inla(y ~ -1 + f(i, model = "rw1", constr = TRUE),
                    data = data.frame(y = y, i = 1:N), control.compute = list(config = TRUE))
  chk <- ng.check(LGM, compute.fixed = FALSE)
  expect_gt(chk$components$i$s0, 10)             # strong non-Gaussian signal
  expect_lt(chk$components$i$p.value, 0.01)
  top <- order(chk$components$i$d, decreasing = TRUE)[1:2]
  expect_true(all(abs(top - c(79, 49)) <= 1))    # jump locations
})
