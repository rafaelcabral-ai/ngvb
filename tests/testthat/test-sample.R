## ngvb_sample() / bayes.factor() / summary.ngvb.samples()

test_that("summary.ngvb.samples: 'all' shows both tables, skipping absent ones", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  set.seed(1); n <- 25
  d <- data.frame(y = rnorm(n), i = 1:n)

  ## no fixed effects: 'fixed' has nothing to pool
  LGM  <- INLA::inla(y ~ -1 + f(i, model = "iid"), data = d,
                     control.compute = list(config = TRUE))
  LnGM <- suppressWarnings(suppressMessages(ngvb(LGM, iter = 3, verbose = FALSE)))
  samples <- suppressWarnings(suppressMessages(ngvb_sample(LnGM, n.samples = 8, verbose = FALSE)))

  ## default what = "all"
  expect_output(res.all <- summary(samples), "Fixed effects")
  expect_output(summary(samples), "Hyperparameters")
  expect_null(res.all$fixed)
  expect_true(is.data.frame(res.all$hyperpar))
  expect_true(all(c("mean", "sd") %in% names(res.all$hyperpar)))

  ## what = "fixed" alone still degrades gracefully (not silent)
  expect_output(res.fixed <- summary(samples, what = "fixed"), "no fixed effects")
  expect_null(res.fixed)

  ## what = "hyperpar" alone matches the "all" pooling
  res.hp <- summary(samples, what = "hyperpar")
  expect_equal(res.hp, res.all$hyperpar)
})

test_that("summary.ngvb.samples: 'all' shows real fixed-effect estimates when present", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  set.seed(2); n <- 30
  x <- rnorm(n)
  d <- data.frame(y = 2 + 0.5 * x + rnorm(n), x = x, i = 1:n)

  LGM  <- INLA::inla(y ~ 1 + x + f(i, model = "iid"), data = d,
                     control.compute = list(config = TRUE))
  LnGM <- suppressWarnings(suppressMessages(ngvb(LGM, iter = 3, verbose = FALSE)))
  samples <- suppressWarnings(suppressMessages(ngvb_sample(LnGM, n.samples = 8, verbose = FALSE)))

  res <- summary(samples, what = "all")
  expect_true(is.data.frame(res$fixed))
  expect_true(all(c("(Intercept)", "x") %in% rownames(res$fixed)))
  expect_true(is.data.frame(res$hyperpar))
})

test_that("summary.ngvb.samples: the 'what' choices are discoverable via args()", {
  formals.what <- eval(formals(summary.ngvb.samples)$what)
  expect_equal(formals.what, c("all", "fixed", "hyperpar"))
})
