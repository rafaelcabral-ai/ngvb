## Additional operator + method coverage.

test_that("SPDE operator equals tau(kappa^4 M0 + 2 kappa^2 M1 + M2)", {
  skip_if_not_installed("fmesher")
  skip_if_not_installed("INLA")
  set.seed(1); loc <- matrix(runif(40), 20, 2)
  mesh <- fmesher::fm_mesh_2d(loc, max.edge = c(0.3, 0.6), cutoff = 0.1)
  spde <- INLA::inla.spde2.matern(mesh)
  op <- ngvb_operator("spde", spde = spde)
  tau <- 1.5; k2 <- 4
  Q <- as.matrix(ngvb_precision(op, theta = c(log(tau), log(k2)), V = op$h))
  M0 <- spde$param.inla$M0; M1 <- spde$param.inla$M1; M2 <- spde$param.inla$M2
  Q_ref <- as.matrix(tau * (k2^2 * M0 + 2 * k2 * M1 + M2))
  expect_equal(Q, Q_ref, tolerance = 1e-8, ignore_attr = TRUE)
})

test_that("intrinsic CAR builds an edge-incidence operator (rank deficiency 1)", {
  W <- matrix(0, 4, 4)
  edges <- rbind(c(1, 2), c(2, 3), c(3, 4))
  for (e in seq_len(nrow(edges))) { W[edges[e, 1], edges[e, 2]] <- 1; W[edges[e, 2], edges[e, 1]] <- 1 }
  op <- ngvb_operator("car", W = W, intrinsic = TRUE)
  expect_equal(op$rankdef, 1L)
  Q <- as.matrix(ngvb_precision(op, theta = 0, V = op$h))
  R <- diag(rowSums(W)) - W            # ICAR structure matrix
  expect_equal(Q, R, tolerance = 1e-8, ignore_attr = TRUE)
})

test_that("SVI method also recovers planted outliers", {
  skip_on_cran(); skip_if_not_installed("INLA")
  set.seed(42); N <- 100
  innov <- rnorm(N, sd = 0.7); innov[c(30, 70)] <- c(7, -7)
  x <- numeric(N); x[1] <- innov[1]; for (i in 2:N) x[i] <- 0.5 * x[i-1] + innov[i]
  y <- x + rnorm(N, sd = 0.3)
  LGM <- INLA::inla(y ~ -1 + f(i, model = "ar1"), data = data.frame(y = y, i = 1:N),
                    control.compute = list(config = TRUE))
  res <- ngvb(LGM, method = "SVI", iter = 6, verbose = FALSE)
  expect_setequal(order(res$V$i, decreasing = TRUE)[1:2], c(30, 70))
})

test_that("ng.check returns a fixed-effect sensitivity matrix without error", {
  skip_on_cran(); skip_if_not_installed("INLA")
  set.seed(1); ng <- 25; nrep <- 5; gi <- rep(1:ng, each = nrep)
  b0 <- rnorm(ng); b0[c(5, 20)] <- c(5, -5); tt <- rnorm(ng * nrep)
  yy <- 2 + b0[gi] + 0.5 * tt + rnorm(ng * nrep, sd = 0.3)
  LGM <- INLA::inla(y ~ 1 + tt + f(g, model = "iid"),
                    data = data.frame(y = yy, g = gi, tt = tt), control.compute = list(config = TRUE))
  chk <- ng.check(LGM, compute.fixed = TRUE)
  expect_true(is.matrix(chk$sens.fixed))
  expect_equal(rownames(chk$sens.fixed), "g")
  expect_true(all(is.finite(chk$sens.fixed)))
})
