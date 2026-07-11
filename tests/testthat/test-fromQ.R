## op_from_Q: canonical signed-incidence factorization Q = D^T D for M-matrices.

test_that("from_Q reproduces the rw1 structure matrix exactly (intrinsic, rankdef 1)", {
  n  <- 12
  Q  <- ngvb_precision(ngvb_operator("rw1", n = n))     # tau = exp(4) at theta.initial
  op <- ngvb_operator("from_Q", Q = Q)
  expect_equal(as.matrix(Matrix::crossprod(op$Dfunc(0))), as.matrix(Q),
               tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(op$rankdef, 1L)
  expect_equal(op$n, n)
})

test_that("from_Q handles a proper CAR structure (full rank, edge + anchor rows)", {
  ## 4-cycle with diagonal excess: Q = diag(3) - A (3 > degree 2, diagonally dominant)
  A <- Matrix::sparseMatrix(i = c(1, 2, 3, 4), j = c(2, 3, 4, 1), x = 1, dims = c(4, 4))
  A <- A + Matrix::t(A)
  Q <- Matrix::Diagonal(4, 3) - A
  op <- ngvb_operator("from_Q", Q = Q)
  expect_equal(op$rankdef, 0L)
  expect_equal(as.matrix(Matrix::crossprod(op$Dfunc(0))), as.matrix(Q),
               tolerance = 1e-12, ignore_attr = TRUE)
  ## 4 edges + 4 anchors
  expect_equal(nrow(op$Dfunc(0)), 8L)
})

test_that("from_Q on a diagonal Q reduces to iid-style anchors only", {
  Q  <- Matrix::Diagonal(5, x = c(1, 2, 3, 4, 5))
  op <- ngvb_operator("from_Q", Q = Q)
  D  <- op$Dfunc(0)
  expect_equal(dim(D), c(5L, 5L))
  expect_equal(op$rankdef, 0L)
  expect_equal(as.matrix(Matrix::crossprod(D)), as.matrix(Q),
               tolerance = 1e-12, ignore_attr = TRUE)
})

test_that("from_Q counts rank deficiency per unanchored connected component", {
  ## two disconnected paths: one intrinsic (zero row sums), one anchored
  D1 <- Matrix::bandSparse(2, 3, k = c(0, 1), diagonals = list(rep(-1, 2), rep(1, 2)))
  Q1 <- Matrix::crossprod(D1)                        # intrinsic path on nodes 1:3
  Q2 <- Matrix::Diagonal(2, 2) - Matrix::sparseMatrix(i = c(1, 2), j = c(2, 1), x = 1)
  Q  <- Matrix::bdiag(Q1, Q2)                        # anchored pair on nodes 4:5
  op <- ngvb_operator("from_Q", Q = Q)
  expect_equal(op$rankdef, 1L)
})

test_that("from_Q refuses positive off-diagonals (rw2-class) and negative row sums", {
  Qrw2 <- ngvb_precision(ngvb_operator("rw2", n = 10))
  expect_error(ngvb_operator("from_Q", Q = Qrw2), "positive off-diagonals")

  Qbad <- Matrix::Diagonal(3, 1) - Matrix::sparseMatrix(i = c(1, 2), j = c(2, 1), x = 2,
                                                        dims = c(3, 3))
  expect_error(ngvb_operator("from_Q", Q = Qbad), "row sums|non-positive off-diagonal")
})

test_that("from_Q on the ICAR structure matches the car operator's precision", {
  ## a small lattice adjacency
  A <- Matrix::sparseMatrix(i = c(1, 1, 2, 3, 2, 4), j = c(2, 3, 4, 4, 3, 1), x = 1,
                            dims = c(4, 4))
  A <- methods::as((A + Matrix::t(A)) > 0, "CsparseMatrix") * 1
  op.car <- ngvb_operator("car", W = A, intrinsic = TRUE)
  Q      <- ngvb_precision(op.car)
  op.q   <- ngvb_operator("from_Q", Q = Q)
  expect_equal(as.matrix(ngvb_precision(op.q, theta = 0)),
               as.matrix(Q), tolerance = 1e-10, ignore_attr = TRUE)
  expect_equal(op.q$rankdef, op.car$rankdef)
})

test_that("ngvb() runs end to end on a from_Q component", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  set.seed(4); n <- 30
  x <- cumsum(rnorm(n)); x[15] <- x[15] + 5
  d <- data.frame(y = x + rnorm(n, sd = 0.3), i = 1:n)
  LGM <- INLA::inla(y ~ -1 + f(i, model = "rw1", scale.model = FALSE),
                    data = d, control.compute = list(config = TRUE))
  ## rw1 structure at unit precision: first-difference incidence crossprod
  D1 <- Matrix::bandSparse(n - 1, n, k = c(0, 1),
                           diagonals = list(rep(-1, n - 1), rep(1, n - 1)))
  Q  <- Matrix::crossprod(D1)
  op <- ngvb_operator("from_Q", Q = Q)
  fitq <- suppressWarnings(suppressMessages(
    ngvb(LGM, components = list(i = op), iter = 3, verbose = FALSE)))
  expect_true(is.finite(fitq$eta[["i"]]))
  expect_true(all(vapply(fitq$V, function(v) all(is.finite(v)), TRUE)))
})
