## ---------------------------------------------------------------------------
## Generalized inverse Gaussian (GIG) moments, parameterized as in the ngvb
## papers: density proportional to  x^(p-1) exp(-0.5 (a x + b / x)),  x > 0.
## Ratios of Bessel functions are computed with expon.scaled = TRUE so the
## exp(z) factors cancel and large arguments do not overflow.
## ---------------------------------------------------------------------------

#' @keywords internal
GIGM1 <- function(p, a, b) {            # E[X]
  z <- sqrt(a * b)
  sqrt(b / a) * besselK(z, p + 1, expon.scaled = TRUE) / besselK(z, p, expon.scaled = TRUE)
}

#' @keywords internal
GIGM2 <- function(p, a, b) {            # E[X^2]
  z <- sqrt(a * b)
  (b / a) * besselK(z, p + 2, expon.scaled = TRUE) / besselK(z, p, expon.scaled = TRUE)
}

#' @keywords internal
GIGMm1 <- function(p, a, b) {           # E[1/X]
  z <- sqrt(a * b)
  sqrt(a / b) * besselK(z, p + 1, expon.scaled = TRUE) / besselK(z, p, expon.scaled = TRUE) - 2 * p / b
}

#' @keywords internal
GIGmode <- function(p, a, b) {
  ((p - 1) + sqrt((p - 1)^2 + a * b)) / a
}

## Sampling from a single GIG (scalar p, a, b). GIGrvg uses the (lambda, chi, psi)
## parameterization with density prop. x^(lambda-1) exp(-0.5 (chi/x + psi x)),
## so a -> psi, b -> chi, p -> lambda.
#' @keywords internal
rGIG <- function(n, p, a, b) GIGrvg::rgig(n, lambda = p, chi = b, psi = a)

## Monte-Carlo GIG moment of given order. Used for the eta full-conditional,
## whose order -N/2 + 1 makes the Bessel-ratio moments numerically unstable.
#' @keywords internal
mGIG <- function(p, a, b, order = 1L, n = 2000L) {
  mean(rGIG(n, p, a, b)^order)
}

## Full q(eta) summary from a SINGLE draw (mean/inverse-mean/sd/median/CI all
## computed from the same n.sampling samples, rather than mGIG's independent
## draws per moment -- less MC noise between the two moments that feed the
## next V-update, and the sd/median/CI come for free from the same draws).
## q05/q95 give a 90% credible interval (not the posterior mean +/- sd).
#' @keywords internal
gig_moments <- function(p, a, b, n = 2000L) {
  x <- rGIG(n, p, a, b)
  list(mean = mean(x), median = stats::median(x), inv_mean = mean(1 / x),
       sd = stats::sd(x),
       q05 = stats::quantile(x, 0.05, names = FALSE),
       q95 = stats::quantile(x, 0.95, names = FALSE),
       samples = x, p = p, a = a, b = b)
}

## Weighted quantile on a discretized (grid, weight) density -- used to give
## scvi_update() the same median/CI summary as gig_moments(), from its quadrature.
#' @keywords internal
weighted_quantile <- function(grid, w, probs) {
  o <- order(grid); grid <- grid[o]; w <- w[o]
  cw <- cumsum(w) / sum(w)
  vapply(probs, function(p) grid[which(cw >= p)[1L]], numeric(1))
}

## Per-index posterior summary of the mixing variables V_i ~ GIG(-1, a_V, b_V_i),
## conditional on the fitted eta (a_V is shared, b_V varies by index). The mean
## is exact (GIGM1, vectorized over b_V); GIGrvg::rgig has no closed-form
## quantile and SEGFAULTS if given a vector chi, so median/q05/q95 are read off
## one Monte-Carlo sample PER INDEX (looped; run once at convergence, not per
## VB iteration, so the cost is negligible).
#' @keywords internal
ngvb_V_summary <- function(a_V, b_V, n = 2000L) {
  mean_v <- GIGM1(-1, a_V, b_V)
  qs <- vapply(b_V, function(b_i) {
    x <- rGIG(n, -1, a_V, b_i)
    c(median = stats::median(x), q05 = stats::quantile(x, 0.05, names = FALSE),
      q95 = stats::quantile(x, 0.95, names = FALSE))
  }, numeric(3))
  list(mean = mean_v, median = qs["median", ], q05 = qs["q05", ], q95 = qs["q95", ])
}

## ---- SCVI: the collapsed (eta integrated out) posterior of eta ------------
## log q(eta) up to a constant (Cabral, Bolin & Rue 2022, Theorem 2), with
## theta = alpha.eta the exponential PC-prior rate and d_i = E[(Dx)_i^2].
#' @keywords internal
eta.prior.f <- function(d, h, theta, N) {
  function(eta) {
    vapply(eta, function(e) {
      value <- sqrt((1 / e) * (d + h^2 / e))
      sum(h) / e - theta * e - N / 2 * log(e) +
        sum(log(besselK(value, -1, expon.scaled = TRUE)) - value) -
        0.5 * sum(log(d * e + h^2))
    }, numeric(1))
  }
}

## SCVI update by DETERMINISTIC quadrature of the 1-D collapsed eta posterior.
## Returns E[eta] and the reciprocal mixing moments E[1/V_i] by integrating on an
## adaptive log-grid around the posterior mode -- no sampling, so the result is
## stable and reproducible (the earlier inverse-CDF spline sampler could overshoot
## on flat CDF regions and return degenerate eta values).
#' @keywords internal
scvi_update <- function(d, h, alpha, N, cap = 500, ngrid = 400L) {
  logf <- eta.prior.f(d, h, alpha, N)
  mode <- tryCatch(stats::optimize(logf, c(1e-4, cap), maximum = TRUE)$maximum,
                   error = function(e) 1)
  lo   <- max(1e-4, mode / 50); hi <- min(cap, max(mode * 50, lo * 10))
  grid <- exp(seq(log(lo), log(hi), length.out = ngrid))
  lp   <- logf(grid); lp[!is.finite(lp)] <- -Inf
  deta <- numeric(ngrid)
  deta[1]              <- grid[2] - grid[1]
  deta[ngrid]          <- grid[ngrid] - grid[ngrid - 1]
  deta[2:(ngrid - 1L)] <- (grid[3:ngrid] - grid[1:(ngrid - 2L)]) / 2
  w <- exp(lp - max(lp)) * deta; w <- w / sum(w)
  eta.mean <- sum(w * grid)
  list(eta    = eta.mean,
       median = weighted_quantile(grid, w, 0.5),
       EVm1   = vapply(seq_along(h),
                       function(i) sum(w * GIGMm1(-1, 1 / grid, d[i] + h[i]^2 / grid)),
                       numeric(1)),
       sd     = sqrt(max(0, sum(w * grid^2) - eta.mean^2)),
       q05    = weighted_quantile(grid, w, 0.05),
       q95    = weighted_quantile(grid, w, 0.95),
       grid = grid, weights = w)
}
