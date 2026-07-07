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
  list(eta  = sum(w * grid),
       EVm1 = vapply(seq_along(h),
                     function(i) sum(w * GIGMm1(-1, 1 / grid, d[i] + h[i]^2 / grid)),
                     numeric(1)))
}
