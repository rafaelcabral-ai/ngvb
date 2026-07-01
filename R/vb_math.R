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

## Inverse-CDF sampler for a 1-D log-density on (supp.min, supp.max).
#' @keywords internal
sampler.inverseCDF <- function(logpdf, supp.min = 0, supp.max = 100,
                               supp.points = 5000, n.samples = 1000) {
  mode  <- stats::optimize(logpdf, interval = c(supp.min, supp.max), maximum = TRUE)$maximum
  lpm   <- logpdf(mode)
  pdf   <- function(x) exp(logpdf(x) - lpm)
  imin <- 1; while (pdf(mode / 1.3^imin) > 1e-7 && mode / 1.3^imin > supp.min) imin <- imin + 1
  imax <- 1; while (pdf(mode * 1.3^imax) > 1e-7 && mode * 1.3^imax < supp.max) imax <- imax + 1
  x     <- seq(max(supp.min + 1e-8, mode / 1.3^imin), min(supp.max, mode * 1.3^imax),
               length.out = supp.points)
  dens  <- pdf(x); cdf <- cumsum(dens); cdf <- cdf / cdf[length(cdf)]
  keep  <- !duplicated(cdf) & cdf < 0.9999999
  stats::spline(x = cdf[keep], y = x[keep], xout = stats::runif(n.samples))$y
}
