## ---------------------------------------------------------------------------
## Variational-Bayes loop for a latent non-Gaussian model (single component,
## structured VI / SVI). Repeatedly fits the conditionally-Gaussian LGM for a
## fixed mixing vector V (via INLA), reads the discrepancy d_i = E[(Dx)_i^2],
## and updates the variational factors q(V) and q(eta). See Cabral, Bolin & Rue
## (2022), Theorem 1.
## ---------------------------------------------------------------------------

#' Fit a latent non-Gaussian model by structured variational inference.
#'
#' @param inla.fit.V Function mapping a list of mixing vectors `V` to an INLA
#'   fit of the LGM conditioned on `V` (must use `control.compute=list(config=TRUE)`).
#' @param op Operator descriptor (provides `Dfunc`, `h`).
#' @param comp.name Name of the latent component (the `f(<name>, ...)` index).
#' @param alpha.eta Rate of the exponential PC prior on the non-Gaussianity
#'   parameter eta.
#' @param iter Maximum number of VB iterations.
#' @param stop.rel.change Stop when the relative change in E[eta] drops below this.
#' @param n.sampling Monte-Carlo samples for the eta moments.
#' @param verbose Print per-iteration progress.
#' @return A list with the final INLA fit, the mixing vector `V`, `eta`, the
#'   discrepancy `d`, and the `eta` trajectory.
#' @export
ngvb_svi <- function(inla.fit.V, op, comp.name,
                     alpha.eta = 1, iter = 10, stop.rel.change = 1e-3,
                     n.sampling = 2000, verbose = TRUE) {
  h <- op$h
  N <- length(h)

  V        <- list(h)            # start at the Gaussian model V = h
  eta      <- 0.5
  Eetam1   <- 1 / eta
  eta.path <- eta

  fit <- inla.fit.V(V)

  for (i in seq_len(iter)) {
    d <- compute_d(fit, op, comp.name)

    ## q(V_i) ~ GIG(-1, E[1/eta], d_i + h_i^2 E[1/eta])
    a_V  <- Eetam1
    b_V  <- d + h^2 * Eetam1
    EV   <- GIGM1(-1, a_V, b_V)
    EVm1 <- GIGMm1(-1, a_V, b_V)

    ## q(eta) ~ GIG(-N/2 + 1, 2 alpha, sum(EV - 2h + h^2 E[1/V]))
    p_eta  <- -N / 2 + 1
    a_eta  <- 2 * alpha.eta
    b_eta  <- sum(EV - 2 * h + h^2 * EVm1)
    Eeta   <- mGIG(p_eta, a_eta, b_eta, order =  1L, n = n.sampling)
    Eetam1 <- mGIG(p_eta, a_eta, b_eta, order = -1L, n = n.sampling)

    V[[1]] <- 1 / EVm1            # so that 1/V = E[1/V_i], the precision weight
    eta    <- Eeta
    eta.path <- c(eta.path, eta)

    fit <- inla.fit.V(V)

    rel <- abs(eta.path[i + 1] - eta.path[i]) / eta.path[i]
    if (verbose) cat(sprintf("  iter %2d:  E[eta] = %.4f   (rel change %.4f)\n", i, eta, rel))
    if (is.finite(rel) && rel < stop.rel.change) {
      if (verbose) cat("  converged\n"); break
    }
  }

  list(fit = fit, V = V[[1]], eta = eta, h = h, d = d, eta.path = eta.path)
}
