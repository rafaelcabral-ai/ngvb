## ---------------------------------------------------------------------------
## Discrepancy statistic d_i = E_q( [D(theta) x]_i^2 ), the message sent from the
## (x, theta) fit to the mixing-variable update. Computed from INLA's Gaussian
## mixture over hyperparameter configurations (config = TRUE):
##   d_i = sum_j w_j * ( (D mu_j)_i^2 + [D Sigma_j D^T]_ii ),
## with w_j the posterior weights, mu_j / Sigma_j the per-config latent
## mean / covariance of the component block (Predictor block removed; see
## the compact-mode indexing note).
## ---------------------------------------------------------------------------

#' @keywords internal
ngvb_component_index <- function(fit, comp.name) {
  ct <- fit$misc$configs$contents
  predlen <- if ("Predictor" %in% ct$tag) ct$length[match("Predictor", ct$tag)] else 0L
  k <- which(ct$tag == comp.name)
  if (length(k) != 1L) stop("ngvb2: component '", comp.name, "' not found in configs")
  start <- ct$start[k] - predlen
  start:(start + ct$length[k] - 1L)
}

#' @keywords internal
ngvb_theta_positions <- function(fit, comp.name) {
  ## positions of this component's hyperparameters within the theta vector
  grep(paste0(" for ", comp.name, "$"), rownames(fit$summary.hyperpar))
}

#' @keywords internal
compute_d <- function(fit, op, comp.name) {
  Dfunc     <- op$Dfunc
  h         <- op$h
  sel       <- ngvb_component_index(fit, comp.name)
  theta.pos <- ngvb_theta_positions(fit, comp.name)
  cfgs      <- fit$misc$configs$config
  nconfigs  <- length(cfgs)

  dmat <- matrix(NA_real_, nconfigs, length(h))
  ll   <- numeric(nconfigs)
  for (j in seq_len(nconfigs)) {
    cf    <- cfgs[[j]]
    m     <- cf$improved.mean[sel]
    Sigma <- as.matrix(cf$Qinv[sel, sel, drop = FALSE])
    Sigma <- Sigma + t(Sigma) - diag(diag(Sigma))      # symmetrize (Qinv stored upper-tri)
    Di    <- Dfunc(cf$theta[theta.pos])
    Dm    <- as.numeric(Di %*% m)
    dmat[j, ] <- Dm^2 + diag(as.matrix(Di %*% Sigma %*% Matrix::t(Di)))
    ll[j] <- cf$log.posterior
  }
  w <- exp(ll - max(ll)); w <- w / sum(w)
  as.numeric(colSums(dmat * w))
}
