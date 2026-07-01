## ---------------------------------------------------------------------------
## ng.check(): diagnostic for the latent Gaussian assumption (Cabral, Bolin &
## Rue, JRSS-B 2025). For each component it computes the Bayes-factor sensitivity
## s0 = sum_i d_i and (Gaussian response) a reference distribution / p-value,
## from INLA's Gaussian mixture over hyperparameter configurations.
##
## Gaussian-response closed form (Eq. 20, Prop. 1), per config k with weight w_k:
##   b      = D mu_k
##   Gamma  = -D Sigma_k D^T + diag(h)
##   d_i    = (b_i^4 + 3 Gamma_ii^2 - 6 b_i^2 Gamma_ii) / (8 h_i^3)
##   s0     = sum_i d_i
##   Var(s0)= (3/8) (h^-3)^T (Gamma .^4) (h^-3)        [reference variance, mean 0]
##   p      = Phi(-s0 / sqrt(Var(s0)))
## Operators (D, h) come from the same registry/auto-detection as ngvb().
## ---------------------------------------------------------------------------

#' Sensitivity of fixed effects to the non-Gaussianity parameter (Theorem 2).
#' @keywords internal
ng_sens_fixed <- function(b, gii, s12, u1, h) {
  ## b = E[(Dx)_i], gii = Var[(Dx)_i] (= Sigma of Dx diag), s12 = Cov(z_fixed,(Dx)_i),
  ## u1 = E[z_fixed]; returns the per-index sensitivity, summed by the caller.
  EDx2  <- gii + b^2
  EDx4  <- 3 * gii^2 + 6 * gii * b^2 + b^4
  EzDx2 <- 2 * s12 * b + u1 * (gii + b^2)
  EzDx4 <- 4 * s12 * (3 * gii * b + b^3) + u1 * (3 * gii^2 + 6 * gii * b^2 + b^4)
  (1 / (8 * h^3)) * (3 * h^2 * u1 - 6 * h * EzDx2 + EzDx4) -
    u1 * (1 / (8 * h^3)) * (3 * h^2 - 6 * h * EDx2 + EDx4)
}

#' Check the latent Gaussian assumption of an INLA fit.
#'
#' @param fit An `inla` object fitted with `control.compute = list(config = TRUE)`.
#' @param selection Optional named list of components to check (default: all random).
#' @param components Optional operator overrides (e.g. SPDE), as in [ngvb()].
#' @param compute.fixed If `TRUE`, also return the sensitivity of each fixed effect
#'   to each component's non-Gaussianity parameter.
#' @return An object of class `ngvb.check`: per component the BF sensitivity `s0`,
#'   the per-index contributions `d`, and (Gaussian response) the reference SD and
#'   p-value; plus `sens.fixed` if requested.
#' @export
ng.check <- function(fit, selection = NULL, components = NULL, compute.fixed = TRUE) {
  if (is.null(fit$misc$configs))
    stop("ngvb2: refit the LGM with control.compute = list(config = TRUE).")
  if (is.null(selection))
    selection <- lapply(fit$summary.random, function(x) seq_len(nrow(x)))
  comp.names <- names(selection)
  ops <- stats::setNames(
    lapply(comp.names, function(cn) ngvb_detect_operator(fit, cn, user.op = components[[cn]])),
    comp.names)

  cfgs    <- fit$misc$configs$config
  nconfig <- length(cfgs)
  w       <- vapply(cfgs, function(c) c$log.posterior, 0)
  w       <- exp(w - max(w)); w <- w / sum(w)
  family  <- fit$.args$family[1]
  gaussian <- identical(family, "gaussian")

  ## fixed-effect bookkeeping (compact-mode: subtract all predictor-block lengths)
  ct       <- fit$misc$configs$contents
  pred.len <- ngvb_predictor_length(ct)
  fixed.names <- rownames(fit$summary.fixed)
  fixed.pos   <- if (length(fixed.names))
    vapply(fixed.names, function(nm) ct$start[match(nm, ct$tag)] - pred.len, 0L) else integer(0)
  if (length(fixed.pos) == 0L) compute.fixed <- FALSE

  per.comp <- list()
  sens.fixed <- if (compute.fixed)
    matrix(NA_real_, length(comp.names), length(fixed.names),
           dimnames = list(comp.names, fixed.names)) else NULL

  for (cn in comp.names) {
    op <- ops[[cn]]; h <- op$h; Dfunc <- op$Dfunc
    sel  <- ngvb_component_index(fit, cn)
    tpos <- ngvb_theta_positions(fit, cn)

    di.mat  <- matrix(0, nconfig, length(h))
    s0.k    <- numeric(nconfig)
    varref.k <- numeric(nconfig)
    sfix.k  <- if (compute.fixed) matrix(0, length(fixed.names), nconfig) else NULL

    for (k in seq_len(nconfig)) {
      cf    <- cfgs[[k]]
      D     <- Dfunc(cf$theta[tpos])
      mu    <- cf$improved.mean[sel]
      Sigma <- as.matrix(cf$Qinv[sel, sel, drop = FALSE])
      Sigma <- Sigma + t(Sigma) - diag(diag(Sigma))
      b     <- as.numeric(D %*% mu)
      DSDt  <- as.matrix(D %*% Sigma %*% Matrix::t(D))
      Gamma <- -DSDt + diag(h)
      gii   <- diag(Gamma)
      di.mat[k, ] <- (b^4 + 3 * gii^2 - 6 * b^2 * gii) / (8 * h^3)
      s0.k[k]     <- sum(di.mat[k, ])
      if (gaussian)
        varref.k[k] <- (3 / 8) * as.numeric(t(h^(-3)) %*% (Gamma^4) %*% (h^(-3)))
      if (compute.fixed) {
        ## cross-covariance Cov((Dx)_i, z_fixed) = D %*% Cov(x, z_fixed);
        ## the (random, fixed) block of the full config covariance (upper-stored).
        cross <- as.matrix(cf$Qinv[sel, fixed.pos, drop = FALSE])
        Sxz   <- as.matrix(D %*% cross)
        muz   <- cf$improved.mean[fixed.pos]
        for (j in seq_along(fixed.pos))
          sfix.k[j, k] <- sum(ng_sens_fixed(b, gii, Sxz[, j], muz[j], h))
      }
    }
    res <- list(s0 = sum(w * s0.k), d = as.numeric(colSums(w * di.mat)))
    if (gaussian) {
      res$var.ref <- sum(w * varref.k)
      res$sd.ref  <- sqrt(res$var.ref)
      res$p.value <- stats::pnorm(-res$s0 / res$sd.ref)
    }
    per.comp[[cn]] <- res
    if (compute.fixed) sens.fixed[cn, ] <- as.numeric(sfix.k %*% w)
  }

  out <- list(components = per.comp, sens.fixed = sens.fixed,
              gaussian = gaussian, selection = selection)
  class(out) <- "ngvb.check"
  out
}

#' @export
print.ngvb.check <- function(x, ...) {
  cat("ngvb2 latent-Gaussianity check\n")
  tab <- do.call(rbind, lapply(names(x$components), function(cn) {
    c <- x$components[[cn]]
    data.frame(component = cn, s0 = round(c$s0, 3),
               sd.ref = if (!is.null(c$sd.ref)) round(c$sd.ref, 3) else NA,
               p.value = if (!is.null(c$p.value)) signif(c$p.value, 3) else NA)
  }))
  print(tab, row.names = FALSE)
  if (!is.null(x$sens.fixed)) {
    cat("\nSensitivity of fixed effects to non-Gaussianity:\n")
    print(round(x$sens.fixed, 4))
  }
  invisible(x)
}
