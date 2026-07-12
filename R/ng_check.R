## ---------------------------------------------------------------------------
## ng.check(): diagnostic for the latent Gaussian assumption (Cabral, Bolin &
## Rue, JRSS-B 2025). For each component it computes the Bayes-factor sensitivity
## s0 = sum_i d_i and (Gaussian response) a reference distribution / p-value,
## from INLA's Gaussian mixture over hyperparameter configurations.
##
## Gaussian-response closed form (Eq. 20, Prop. 1), evaluated at the posterior
## MODE gamma-hat of the hyperparameters (the config of largest weight; the
## reference variance is derived conditional on gamma-hat, so it is not a
## mixture over configs):
##   b      = D mu
##   Gamma  = -D Sigma D^T + diag(h)
##   d_i    = (b_i^4 + 3 Gamma_ii^2 - 6 b_i^2 Gamma_ii) / (8 h_i^3)
##   s0     = sum_i d_i
##   Var(s0)= (3/8) (h^-3)^T (Gamma .^4) (h^-3)        [reference variance, mean 0]
##   p      = Phi(-s0 / sqrt(Var(s0)))
## `s0.mixture` / `d.mixture` additionally give the hyperparameter-weighted
## average. Operators (D, h) come from the same registry/auto-detection as ngvb().
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
#' @param compute.random If `TRUE`, also return, for each checked component, the
#'   sensitivity of *its own* posterior mean at every node to its own
#'   non-Gaussianity -- i.e. how much each node's fitted value would move under
#'   [ngvb()], estimated from the base (Gaussian) fit alone, without actually
#'   fitting the non-Gaussian model. This is the same Bayes-factor-sensitivity
#'   theory used for `sens.fixed` (Theorem 2 of Cabral, Bolin & Rue, JRSS-B
#'   2025), applied to the component's own field instead of a fixed effect: it
#'   reuses the cross-covariance already computed for `d`, so it costs one
#'   extra matrix multiply per hyperparameter configuration, not a refit.
#'   Useful as a preview -- e.g. mapped over an SPDE mesh -- of where a
#'   subsequent `ngvb()` fit is expected to change the predictions the most,
#'   before actually running it.
#' @param plot If `TRUE` (default), draw the diagnostic plots (see [plot.ngvb.check()]):
#'   per-index Bayes-factor sensitivity, and the observed overall sensitivity against
#'   its Gaussian reference distribution.
#' @return An object of class `ngvb.check`. Per component, evaluated at the
#'   hyperparameter posterior mode \eqn{\hat\gamma}: the BF sensitivity `s0`, the
#'   per-index contributions `d`, and (Gaussian response) the reference SD `sd.ref`
#'   and `p.value`; the hyperparameter-mixture averages are also kept as
#'   `s0.mixture` / `d.mixture`. Plus `sens.fixed` if `compute.fixed = TRUE`, and
#'   `sens.random` (a named list, one vector per checked component, aligned with
#'   `selection`) if `compute.random = TRUE`.
#' @seealso [ngvb()]
#' @examples
#' \donttest{
#' if (requireNamespace("INLA", quietly = TRUE)) {
#'   set.seed(1); n <- 100
#'   x <- cumsum(rnorm(n, sd = 0.3)); x[50:n] <- x[50:n] + 6
#'   y <- x + rnorm(n, sd = 0.4)
#'   LGM <- INLA::inla(y ~ -1 + f(i, model = "rw1", constr = TRUE),
#'                     data = data.frame(y = y, i = 1:n),
#'                     control.compute = list(config = TRUE))
#'   ng.check(LGM)      # small p-value flags departure from latent Gaussianity
#' }
#' }
#' @export
ng.check <- function(fit, selection = NULL, components = NULL, compute.fixed = TRUE,
                     compute.random = FALSE, plot = TRUE) {
  .need_inla()
  if (is.null(fit$misc$configs))
    stop("ngvb: refit the LGM with control.compute = list(config = TRUE).")
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
  sens.random <- if (compute.random)
    stats::setNames(vector("list", length(comp.names)), comp.names) else NULL

  for (cn in comp.names) {
    op <- ops[[cn]]; h <- op$h; Dfunc <- op$Dfunc
    sel  <- ngvb_component_index(fit, cn)
    tpos <- ngvb_theta_positions(fit, cn)

    di.mat  <- matrix(0, nconfig, length(h))
    s0.k    <- numeric(nconfig)
    varref.k <- numeric(nconfig)
    sfix.k  <- if (compute.fixed) matrix(0, length(fixed.names), nconfig) else NULL
    srand.k <- if (compute.random) matrix(0, length(sel), nconfig) else NULL

    for (k in seq_len(nconfig)) {
      cf    <- cfgs[[k]]
      D     <- Dfunc(cf$theta[tpos])
      mu    <- cf$improved.mean[sel]
      Sigma <- as.matrix(cf$Qinv[sel, sel, drop = FALSE])
      Sigma <- Sigma + t(Sigma) - diag(diag(Sigma))
      ## Cov((Dx)_i, x_j) = D %*% Cov(x, x_j); shared by Gamma below and, when
      ## requested, by sens.random (Cov((Dx)_i, x_j) for x_j in this SAME
      ## component -- the fixed-effect cross-covariance below is the same
      ## construction against a different set of columns).
      DS    <- D %*% Sigma
      b     <- as.numeric(D %*% mu)
      DSDt  <- as.matrix(DS %*% Matrix::t(D))
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
      if (compute.random) {
        ## Sensitivity of this component's OWN posterior mean at node j to its
        ## OWN non-Gaussianity: same Theorem-2 construction as sens.fixed,
        ## with the field's own nodes standing in for the "fixed effect".
        DSm <- as.matrix(DS)
        for (j in seq_along(sel))
          srand.k[j, k] <- sum(ng_sens_fixed(b, gii, DSm[, j], mu[j], h))
      }
    }
    ## The diagnostic s0(y, gamma-hat) and its reference variance (Prop. 1) are
    ## defined at the posterior MODE of the hyperparameters -- the reference
    ## variance is derived conditional on gamma-hat, so it must not be mixed
    ## across configs. Report the mode-config value as the headline (matching the
    ## reference implementation's s0.mode / var.ref.mode) and keep the
    ## hyperparameter-mixture average as a secondary field.
    km  <- which.max(w)                                   # gamma-hat
    res <- list(s0 = s0.k[km], d = di.mat[km, ],
                s0.mixture = sum(w * s0.k),
                d.mixture  = as.numeric(colSums(w * di.mat)))
    if (gaussian) {
      res$var.ref <- varref.k[km]
      res$sd.ref  <- sqrt(res$var.ref)
      res$p.value <- stats::pnorm(-res$s0 / res$sd.ref)
    }
    per.comp[[cn]] <- res
    if (compute.fixed) sens.fixed[cn, ] <- as.numeric(sfix.k %*% w)
    if (compute.random) sens.random[[cn]] <- as.numeric(srand.k %*% w)
  }

  out <- list(components = per.comp, sens.fixed = sens.fixed, sens.random = sens.random,
              gaussian = gaussian, selection = selection)
  class(out) <- "ngvb.check"
  if (isTRUE(plot)) print(plot(out))
  invisible(out)
}

#' Diagnostic plots for a latent-Gaussianity check.
#'
#' For each checked component, draws (left) the per-index Bayes-factor sensitivity
#' \eqn{d_i(y)} -- spikes locate where the Gaussian assumption is least adequate --
#' and (right, for a Gaussian response) the observed overall sensitivity
#' \eqn{s_0=\sum_i d_i} against its Gaussian reference distribution; an observed value
#' far in the tail (small p-value) signals latent non-Gaussianity.
#'
#' @param x An `ngvb.check` object from [ng.check()].
#' @param ... Ignored.
#' @return A `patchwork` object combining the per-component diagnostic panels;
#'   called mainly for the plot it draws.
#' @method plot ngvb.check
#' @export
plot.ngvb.check <- function(x, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("patchwork", quietly = TRUE))
    stop("plot.ngvb.check() needs the 'ggplot2' and 'patchwork' packages.")
  accent <- "#C64A2E"
  panels <- list()
  for (cn in names(x$components)) {
    cc <- x$components[[cn]]; d <- cc$d
    j  <- which.max(abs(d))
    df <- data.frame(idx = seq_along(d), d = d, hit = seq_along(d) == j)
    pD <- ggplot2::ggplot(df, ggplot2::aes(x = .data$idx, y = .data$d)) +
      ggplot2::geom_segment(ggplot2::aes(xend = .data$idx, yend = 0),
                            colour = "grey60") +
      ggplot2::geom_point(ggplot2::aes(colour = .data$hit, size = .data$hit)) +
      ggplot2::scale_colour_manual(values = c(`FALSE` = "grey35", `TRUE` = accent),
                                   guide = "none") +
      ggplot2::scale_size_manual(values = c(`FALSE` = 1, `TRUE` = 2.4), guide = "none") +
      ggplot2::labs(x = "index i", y = expression(d[i](y)),
                    title = paste0("BF sensitivity per index: ", cn)) +
      ngvb_theme()

    if (!is.null(cc$sd.ref) && is.finite(cc$sd.ref) && cc$sd.ref > 0) {
      lim <- max(4 * cc$sd.ref, abs(cc$s0) * 1.15)
      xs  <- seq(-lim, lim, length.out = 400)
      rd  <- data.frame(s = xs, dens = stats::dnorm(xs, 0, cc$sd.ref))
      pR <- ggplot2::ggplot(rd, ggplot2::aes(x = .data$s, y = .data$dens)) +
        ggplot2::geom_area(fill = "grey85", colour = "grey55") +
        ggplot2::geom_vline(xintercept = cc$s0, colour = accent, linewidth = 1) +
        ggplot2::labs(x = expression(s[0](y)), y = "reference density",
                      title = paste0("Observed vs reference: ", cn),
                      subtitle = sprintf("observed s0 = %.3g  (p = %.3f)",
                                         cc$s0, cc$p.value)) +
        ngvb_theme()
    } else {
      pR <- ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0, y = 0,
                          label = "reference only for\na Gaussian response",
                          colour = "grey35") +
        ggplot2::theme_void()
    }
    panels <- c(panels, list(pD, pR))
  }
  patchwork::wrap_plots(panels, ncol = 2, byrow = TRUE)
}

#' @export
print.ngvb.check <- function(x, ...) {
  cat("ngvb latent-Gaussianity check (at the hyperparameter posterior mode)\n")
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
