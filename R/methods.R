## S3 methods for the object returned by ngvb().
##
## Terminology used throughout (see also the vignette): for a mixing variable
## V_i, "mean" = E[V_i] (exact, from the GIG moment formula); "median"/"90% CI"
## = read off one Monte-Carlo sample of q(V_i), conditional on the fitted eta
## (see ngvb_V_summary()). Neither is the same as `x$V`, which is 1/E[1/V_i]
## -- the precision-consistent value the algorithm actually plugs into INLA's
## Q(V) at convergence, pulled below the mean by Jensen's inequality for the
## right-skewed q(V). All three are reported so the two roles (what INLA was
## refit with vs. what the posterior of V_i actually looks like) aren't conflated.

#' @export
print.ngvb <- function(x, ...) {
  cat("Latent non-Gaussian model fit (ngvb)\n")
  cat("  components:", paste(sprintf("%s [%s]", x$comp.names,
                                     vapply(x$ops, `[[`, "", "type")), collapse = ", "), "\n")
  cat("  iterations:", x$iterations, "\n")
  cat("  eta -- mean (median) [90% CI]:\n")
  for (cn in x$comp.names) {
    q <- x$eta.q[[cn]]
    cat(sprintf("    %-12s %.3f (%.3f) [%.3f, %.3f]\n",
                cn, x$eta[[cn]], q$median, q$q05, q$q95))
  }
  invisible(x)
}

#' Summary of an ngvb fit: the underlying INLA summaries plus, per component, the
#' non-Gaussianity parameter and the most strongly down-weighted (outlying) indices.
#' @param object An `ngvb` object.
#' @param n.flag Number of top mixing-variable indices to report per component.
#' @param ... Ignored.
#' @return The `ngvb` object, invisibly; called for the summary it prints.
#' @method summary ngvb
#' @export
summary.ngvb <- function(object, n.flag = 5, ...) {
  cat("Latent non-Gaussian model (ngvb)\n\n")
  cat("Fixed effects:\n"); print(round(object$fit$summary.fixed, 4))
  cat("\nNon-Gaussianity -- eta: mean (median) [90% CI]:\n")
  for (cn in object$comp.names) {
    Vm <- object$V.summary[[cn]]$mean; h <- object$h[[cn]]; q <- object$eta.q[[cn]]
    ord <- order(Vm / h, decreasing = TRUE)[seq_len(min(n.flag, length(Vm)))]
    cat(sprintf(paste0("  %-12s [%s]: eta = %.3f (%.3f) [%.3f, %.3f];\n",
                       "                 most non-Gaussian indices (E[V]/h): %s\n"),
                cn, object$ops[[cn]]$type, object$eta[[cn]], q$median, q$q05, q$q95,
                paste(sprintf("%d (%.1f)", ord, (Vm / h)[ord]), collapse = ", ")))
  }
  invisible(object)
}

#' Posterior density of a component's non-Gaussianity parameter eta.
#'
#' For `method = "SVI"` this is a kernel density estimate over the GIG
#' variational factor's Monte-Carlo sample; for `"SCVI"` it is the weighted
#' quadrature grid underlying `scvi_update()`. Either way it is `q(eta)`, the
#' actual (approximate) posterior -- not just its mean.
#' @param x An `ngvb` object.
#' @param component Component name (default: the first).
#' @param ... Passed to [stats::density()].
#' @return A `density` object (see [stats::density()]); `plot()` it directly.
#' @export
density.ngvb <- function(x, component = x$comp.names[1L], ...) {
  q <- x$eta.q[[component]]
  if (!is.null(q$samples)) return(suppressWarnings(stats::density(q$samples, ...)))
  ## trim the quadrature grid to where the mass actually is -- its outer edge
  ## is set generously wide (scvi_update()'s search range), which would
  ## otherwise stretch the plotted x-axis far past where q(eta) has any mass.
  w  <- q$weights / sum(q$weights)
  cw <- cumsum(w[order(q$grid)])
  keep <- q$grid >= q$grid[order(q$grid)][which(cw >= 5e-4)[1L]] &
          q$grid <= q$grid[order(q$grid)][which(cw >= 1 - 5e-4)[1L]]
  suppressWarnings(stats::density(q$grid[keep], weights = w[keep] / sum(w[keep]), ...))
}

#' Diagnostic plots for an ngvb fit.
#'
#' For each component two panels are drawn: the mixing weights V_i/h_i (posterior
#' mean with a 90 percent credible interval, and the plug-in value actually refit
#' with INLA), and the posterior density of eta (with its mean and median marked).
#' A final panel shows how eta moved across the variational iterations. The plots
#' are ggplot2 objects assembled with patchwork.
#' @param x An `ngvb` object.
#' @param ... Ignored.
#' @return (Invisibly) the `ngvb` object; called for the plot it draws.
#' @export
plot.ngvb <- function(x, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("patchwork", quietly = TRUE))
    stop("plot.ngvb() needs the 'ggplot2' and 'patchwork' packages.")
  accent <- "#C64A2E"
  panels <- list()

  for (cn in x$comp.names) {
    Vs <- x$V.summary[[cn]]; h <- x$h[[cn]]; idx <- seq_along(h)
    df <- data.frame(idx = idx, mean = Vs$mean / h, lo = Vs$q05 / h,
                     hi = Vs$q95 / h, plug = x$V[[cn]] / h)
    pV <- ggplot2::ggplot(df, ggplot2::aes(x = .data$idx)) +
      ggplot2::geom_hline(yintercept = 1, linetype = 2, colour = "grey45") +
      ggplot2::geom_linerange(ggplot2::aes(ymin = .data$lo, ymax = .data$hi),
                              colour = "grey70") +
      ggplot2::geom_point(ggplot2::aes(y = .data$mean), size = 1.1, colour = "grey20") +
      ggplot2::geom_point(ggplot2::aes(y = .data$plug), shape = 4, size = 1.3,
                          colour = accent) +
      ggplot2::labs(x = "index", y = expression(V[i] / h[i]),
                    title = paste0("Mixing weights: ", cn),
                    subtitle = "point = mean, bar = 90% CI, x = INLA plug-in") +
      ngvb_theme()

    dens <- density.ngvb(x, cn)
    dd   <- data.frame(eta = dens$x, dens = dens$y)
    q    <- x$eta.q[[cn]]
    pE <- ggplot2::ggplot(dd, ggplot2::aes(x = .data$eta, y = .data$dens)) +
      ggplot2::geom_area(fill = "grey85", colour = "grey55") +
      ggplot2::geom_vline(xintercept = x$eta[[cn]], colour = accent, linewidth = 0.9) +
      ggplot2::geom_vline(xintercept = q$median, colour = accent, linetype = 2) +
      ggplot2::labs(x = expression(eta), y = "posterior density",
                    title = paste0("Non-Gaussianity: ", cn),
                    subtitle = sprintf("mean = %.3f (solid), median = %.3f (dashed)",
                                       x$eta[[cn]], q$median)) +
      ngvb_theme()
    panels <- c(panels, list(pV, pE))
  }

  eh <- as.data.frame(x$eta.hist); eh$iter <- seq_len(nrow(eh)) - 1L
  long <- stats::reshape(eh, varying = x$comp.names, v.names = "eta",
                         timevar = "component", times = x$comp.names,
                         direction = "long")
  pT <- ggplot2::ggplot(long, ggplot2::aes(x = .data$iter, y = .data$eta,
                                           colour = .data$component)) +
    ggplot2::geom_line() + ggplot2::geom_point(size = 1) +
    ggplot2::scale_colour_brewer(palette = "Dark2", name = NULL) +
    ggplot2::labs(x = "iteration", y = expression(E * group("[", eta, "]")),
                  title = "Convergence of E[eta]") +
    ngvb_theme()

  patchwork::wrap_plots(c(panels, list(pT)), ncol = 2, byrow = TRUE)
}

#' Fitted values of the underlying INLA fit.
#' @param object An `ngvb` object.
#' @param ... Ignored.
#' @export
fitted.ngvb <- function(object, ...) object$fit$summary.fitted.values
