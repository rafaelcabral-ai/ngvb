## S3 methods for the object returned by ngvb().

#' @export
print.ngvb <- function(x, ...) {
  cat("Latent non-Gaussian model fit (ngvb2)\n")
  cat("  components:", paste(sprintf("%s [%s]", x$comp.names,
                                     vapply(x$ops, `[[`, "", "type")), collapse = ", "), "\n")
  cat("  iterations:", x$iterations, "\n")
  cat("  E[eta]    :", paste(sprintf("%s=%.3f", names(x$eta), x$eta), collapse = ", "), "\n")
  invisible(x)
}

#' Summary of an ngvb fit: the underlying INLA summaries plus, per component, the
#' non-Gaussianity parameter and the most strongly down-weighted (outlying) indices.
#' @param object An `ngvb` object.
#' @param n.flag Number of top mixing-variable indices to report per component.
#' @param ... Ignored.
#' @method summary ngvb
#' @export
summary.ngvb <- function(object, n.flag = 5, ...) {
  cat("Latent non-Gaussian model (ngvb2)\n\n")
  cat("Fixed effects:\n"); print(round(object$fit$summary.fixed, 4))
  cat("\nNon-Gaussianity:\n")
  for (cn in object$comp.names) {
    V <- object$V[[cn]]; h <- object$h[[cn]]
    ord <- order(V / h, decreasing = TRUE)[seq_len(min(n.flag, length(V)))]
    cat(sprintf("  %-12s [%s]: E[eta] = %.3f;  most non-Gaussian indices (V/h): %s\n",
                cn, object$ops[[cn]]$type, object$eta[[cn]],
                paste(sprintf("%d (%.1f)", ord, (V / h)[ord]), collapse = ", ")))
  }
  invisible(object)
}

#' Diagnostic plots for an ngvb fit: mixing variables per component and the
#' evolution of the non-Gaussianity parameters.
#' @param x An `ngvb` object.
#' @param ... Passed to `plot`.
#' @export
plot.ngvb <- function(x, ...) {
  ncomp <- length(x$comp.names)
  op <- graphics::par(mfrow = c(1, ncomp + 1)); on.exit(graphics::par(op))
  for (cn in x$comp.names) {
    plot(x$V[[cn]] / x$h[[cn]], type = "h", xlab = "index", ylab = "V / h",
         main = paste("mixing weights:", cn), ...)
    graphics::abline(h = 1, lty = 2, col = "grey")
  }
  matplot(0:(nrow(x$eta.hist) - 1), x$eta.hist, type = "b", pch = 1,
          xlab = "iteration", ylab = expression(E * "[" * eta * "]"),
          main = "non-Gaussianity")
  graphics::legend("topleft", legend = x$comp.names, col = seq_len(ncomp), lty = 1, bty = "n")
  invisible(x)
}

#' Fitted values of the underlying INLA fit.
#' @param object An `ngvb` object.
#' @param ... Ignored.
#' @export
fitted.ngvb <- function(object, ...) object$fit$summary.fitted.values
