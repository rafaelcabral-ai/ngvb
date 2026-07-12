## ---------------------------------------------------------------------------
## Posterior-over-V sampling and a proper Bayes factor.
##
## The variational fit gives q(V), an approximation to the posterior of the
## mixing variables. Sampling V from q(V) and refitting the conditionally-Gaussian
## LGM with INLA at each draw yields a set of ordinary `inla` fits that together
## represent the latent non-Gaussian model with V *integrated out*, rather than
## held at a single plug-in value. From those fits we can form a marginal-
## likelihood Bayes factor against the Gaussian model (V = h), and pool any
## posterior summary INLA reports.
##
## Marginal likelihood of the LnGM:
##   Z_LnGM = INT m(V) p(V | eta) dV,   m(V) = p(y | V) = exp(inla mlik at V),
## estimated by importance sampling with q(V) as proposal,
##   Z_LnGM ~ mean_m  m(V_m) p(V_m | eta) / q(V_m),   V_m ~ q(V).
## The Gaussian model has V = h, so Z_LGM = m(h). Both marginal likelihoods come
## from the *same* engine, so the Bayes factor Z_LnGM / Z_LGM is internally
## consistent. eta is held at its variational point estimate.
## ---------------------------------------------------------------------------

## log density of the GIG variational factor q(V_i) = GIG(lambda=-1, chi=b, psi=a).
#' @keywords internal
.ldgig <- function(x, lambda, chi, psi) {
  z <- sqrt(chi * psi)
  logK <- log(besselK(z, abs(lambda), expon.scaled = TRUE)) - z   # log K_lambda(z)
  0.5 * lambda * (log(psi) - log(chi)) - log(2) - logK +
    (lambda - 1) * log(x) - 0.5 * (chi / x + psi * x)
}

## log density of the inverse-Gaussian prior p(V_i | eta) = IG(mean = h, shape = h^2/eta).
#' @keywords internal
.ldig <- function(x, h, eta) {
  lambda <- h^2 / eta
  0.5 * (log(lambda) - log(2 * pi) - 3 * log(x)) - lambda * (x - h)^2 / (2 * h^2 * x)
}

#' Sample the mixing variables and refit INLA at each draw.
#'
#' Draws `n.samples` values of the mixing vector `V` from the variational
#' posterior `q(V)` of a fitted [ngvb()] object, and refits the
#' conditionally-Gaussian LGM with INLA at each draw. The returned object holds
#' the list of `inla` fits together with the importance weights needed to treat
#' `V` as integrated out (see [bayes.factor()], [summary.ngvb.samples()]).
#'
#' @param object An `ngvb` object.
#' @param n.samples Number of posterior draws of `V` (each is one `inla` refit).
#' @param seed Optional integer seed for reproducibility.
#' @param verbose Show a progress bar.
#' @return An object of class `ngvb.samples`: a list with the `inla` `fits`, the
#'   drawn `V`, per-draw log marginal likelihoods `logm` and log importance
#'   weights `logw`, and the Gaussian-model log marginal likelihood `logm.LGM`.
#' @seealso [bayes.factor()]
#' @export
ngvb_sample <- function(object, n.samples = 50, seed = NULL,
                        verbose = interactive()) {
  .need_inla()
  stopifnot(inherits(object, "ngvb"))
  if (is.null(object$inla.fit.V))
    stop("ngvb: this ngvb object predates sampling support; refit with ngvb().")
  if (!is.null(seed)) set.seed(seed)
  cn.all <- object$comp.names
  h  <- object$h; d <- object$d; eta <- object$eta
  ## GIG parameters of q(V_i) per component (a = 1/eta, b = d_i + h_i^2/eta)
  aV <- stats::setNames(lapply(cn.all, function(cn) rep(1 / eta[[cn]], length(h[[cn]]))), cn.all)
  bV <- stats::setNames(lapply(cn.all, function(cn) d[[cn]] + h[[cn]]^2 / eta[[cn]]), cn.all)

  on.exit(suppressWarnings(rm(list = paste0("ngvb.model.", cn.all), envir = globalenv())),
          add = TRUE)

  ## Gaussian model (V = h): marginal likelihood and fit
  fit.LGM  <- object$inla.fit.V(h)
  logm.LGM <- fit.LGM$mlik[1, 1]

  fits <- vector("list", n.samples)
  logm <- numeric(n.samples); logw <- numeric(n.samples)
  if (verbose) pb <- utils::txtProgressBar(min = 0, max = n.samples, style = 3)
  for (m in seq_len(n.samples)) {
    ## Draw V_i ~ GIG(-1, chi = b_i, psi = a_i) PER INDEX -- rgig does not
    ## vectorise over its parameters (see .rgig_vec / ngvb_V_summary).
    Vm <- stats::setNames(lapply(cn.all, function(cn) .rgig_vec(aV[[cn]], bV[[cn]])), cn.all)
    fm <- object$inla.fit.V(Vm)
    fits[[m]] <- fm
    logm[m] <- fm$mlik[1, 1]
    ## log importance weight: sum over components/indices of log p(V|eta) - log q(V)
    logw[m] <- sum(vapply(cn.all, function(cn)
      sum(.ldig(Vm[[cn]], h[[cn]], eta[[cn]]) -
          .ldgig(Vm[[cn]], -1, bV[[cn]], aV[[cn]])), numeric(1)))
    if (verbose) utils::setTxtProgressBar(pb, m)
  }
  if (verbose) close(pb)

  out <- list(fits = fits, V.samples = NULL, logm = logm, logw = logw,
              logm.LGM = logm.LGM, fit.LGM = fit.LGM,
              comp.names = cn.all, n.samples = n.samples)
  class(out) <- "ngvb.samples"
  out
}

#' @keywords internal
.logmeanexp <- function(x) { mx <- max(x); mx + log(mean(exp(x - mx))) }

#' Bayes factor of a latent non-Gaussian model against its Gaussian counterpart.
#'
#' Returns the marginal-likelihood Bayes factor \eqn{Z_{\mathrm{LnGM}}/Z_{\mathrm{LGM}}}
#' with the mixing variables `V` integrated out by importance sampling from
#' `q(V)`. A value above one favours the non-Gaussian model; the usual Jeffreys
#' reading is that \eqn{\log_{10}} above 0.5, 1, 2 is substantial, strong, decisive.
#'
#' @param x An `ngvb` object (sampled internally) or an `ngvb.samples` object.
#' @param ... For the `ngvb` method, passed to [ngvb_sample()] (e.g. `n.samples`, `seed`).
#' @return A list with `BF`, `log10BF`, `logBF`, and the effective sample size
#'   `ess` of the importance weights (a low `ess` means more draws are needed).
#' @export
bayes.factor <- function(x, ...) UseMethod("bayes.factor")

#' @rdname bayes.factor
#' @export
bayes.factor.ngvb <- function(x, ...) bayes.factor(ngvb_sample(x, ...))

#' @rdname bayes.factor
#' @export
bayes.factor.ngvb.samples <- function(x, ...) {
  lg   <- x$logm + x$logw                       # log of the per-draw IS terms for Z
  logZ <- .logmeanexp(lg)                        # log Z_LnGM
  logBF <- logZ - x$logm.LGM
  w <- exp(lg - max(lg)); w <- w / sum(w)        # ESS of the marginal-likelihood weights
  ess <- 1 / sum(w^2)
  out <- list(BF = exp(logBF), log10BF = logBF / log(10), logBF = logBF,
              ess = ess, n.samples = x$n.samples)
  class(out) <- "ngvb.bf"
  out
}

#' @export
print.ngvb.bf <- function(x, ...) {
  strength <- if (x$log10BF > 2) "decisive" else if (x$log10BF > 1) "strong" else
              if (x$log10BF > 0.5) "substantial" else if (x$log10BF > -0.5) "weak/none" else
              "favours the Gaussian model"
  cat(sprintf("Bayes factor (non-Gaussian vs Gaussian): %.3g\n", x$BF))
  cat(sprintf("  log10 BF = %.2f  (%s)\n", x$log10BF, strength))
  cat(sprintf("  weight ESS = %.1f of %d draws\n", x$ess, x$n.samples))
  invisible(x)
}

#' @export
print.ngvb.samples <- function(x, ...) {
  bf <- bayes.factor(x)
  cat("ngvb posterior-V samples\n")
  cat("  draws        :", x$n.samples, "\n")
  cat("  components   :", paste(x$comp.names, collapse = ", "), "\n")
  cat(sprintf("  Bayes factor : %.3g   (log10 = %.2f)\n", bf$BF, bf$log10BF))
  cat(sprintf("  weight ESS   : %.1f of %d draws\n", bf$ess, x$n.samples))
  invisible(x)
}

## Pool one INLA summary table (`"summary.fixed"` or `"summary.hyperpar"`)
## across the sampled fits with the normalized importance weights `w`.
## Returns NULL if that table is absent/empty in every draw (e.g. "fixed" on
## a model with no fixed effects) rather than erroring.
#' @keywords internal
.pool_ngvb_summary <- function(object, w, key) {
  tabs <- lapply(object$fits, function(f) f[[key]])
  ok <- !vapply(tabs, is.null, logical(1)) & vapply(tabs, function(t) nrow(t) > 0, logical(1))
  if (!any(ok)) return(NULL)
  rn <- rownames(tabs[[which(ok)[1]]])
  mean.mat <- vapply(tabs[ok], function(t) t[rn, "mean"], numeric(length(rn)))
  sd.mat   <- vapply(tabs[ok], function(t) t[rn, "sd"],   numeric(length(rn)))
  ww <- w[ok] / sum(w[ok])
  pooled.mean <- as.numeric(mean.mat %*% ww)
  ## law of total variance across the mixture of fits
  pooled.var  <- as.numeric((sd.mat^2 + mean.mat^2) %*% ww) - pooled.mean^2
  data.frame(mean = round(pooled.mean, 4), sd = round(sqrt(pmax(0, pooled.var)), 4),
             row.names = rn)
}

## Pool summary.random$<comp.name> across the sampled fits (same importance-
## weighted mean/sd as .pool_ngvb_summary(), but summary.random is a *list* of
## per-component tables keyed by node ID rather than one shared table).
#' @keywords internal
.pool_ngvb_random <- function(object, w, comp.name) {
  tabs <- lapply(object$fits, function(f) f$summary.random[[comp.name]])
  ok <- !vapply(tabs, is.null, logical(1)) & vapply(tabs, function(t) nrow(t) > 0, logical(1))
  if (!any(ok)) return(NULL)
  rn       <- tabs[[which(ok)[1]]]$ID
  mean.mat <- vapply(tabs[ok], function(t) t$mean, numeric(length(rn)))
  sd.mat   <- vapply(tabs[ok], function(t) t$sd,   numeric(length(rn)))
  ww <- w[ok] / sum(w[ok])
  pooled.mean <- as.numeric(mean.mat %*% ww)
  ## law of total variance across the mixture of fits
  pooled.var  <- as.numeric((sd.mat^2 + mean.mat^2) %*% ww) - pooled.mean^2
  data.frame(ID = rn, mean = round(pooled.mean, 4), sd = round(sqrt(pmax(0, pooled.var)), 4))
}

#' Importance-weighted posterior summary of the sampled fits.
#'
#' Pools a chosen INLA summary across the sampled fits with the importance
#' weights, giving posterior means and standard deviations for the latent
#' non-Gaussian model (with `V` integrated out), plus the Bayes factor.
#'
#' @param object An `ngvb.samples` object.
#' @param what One of `"all"` (default), `"fixed"`, `"hyperpar"`, or `"random"`
#'   -- which INLA summary table(s) to pool. `"random"` pools `summary.random`
#'   for every component (the latent field itself: one table per component,
#'   keyed by node `ID`, e.g. mesh node or time index). `"all"` shows and
#'   returns all three, silently skipping any that don't apply to this model
#'   (e.g. `"fixed"` on a fixed-effect-free model); random-effect tables are
#'   printed as a head (all rows are still returned). Run
#'   `args(summary.ngvb.samples)` or `?summary.ngvb.samples` to see this list
#'   again.
#' @param verbose If `TRUE` (default), print the report as a side effect (as
#'   `ng.check()` does for its plot). Set `FALSE` to only get the return
#'   value back, e.g. when pooling `what = "random"` into a plot without the
#'   table being echoed.
#' @param ... Ignored.
#' @return Invisibly: a single data frame of importance-weighted posterior
#'   means/sds for `what = "fixed"` or `"hyperpar"`; a named list of one such
#'   data frame per component (keyed by node `ID`) for `what = "random"`; for
#'   `"all"`, a list `list(fixed = , hyperpar = , random = )` with any
#'   inapplicable element `NULL`. Printed as a side effect when `verbose = TRUE`.
#' @method summary ngvb.samples
#' @export
summary.ngvb.samples <- function(object, what = c("all", "fixed", "hyperpar", "random"),
                                 verbose = TRUE, ...) {
  what <- match.arg(what)
  bf <- bayes.factor(object)
  if (verbose) {
    cat("Latent non-Gaussian model, V integrated out over", object$n.samples, "draws\n")
    cat(sprintf("Bayes factor vs Gaussian model: %.3g (log10 = %.2f), weight ESS %.1f\n\n",
                bf$BF, bf$log10BF, bf$ess))
  }
  ## normalized importance weights (marginal-likelihood weighted)
  lw <- object$logm + object$logw; w <- exp(lw - max(lw)); w <- w / sum(w)

  show_one <- function(label, key) {
    res <- .pool_ngvb_summary(object, w, key)
    if (verbose) {
      cat(label, ":\n", sep = "")
      if (is.null(res)) cat("  (none)\n") else print(res)
      cat("\n")
    }
    res
  }

  show_random <- function(full) {
    res <- stats::setNames(
      lapply(object$comp.names, function(cn) .pool_ngvb_random(object, w, cn)),
      object$comp.names)
    if (verbose) {
      cat("Random effects:\n")
      for (cn in object$comp.names) {
        tab <- res[[cn]]
        if (is.null(tab)) { cat("  $", cn, ": (none)\n", sep = ""); next }
        cat("  $", cn, " (", nrow(tab), " nodes)\n", sep = "")
        print(if (full) tab else utils::head(tab, 6))
        if (!full && nrow(tab) > 6)
          cat("  ... ", nrow(tab) - 6,
              " more row(s); use what = \"random\" to print in full\n", sep = "")
      }
      cat("\n")
    }
    res
  }

  if (what == "all") {
    res <- list(fixed    = show_one("Fixed effects", "summary.fixed"),
               hyperpar = show_one("Hyperparameters", "summary.hyperpar"),
               random   = show_random(full = FALSE))
    return(invisible(res))
  }

  if (what == "random") return(invisible(show_random(full = TRUE)))

  key <- if (what == "fixed") "summary.fixed" else "summary.hyperpar"
  res <- .pool_ngvb_summary(object, w, key)
  if (verbose) {
    if (is.null(res)) cat("(no", what, "effects)\n") else print(res)
  }
  invisible(res)
}
