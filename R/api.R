## ---------------------------------------------------------------------------
## User-facing API: ngvb(fit). Auto-detects each selected component, rebuilds the
## INLA formula so those components use the ngvb engine (conditioned on V), and
## runs the multi-component structured-VI loop (block-diagonal D; one joint
## inla() call per iteration -- the fix for the old single-component loop bug).
## ---------------------------------------------------------------------------

#' Rewrite an INLA formula so selected f(<name>, ...) terms use the ngvb engine.
#' The engine objects live in `model.env` as `ngvb.model.<name>`; model-specific
#' args (hyper, the original model) are dropped, positional args (index, optional
#' weight) preserved, and `constr = TRUE` added for intrinsic components.
#' @keywords internal
ngvb_swap_formula <- function(formula, comp.names, rankdef.map, model.env) {
  walk <- function(e) {
    if (!is.call(e)) return(e)
    if (identical(e[[1L]], as.name("f"))) {
      idx <- tryCatch(as.character(e[[2L]]), error = function(...) "")
      if (length(idx) == 1L && idx %in% comp.names) {
        al <- as.list(e)[-1L]
        nm <- names(al); if (is.null(nm)) nm <- rep("", length(al))
        positional <- al[nm == ""]
        keep <- positional[seq_len(min(2L, length(positional)))]   # index (+ optional weight)
        newargs <- c(keep, list(model = as.name(paste0("ngvb.model.", idx))))
        if (isTRUE(rankdef.map[[idx]] > 0L)) newargs <- c(newargs, list(constr = TRUE))
        return(as.call(c(list(as.name("f")), newargs)))
      }
      return(e)
    }
    for (i in seq_along(e)) e[[i]] <- walk(e[[i]])
    e
  }
  formula[[length(formula)]] <- walk(formula[[length(formula)]])   # rewrite the RHS
  environment(formula) <- model.env
  formula
}

#' Build the conditional-on-V INLA fitter for a set of components.
#' @keywords internal
ngvb_make_fit_V <- function(fit, ops, comp.names) {
  arglist <- fit$.args
  arglist$control.compute$config <- TRUE
  rankdef.map <- stats::setNames(lapply(ops, `[[`, "rankdef"), comp.names)
  ## INLA resolves the formula's `model =` argument (and re-evaluates data) in the
  ## global environment, so the engine objects must live there (cleaned up by ngvb()).
  model.env   <- globalenv()
  arglist$formula <- ngvb_swap_formula(arglist$formula, comp.names, rankdef.map, model.env)
  function(V) {
    for (cn in comp.names)
      assign(paste0("ngvb.model.", cn), ngvb_rgeneric(ops[[cn]], V = V[[cn]]), envir = model.env)
    do.call(INLA::inla, arglist)
  }
}

#' Multi-component structured variational inference loop.
#' @keywords internal
ngvb_vb <- function(inla.fit.V, ops, comp.names, method = c("SVI", "SCVI"),
                    alpha.eta = 2, identify.scale = TRUE, iter = 10,
                    stop.rel.change = 1e-3, n.sampling = 2000, verbose = TRUE) {
  method <- match.arg(method)
  ncomp  <- length(comp.names)
  if (length(alpha.eta) == 1L) alpha.eta <- rep(alpha.eta, ncomp)
  h      <- stats::setNames(lapply(ops, `[[`, "h"), comp.names)
  V      <- h
  eta    <- stats::setNames(rep(0.5, ncomp), comp.names)
  Eetam1 <- stats::setNames(1 / eta, comp.names)
  eta.hist <- matrix(eta, nrow = 1, dimnames = list(NULL, comp.names))
  eta.q  <- stats::setNames(vector("list", ncomp), comp.names)   # q(eta) summary, last iteration

  fit <- inla.fit.V(V)
  d   <- NULL
  converged <- FALSE
  ## progress bar only in an interactive session; its carriage returns clutter
  ## knitr / script logs, whereas the convergence summary below stays on `verbose`.
  show.pb <- isTRUE(verbose) && interactive()
  if (show.pb) pb <- utils::txtProgressBar(min = 0, max = iter, style = 3)
  for (it in seq_len(iter)) {
    ## d_i = E[(Dx)_i^2] >= 0 in exact arithmetic; clamp tiny negative round-off
    ## (seen on near-Gaussian components) so the GIG scale b_V = d + h^2/eta stays
    ## positive and the variational moments are well defined.
    d <- stats::setNames(lapply(comp.names, function(cn)
      pmax(compute_d(fit, ops[[cn]], cn), 0)), comp.names)
    ## Scale-identification constraint. The conditional precision
    ## Q = tau * D0' diag(1/V) D0 is invariant under (tau, V) -> (c*tau, c*V),
    ## an exact flat ridge broken only softly by the V-prior. Re-anchoring each
    ## component's discrepancy so sum(d) = sum(h) pins that scale to tau, leaving
    ## eta to respond only to the RELATIVE pattern of the increments (genuine
    ## outliers), not their absolute level -- which is tau's job. This removes the
    ## tau<->eta runaway (e.g. intrinsic CAR) without touching well-scaled fits.
    if (identify.scale)
      d <- stats::setNames(lapply(comp.names, function(cn) {
        s <- sum(d[[cn]])
        if (is.finite(s) && s > 0) d[[cn]] * sum(h[[cn]]) / s else d[[cn]]
      }), comp.names)
    for (k in seq_len(ncomp)) {
      cn <- comp.names[k]; hk <- h[[cn]]; Nk <- length(hk); dk <- d[[cn]]
      if (method == "SVI") {
        ## structured VI: closed-form GIG variational factors (Theorem 1)
        a_V  <- Eetam1[cn]; b_V <- dk + hk^2 * Eetam1[cn]
        EV   <- GIGM1(-1, a_V, b_V); EVm1 <- GIGMm1(-1, a_V, b_V)
        p_e  <- -Nk / 2 + 1; a_e <- 2 * alpha.eta[k]; b_e <- sum(EV - 2 * hk + hk^2 * EVm1)
        gm   <- gig_moments(p_e, a_e, b_e, n = n.sampling)   # one draw, both moments + sd/CI
        eta[cn]    <- gm$mean
        Eetam1[cn] <- gm$inv_mean
        V[[cn]]    <- 1 / EVm1
        eta.q[[cn]] <- gm
      } else {
        ## structured & collapsed VI (Theorem 2): E[eta] and E[1/V_i] by
        ## deterministic quadrature of the collapsed eta posterior (stable).
        upd     <- scvi_update(dk, hk, alpha.eta[k], Nk)
        V[[cn]] <- 1 / upd$EVm1
        eta[cn] <- upd$eta
        eta.q[[cn]] <- upd[c("eta", "median", "sd", "q05", "q95", "grid", "weights")]
      }
    }
    eta.hist <- rbind(eta.hist, eta)
    fit <- inla.fit.V(V)
    rel <- max(abs(eta - eta.hist[it, ]) / eta.hist[it, ])
    if (show.pb) utils::setTxtProgressBar(pb, it)
    if (is.finite(rel) && rel < stop.rel.change) { converged <- TRUE; break }
  }
  if (show.pb) { utils::setTxtProgressBar(pb, iter); close(pb) }
  if (verbose) {
    cat(sprintf("ngvb: %s after %d iteration(s);  E[eta] = %s\n",
                if (converged) "converged" else "reached the iteration limit",
                nrow(eta.hist) - 1L, paste(sprintf("%.3f", eta), collapse = ", ")))
  }
  ngvb_check_degeneracy(V, h, comp.names, alpha.eta, verbose)
  ## Per-index V_i summary (mean exact; median/90% CI by one-time sampling),
  ## conditional on the final eta point estimate -- NOT the same quantity as
  ## `V` itself, which is 1/E[1/V_i] (the precision-consistent plug-in the
  ## algorithm actually re-fits INLA with, pulled below the mean by Jensen's
  ## inequality for a right-skewed q(V)). See plot.ngvb()/summary.ngvb().
  V.summary <- stats::setNames(lapply(comp.names, function(cn) {
    a_V <- 1 / eta[[cn]]; b_V <- d[[cn]] + h[[cn]]^2 * a_V
    ngvb_V_summary(a_V, b_V, n = n.sampling)
  }), comp.names)
  out <- list(fit = fit, V = V, V.summary = V.summary, eta = eta, eta.q = eta.q,
              h = h, d = d, eta.hist = eta.hist, alpha.eta = alpha.eta,
              ops = ops, comp.names = comp.names, iterations = nrow(eta.hist) - 1L,
              inla.fit.V = inla.fit.V)   # reused by ngvb_sample() to refit at drawn V
  class(out) <- "ngvb"
  out
}

#' Warn when a component's mixing weights have shrunk nearly uniformly across
#' (almost) every index, rather than concentrating on a few outliers.
#'
#' SVI and SCVI target the same fixed point and, when a component's per-index
#' discrepancies d_i are small but roughly HOMOGENEOUS (no real outlier/inlier
#' split), that shared fixed point can drift to a large, mostly prior- and
#' N-driven eta rather than one reflecting genuine local non-Gaussianity: tau
#' and eta become weakly identified against each other and reinforce each
#' other iteration over iteration. A uniform V/h << 1 across almost the whole
#' component is the signature of that regime, not of real outlier detection.
#' @keywords internal
ngvb_check_degeneracy <- function(V, h, comp.names, alpha.eta, verbose,
                                  frac.threshold = 0.8, ratio.threshold = 0.5,
                                  inflate.threshold = 2) {
  for (k in seq_along(comp.names)) {
    cn <- comp.names[k]
    r  <- V[[cn]] / h[[cn]]
    frac.low <- mean(r < ratio.threshold)
    ## Degenerate = almost everything shrunk AND nothing inflated. A genuine
    ## heavy-tailed fit also shrinks most indices below h, but it pays for that
    ## with a few strongly inflated ones (V/h large); the presence of any real
    ## outlier (max V/h above inflate.threshold) rules degeneracy out.
    if (frac.low > frac.threshold && max(r) < inflate.threshold) {
      msg <- sprintf(paste0(
        "ngvb: component '%s' -- %.0f%% of indices show V/h < %.2g (mostly-uniform ",
        "shrinkage, not a few outliers). This is the signature of a weakly-",
        "identified fit (eta and the component's precision compete to explain ",
        "the same homogeneous discrepancy) rather than genuine non-Gaussianity. ",
        "Consider a stronger alpha.eta (currently %.3g) to keep the fit out of ",
        "this regime, and inspect plot(fit)'s per-index panel for '%s'."),
        cn, 100 * frac.low, ratio.threshold, alpha.eta[k], cn)
      if (verbose) warning(msg, call. = FALSE, immediate. = TRUE) else warning(msg, call. = FALSE)
    }
  }
  invisible(NULL)
}

#' Fit a latent non-Gaussian model from a fitted INLA (LGM) object.
#'
#' @param fit An `inla` object (fitted with `control.compute = list(config = TRUE)`).
#' @param selection Optional named list of component indices to extend to
#'   non-Gaussianity (default: all random effects).
#' @param components Optional named list of operator descriptors overriding
#'   auto-detection (required for SPDE: `list(s = ngvb_operator("spde", spde = spde))`).
#' @param method Variational algorithm: `"SVI"` (structured, mean-field; the
#'   default -- more reliable) or `"SCVI"` (structured & collapsed; reaches a
#'   fixed point in far fewer iterations, but its collapsed marginal can
#'   over-shrink a weak-but-genuine effect all the way to Gaussian, `eta ~ 0`,
#'   where SVI holds a small positive value). The two agree when the signal is
#'   strong or clearly absent and disagree in the weak-signal regime; prefer the
#'   default `"SVI"` unless you specifically need SCVI's speed and have checked
#'   the two give the same answer on your problem. Regardless of method, when a
#'   component's discrepancies are small but roughly homogeneous across all its
#'   indices (no real outlier/inlier split), the fit can drift toward a large,
#'   mostly prior-driven eta rather than genuine non-Gaussianity -- see
#'   `alpha.eta` and the degeneracy warning this function may emit.
#' @param alpha.eta Exponential-PC-prior rate(s) on the non-Gaussianity
#'   parameter(s): larger values shrink harder toward the Gaussian model
#'   (`eta = 0`). If a fit triggers the degeneracy warning, raise this (e.g.
#'   by 2-5x) before trusting the result.
#' @param identify.scale Enforce the scale-identifiability constraint (default
#'   `TRUE`). The conditional precision is invariant under
#'   \eqn{(\tau, \mathbf V) \mapsto (c\tau, c\mathbf V)}, an exact flat ridge in
#'   the likelihood that the VB loop can otherwise walk -- letting the mixing
#'   variables absorb the overall scale and inflating `eta` spuriously (most
#'   visibly for intrinsic models such as `besag`/`rw`). With `identify.scale`
#'   on, each component's discrepancies are re-anchored so their total matches
#'   the Gaussian total each iteration, pinning the scale to `tau` and leaving
#'   `eta` to respond only to *relative* departures (genuine outliers). Turn off
#'   only to reproduce the unconstrained updates.
#' @param iter,stop.rel.change,n.sampling,verbose VB controls.
#' @return An object of class `ngvb` with the final INLA `fit`, the mixing
#'   vectors `V`, the non-Gaussianity parameters `eta`, and their trajectory.
#' @seealso [ng.check()], [ngvb_operator()]
#' @examples
#' \donttest{
#' if (requireNamespace("INLA", quietly = TRUE)) {
#'   data(jumpts)   # a series with two abrupt jumps -- see the package vignette
#'   LGM  <- INLA::inla(y ~ -1 + f(x, model = "rw1"), data = jumpts,
#'                      control.compute = list(config = TRUE))
#'   LnGM <- ngvb(LGM, iter = 10)     # non-Gaussian extension, model auto-detected
#'   summary(LnGM)
#' }
#' }
#' @export
ngvb <- function(fit, selection = NULL, components = NULL, method = c("SVI", "SCVI"),
                 alpha.eta = 2, identify.scale = TRUE, iter = 30,
                 stop.rel.change = 1e-3, n.sampling = 2000, verbose = TRUE) {
  .need_inla()
  method <- match.arg(method)
  if (is.null(fit$misc$configs))
    stop("ngvb2: refit the LGM with control.compute = list(config = TRUE).")
  if (is.null(selection))
    selection <- lapply(fit$summary.random, function(x) seq_len(nrow(x)))
  comp.names <- names(selection)
  if (length(comp.names) == 0L) stop("ngvb2: no random components found in the fit.")
  on.exit(suppressWarnings(rm(list = paste0("ngvb.model.", comp.names), envir = globalenv())),
          add = TRUE)

  ops <- stats::setNames(
    lapply(comp.names, function(cn) ngvb_detect_operator(fit, cn, user.op = components[[cn]])),
    comp.names)
  if (verbose)
    cat("Components:", paste(sprintf("%s [%s]", comp.names, vapply(ops, `[[`, "", "type")),
                             collapse = ", "), "\n")

  inla.fit.V <- ngvb_make_fit_V(fit, ops, comp.names)
  ngvb_vb(inla.fit.V, ops, comp.names, method = method, alpha.eta = alpha.eta,
          identify.scale = identify.scale, iter = iter,
          stop.rel.change = stop.rel.change,
          n.sampling = n.sampling, verbose = verbose)
}
