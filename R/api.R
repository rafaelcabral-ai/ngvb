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
ngvb_vb <- function(inla.fit.V, ops, comp.names, method = c("SCVI", "SVI"),
                    alpha.eta = 1, iter = 10, stop.rel.change = 1e-3,
                    n.sampling = 2000, verbose = TRUE) {
  method <- match.arg(method)
  ncomp  <- length(comp.names)
  if (length(alpha.eta) == 1L) alpha.eta <- rep(alpha.eta, ncomp)
  h      <- stats::setNames(lapply(ops, `[[`, "h"), comp.names)
  V      <- h
  eta    <- stats::setNames(rep(0.5, ncomp), comp.names)
  Eetam1 <- stats::setNames(1 / eta, comp.names)
  eta.hist <- matrix(eta, nrow = 1, dimnames = list(NULL, comp.names))

  fit <- inla.fit.V(V)
  d   <- NULL
  for (it in seq_len(iter)) {
    d <- stats::setNames(lapply(comp.names, function(cn) compute_d(fit, ops[[cn]], cn)), comp.names)
    for (k in seq_len(ncomp)) {
      cn <- comp.names[k]; hk <- h[[cn]]; Nk <- length(hk); dk <- d[[cn]]
      if (method == "SVI") {
        ## structured VI: closed-form GIG variational factors (Theorem 1)
        a_V  <- Eetam1[cn]; b_V <- dk + hk^2 * Eetam1[cn]
        EV   <- GIGM1(-1, a_V, b_V); EVm1 <- GIGMm1(-1, a_V, b_V)
        p_e  <- -Nk / 2 + 1; a_e <- 2 * alpha.eta[k]; b_e <- sum(EV - 2 * hk + hk^2 * EVm1)
        eta[cn]    <- mGIG(p_e, a_e, b_e, order =  1L, n = n.sampling)
        Eetam1[cn] <- mGIG(p_e, a_e, b_e, order = -1L, n = n.sampling)
        V[[cn]]    <- 1 / EVm1
      } else {
        ## structured & collapsed VI (Theorem 2): sample eta from its collapsed
        ## posterior, then E[1/V_i] = mean over eta-samples of the closed-form
        ## GIG(-1, 1/eta, d_i + h_i^2/eta) reciprocal moment (no GIG sampling).
        es   <- sampler.inverseCDF(eta.prior.f(dk, hk, alpha.eta[k], Nk),
                                   supp.min = 0, supp.max = 100, n.samples = n.sampling)
        es   <- es[is.finite(es) & es > 0]
        EVm1 <- vapply(seq_along(hk),
                       function(i) mean(GIGMm1(-1, 1 / es, dk[i] + hk[i]^2 / es)), numeric(1))
        V[[cn]] <- 1 / EVm1
        eta[cn] <- mean(es)
      }
    }
    eta.hist <- rbind(eta.hist, eta)
    fit <- inla.fit.V(V)
    rel <- max(abs(eta - eta.hist[it, ]) / eta.hist[it, ])
    if (verbose) cat(sprintf("  iter %2d:  E[eta] = %s   (max rel change %.4f)\n",
                             it, paste(sprintf("%.3f", eta), collapse = ", "), rel))
    if (is.finite(rel) && rel < stop.rel.change) { if (verbose) cat("  converged\n"); break }
  }
  out <- list(fit = fit, V = V, eta = eta, h = h, d = d, eta.hist = eta.hist,
              ops = ops, comp.names = comp.names, iterations = nrow(eta.hist) - 1L)
  class(out) <- "ngvb"
  out
}

#' Fit a latent non-Gaussian model from a fitted INLA (LGM) object.
#'
#' @param fit An `inla` object (fitted with `control.compute = list(config = TRUE)`).
#' @param selection Optional named list of component indices to extend to
#'   non-Gaussianity (default: all random effects).
#' @param components Optional named list of operator descriptors overriding
#'   auto-detection (required for SPDE: `list(s = ngvb_operator("spde", spde = spde))`).
#' @param method Variational algorithm: `"SCVI"` (structured & collapsed, the default)
#'   or `"SVI"` (structured).
#' @param alpha.eta Exponential-PC-prior rate(s) on the non-Gaussianity parameter(s).
#' @param iter,stop.rel.change,n.sampling,verbose VB controls.
#' @return A list with the final `fit`, mixing vectors `V`, `eta`, etc.
#' @export
ngvb <- function(fit, selection = NULL, components = NULL, method = c("SCVI", "SVI"),
                 alpha.eta = 1, iter = 10, stop.rel.change = 1e-3,
                 n.sampling = 2000, verbose = TRUE) {
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
          iter = iter, stop.rel.change = stop.rel.change,
          n.sampling = n.sampling, verbose = verbose)
}
