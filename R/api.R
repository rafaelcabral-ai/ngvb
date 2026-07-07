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
  converged <- FALSE
  if (verbose) pb <- utils::txtProgressBar(min = 0, max = iter, style = 3)
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
        ## structured & collapsed VI (Theorem 2): E[eta] and E[1/V_i] by
        ## deterministic quadrature of the collapsed eta posterior (stable).
        upd     <- scvi_update(dk, hk, alpha.eta[k], Nk)
        V[[cn]] <- 1 / upd$EVm1
        eta[cn] <- upd$eta
      }
    }
    eta.hist <- rbind(eta.hist, eta)
    fit <- inla.fit.V(V)
    rel <- max(abs(eta - eta.hist[it, ]) / eta.hist[it, ])
    if (verbose) utils::setTxtProgressBar(pb, it)
    if (is.finite(rel) && rel < stop.rel.change) { converged <- TRUE; break }
  }
  if (verbose) {
    utils::setTxtProgressBar(pb, iter); close(pb)
    cat(sprintf("ngvb: %s after %d iteration(s);  E[eta] = %s\n",
                if (converged) "converged" else "reached the iteration limit",
                nrow(eta.hist) - 1L, paste(sprintf("%.3f", eta), collapse = ", ")))
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
#' @param method Variational algorithm: `"SVI"` (structured, the default -- robust)
#'   or `"SCVI"` (structured & collapsed; more accurate when non-Gaussianity is
#'   clearly present, but its collapsed eta-posterior can drift toward N/2 when the
#'   non-Gaussian signal is weak).
#' @param alpha.eta Exponential-PC-prior rate(s) on the non-Gaussianity parameter(s).
#' @param iter,stop.rel.change,n.sampling,verbose VB controls.
#' @return An object of class `ngvb` with the final INLA `fit`, the mixing
#'   vectors `V`, the non-Gaussianity parameters `eta`, and their trajectory.
#' @seealso [ng.check()], [ngvb_operator()]
#' @examples
#' \donttest{
#' if (requireNamespace("INLA", quietly = TRUE)) {
#'   set.seed(1); n <- 100
#'   x <- cumsum(rnorm(n, sd = 0.3)); x[50:n] <- x[50:n] + 6      # a level shift
#'   y <- x + rnorm(n, sd = 0.4)
#'   LGM  <- INLA::inla(y ~ -1 + f(i, model = "rw1", constr = TRUE),
#'                      data = data.frame(y = y, i = 1:n),
#'                      control.compute = list(config = TRUE))
#'   LnGM <- ngvb(LGM, iter = 5)     # non-Gaussian extension, model auto-detected
#'   summary(LnGM)
#' }
#' }
#' @export
ngvb <- function(fit, selection = NULL, components = NULL, method = c("SVI", "SCVI"),
                 alpha.eta = 1, iter = 20, stop.rel.change = 1e-3,
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
