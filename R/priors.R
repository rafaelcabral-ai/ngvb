## ---------------------------------------------------------------------------
## Carrying a user's f()-term precision prior into the ngvb engine.
##
## ngvb() rewrites each selected f(<name>, model = ..., hyper = ...) term to use
## the rgeneric engine, which drops the original model-specific arguments. The
## structure (model/graph/scale.model) MUST be dropped -- the engine rebuilds it
## -- but the *precision prior* the user set in `hyper` should be honored where we
## can. The engine's precision hyperparameter (theta[1] = log tau) takes any log
## prior, so we translate the common families (pc.prec, loggamma) into the
## operator's `logprior`; anything else (an unknown family, a fixed precision, or
## a prior on a secondary hyperparameter such as an AR1 correlation) is dropped
## with a warning and the operator's default PC prior is used instead.
## ---------------------------------------------------------------------------

## log prior of theta = log(tau) for a Gamma(shape = a, rate = b) precision
## ("loggamma" in INLA), including the d tau / d theta = tau Jacobian:
##   p(theta) = Gamma(exp(theta); a, b) * exp(theta)
#' @keywords internal
.loggamma_prec_logprior <- function(a, b) {
  function(lprec) a * log(b) - lgamma(a) + a * lprec - b * exp(lprec)
}

## Translate one INLA precision-hyper spec (list(prior=, param=, fixed=)) into a
## log-prior closure on theta = log(tau). Returns a list with the closure `lp`
## (NULL if unmappable), the `family` name, and a `mappable` flag.
#' @keywords internal
.prec_logprior_from_hyper <- function(spec) {
  if (is.null(spec) || !is.list(spec)) return(NULL)
  if (isTRUE(spec$fixed)) return(list(lp = NULL, family = "fixed", mappable = FALSE))
  prior <- spec$prior; param <- spec$param
  if (is.null(prior)) return(NULL)
  fam <- tolower(prior)
  ok2 <- is.numeric(param) && length(param) >= 2L
  if (fam %in% c("pc.prec", "pcprec")) {
    if (!ok2) return(list(lp = NULL, family = prior, mappable = FALSE))
    return(list(lp = .pc_prec_logprior(param[1L], param[2L]),
                family = "pc.prec", mappable = TRUE))
  }
  if (fam %in% c("loggamma", "log.gamma")) {
    if (!ok2) return(list(lp = NULL, family = prior, mappable = FALSE))
    return(list(lp = .loggamma_prec_logprior(param[1L], param[2L]),
                family = "loggamma", mappable = TRUE))
  }
  list(lp = NULL, family = prior, mappable = FALSE)
}

## Extract the (evaluated) `hyper` argument of each selected f(<name>, ...) term
## from a model formula. Returns a named list (one entry per comp.name, NULL when
## the term has no `hyper`).
#' @keywords internal
ngvb_extract_f_hyper <- function(formula, comp.names) {
  env <- environment(formula); if (is.null(env)) env <- parent.frame()
  out <- stats::setNames(vector("list", length(comp.names)), comp.names)
  walk <- function(e) {
    if (!is.call(e)) return(invisible())
    if (identical(e[[1L]], as.name("f"))) {
      al <- as.list(e)[-1L]
      nm <- names(al); if (is.null(nm)) nm <- rep("", length(al))
      pos <- al[nm == ""]
      idx <- if (length(pos)) tryCatch(as.character(pos[[1L]]), error = function(...) "") else ""
      if (length(idx) == 1L && idx %in% comp.names && "hyper" %in% nm)
        out[[idx]] <<- tryCatch(eval(al[["hyper"]], env), error = function(...) NULL)
    }
    for (i in seq_along(e)) if (is.call(e[[i]])) walk(e[[i]])
    invisible()
  }
  walk(formula[[length(formula)]])
  out
}

## Rebuild op$logprior so its precision piece is the user's mapped prior, leaving
## any secondary-hyperparameter prior intact. Warns (and falls back to the default)
## when the prior cannot be carried over.
#' @keywords internal
ngvb_apply_user_prec_prior <- function(op, hyper, comp.name, verbose = TRUE) {
  if (is.null(hyper) || !length(hyper)) return(op)
  warn <- function(msg)
    if (isTRUE(verbose)) warning(msg, call. = FALSE, immediate. = TRUE)
    else warning(msg, call. = FALSE)

  prec.names <- c("prec", "theta1", "theta")
  which.prec <- intersect(prec.names, names(hyper))
  spec       <- if (length(which.prec)) hyper[[which.prec[1L]]] else NULL
  secondary  <- setdiff(names(hyper), which.prec)

  if (!is.null(spec)) {
    map <- .prec_logprior_from_hyper(spec)
    if (is.null(op$prec.logprior)) {
      warn(sprintf(paste0("ngvb2: a precision prior was set on '%s' in f(), but the '%s' ",
                          "engine is not parameterized by a single precision, so it cannot be ",
                          "carried over; using the operator's default prior."),
                   comp.name, op$type))
    } else if (isTRUE(map$mappable)) {
      base <- op$logprior; old <- op$prec.logprior; new <- map$lp
      op$logprior     <- function(theta) base(theta) - old(theta[1L]) + new(theta[1L])
      op$prec.logprior <- new
      if (isTRUE(verbose))
        message(sprintf("ngvb2: carried the %s precision prior from f() into component '%s'.",
                        map$family, comp.name))
    } else {
      reason <- if (!is.null(map) && identical(map$family, "fixed"))
        "a fixed precision is not supported by the ngvb engine"
      else sprintf("the '%s' prior family is not one ngvb can map (only pc.prec and loggamma)",
                   if (is.null(map)) "unknown" else map$family)
      warn(sprintf(paste0("ngvb2: the precision prior on '%s' was dropped -- %s. Using the ",
                          "default PC prior; set it explicitly with components = list(%s = ",
                          "ngvb_operator(..., pc.prec = c(U = , alpha = )))."),
                   comp.name, reason, comp.name))
    }
  }
  if (length(secondary))
    warn(sprintf(paste0("ngvb2: prior(s) on %s for component '%s' were not carried over ",
                        "(ngvb maps only the precision prior); the engine default is used."),
                 paste(sprintf("'%s'", secondary), collapse = ", "), comp.name))
  op
}
