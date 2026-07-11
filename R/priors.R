## ---------------------------------------------------------------------------
## Carrying the LGM's precision priors into the ngvb engine.
##
## ngvb() rewrites each selected f() term to use the rgeneric engine, which
## drops the original model-specific arguments. The structure (model/graph/
## scale.model) MUST be dropped -- the engine rebuilds it -- but the *precision
## prior* of the original fit is honored where we can, so the LnGM is the exact
## non-Gaussian extension of the LGM the user fitted.
##
## The prior is read from fit$all.hyper$random -- INLA's own normalized record
## of the prior it ACTUALLY used for every hyperparameter (user-set or default,
## params fully resolved). It is deliberately NOT parsed from the model
## formula: re-evaluating f(..., hyper = ) in environment(formula) breaks
## silently when ngvb() is called from inside a function (the argument no
## longer resolves in the stored environment), which used to produce
## default-prior fits with no indication anything was wrong.
##
## The engine's precision hyperparameter (theta[1] = log tau) takes any log
## prior; we translate the families INLA uses most: pc.prec, loggamma, and
## normal (Gaussian on log precision). Anything else (an unknown family or a
## fixed precision) falls back to the operator's default PC prior WITH A
## WARNING -- never silently.
## ---------------------------------------------------------------------------

## log prior of theta = log(tau) for a Gamma(shape = a, rate = b) precision
## ("loggamma" in INLA), including the d tau / d theta = tau Jacobian:
##   p(theta) = Gamma(exp(theta); a, b) * exp(theta)
#' @keywords internal
.loggamma_prec_logprior <- function(a, b) {
  function(lprec) a * log(b) - lgamma(a) + a * lprec - b * exp(lprec)
}

## log prior of theta = log(tau) for INLA's "normal"/"gaussian" prior, which is
## a Gaussian placed directly on the internal scale theta with
## param = (mean, precision). No Jacobian: it already is the density of theta.
#' @keywords internal
.normal_prec_logprior <- function(mean, prec) {
  function(lprec) stats::dnorm(lprec, mean = mean, sd = 1 / sqrt(prec), log = TRUE)
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
  if (fam %in% c("normal", "gaussian")) {
    if (!ok2 || param[2L] <= 0) return(list(lp = NULL, family = prior, mappable = FALSE))
    return(list(lp = .normal_prec_logprior(param[1L], param[2L]),
                family = "normal", mappable = TRUE))
  }
  list(lp = NULL, family = prior, mappable = FALSE)
}

## Read each selected component's hyperparameter record from
## fit$all.hyper$random, whose entries are parallel to names(fit$summary.random).
## Returns a named list (NULL entry when the component has no record -- which
## ngvb_apply_user_prec_prior treats as a loud fallback, not a silent one).
#' @keywords internal
ngvb_hyper_from_fit <- function(fit, comp.names) {
  cn.all <- names(fit$summary.random)
  ah     <- fit$all.hyper$random
  out    <- stats::setNames(vector("list", length(comp.names)), comp.names)
  if (is.null(ah)) return(out)
  for (cn in comp.names) {
    k <- match(cn, cn.all)
    if (!is.na(k) && k <= length(ah)) out[[cn]] <- ah[[k]]$hyper
  }
  out
}

## Locate the precision entry in a hyper list. Supports both INLA's all.hyper
## layout (entries carry name/short.name) and a plain user-style
## list(prec = list(...)). Returns the spec, the display labels of the FREE
## (non-fixed) secondary hyperparameters, and whether the input came from
## all.hyper (which determines message vs warning for secondaries).
#' @keywords internal
.find_prec_entry <- function(hyper) {
  nm <- names(hyper)
  if (is.null(nm)) nm <- rep("", length(hyper))
  short <- vapply(hyper, function(hh)
    if (is.list(hh) && !is.null(hh$short.name)) as.character(hh$short.name) else NA_character_,
    NA_character_)
  from.all.hyper <- any(!is.na(short))

  k <- which(short == "prec")
  if (!length(k)) {
    lname <- vapply(hyper, function(hh)
      if (is.list(hh) && !is.null(hh$name)) tolower(as.character(hh$name)) else "", "")
    k <- which(grepl("precision", lname, fixed = TRUE))
  }
  if (!length(k)) k <- which(nm %in% c("prec", "theta1", "theta"))
  if (!length(k))
    return(list(spec = NULL, secondary = nm, from.all.hyper = from.all.hyper))
  k <- k[1L]

  free <- vapply(hyper, function(hh) !isTRUE(hh$fixed), TRUE)
  lab  <- ifelse(is.na(short), nm, short)
  list(spec = hyper[[k]],
       secondary = lab[free & seq_along(hyper) != k],
       from.all.hyper = from.all.hyper)
}

## Rebuild op$logprior so its precision piece is the fit's prior, leaving any
## secondary-hyperparameter prior intact. Every fallback is loud: a warning for
## anything that cannot be carried, a message for what was carried.
#' @keywords internal
ngvb_apply_user_prec_prior <- function(op, hyper, comp.name, verbose = TRUE) {
  warn <- function(msg)
    if (isTRUE(verbose)) warning(msg, call. = FALSE, immediate. = TRUE)
    else warning(msg, call. = FALSE)

  if (is.null(hyper) || !length(hyper)) {
    warn(sprintf(paste0("ngvb: no hyperparameter record found for '%s' in the fit ",
                        "(fit$all.hyper); using the '%s' operator's default prior."),
                 comp.name, op$type))
    return(op)
  }

  found <- .find_prec_entry(hyper)
  spec  <- found$spec

  if (is.null(spec)) {
    ## In a user-style list, having no precision entry just means there is
    ## nothing to carry (e.g. only a correlation prior was set) -- the
    ## secondary handling below still reports it. In an all.hyper record a
    ## missing precision entry is anomalous and deserves a warning.
    if (isTRUE(found$from.all.hyper))
      warn(sprintf(paste0("ngvb: could not identify a precision hyperparameter for '%s'; ",
                          "using the '%s' operator's default prior."), comp.name, op$type))
  } else {
    map <- .prec_logprior_from_hyper(spec)
    if (is.null(op$prec.logprior)) {
      warn(sprintf(paste0("ngvb: a precision prior is set on '%s' in the fit, but the '%s' ",
                          "engine is not parameterized by a single precision, so it cannot be ",
                          "carried over; using the operator's default prior."),
                   comp.name, op$type))
    } else if (isTRUE(map$mappable)) {
      base <- op$logprior; old <- op$prec.logprior; new <- map$lp
      op$logprior      <- function(theta) base(theta) - old(theta[1L]) + new(theta[1L])
      op$prec.logprior <- new
      if (isTRUE(verbose))
        message(sprintf("ngvb: carried the %s precision prior from the INLA fit into component '%s'.",
                        map$family, comp.name))
    } else {
      reason <- if (identical(map$family, "fixed"))
        "a fixed precision is not supported by the ngvb engine"
      else sprintf("the '%s' prior family is not one ngvb can map (pc.prec, loggamma, normal)",
                   if (is.null(map)) "unknown" else map$family)
      warn(sprintf(paste0("ngvb: the precision prior on '%s' was dropped -- %s. Using the ",
                          "default PC prior; set it explicitly with components = list(%s = ",
                          "ngvb_operator(..., pc.prec = c(U = , alpha = )))."),
                   comp.name, reason, comp.name))
    }
  }

  if (length(found$secondary)) {
    txt <- paste(sprintf("'%s'", found$secondary), collapse = ", ")
    if (isTRUE(found$from.all.hyper)) {
      ## all.hyper always records every hyperparameter (usually at its INLA
      ## default), so this is routine, not a user error: inform, don't warn.
      if (isTRUE(verbose))
        message(sprintf(paste0("ngvb: component '%s' uses the engine's default prior(s) for %s ",
                               "(only the precision prior is carried over)."), comp.name, txt))
    } else {
      warn(sprintf(paste0("ngvb: prior(s) on %s for component '%s' were not carried over ",
                          "(ngvb maps only the precision prior); the engine default is used."),
                   txt, comp.name))
    }
  }
  op
}
