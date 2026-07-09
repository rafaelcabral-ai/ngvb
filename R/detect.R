## ---------------------------------------------------------------------------
## Auto-detection: map a fitted INLA component to an ngvb operator.
## Model TYPE comes from fit$model.random (parallel to names(fit$summary.random));
## the operator is then rebuilt from the fit where possible:
##   iid / rw1 / rw2 / ar1  -> analytic (dimension from summary.random)
##   besag (ICAR)           -> adjacency recovered from config$Qprior
##   SPDE2                  -> NOT reconstructable from the fit; user passes it.
## ---------------------------------------------------------------------------

.ngvb_model_map <- c(
  "IID model"           = "iid",
  "RW1 model"           = "rw1",
  "RW2 model"           = "rw2",
  "AR1 model"           = "ar1",
  "Besags ICAR model"   = "car_icar",
  "SPDE2 model"         = "spde",
  "Seasonal model"      = "seasonal"
  ## Note: "Generic0 model" is intentionally NOT auto-detected. generic0 is
  ## available as a deliberate manual operator, ngvb_operator("generic0", C = ),
  ## passed via components = -- see its caveats (the Cmatrix must be positive
  ## semi-definite, and an intrinsic C needs its own null-space constraints).
)

#' Find the `model =` expression of the f(<comp>, ...) term in a formula.
#' @keywords internal
ngvb_find_f_model <- function(formula, comp.name) {
  found <- NULL
  walk <- function(e) {
    if (is.call(e)) {
      if (identical(e[[1L]], as.name("f"))) {
        idx <- tryCatch(as.character(e[[2L]]), error = function(...) "")
        al  <- as.list(e)
        if (length(idx) == 1L && idx == comp.name && "model" %in% names(al))
          found <<- al[["model"]]
      }
      for (i in seq_along(e)) walk(e[[i]])
    }
  }
  walk(formula[[length(formula)]])
  found
}

#' Find and evaluate a named argument of the f(<comp>, ...) term in a formula.
#' Used to recover, e.g., `season.length` or `Cmatrix` from the fitted call.
#' @keywords internal
ngvb_find_f_arg <- function(fit, comp.name, argname) {
  frm <- fit$.args$formula
  found <- NULL
  walk <- function(e) {
    if (is.call(e)) {
      if (identical(e[[1L]], as.name("f"))) {
        al <- as.list(e)
        idx <- tryCatch(as.character(al[[2L]]), error = function(...) "")
        if (length(idx) == 1L && idx == comp.name && argname %in% names(al))
          found <<- al[[argname]]
      }
      for (i in seq_along(e)) walk(e[[i]])
    }
  }
  walk(frm[[length(frm)]])
  if (is.null(found)) return(NULL)
  tryCatch(eval(found, environment(frm)), error = function(e) NULL)
}

#' Recover the inla.spde2 object referenced by a fitted SPDE component.
#' The object is not stored in the fit, but the formula still references it in
#' its environment, so we evaluate the `model =` symbol there.
#' @keywords internal
ngvb_recover_spde <- function(fit, comp.name) {
  frm   <- fit$.args$formula
  mexpr <- ngvb_find_f_model(frm, comp.name)
  if (is.null(mexpr))
    stop("ngvb2: could not locate the model of SPDE component '", comp.name, "' in the formula.")
  obj <- tryCatch(eval(mexpr, environment(frm)), error = function(e) NULL)
  if (!inherits(obj, "inla.spde2"))
    stop("ngvb2: could not recover the inla.spde2 object for '", comp.name,
         "' (it is not retained in the fit). Supply it via components = list(",
         comp.name, " = ngvb_operator('spde', spde = <your spde>)).")
  obj
}

#' Recover a binary adjacency from a fitted besag/ICAR component's prior precision.
#' config$Qprior block = tau * (diag(nnbs) - A); off-diagonal non-zeros mark neighbours.
#' @keywords internal
ngvb_recover_adjacency <- function(fit, comp.name) {
  cf  <- fit$misc$configs$config[[1L]]
  sel <- ngvb_component_index(fit, comp.name)
  Qp  <- as.matrix(cf$Qprior[sel, sel, drop = FALSE])
  Qp  <- Qp + t(Qp) - diag(diag(Qp))                 # symmetrize (stored upper-tri)
  nb  <- (abs(Qp) > 1e-8) & (row(Qp) != col(Qp))
  W   <- matrix(0, nrow(Qp), ncol(Qp)); W[nb] <- 1
  methods::as(Matrix::Matrix(W, sparse = TRUE), "CsparseMatrix")
}

#' Detect (or accept) the ngvb operator for one fitted component.
#' @keywords internal
ngvb_detect_operator <- function(fit, comp.name, user.op = NULL) {
  if (!is.null(user.op)) return(user.op)
  cn.all <- names(fit$summary.random)
  k  <- match(comp.name, cn.all)
  if (is.na(k)) stop("ngvb2: component '", comp.name, "' not found among random effects.")
  mr  <- fit$model.random[k]
  typ <- unname(.ngvb_model_map[mr])
  if (is.na(typ))
    stop("ngvb2: cannot auto-detect operator for '", comp.name,
         "' (model.random = '", mr, "'). Supply it via components = list(",
         comp.name, " = ngvb_operator(...)).")
  n <- nrow(fit$summary.random[[comp.name]])
  switch(typ,
    iid = ngvb_operator("iid", n = n),
    rw1 = ngvb_operator("rw1", n = n),
    rw2 = ngvb_operator("rw2", n = n),
    ar1 = ngvb_operator("ar1", n = n),
    car_icar = ngvb_operator("car", W = ngvb_recover_adjacency(fit, comp.name), intrinsic = TRUE),
    spde = ngvb_operator("spde", spde = ngvb_recover_spde(fit, comp.name)),
    seasonal = {
      s <- ngvb_find_f_arg(fit, comp.name, "season.length")
      if (is.null(s)) stop("ngvb2: could not recover season.length for '", comp.name,
                           "'. Supply components = list(", comp.name,
                           " = ngvb_operator('seasonal', n = ", n, ", season = <s>)).")
      ngvb_operator("seasonal", n = n, season = as.integer(s))
    }
  )
}
