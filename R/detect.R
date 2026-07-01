## ---------------------------------------------------------------------------
## Auto-detection: map a fitted INLA component to an ngvb operator.
## Model TYPE comes from fit$model.random (parallel to names(fit$summary.random));
## the operator is then rebuilt from the fit where possible:
##   iid / rw1 / rw2 / ar1  -> analytic (dimension from summary.random)
##   besag (ICAR)           -> adjacency recovered from config$Qprior
##   SPDE2                  -> NOT reconstructable from the fit; user passes it.
## ---------------------------------------------------------------------------

.ngvb_model_map <- c(
  "IID model"          = "iid",
  "RW1 model"          = "rw1",
  "RW2 model"          = "rw2",
  "AR1 model"          = "ar1",
  "Besags ICAR model"  = "car_icar",
  "SPDE2 model"        = "spde"
)

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
    spde = stop("ngvb2: SPDE component '", comp.name,
                "' cannot be rebuilt from the fit (mesh not retained). Supply it via ",
                "components = list(", comp.name, " = ngvb_operator('spde', spde = <your spde>)).")
  )
}
