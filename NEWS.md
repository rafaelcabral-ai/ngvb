# ngvb2 0.0.0.9000

Initial rebuild of `ngvb` on top of R-INLA.

## Engine
* A single unified `rgeneric` engine implements the conditional precision
  `Q(theta, V) = D(theta)^T diag(1/V) D(theta)` for every model, and **owns the
  V-dependent normalizing constant** — fixing the bias of INLA's `generic0`,
  which drops `0.5*log|Cmatrix|` (and hence the per-iteration `-0.5*sum(log V)`
  term). Operators may supply an analytic `lognc` (SAR/CAR use the exact
  eigenvalue log-determinant, stable near the singular boundary).

## Models
* Operators: `iid`, `rw1`, `rw2`, `ar1`, `sar`, `car` (proper + intrinsic ICAR),
  `spde`/Matern, and `ou` (Ornstein–Uhlenbeck, irregular times). Each is
  validated against the analytic precision or native INLA.

## User interface
* `ngvb(fit)` extends a fitted latent Gaussian model to non-Gaussian: it
  auto-detects each component (`fit$model.random`), rewrites the formula, and
  runs a multi-component (block-diagonal) variational-Bayes loop — SCVI by
  default, `method = "SVI"` available. SPDE components are supplied via
  `components = list(s = ngvb_operator("spde", spde = spde))`.
* `ng.check(fit)` computes the Bayes-factor sensitivity diagnostic of Cabral,
  Bolin & Rue (JRSS-B 2025); reproduces the reference implementation exactly and
  additionally handles `inla.stack` (SPDE) fits.
* `print`/`summary`/`plot`/`fitted` methods for the fitted object.

## Not yet implemented
* A `cgeneric` fast engine (the `rgeneric` engine covers all models meanwhile).
* AR(p) for p > 1, separable space-time, BYM2.
