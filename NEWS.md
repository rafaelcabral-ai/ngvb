# ngvb2 0.1.0

First release: a rebuild of `ngvb` on top of R-INLA. INLA is a `Suggests`
dependency (installed from its own repository) so the package installs and
checks cleanly on CRAN.

## Engine
* A single unified `rgeneric` engine implements the conditional precision
  `Q(theta, V) = D(theta)^T diag(1/V) D(theta)` for every model, and **owns the
  V-dependent normalizing constant** — fixing the bias of INLA's `generic0`,
  which drops `0.5*log|Cmatrix|` (and hence the per-iteration `-0.5*sum(log V)`
  term). Operators may supply an analytic `lognc` (SAR/CAR use the exact
  eigenvalue log-determinant, stable near the singular boundary).

## Models
* Operators: `iid`, `rw1`, `rw2`, `ar1`, `sar`, `car` (proper + intrinsic ICAR),
  `spde`/Matern, `ou` (Ornstein–Uhlenbeck, irregular times), `seasonal`, and
  `generic0` (any user-supplied structure matrix `Cmatrix`, factored via its
  eigendecomposition; handles proper and intrinsic C). Each is validated against
  the analytic precision or native INLA. `seasonal` is auto-detected from a
  fitted `inla` object; `generic0` is a deliberate manual operator
  (`ngvb_operator("generic0", C = )`, passed via `components =`) because its
  `Cmatrix` must be positive semi-definite and an intrinsic `C` needs its own
  null-space constraints.
* `ngvb_sample()` draws `V` from its variational posterior `q(V)` and refits INLA
  at each draw (`ngvb.samples`); `bayes.factor()` returns the marginal-likelihood
  Bayes factor of the LnGM vs the Gaussian LGM with `V` integrated out (by
  importance sampling), plus the weight ESS. `summary.ngvb.samples()` pools any
  INLA summary (fixed effects, hyperparameters) over the draws.

## User interface
* `ngvb(fit)` extends a fitted latent Gaussian model to non-Gaussian: it
  auto-detects each component (`fit$model.random`), rewrites the formula, and
  runs a multi-component (block-diagonal) variational-Bayes loop — SVI (the
  more reliable structured mean-field) by default; `method = "SCVI"` (collapsed)
  is faster but can over-shrink a weak-but-genuine effect to `eta ~ 0`, so it is
  opt-in. SPDE components are supplied via
  `components = list(s = ngvb_operator("spde", spde = spde))`.
* **The precision prior set on an `f()` term is carried into the fit.** When the
  original `f(<name>, ..., hyper = list(prec = ...))` uses a `pc.prec` or
  `loggamma` prior, `ngvb()` translates it onto the engine's precision
  hyperparameter (so your prior persists in `LnGM$fit`) and reports that it did.
  Any other prior family, a fixed precision, or a prior on a secondary
  hyperparameter (e.g. an AR1 correlation) cannot be mapped: it is dropped with a
  warning and the operator's default PC prior is used — set it explicitly with
  `components = list(<name> = ngvb_operator(..., pc.prec = c(U = , alpha = )))`.
  (Structure arguments such as `graph`/`scale.model` are always rebuilt by the
  engine and do not carry over.)
* `ng.check(fit)` computes the Bayes-factor sensitivity diagnostic of Cabral,
  Bolin & Rue (JRSS-B 2025); reproduces the reference implementation exactly and
  additionally handles `inla.stack` (SPDE) fits.
* `print`/`summary`/`plot`/`fitted`/`density` methods for the fitted object;
  `eta` is reported with a Monte-Carlo/quadrature 90% interval (`$eta.q`) and
  each `V_i` with a mean and 90% CI (`$V.summary`), not just point estimates.
  `plot.ngvb()` and `plot.ngvb.check()` are ggplot2/patchwork figures.
* `ng.check()` reports `s0` and its reference at the hyperparameter posterior
  mode (matching the reference implementation's `s0.mode`/`var.ref.mode`); the
  hyperparameter-mixture averages are kept as `s0.mixture`/`d.mixture`.
* `areal.plot()` and `geo.plot()` map areal and geostatistical data with
  ggplot2/sf (static, reproducible replacements for the old leaflet helpers).
* **Scale identifiability (`identify.scale = TRUE`, default).** The conditional
  precision is invariant under `(tau, V) -> (c*tau, c*V)`, an exact flat ridge
  the VB loop can walk on weak/homogeneous signal — letting `V` absorb the scale
  and inflating `eta` spuriously (intrinsic `besag`/`rw` are the textbook case).
  Each iteration the constraint re-anchors a component's discrepancies so their
  total matches the Gaussian total, pinning the scale to `tau` and leaving `eta`
  to respond only to *relative* departures (genuine outliers). Removes the
  runaway (e.g. besag `eta` 3.4 -> 0.1, matching a field with ~0 non-Gaussian
  edges) and leaves well-scaled fits essentially unchanged. Turn off to recover
  the unconstrained updates.
* `ngvb()` warns when a component's mixing weights shrink nearly uniformly
  across (almost) every index *and none is inflated* — the residual signature of
  `eta` and the component's own precision becoming weakly identified (a genuine
  heavy-tailed fit also shrinks most indices, but pays with a few strongly
  inflated ones, so it is not flagged). Raise `alpha.eta` (default 2) in response.

## Not yet implemented
* A `cgeneric` fast engine (the `rgeneric` engine covers all models meanwhile).
* AR(p) for p > 1, separable space-time, BYM2.
