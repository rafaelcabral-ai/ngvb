# Fit a latent non-Gaussian model from a fitted INLA (LGM) object.

Fit a latent non-Gaussian model from a fitted INLA (LGM) object.

## Usage

``` r
ngvb(
  fit,
  selection = NULL,
  components = NULL,
  method = c("SVI", "SCVI"),
  alpha.eta = 2,
  identify.scale = TRUE,
  iter = 30,
  stop.rel.change = 0.001,
  n.sampling = 2000,
  verbose = TRUE
)
```

## Arguments

- fit:

  An \`inla\` object (fitted with \`control.compute = list(config =
  TRUE)\`).

- selection:

  Optional named list of component indices to extend to non-Gaussianity
  (default: all random effects).

- components:

  Optional named list of operator descriptors overriding auto-detection
  (required for SPDE: \`list(s = ngvb_operator("spde", spde = spde))\`).
  A component listed here is taken as the full prior specification
  (including its \`pc.prec\`), so its \`f()\`-term \`hyper\` is not
  consulted.

- method:

  Variational algorithm: \`"SVI"\` (structured, mean-field; the default
  – more reliable) or \`"SCVI"\` (structured & collapsed; reaches a
  fixed point in far fewer iterations, but its collapsed marginal can
  over-shrink a weak-but-genuine effect all the way to Gaussian, \`eta ~
  0\`, where SVI holds a small positive value). The two agree when the
  signal is strong or clearly absent and disagree in the weak-signal
  regime; prefer the default \`"SVI"\` unless you specifically need
  SCVI's speed and have checked the two give the same answer on your
  problem. Regardless of method, when a component's discrepancies are
  small but roughly homogeneous across all its indices (no real
  outlier/inlier split), the fit can drift toward a large, mostly
  prior-driven eta rather than genuine non-Gaussianity – see
  \`alpha.eta\` and the degeneracy warning this function may emit.

- alpha.eta:

  Exponential-PC-prior rate(s) on the non-Gaussianity parameter(s):
  larger values shrink harder toward the Gaussian model (\`eta = 0\`).
  If a fit triggers the degeneracy warning, raise this (e.g. by 2-5x)
  before trusting the result.

- identify.scale:

  Enforce the scale-identifiability constraint (default \`TRUE\`). The
  conditional precision is invariant under \\(\tau, \mathbf V) \mapsto
  (c\tau, c\mathbf V)\\, an exact flat ridge in the likelihood that the
  VB loop can otherwise walk – letting the mixing variables absorb the
  overall scale and inflating \`eta\` spuriously (most visibly for
  intrinsic models such as \`besag\`/\`rw\`). With \`identify.scale\`
  on, each component's discrepancies are re-anchored so their total
  matches the Gaussian total each iteration, pinning the scale to
  \`tau\` and leaving \`eta\` to respond only to \*relative\* departures
  (genuine outliers). Turn off only to reproduce the unconstrained
  updates.

- iter, stop.rel.change, n.sampling, verbose:

  VB controls.

## Value

An object of class \`ngvb\` with the final INLA \`fit\`, the mixing
vectors \`V\`, the non-Gaussianity parameters \`eta\`, and their
trajectory.

## Priors

The ngvb engine rebuilds each selected component's structure, so
\`model\`/\`graph\`/\`scale.model\` on the original \`f()\` term are
replaced. A precision prior set there (\`hyper = list(prec = ...)\`)
\*is\* honored when it is a \`pc.prec\` or \`loggamma\` prior — it is
translated onto the engine's precision hyperparameter and a message
reports it. Any other family, a fixed precision, or a prior on a
secondary hyperparameter (e.g. an AR1 correlation) is dropped with a
warning and the operator's default PC prior is used; to control it, pass
the operator via \`components\` with an explicit \`pc.prec = c(U = ,
alpha = )\` (the PC prior sets P(sigma \> U) = alpha, with sigma =
1/sqrt(precision)).

## See also

\[ng.check()\], \[ngvb_operator()\]

## Examples

``` r
# \donttest{
if (requireNamespace("INLA", quietly = TRUE)) {
  data(jumpts)   # a series with two abrupt jumps -- see the package vignette
  LGM  <- INLA::inla(y ~ -1 + f(x, model = "rw1"), data = jumpts,
                     control.compute = list(config = TRUE))
  LnGM <- ngvb(LGM, iter = 10)     # non-Gaussian extension, model auto-detected
  summary(LnGM)
}
#> Components: x [rw1] 
#> ngvb: carried the loggamma precision prior from the INLA fit into component 'x'.
#> ngvb: reached the iteration limit after 10 iteration(s);  E[eta] = 2.442
#> Latent non-Gaussian model (ngvb)
#> 
#> Fixed effects:
#> data frame with 0 columns and 0 rows
#> 
#> Non-Gaussianity -- eta: mean (median) [90% CI]:
#>   x            [rw1]: eta = 2.442 (2.408) [1.966, 3.035];
#>                  most non-Gaussian indices (E[V]/h): 40 (10.6), 20 (6.9), 16 (1.2), 14 (1.2), 90 (1.1)
# }
```
