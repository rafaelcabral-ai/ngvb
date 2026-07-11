# Check the latent Gaussian assumption of an INLA fit.

Check the latent Gaussian assumption of an INLA fit.

## Usage

``` r
ng.check(
  fit,
  selection = NULL,
  components = NULL,
  compute.fixed = TRUE,
  plot = TRUE
)
```

## Arguments

- fit:

  An \`inla\` object fitted with \`control.compute = list(config =
  TRUE)\`.

- selection:

  Optional named list of components to check (default: all random).

- components:

  Optional operator overrides (e.g. SPDE), as in \[ngvb()\].

- compute.fixed:

  If \`TRUE\`, also return the sensitivity of each fixed effect to each
  component's non-Gaussianity parameter.

- plot:

  If \`TRUE\` (default), draw the diagnostic plots (see
  \[plot.ngvb.check()\]): per-index Bayes-factor sensitivity, and the
  observed overall sensitivity against its Gaussian reference
  distribution.

## Value

An object of class \`ngvb.check\`. Per component, evaluated at the
hyperparameter posterior mode \\\hat\gamma\\: the BF sensitivity \`s0\`,
the per-index contributions \`d\`, and (Gaussian response) the reference
SD \`sd.ref\` and \`p.value\`; the hyperparameter-mixture averages are
also kept as \`s0.mixture\` / \`d.mixture\`. Plus \`sens.fixed\` if
requested.

## See also

\[ngvb()\]

## Examples

``` r
# \donttest{
if (requireNamespace("INLA", quietly = TRUE)) {
  set.seed(1); n <- 100
  x <- cumsum(rnorm(n, sd = 0.3)); x[50:n] <- x[50:n] + 6
  y <- x + rnorm(n, sd = 0.4)
  LGM <- INLA::inla(y ~ -1 + f(i, model = "rw1", constr = TRUE),
                    data = data.frame(y = y, i = 1:n),
                    control.compute = list(config = TRUE))
  ng.check(LGM)      # small p-value flags departure from latent Gaussianity
}

# }
```
