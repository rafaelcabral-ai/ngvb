# Importance-weighted posterior summary of the sampled fits.

Pools a chosen INLA summary across the sampled fits with the importance
weights, giving posterior means and standard deviations for the latent
non-Gaussian model (with \`V\` integrated out), plus the Bayes factor.

## Usage

``` r
# S3 method for class 'ngvb.samples'
summary(
  object,
  what = c("all", "fixed", "hyperpar", "random"),
  verbose = TRUE,
  ...
)
```

## Arguments

- object:

  An \`ngvb.samples\` object.

- what:

  One of \`"all"\` (default), \`"fixed"\`, \`"hyperpar"\`, or
  \`"random"\` – which INLA summary table(s) to pool. \`"random"\` pools
  \`summary.random\` for every component (the latent field itself: one
  table per component, keyed by node \`ID\`, e.g. mesh node or time
  index). \`"all"\` shows and returns all three, silently skipping any
  that don't apply to this model (e.g. \`"fixed"\` on a
  fixed-effect-free model); random-effect tables are printed as a head
  (all rows are still returned). Run \`args(summary.ngvb.samples)\` or
  \`?summary.ngvb.samples\` to see this list again.

- verbose:

  If \`TRUE\` (default), print the report as a side effect (as
  \`ng.check()\` does for its plot). Set \`FALSE\` to only get the
  return value back, e.g. when pooling \`what = "random"\` into a plot
  without the table being echoed.

- ...:

  Ignored.

## Value

Invisibly: a single data frame of importance-weighted posterior
means/sds for \`what = "fixed"\` or \`"hyperpar"\`; a named list of one
such data frame per component (keyed by node \`ID\`) for \`what =
"random"\`; for \`"all"\`, a list \`list(fixed = , hyperpar = , random =
)\` with any inapplicable element \`NULL\`. Printed as a side effect
when \`verbose = TRUE\`.
