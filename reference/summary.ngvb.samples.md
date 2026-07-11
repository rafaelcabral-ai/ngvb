# Importance-weighted posterior summary of the sampled fits.

Pools a chosen INLA summary across the sampled fits with the importance
weights, giving posterior means and standard deviations for the latent
non-Gaussian model (with \`V\` integrated out), plus the Bayes factor.

## Usage

``` r
# S3 method for class 'ngvb.samples'
summary(object, what = c("all", "fixed", "hyperpar"), ...)
```

## Arguments

- object:

  An \`ngvb.samples\` object.

- what:

  One of \`"all"\` (default), \`"fixed"\`, or \`"hyperpar"\` – which
  INLA summary table(s) to pool. \`"all"\` shows and returns both,
  silently skipping either one that doesn't apply to this model (e.g.
  \`"fixed"\` on a fixed-effect-free model). Run
  \`args(summary.ngvb.samples)\` or \`?summary.ngvb.samples\` to see
  this list again.

- ...:

  Ignored.

## Value

Called for the summary it prints. Invisibly: a single data frame of
importance-weighted posterior means/sds for \`what = "fixed"\` or
\`"hyperpar"\`; for \`"all"\`, a list \`list(fixed = , hyperpar = )\`
with either element \`NULL\` if that table doesn't apply.
