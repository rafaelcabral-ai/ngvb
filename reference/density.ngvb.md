# Posterior density of a component's non-Gaussianity parameter eta.

For \`method = "SVI"\` this is a kernel density estimate over the GIG
variational factor's Monte-Carlo sample; for \`"SCVI"\` it is the
weighted quadrature grid underlying \`scvi_update()\`. Either way it is
\`q(eta)\`, the actual (approximate) posterior – not just its mean.

## Usage

``` r
# S3 method for class 'ngvb'
density(x, component = x$comp.names[1L], ...)
```

## Arguments

- x:

  An \`ngvb\` object.

- component:

  Component name (default: the first).

- ...:

  Passed to \[stats::density()\].

## Value

A \`density\` object (see \[stats::density()\]); \`plot()\` it directly.
