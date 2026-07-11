# Bayes factor of a latent non-Gaussian model against its Gaussian counterpart.

Returns the marginal-likelihood Bayes factor
\\Z\_{\mathrm{LnGM}}/Z\_{\mathrm{LGM}}\\ with the mixing variables \`V\`
integrated out by importance sampling from \`q(V)\`. A value above one
favours the non-Gaussian model; the usual Jeffreys reading is that
\\\log\_{10}\\ above 0.5, 1, 2 is substantial, strong, decisive.

## Usage

``` r
bayes.factor(x, ...)

# S3 method for class 'ngvb'
bayes.factor(x, ...)

# S3 method for class 'ngvb.samples'
bayes.factor(x, ...)
```

## Arguments

- x:

  An \`ngvb\` object (sampled internally) or an \`ngvb.samples\` object.

- ...:

  For the \`ngvb\` method, passed to \[ngvb_sample()\] (e.g.
  \`n.samples\`, \`seed\`).

## Value

A list with \`BF\`, \`log10BF\`, \`logBF\`, and the effective sample
size \`ess\` of the importance weights (a low \`ess\` means more draws
are needed).
