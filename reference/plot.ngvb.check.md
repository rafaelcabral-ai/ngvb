# Diagnostic plots for a latent-Gaussianity check.

For each checked component, draws (left) the per-index Bayes-factor
sensitivity \\d_i(y)\\ – spikes locate where the Gaussian assumption is
least adequate – and (right, for a Gaussian response) the observed
overall sensitivity \\s_0=\sum_i d_i\\ against its Gaussian reference
distribution; an observed value far in the tail (small p-value) signals
latent non-Gaussianity.

## Usage

``` r
# S3 method for class 'ngvb.check'
plot(x, ...)
```

## Arguments

- x:

  An \`ngvb.check\` object from \[ng.check()\].

- ...:

  Ignored.

## Value

A \`patchwork\` object combining the per-component diagnostic panels;
called mainly for the plot it draws.
