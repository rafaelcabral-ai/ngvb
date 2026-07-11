# Diagnostic plots for an ngvb fit.

For each component two panels are drawn: the mixing weights V_i/h_i
(posterior mean with a 90 percent credible interval, and the plug-in
value actually refit with INLA), and the posterior density of eta (with
its mean and median marked). A final panel shows how eta moved across
the variational iterations. The plots are ggplot2 objects assembled with
patchwork.

## Usage

``` r
# S3 method for class 'ngvb'
plot(x, ...)
```

## Arguments

- x:

  An \`ngvb\` object.

- ...:

  Ignored.

## Value

(Invisibly) the \`ngvb\` object; called for the plot it draws.
