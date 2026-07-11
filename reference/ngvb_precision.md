# Conditional precision matrix Q(theta, V) for an operator (main-process helper).

Used for validation and for computing the discrepancy statistic d_i.

## Usage

``` r
ngvb_precision(op, theta = op$theta.initial, V = op$h)
```

## Arguments

- op:

  Operator descriptor.

- theta:

  Hyperparameters (internal scale). Defaults to \`op\$theta.initial\`.

- V:

  Mixing vector. Defaults to \`op\$h\` (Gaussian model).

## Value

A sparse precision matrix \`Q(theta, V) = D(theta)^T diag(1/V)
D(theta)\`.
