# Build an INLA rgeneric model for an ngvb operator at a fixed V.

Build an INLA rgeneric model for an ngvb operator at a fixed V.

## Usage

``` r
ngvb_rgeneric(op, V = op$h)
```

## Arguments

- op:

  An operator descriptor from \[ngvb_operator()\].

- V:

  Numeric vector of mixing variables (length = nrow of the latent
  field). The precision uses \`diag(1/V)\`; pass \`op\$h\` for the
  Gaussian model.

## Value

An object usable as \`f(idx, model = \<this\>)\` in an INLA formula.
