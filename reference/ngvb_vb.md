# Multi-component structured variational inference loop.

Multi-component structured variational inference loop.

## Usage

``` r
ngvb_vb(
  inla.fit.V,
  ops,
  comp.names,
  method = c("SVI", "SCVI"),
  alpha.eta = 2,
  identify.scale = TRUE,
  iter = 10,
  stop.rel.change = 0.001,
  n.sampling = 2000,
  verbose = TRUE
)
```
