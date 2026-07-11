# Recover the inla.spde2 object referenced by a fitted SPDE component. The object is not stored in the fit, but the formula still references it in its environment, so we evaluate the \`model =\` symbol there.

Recover the inla.spde2 object referenced by a fitted SPDE component. The
object is not stored in the fit, but the formula still references it in
its environment, so we evaluate the \`model =\` symbol there.

## Usage

``` r
ngvb_recover_spde(fit, comp.name)
```
