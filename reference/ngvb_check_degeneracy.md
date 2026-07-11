# Warn when a component's mixing weights have shrunk nearly uniformly across (almost) every index, rather than concentrating on a few outliers.

SVI and SCVI target the same fixed point and, when a component's
per-index discrepancies d_i are small but roughly HOMOGENEOUS (no real
outlier/inlier split), that shared fixed point can drift to a large,
mostly prior- and N-driven eta rather than one reflecting genuine local
non-Gaussianity: tau and eta become weakly identified against each other
and reinforce each other iteration over iteration. A uniform V/h \<\< 1
across almost the whole component is the signature of that regime, not
of real outlier detection.

## Usage

``` r
ngvb_check_degeneracy(
  V,
  h,
  comp.names,
  alpha.eta,
  verbose,
  frac.threshold = 0.8,
  ratio.threshold = 0.5,
  inflate.threshold = 2
)
```
