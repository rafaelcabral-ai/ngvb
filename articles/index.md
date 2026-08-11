# Articles

### Getting started

- [Get started with
  ngvb](https://rafaelcabral-ai.github.io/ngvb/articles/ngvb.md):

  What a latent non-Gaussian model is, which INLA components are
  supported, and a first end-to-end example on a time series with two
  jumps.

### Worked examples

One data type per article. Each is self-contained and runs against a
real INLA installation.

- [Geostatistical data: a Matern (SPDE)
  field](https://rafaelcabral-ai.github.io/ngvb/articles/geostatistical.md):

  A Matern field fitted to Pacific Northwest pressure readings, with a
  map of where the spatial predictions will move before the non-Gaussian
  model is ever fitted.

- [Areal data: an intrinsic model
  (besag)](https://rafaelcabral-ai.github.io/ngvb/articles/areal.md):

  An intrinsic conditional autoregression for burglary rates across
  Columbus, Ohio, and what the diagnostic says when a field really is
  Gaussian.

- [Random intercepts and
  slopes](https://rafaelcabral-ai.github.io/ngvb/articles/longitudinal.md):

  Subject-level random effects in a longitudinal growth study, and what
  changes when the random effects are allowed to be non-Gaussian.

- [Custom models: your own precision
  matrix](https://rafaelcabral-ai.github.io/ngvb/articles/custom-models.md):

  Building a simultaneous autoregression by hand with ngvb_custom(),
  argument by argument, and borrowing ready-made dependency matrices
  from the ngme2 package.
