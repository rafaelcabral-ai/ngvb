# Package index

## Diagnose and fit

The two main functions. Both take an ordinary fitted `inla` object.

- [`ng.check()`](https://rafaelcabral-ai.github.io/ngvb/reference/ng.check.md)
  : Check the latent Gaussian assumption of an INLA fit.
- [`ngvb()`](https://rafaelcabral-ai.github.io/ngvb/reference/ngvb.md) :
  Fit a latent non-Gaussian model from a fitted INLA (LGM) object.

## Model comparison

Integrate the mixing variables out and compare against the Gaussian
model.

- [`ngvb_sample()`](https://rafaelcabral-ai.github.io/ngvb/reference/ngvb_sample.md)
  : Sample the mixing variables and refit INLA at each draw.
- [`bayes.factor()`](https://rafaelcabral-ai.github.io/ngvb/reference/bayes.factor.md)
  : Bayes factor of a latent non-Gaussian model against its Gaussian
  counterpart.

## Operators and custom models

The dependency-matrix building blocks, and how to add your own.

- [`ngvb_operator()`](https://rafaelcabral-ai.github.io/ngvb/reference/ngvb_operator.md)
  : Construct an ngvb operator descriptor.
- [`ngvb_custom()`](https://rafaelcabral-ai.github.io/ngvb/reference/ngvb_custom.md)
  : Define a custom latent non-Gaussian model from a precision "square
  root".
- [`ngvb_rgeneric()`](https://rafaelcabral-ai.github.io/ngvb/reference/ngvb_rgeneric.md)
  : Build an INLA rgeneric model for an ngvb operator at a fixed V.
- [`ngvb_precision()`](https://rafaelcabral-ai.github.io/ngvb/reference/ngvb_precision.md)
  : Conditional precision matrix Q(theta, V) for an operator
  (main-process helper).

## Plotting

Diagnostic plots and maps for areal and geostatistical data.

- [`plot(`*`<ngvb>`*`)`](https://rafaelcabral-ai.github.io/ngvb/reference/plot.ngvb.md)
  : Diagnostic plots for an ngvb fit.
- [`plot(`*`<ngvb.check>`*`)`](https://rafaelcabral-ai.github.io/ngvb/reference/plot.ngvb.check.md)
  : Diagnostic plots for a latent-Gaussianity check.
- [`areal.plot()`](https://rafaelcabral-ai.github.io/ngvb/reference/areal.plot.md)
  : Plot areal (regional) data on a map.
- [`geo.plot()`](https://rafaelcabral-ai.github.io/ngvb/reference/geo.plot.md)
  : Plot geostatistical (point-referenced) data, optionally over a mesh.

## Methods

S3 methods for the fitted objects.

- [`summary(`*`<ngvb>`*`)`](https://rafaelcabral-ai.github.io/ngvb/reference/summary.ngvb.md)
  : Summary of an ngvb fit: the underlying INLA summaries plus, per
  component, the non-Gaussianity parameter and the most strongly
  down-weighted (outlying) indices.
- [`density(`*`<ngvb>`*`)`](https://rafaelcabral-ai.github.io/ngvb/reference/density.ngvb.md)
  : Posterior density of a component's non-Gaussianity parameter eta.
- [`fitted(`*`<ngvb>`*`)`](https://rafaelcabral-ai.github.io/ngvb/reference/fitted.ngvb.md)
  : Fitted values of the underlying INLA fit.
- [`summary(`*`<ngvb.samples>`*`)`](https://rafaelcabral-ai.github.io/ngvb/reference/summary.ngvb.samples.md)
  : Importance-weighted posterior summary of the sampled fits.

## Data

- [`jumpts`](https://rafaelcabral-ai.github.io/ngvb/reference/jumpts.md)
  : Simulated time series with two jumps
- [`Orthodont`](https://rafaelcabral-ai.github.io/ngvb/reference/Orthodont.md)
  : Orthodontic growth data
