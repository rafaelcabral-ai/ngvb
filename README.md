# ngvb2

<!-- badges: start -->
[![R-CMD-check](https://github.com/rafaelcabral-ai/ngvb2/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/rafaelcabral-ai/ngvb2/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

**Latent non-Gaussian modelling with R-INLA and variational Bayes.** A rebuild of
[`ngvb`](https://github.com/rafaelcabral96/ngvb).

`ngvb2` takes an ordinary R-INLA latent *Gaussian* model, **checks** whether the Gaussian
assumption is adequate, and **extends** it to a latent *non-Gaussian* model, each in one line.

```r
LGM <- inla(y ~ f(s, model = "rw1"), data = d, control.compute = list(config = TRUE))

ng.check(LGM)       # is the latent Gaussian assumption adequate, and where not?
LnGM <- ngvb(LGM)   # fit the non-Gaussian extension
bayes.factor(LnGM)  # how much better is it, with V integrated out?
```

## Installation

`ngvb2` builds on [INLA](https://www.r-inla.org/), which is not on CRAN, so install it first:

```r
install.packages("INLA",
  repos = c(getOption("repos"), INLA = "https://inla.r-inla-download.org/R/stable"),
  dep = TRUE)

# install.packages("remotes")
remotes::install_github("rafaelcabral-ai/ngvb2")
```

## The idea

A latent Gaussian field has precision `Q = Dᵀ D` for a model-specific *dependency matrix* `D`
(the difference operator for a random walk, `I − ρW` for a SAR, `κ²C + G` for a Matérn/SPDE, and
so on). The non-Gaussian extension replaces the Gaussian driving noise, which is equivalent to
conditioning on **mixing variables** `V`:

```
Q(θ, V) = D(θ)ᵀ diag(1/V) D(θ)
```

Where `Vᵢ = hᵢ` the model is Gaussian; where `Vᵢ` is inflated the increment is allowed to be
unusually large. A single non-Gaussianity parameter `η ≥ 0` controls how far `V` may stray, and
its exponential prior shrinks back to the Gaussian model unless the data pull away. `ngvb2` fits
this with a variational-Bayes loop that alternates an INLA fit (for fixed `V`) with closed-form
updates of `V` and `η`. One unified engine implements `Q(θ, V)` for every model and, unlike
INLA's `generic0`, correctly owns the `V`-dependent normalizing constant.

## Models

`ngvb(fit)` and `ng.check(fit)` auto-detect the component type from the fitted object and
rebuild the operator. No arguments are needed for the auto-detected models:

| Model | `f(..., model = )` | How |
|---|---|---|
| i.i.d. random effects | `"iid"` | auto |
| Random walks | `"rw1"`, `"rw2"` | auto |
| Autoregressive order 1 | `"ar1"` | auto |
| Intrinsic CAR (areal) | `"besag"` | auto (graph recovered) |
| Matérn / SPDE | `inla.spde2.matern(mesh)` | auto (read from the fit) |
| Seasonal | `"seasonal"` | auto |
| Structure matrix | `"generic0"` | `ngvb_operator("generic0", C = )` |
| SAR, proper CAR, OU, AR(p) | — | `ngvb_operator()` / `ngvb_custom()` |
| Your own precision | — | `ngvb_custom(D, h, ...)` |

Additive models just work: each `f()` term gets its own mixing variables and non-Gaussianity
parameter, all fit jointly.

A `pc.prec` or `loggamma` **precision prior** you set on a component (`f(s, …, hyper = list(prec = …))`)
is carried into the non-Gaussian fit; other prior families are dropped with a warning, in which case
set the prior explicitly via `components = list(s = ngvb_operator(…, pc.prec = c(U = , alpha = )))`.

## Documentation

See the [Get started](https://rafaelcabral-ai.github.io/ngvb2/articles/ngvb2.html) article for
worked examples across time series, areal, geostatistical and custom models, and the
[reference](https://rafaelcabral-ai.github.io/ngvb2/reference/index.html) for the full API.

## References

- Cabral, Bolin & Rue (2024). *Fitting latent non-Gaussian models using variational Bayes and
  Laplace approximations.* (the `ngvb` method)
- Cabral, Bolin & Rue (2025). *Robustness, model checking, and hierarchical models.* JRSS-B
  87(3):632–652. (the `ng.check` diagnostic)
