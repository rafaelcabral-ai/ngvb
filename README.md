# ngvb

<!-- badges: start -->
[![R-CMD-check](https://github.com/rafaelcabral-ai/ngvb/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/rafaelcabral-ai/ngvb/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/rafaelcabral-ai/ngvb/actions/workflows/pkgdown.yaml/badge.svg)](https://rafaelcabral-ai.github.io/ngvb/)
[![License: GPL v3](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![JASA 2024](https://img.shields.io/badge/JASA%202024-VB%20%2B%20Laplace-b31b1b.svg)](https://doi.org/10.1080/01621459.2023.2296704)
[![JRSS-B 2025](https://img.shields.io/badge/JRSS--B%202025-Model%20checking-b31b1b.svg)](https://doi.org/10.1093/jrsssb/qkae107)
<!-- badges: end -->

**ngvb** (non-Gaussian variational Bayes) takes an ordinary R-INLA latent *Gaussian* model (LGM), **checks** whether the latent Gaussian
assumption is adequate, and **extends** it to a latent *non-Gaussian* model (LnGM), each in one line.

```r
LGM <- inla(y ~ f(s, model = "rw1"), data = d, control.compute = list(config = TRUE))

ng.check(LGM)       # is the latent Gaussian assumption adequate, and where not?
LnGM <- ngvb(LGM)   # fit the non-Gaussian extension

samples <- ngvb_sample(LnGM, n.samples = 30) # sample the non-gaussian mixing variables V from their variational posterior and fit several R-INLA models
bayes.factor(samples)               # how much better the LnGM is, with V integrated out?
summary(samples)                    # importance-weighted summairies of fixed effects + hyperparameters of the LnGM
```

## Installation

`ngvb` builds on [INLA](https://www.r-inla.org/), which is not on CRAN, so install it first:

```r
install.packages("INLA",
  repos = c(getOption("repos"), INLA = "https://inla.r-inla-download.org/R/stable"),
  dep = TRUE)

# install.packages("remotes")
remotes::install_github("rafaelcabral-ai/ngvb")
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
its exponential prior shrinks back to the Gaussian model unless the data pull away. `ngvb` fits
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
| Matérn / SPDE | `inla.spde2.pcmatern(mesh, prior.range =, prior.sigma =)` | auto (read from the fit) |
| Seasonal | `"seasonal"` | auto |
| Structure matrix (any PSD `Cmatrix`) | `"generic0"` | `ngvb_operator("generic0", C = )` |
| Any CAR-type precision `Q` (non-positive off-diagonals) | — | `ngvb_operator("from_Q", Q = )` |
| SAR, proper CAR, OU | — | `ngvb_operator()` / `ngvb_custom()` |
| Your own precision, or borrowed from [ngme2](https://davidbolin.github.io/ngme2/) | — | `ngvb_custom(D, h, ...)` |

Additive models just work: each `f()` term gets its own mixing variables and non-Gaussianity
parameter, all fit jointly. `from_Q` covers the whole conditional-autoregression class (i.i.d.,
random walk, ICAR, proper CAR) with one canonical dependency-matrix factorization, so it's the
right fallback for a CAR-type component that isn't in the table above; models with positive
off-diagonal precision entries (RW2, AR(*p* > 1), SPDE) need their dedicated operator instead.

A `pc.prec`, `loggamma`, or `normal` **precision prior** you set on a component (`f(s, …, hyper = list(prec = …))`)
is carried into the non-Gaussian fit — including INLA's own default prior if you didn't set one —
so the LnGM is the exact non-Gaussian extension of the LGM you fitted; other prior families are
dropped with a warning, in which case set the prior explicitly via
`components = list(s = ngvb_operator(…, pc.prec = c(U = , alpha = )))`.

## Documentation

See the [Get started](https://rafaelcabral-ai.github.io/ngvb/articles/ngvb.html) article for
worked examples across time series, areal, geostatistical and custom models, and the
[reference](https://rafaelcabral-ai.github.io/ngvb/reference/index.html) for the full API.

## References

- Cabral, Bolin & Rue (2024). [*Fitting latent non-Gaussian models using variational Bayes and
  Laplace approximations.*](https://doi.org/10.1080/01621459.2023.2296704) JASA 119(548):2983–2995.
  (the `ngvb` method)
- Cabral, Bolin & Rue (2025). [*Robustness, model checking, and hierarchical models.*](https://doi.org/10.1093/jrsssb/qkae107)
  JRSS-B 87(3):632–652. (the `ng.check` diagnostic)
