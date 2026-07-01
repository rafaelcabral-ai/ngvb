# ngvb2

**Latent non-Gaussian models via R-INLA and variational Bayes** — a rebuild of
[`ngvb`](https://github.com/rafaelcabral96/ngvb).

`ngvb2` lets you take an ordinary R-INLA latent *Gaussian* model, **check** whether the
Gaussian assumption is adequate, and **extend** it to a latent *non-Gaussian* model — each in
one line.

```r
LGM  <- inla(y ~ f(s, model = "rw1"), data = d, control.compute = list(config = TRUE))

ng.check(LGM)   # is the latent Gaussian assumption adequate?
LnGM <- ngvb(LGM)   # fit the non-Gaussian extension
```

## The idea

A latent Gaussian field has precision `Q = Dᵀ D` for a model-specific *dependency matrix* `D`
(the difference operator for a random walk, `I − ρW` for a SAR, `κ²C + G` for a Matérn/SPDE, …).
The non-Gaussian extension replaces the driving noise, which is equivalent to conditioning on
**mixing variables** `V`:

```
Q(θ, V) = D(θ)ᵀ diag(1/V) D(θ)
```

`ngvb2` fits this with a variational-Bayes loop that alternates an INLA fit (for fixed `V`) with
closed-form updates of `V` and the non-Gaussianity parameter `η`. A single unified `rgeneric`
engine implements `Q(θ, V)` for every model and — unlike INLA's `generic0` — correctly owns the
`V`-dependent normalizing constant.

## Models

`ngvb(fit)` and `ng.check(fit)` auto-detect the component type from the fitted object
(`fit$model.random`) and rebuild the operator:

| Model | `f(..., model = )` | Auto from fit? |
|---|---|---|
| i.i.d. random effects | `"iid"` | yes |
| Random walks | `"rw1"`, `"rw2"` | yes |
| Autoregressive | `"ar1"` | yes |
| CAR / ICAR (areal) | `"besag"`, `"besagproper"` | yes (graph recovered) |
| SAR (areal) | *(not native to INLA)* | via `ngvb_operator("sar", W)` |
| Matérn / SPDE | `inla.spde2.matern(mesh)` | pass the mesh: `components = list(s = ngvb_operator("spde", spde = spde))` |

Only SPDE needs an object passed (the mesh isn't retained in the fit); everything else is a
true one-liner on the existing fit.

## Multiple components

Additive models just work — each `f()` term gets its own mixing variables and non-Gaussianity
parameter, all fit jointly:

```r
LGM  <- inla(y ~ f(time, model = "ar1") + f(subj, model = "iid"), data = d,
             control.compute = list(config = TRUE))
LnGM <- ngvb(LGM)          # extends both components; each gets its own eta
summary(LnGM)
```

## Installation

Requires [R-INLA](https://www.r-inla.org/):

```r
install.packages("INLA", repos = c(getOption("repos"), INLA = "https://inla.r-inla-download.org/R/stable"))
# then, from the ngvb2 source directory:
# devtools::install("ngvb2")
```

## References

- Cabral, Bolin & Rue (2024). *Fitting latent non-Gaussian models using variational Bayes and
  Laplace approximations.* (the `ngvb` method)
- Cabral, Bolin & Rue (2025). *Robustness, model checking, and hierarchical models.* JRSS-B
  87(3):632–652. (the `ng.check` diagnostic)

## Status

Rebuild in progress. Engine + operators (iid, rw1, rw2, ar1, SAR, CAR, SPDE), `ngvb()`,
`ng.check()`, and a passing test suite are in place; a `cgeneric` fast engine and additional
models (AR(p), Ornstein–Uhlenbeck, separable space-time) are planned.
