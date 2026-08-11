# Get started with ngvb

``` r
library(ngvb)
library(ggplot2)
library(INLA)     # ngvb uses INLA as a backend; attach it to call inla() directly
```

## Overview

A latent Gaussian model (LGM), the workhorse of R-INLA, places a
Gaussian prior on a latent field: a random walk, an autoregression, a
spatial effect, or a set of random effects. That prior is convenient but
rigid, because it assumes Gaussianity everywhere. Real data often break
this assumption in just a few places, such as a sudden jump in a time
series, one outlying subject, a sharp boundary between regions, or a
local hotspot in space. A Gaussian model then has to compromise. It
either over-smooths the interesting feature or inflates the variance
everywhere to accommodate it.

`ngvb` relaxes that single assumption. It keeps everything you already
have, the same INLA model and the same code, and replaces the Gaussian
driving noise with a heavier-tailed one. The result adds flexibility
only where the data ask for it and falls back to the Gaussian model
everywhere else. Two functions do the work, and both take an ordinary
`inla` fit that was run with `control.compute = list(config = TRUE)`:

- `ng.check(fit)` asks whether you should relax the Gaussian assumption,
  and where.
- `ngvb(fit)` fits the non-Gaussian extension.

## The model in one paragraph

Every LGM writes its latent field $`\mathbf x`$ through a model-specific
dependency matrix $`\mathbf D(\boldsymbol\theta)`$ that turns
$`\mathbf x`$ into independent Gaussian increments,
``` math
\mathbf D(\boldsymbol\theta)\,\mathbf x \;\overset{d}{=}\; \boldsymbol\Lambda,
\qquad \Lambda_i \sim N(0,\,h_i).
```
For a random walk $`\Lambda_i=x_{i}-x_{i-1}`$, for an autoregression
$`\Lambda_i=x_i-\rho x_{i-1}`$, and so on. Equivalently, the Gaussian
field has precision
$`\mathbf Q(\boldsymbol\theta)=\mathbf D(\boldsymbol\theta)^\top\mathrm{diag}(1/\mathbf h)\,
\mathbf D(\boldsymbol\theta)`$. The non-Gaussian model keeps the same
$`\mathbf D`$ but lets each increment have its own variance $`V_i`$,
``` math
\Lambda_i \mid V_i \sim N(0,\,V_i),\qquad
  V_i \sim \mathrm{IG}(h_i,\,\eta^{-1}h_i^2),\qquad E[V_i]=h_i .
```
The mixing variables $`\mathbf V`$ carry the local flexibility. Where
$`V_i=h_i`$ the model is Gaussian, and where $`V_i`$ is inflated the
increment is allowed to be unusually large. A single non-Gaussianity
parameter $`\eta\ge 0`$ controls how far $`\mathbf V`$ may stray from
$`\mathbf h`$, and $`\eta=0`$ is exactly the Gaussian model. Conditional
on $`\mathbf V`$ the field is still Gaussian, with precision
``` math
\mathbf Q(\boldsymbol\theta,\mathbf V)=\mathbf D(\boldsymbol\theta)^\top\,
   \mathrm{diag}(1/\mathbf V)\,\mathbf D(\boldsymbol\theta).
```
`ngvb` fits the model with a short variational-Bayes loop that
repeatedly calls INLA with $`\mathbf V`$ fixed and then updates
$`\mathbf V`$ and $`\eta`$. The prior on $`\eta`$ is an exponential (a
penalised-complexity prior with rate `alpha.eta`, default $`2`$) that
shrinks toward the Gaussian model, so non-Gaussianity has to earn its
place. Increase `alpha.eta` for stronger shrinkage.

Two variational algorithms are available. The default, `method = "SVI"`,
is a structured mean-field scheme that updates $`q(\mathbf V)`$ and
$`q(\eta)`$ in turn, and it is the more reliable choice. The
alternative, `method = "SCVI"`, integrates $`\mathbf V`$ out of the
$`\eta`$ update analytically and so converges in fewer iterations, but
its collapsed marginal can over-shrink a weak yet genuine effect all the
way to the Gaussian model.

## Which models and likelihoods are supported

The non-Gaussianity in `ngvb` lives entirely in the latent field, on its
driving noise, and never touches the observation model. So **every INLA
likelihood is available**: Gaussian, Poisson, binomial, negative
binomial, survival families, and the rest all work unchanged. You just
fit your LGM as usual and hand it to
[`ng.check()`](https://rafaelcabral-ai.github.io/ngvb/reference/ng.check.md)
or [`ngvb()`](https://rafaelcabral-ai.github.io/ngvb/reference/ngvb.md).

The latent components split into three groups.

**Detected automatically.** If your fit uses one of these,
[`ngvb()`](https://rafaelcabral-ai.github.io/ngvb/reference/ngvb.md) and
[`ng.check()`](https://rafaelcabral-ai.github.io/ngvb/reference/ng.check.md)
find it with no extra arguments: `iid`, `rw1`, `rw2`, `ar1`, `besag`
(the intrinsic CAR), `spde` / `spde2` (Matern fields), and `seasonal`.
Each of these is a single operator applied to independent noise, so it
has the dependency matrix $`\mathbf D`$ the method needs.

**Added by hand.** Some models are not auto-detected but are still
described by
$`\mathbf D(\boldsymbol\theta)\,\mathbf x \;\overset{d}{=}\; \boldsymbol\Lambda`$,
so you can build them with
[`ngvb_operator()`](https://rafaelcabral-ai.github.io/ngvb/reference/ngvb_operator.md)
which expects the matrix $`\mathbf D`$, the vector $`\mathbf h`$. This
covers the simultaneous autoregression (`sar`) or the Ornstein-Uhlenbeck
process (`ou`). It is also possible to use the matrices and operators of
the `ngme2` package to define the LnGM; both are covered in [Custom
models](https://rafaelcabral-ai.github.io/ngvb/articles/custom-models.md).

**Not compatible.** Some components are not a single operator applied to
one field, and the method does not apply to them. The main cases are the
additive spatial models `bym`, `bym2`, `besag2` (a latent field that is
a sum of parts), the conditionally specified `besagproper`, the copy and
mixture constructions `generic1`, `generic2`, `generic3`, `copy`,
`scopy`, `z`, the measurement-error models `mec`, `meb`, and the
nonlinear-covariate effects `sigm`, `revsigm`, `log1exp`, `logdist`. For
an additive model such as `bym`, the practical fix is to split it into
its `besag` and `iid` parts as two separate
[`f()`](https://rdrr.io/pkg/INLA/man/f.html) terms, which
[`ngvb()`](https://rafaelcabral-ai.github.io/ngvb/reference/ngvb.md)
then extends jointly.

Priors carry over the same way: `pc.prec`, `loggamma`, and `normal`
precision priors on the original
[`f()`](https://rdrr.io/pkg/INLA/man/f.html) term are read from the fit
and reused automatically; anything else falls back to a default prior
with a warning.

## A time series with jumps

`jumpts` is a series that is smooth apart from two abrupt jumps.

``` r
ggplot(jumpts, aes(x, y)) +
  geom_line(colour = "grey70") +
  geom_point(colour = "grey30", size = 1.3) +
  labs(title = "jumpts", x = "time", y = "y") +
  theme_minimal(base_size = 12)
```

![The jumpts series: a smooth trend broken by two abrupt vertical
jumps.](ngvb_files/figure-html/rw1-data-1.png)

A first-order random walk models the series through its increments,
``` math
x_i - x_{i-1} \sim N(0,\,1/\tau),
```
which is the operator $`\mathbf D = \sqrt{\tau}\,\mathbf D_0`$ with
$`\mathbf D_0`$ the first-difference matrix. A Gaussian RW1 forces every
increment to share the same variance, so it cannot tell a jump from
ordinary noise. Fit it with INLA:

``` r
LGM <- inla(y ~ -1 + f(x, model = "rw1"), data = jumpts,
            control.compute = list(config = TRUE))
```

Now check whether the latent Gaussian assumption holds:

``` r
check <- ng.check(LGM, compute.fixed = FALSE)
```

![Diagnostic for the RW1 fit: per-index sensitivity spiking at the two
jumps, and the observed overall sensitivity far in the tail of its
Gaussian reference
distribution.](ngvb_files/figure-html/rw1-check-1.png)

The left panel is the per-index sensitivity $`d_i(y)`$, which measures
how much the evidence would improve if increment $`i`$ were allowed to
be non-Gaussian. It spikes at the two jumps. The right panel compares
the observed overall sensitivity $`s_0=\sum_i d_i`$ (the vertical line)
against the distribution it would have under a truly Gaussian field (the
grey curve). The observed value sits far out in the tail, with
$`p\approx0`$, so the Gaussian RW1 is inadequate and we know exactly
where.

Now we extend it to non-Gaussianity. The model is auto-detected, so this
is a single call:

``` r
LnGM <- ngvb(LGM)
#> Components: x [rw1]
#> ngvb: converged after 24 iteration(s);  E[eta] = 3.748
plot(LnGM)
```

![ngvb fit summary: mixing weights inflated at the two jumps and near
one elsewhere, the posterior of the non-Gaussianity parameter, and its
convergence across iterations.](ngvb_files/figure-html/rw1-fit-1.png)

For each component the plot shows the mixing weights $`V_i/h_i`$ (the
posterior mean with a 90 percent credible interval, and a cross at the
plug-in value $`1/E[1/V_i]`$ that INLA is actually refit with, which
sits a little below the mean because $`q(V_i)`$ is right-skewed). The
final panel tracks $`E[\eta]`$ across the iterations. The mixing
variables are inflated precisely at the two jumps and sit near one
elsewhere, so the model bought local flexibility exactly where it was
needed and stayed Gaussian in between. The non-Gaussianity parameter is
clearly positive.

Beyond the point estimate, you can integrate the mixing variables V out
entirely by sampling from their variational posterior and refitting the
LGM at each draw. This gives the Bayes factor against the Gaussian
model, and, pooled with importance weights, posterior summaries of
anything the underlying inla() fit reports:

``` r
samples <- ngvb_sample(LnGM, n.samples = 30)
bayes.factor(samples)               # LnGM vs. LGM marginal-likelihood ratio
#> Bayes factor (non-Gaussian vs Gaussian): 8.32e+23
#>   log10 BF = 23.92  (decisive)
#>   weight ESS = 2.1 of 30 draws
summary(samples)                    # importance-weighted summaries, V integrated out
#> Latent non-Gaussian model, V integrated out over 30 draws
#> Bayes factor vs Gaussian model: 8.32e+23 (log10 = 23.92), weight ESS 2.1
#> 
#> Fixed effects:
#>   (none)
#> 
#> Hyperparameters:
#>                                              mean        sd
#> Precision for the Gaussian observations 9760.4055 23400.414
#> Theta1 for x                               1.8318     0.456
#> 
#> Random effects:
#>   $x (100 nodes)
#>   ID    mean     sd
#> 1  1 -0.8562 0.1148
#> 2  2 -0.6559 0.1072
#> 3  3 -0.7167 0.1038
#> 4  4 -0.6851 0.1102
#> 5  5 -0.9215 0.1185
#> 6  6 -0.6487 0.1278
#>   ... 94 more row(s); use what = "random" to print in full
```

On the usual scale a $`\log_{10}`$ Bayes factor above $`0.5`$, $`1`$, or
$`2`$ is substantial, strong, or decisive evidence for the non-Gaussian
model. Here it is far into the decisive range. The `ess` field is the
effective number of draws behind the estimate; if it is small, raise
`n.samples`.

## Where to next

The same two calls carry over to every other kind of data:

- [Geostatistical
  data](https://rafaelcabral-ai.github.io/ngvb/articles/geostatistical.md)
  fits a Matern (SPDE) field to pressure readings and maps, before
  [`ngvb()`](https://rafaelcabral-ai.github.io/ngvb/reference/ngvb.md)
  is ever run, where the spatial predictions will move.
- [Areal data](https://rafaelcabral-ai.github.io/ngvb/articles/areal.md)
  applies the intrinsic CAR (`besag`) model to burglary rates across
  Columbus, Ohio.
- [Random intercepts and
  slopes](https://rafaelcabral-ai.github.io/ngvb/articles/longitudinal.md)
  treats subject-level random effects in a longitudinal growth study.
- [Custom
  models](https://rafaelcabral-ai.github.io/ngvb/articles/custom-models.md)
  builds a simultaneous autoregression by hand with
  [`ngvb_custom()`](https://rafaelcabral-ai.github.io/ngvb/reference/ngvb_custom.md),
  and borrows ready-made operators from the `ngme2` package.
