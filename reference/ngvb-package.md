# ngvb: Latent Non-Gaussian Models via INLA and Variational Bayes

"ngvb" stands for non-Gaussian variational Bayes. Checks the latent
Gaussian assumption (\[ng.check()\]) and fits latent non-Gaussian models
(\[ngvb()\]) on top of R-INLA. A single unified rgeneric engine
implements the conditional precision \\Q(\theta, V) = D(\theta)^T
\mathrm{diag}(1/V) D(\theta)\\ and owns the V-dependent normalizing
constant, driven by a model registry (\[ngvb_operator()\]) covering iid,
rw1, rw2, ar1, Ornstein-Uhlenbeck, intrinsic-CAR (besag), SAR,
SPDE/Matern, seasonal and generic structure-matrix models; anything else
is added with \[ngvb_custom()\]. \[ngvb_sample()\] and
\[bayes.factor()\] integrate the mixing variables out for a
marginal-likelihood comparison against the Gaussian model.

## Details

INLA is a Suggests dependency (it is not on CRAN); install it from
\<https://www.r-inla.org/download-install\>.

## See also

Useful links:

- <https://github.com/rafaelcabral-ai/ngvb>

- <https://rafaelcabral-ai.github.io/ngvb/>

- Report bugs at <https://github.com/rafaelcabral-ai/ngvb/issues>

## Author

**Maintainer**: Rafael Cabral <rafael.medeiros.cabral@gmail.com>
\[copyright holder\]
