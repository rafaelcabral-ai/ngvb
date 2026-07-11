# Sample the mixing variables and refit INLA at each draw.

Draws \`n.samples\` values of the mixing vector \`V\` from the
variational posterior \`q(V)\` of a fitted \[ngvb()\] object, and refits
the conditionally-Gaussian LGM with INLA at each draw. The returned
object holds the list of \`inla\` fits together with the importance
weights needed to treat \`V\` as integrated out (see \[bayes.factor()\],
\[summary.ngvb.samples()\]).

## Usage

``` r
ngvb_sample(object, n.samples = 50, seed = NULL, verbose = interactive())
```

## Arguments

- object:

  An \`ngvb\` object.

- n.samples:

  Number of posterior draws of \`V\` (each is one \`inla\` refit).

- seed:

  Optional integer seed for reproducibility.

- verbose:

  Show a progress bar.

## Value

An object of class \`ngvb.samples\`: a list with the \`inla\` \`fits\`,
the drawn \`V\`, per-draw log marginal likelihoods \`logm\` and log
importance weights \`logw\`, and the Gaussian-model log marginal
likelihood \`logm.LGM\`.

## See also

\[bayes.factor()\]
