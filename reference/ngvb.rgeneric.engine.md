# rgeneric model function implementing the ngvb conditional precision.

Not called directly; passed to \`INLA::inla.rgeneric.define()\` by
\[ngvb_rgeneric()\]. Reads \`Dfunc\`, \`Vinv\`, \`rankdef\`, \`ntheta\`,
\`theta.initial\`, \`logprior\`, \`graph.pattern\` from its definition
environment.

## Usage

``` r
ngvb.rgeneric.engine(
  cmd = c("graph", "Q", "mu", "initial", "log.norm.const", "log.prior", "quit"),
  theta = NULL
)
```
