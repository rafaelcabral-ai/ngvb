# Recover a binary adjacency from a fitted besag/ICAR component's prior precision. config\$Qprior block = tau \* (diag(nnbs) - A); off-diagonal non-zeros mark neighbours.

Recover a binary adjacency from a fitted besag/ICAR component's prior
precision. config\$Qprior block = tau \* (diag(nnbs) - A); off-diagonal
non-zeros mark neighbours.

## Usage

``` r
ngvb_recover_adjacency(fit, comp.name)
```
