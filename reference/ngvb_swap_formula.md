# Rewrite an INLA formula so selected f(\<name\>, ...) terms use the ngvb engine. The engine objects live in \`model.env\` as \`ngvb.model.\<name\>\`; model-specific args (hyper, the original model) are dropped, positional args (index, optional weight) preserved, and \`constr = TRUE\` added for intrinsic components.

Rewrite an INLA formula so selected f(\<name\>, ...) terms use the ngvb
engine. The engine objects live in \`model.env\` as
\`ngvb.model.\<name\>\`; model-specific args (hyper, the original model)
are dropped, positional args (index, optional weight) preserved, and
\`constr = TRUE\` added for intrinsic components.

## Usage

``` r
ngvb_swap_formula(formula, comp.names, rankdef.map, model.env)
```
