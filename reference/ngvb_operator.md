# Construct an ngvb operator descriptor.

Construct an ngvb operator descriptor.

## Usage

``` r
ngvb_operator(type, ...)
```

## Arguments

- type:

  Model type: one of \`"iid"\`, \`"rw1"\`, \`"rw2"\`, \`"ar1"\`,
  \`"sar"\`, \`"car"\`, \`"spde"\`, \`"ou"\`, \`"seasonal"\`.

- ...:

  Model-specific arguments (e.g. \`n\` for rw/ar/iid, \`W\` for sar/car,
  \`loc\` for ou). For \`"spde"\`, pass \`spde =
  INLA::inla.spde2.pcmatern(mesh, prior.range = , prior.sigma = )\` (or
  \`mesh\` + \`prior.range\` + \`prior.sigma\` and ngvb builds it); the
  non-Gaussian SPDE uses that PC prior on the practical range and
  marginal SD. Plain \`inla.spde2.matern()\` fields are rejected — only
  the PC-prior parameterization is supported.

## Value

An operator descriptor: a list with \`Dfunc(theta)\`, the constant
vector \`h\`, \`rankdef\`, and prior/graph metadata.

## Examples

``` r
op <- ngvb_operator("rw1", n = 10)
dim(op$Dfunc(0))       # D(theta): the 9 x 10 first-difference operator
#> [1]  9 10
ngvb_precision(op)     # Q = D^T diag(1/h) D  (the RW1 structure matrix)
#> 10 x 10 sparse Matrix of class "dgCMatrix"
#>                                                                            
#>  [1,]  54.59815 -54.59815   .         .         .         .         .      
#>  [2,] -54.59815 109.19630 -54.59815   .         .         .         .      
#>  [3,]   .       -54.59815 109.19630 -54.59815   .         .         .      
#>  [4,]   .         .       -54.59815 109.19630 -54.59815   .         .      
#>  [5,]   .         .         .       -54.59815 109.19630 -54.59815   .      
#>  [6,]   .         .         .         .       -54.59815 109.19630 -54.59815
#>  [7,]   .         .         .         .         .       -54.59815 109.19630
#>  [8,]   .         .         .         .         .         .       -54.59815
#>  [9,]   .         .         .         .         .         .         .      
#> [10,]   .         .         .         .         .         .         .      
#>                                    
#>  [1,]   .         .         .      
#>  [2,]   .         .         .      
#>  [3,]   .         .         .      
#>  [4,]   .         .         .      
#>  [5,]   .         .         .      
#>  [6,]   .         .         .      
#>  [7,] -54.59815   .         .      
#>  [8,] 109.19630 -54.59815   .      
#>  [9,] -54.59815 109.19630 -54.59815
#> [10,]   .       -54.59815  54.59815
```
