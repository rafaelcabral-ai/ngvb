# Simulated time series with two jumps

A 100-point series used to illustrate a non-Gaussian RW1 latent process.
The underlying signal is smooth except for two sudden jumps, which a
Gaussian RW1 over-smooths and a non-Gaussian RW1 accommodates.

## Usage

``` r
jumpts
```

## Format

A data frame with 100 rows and 2 columns:

- x:

  integer time index, 1..100

- y:

  observed value

## Examples

``` r
plot(jumpts)
```
