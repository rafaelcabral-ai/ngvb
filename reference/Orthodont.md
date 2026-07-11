# Orthodontic growth data

Growth measurements on 27 children (16 boys, 11 girls) at ages 8, 10,
12, 14, preprocessed for a linear mixed model with random intercepts and
slopes. Used to illustrate a two-component (random intercept + random
slope) latent non-Gaussian model.

## Usage

``` r
Orthodont
```

## Format

A data frame with 108 rows and 6 columns:

- subject:

  subject id (random-intercept index)

- subject2:

  subject id (random-slope index; distinct name required by INLA)

- Female:

  indicator for female subjects

- time:

  age recentred to 0,1,2,3

- tF:

  Female x time interaction

- value:

  measured distance

## Source

Derived from the `Orthodont` data of Potthoff and Roy (1964), as
distributed with the nlme package.

## Examples

``` r
summary(Orthodont)
#>     subject       Female            time          value             tF       
#>  Min.   : 1   Min.   :0.0000   Min.   :0.00   Min.   :16.50   Min.   :0.000  
#>  1st Qu.: 7   1st Qu.:0.0000   1st Qu.:0.75   1st Qu.:22.00   1st Qu.:0.000  
#>  Median :14   Median :0.0000   Median :1.50   Median :23.75   Median :0.000  
#>  Mean   :14   Mean   :0.4074   Mean   :1.50   Mean   :24.02   Mean   :1.019  
#>  3rd Qu.:21   3rd Qu.:1.0000   3rd Qu.:2.25   3rd Qu.:26.00   3rd Qu.:2.000  
#>  Max.   :27   Max.   :1.0000   Max.   :3.00   Max.   :31.50   Max.   :4.000  
#>     subject2 
#>  Min.   : 1  
#>  1st Qu.: 7  
#>  Median :14  
#>  Mean   :14  
#>  3rd Qu.:21  
#>  Max.   :27  
```
