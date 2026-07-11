# Plot areal (regional) data on a map.

Plot areal (regional) data on a map.

## Usage

``` r
areal.plot(map, data, palette = "YlOrRd", title = "Value")
```

## Arguments

- map:

  An \`sf\` object with one row per region (its geometry column is
  used).

- data:

  Numeric vector, one value per region, in the row order of \`map\`.

- palette:

  A \`RColorBrewer\` sequential palette name; the default \`"YlOrRd"\`
  shades low values yellow and high values red.

- title:

  Legend title.

## Value

A \`ggplot\` object.

## Examples

``` r
# \donttest{
if (requireNamespace("sf", quietly = TRUE)) {
  nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
  areal.plot(nc, nc$BIR74, title = "births")
}

# }
```
