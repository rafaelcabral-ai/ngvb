# Plot geostatistical (point-referenced) data, optionally over a mesh.

Plot geostatistical (point-referenced) data, optionally over a mesh.

## Usage

``` r
geo.plot(data, mesh = NULL, n = 0, palette = "YlOrRd", title = NULL)
```

## Arguments

- data:

  A data frame whose first three columns are longitude, latitude and the
  value to colour by.

- mesh:

  Optional \`inla.mesh\`/\`fmesher\` mesh; its triangulation is drawn
  underneath.

- n:

  Highlight (enlarge) the \`n\` locations with the largest values.

- palette:

  A \`RColorBrewer\` sequential palette name (default \`"YlOrRd"\`).

- title:

  Legend title (defaults to the third column's name).

## Value

A \`ggplot\` object.
