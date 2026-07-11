## ---------------------------------------------------------------------------
## Map helpers for areal and geostatistical data. These are static ggplot2 / sf
## versions of the leaflet helpers from the original ngvb package, kept under the
## same names and roles (areal.plot, geo.plot) but reproducible offline: no tile
## server, no leaflet, so they render inside a vignette and pass R CMD check.
## ---------------------------------------------------------------------------

#' A shared, recessive ggplot theme for the ngvb plots.
#' @keywords internal
ngvb_theme <- function(base = 12) {
  ggplot2::theme_minimal(base_size = base) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = base),
      legend.position = "right")
}

#' Plot areal (regional) data on a map.
#'
#' @param map An `sf` object with one row per region (its geometry column is used).
#' @param data Numeric vector, one value per region, in the row order of `map`.
#' @param palette A `RColorBrewer` sequential palette name; the default `"YlOrRd"`
#'   shades low values yellow and high values red.
#' @param title Legend title.
#' @return A `ggplot` object.
#' @examples
#' \donttest{
#' if (requireNamespace("sf", quietly = TRUE)) {
#'   nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
#'   areal.plot(nc, nc$BIR74, title = "births")
#' }
#' }
#' @export
areal.plot <- function(map, data, palette = "YlOrRd", title = "Value") {
  if (!requireNamespace("sf", quietly = TRUE))
    stop("areal.plot() needs the 'sf' package.")
  map[[".ngvb_value"]] <- as.numeric(data)
  ## Plot in the geometry's own coordinates. Dropping the CRS stops ggplot from
  ## attempting a datum transform, which fails when PROJ cannot build a pipeline
  ## for the map's CRS (e.g. an undefined engineering CRS, as in the Columbus map).
  map <- suppressWarnings(sf::st_set_crs(map, NA))
  ggplot2::ggplot(map) +
    ggplot2::geom_sf(ggplot2::aes(fill = .data[[".ngvb_value"]]),
                     colour = "grey45", linewidth = 0.15) +
    ggplot2::scale_fill_distiller(palette = palette, direction = 1, name = title) +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(legend.position = "right")
}

#' Plot geostatistical (point-referenced) data, optionally over a mesh.
#'
#' @param data A data frame whose first three columns are longitude, latitude and
#'   the value to colour by.
#' @param mesh Optional `inla.mesh`/`fmesher` mesh; its triangulation is drawn underneath.
#' @param n Highlight (enlarge) the `n` locations with the largest values.
#' @param palette A `RColorBrewer` sequential palette name (default `"YlOrRd"`).
#' @param title Legend title (defaults to the third column's name).
#' @return A `ggplot` object.
#' @export
geo.plot <- function(data, mesh = NULL, n = 0, palette = "YlOrRd", title = NULL) {
  if (is.null(title)) title <- colnames(data)[3]
  d <- data.frame(lon = data[[1]], lat = data[[2]], value = data[[3]])
  big <- rep(FALSE, nrow(d))
  if (n > 0) big[order(d$value, decreasing = TRUE)[seq_len(n)]] <- TRUE
  d$big <- big

  p <- ggplot2::ggplot()
  if (!is.null(mesh)) {                                  # triangulation underlay
    tv <- mesh$graph$tv; loc <- mesh$loc
    e  <- rbind(tv[, c(1, 2)], tv[, c(2, 3)], tv[, c(3, 1)])
    segs <- data.frame(x = loc[e[, 1], 1], y = loc[e[, 1], 2],
                       xend = loc[e[, 2], 1], yend = loc[e[, 2], 2])
    p <- p + ggplot2::geom_segment(
      data = segs, ggplot2::aes(x = .data$x, y = .data$y,
                                xend = .data$xend, yend = .data$yend),
      colour = "grey80", linewidth = 0.2)
  }
  p +
    ggplot2::geom_point(
      data = d, ggplot2::aes(x = .data$lon, y = .data$lat,
                             fill = .data$value, size = .data$big),
      shape = 21, colour = "grey20", stroke = 0.3) +
    ggplot2::scale_fill_distiller(palette = palette, direction = 1, name = title) +
    ggplot2::scale_size_manual(values = c(`FALSE` = 1.8, `TRUE` = 3.8), guide = "none") +
    ggplot2::coord_equal() +
    ggplot2::labs(x = "longitude", y = "latitude") +
    ngvb_theme()
}
