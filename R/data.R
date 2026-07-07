#' Simulated time series with two jumps
#'
#' A 100-point series used to illustrate a non-Gaussian RW1 latent process. The
#' underlying signal is smooth except for two sudden jumps, which a Gaussian RW1
#' over-smooths and a non-Gaussian RW1 accommodates.
#'
#' @format A data frame with 100 rows and 2 columns:
#' \describe{
#'   \item{x}{integer time index, 1..100}
#'   \item{y}{observed value}
#' }
#' @examples
#' plot(jumpts)
"jumpts"

#' Orthodontic growth data
#'
#' Growth measurements on 27 children (16 boys, 11 girls) at ages 8, 10, 12, 14,
#' preprocessed for a linear mixed model with random intercepts and slopes. Used
#' to illustrate a two-component (random intercept + random slope) latent
#' non-Gaussian model.
#'
#' @format A data frame with 108 rows and 6 columns:
#' \describe{
#'   \item{subject}{subject id (random-intercept index)}
#'   \item{subject2}{subject id (random-slope index; distinct name required by INLA)}
#'   \item{Female}{indicator for female subjects}
#'   \item{time}{age recentred to 0,1,2,3}
#'   \item{tF}{Female x time interaction}
#'   \item{value}{measured distance}
#' }
#' @source Derived from the \code{Orthodont} data of Potthoff and Roy (1964),
#'   as distributed with the \pkg{nlme} package.
#' @examples
#' summary(Orthodont)
"Orthodont"
