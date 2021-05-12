#' Simulated point process data set
#'
#' A simulated dataset with seven time periods over a unit square with population
#' density given by \code{square_pop}. Primarily used for testing.
#'
#' @format A data frame with 183 rows and 3 variables:
#' \describe{
#'   \item{x}{the x coordinate of the point}
#'   \item{y}{the y coordinate of the point}
#'   \item{t}{the time period of the point}
#' }
"dat"

#' A unit square spatial polygon
#'
#' A SpatialPolygonsDataFrame describing a unit square
"square"

#' Simulated population density across unit square
#'
#' A SpatialPolygonsDataFrame covering a unit square with 11x11 square polygons, each with a
#' different simulated population density. Used to simulate the data in \code{dat}
"square_pop"
