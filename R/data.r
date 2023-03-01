#' UK postcode location and population.
#'
#' A dataset containing the population and the centroid coordinates of all UK postcodes.
# Overwrite this variable if you wish to use a different geographic and demographic data.
#'
#' @format A data frame with 1048575 rows and 7 variables:
#' \describe{
#'   \item{postcode}{}
#'   \item{Total}{}
#'   \item{latitude}{}
#'   \item{longitude}{}
#' }
#' @source \url{https://www.ons.gov.uk}
"postcode.data"

#' Real baseline matrix corresponding to the simulated data in \code{simulation_data}.
#'
#' @format A \code{Matrix} with 10000 rows and 100 columns.
"baseline_for_sim"


#' Simulation data.
#'
#' A dataset containing.
#'
#' @format A data frame with 5028 rows and 9 variables:
#' \describe{
#'   \item{week}
#'   \item{postcode}{}
#'   \item{population}{}
#'   \item{sim}{}
#'   \item{type}{}
#'   \item{latitude}{}
#'   \item{longitude}{}
#'   \item{y}{} 
#'   \item{x}{}  
#' }
#' @source \url{https://www.ons.gov.uk}
"simulation_data"

#' GB boundaries.
#'
#' A dataset containing the boundaries of the Euro constituencies of all Great Britain,
#' Coordinate reference systems (CRS) EPSG code 4326, which uses units of longitude and latitude
#' on the  World Geodetic System 1984 (WGS84) ellipsoid.
#' This information is license under the Open Government Licence v3.0, see \url{https://www.nationalarchives.gov.uk/doc/open-government-licence/version/3/}.
#'  
#'
#' @import sf
#' @format Simple feature collection.
#' @source \url{https://osdatahub.os.uk/}
"GB_region_boundaries"
