#' @useDynLib heatindex, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @title
#' Wet-bulb temperature
#' @description
#' `wetbulb` calculates the thermodynamic (or psychrometric) wet-bulb (or ice-bulb) temperature using the Rankine-Kirchhoff approximations.
#' @author
#' David M. Romps \email{romps@berkeley.edu}
#' @param p The total air pressuire in Pa. This can be a single number, a vector, a matrix, or an array.
#' @param T The absolute air temperature in Kelvin. This can be a single number, a vector, a matrix, or an array.
#' @param rh The relative humidity of the air, with values in the range of 0 to 1, with respect to saturation over liquid water for air temperatures over 273.16 K and with respect to saturation over ice for air temperatures lower than 273.16 K. This can be a single number, a vector, a matrix, or an array.
#' @param psychrometric A logical indicating whether to calculate the thermodynamic wet-bulb temperature (if `FALSE`) or the psychometric (a.k.a., ventilated or aspirated) wet-bulb temperature (if `TRUE`). Default is `FALSE`.
#' @param icebulb A logical indicating whether to calculate the temperature of an ice-bulb (if `TRUE`) or wet-bulb (if `FALSE`). Default is `FALSE`.
#' @param verbose A logical indicating whether or not to print warning messages. Default is `TRUE`.
#' @param lewis The Lewis number for moist air. Default is `0.85`.
#' @return The values of the wet-bulb temperature, in Kelvin, in the same shape as `p`, `T`, and `rh`.
#' @examples
#' wetbulb(1e5,300,0)
#' wetbulb(1e5,301:310,0)
#' wetbulb(1e5,301:310,1:10/10)
#' wetbulb(1:10*1e4,301:310,1:10/10)
#' @references
#' Romps, D. M. (2025). Wet-bulb temperature from pressure, relative humidity, and air temperature. In review.
#' @export
wetbulb <- function(p, T, rh, psychrometric=FALSE, icebulb=FALSE, verbose=TRUE, lewis=0.85) {
  args <- list(p = p, T = T, rh = rh)
  sizes <- sapply(args, length)
  dims <- lapply(args, function(x) if (is.null(dim(x))) length(x) else dim(x))
  maxsize <- max(sizes)
  formatDim <- function(d) if (length(d) == 1) as.character(d) else paste(d, collapse = " x ")
  if (maxsize > 1) {
    idx <- which.max(sizes)
    ref <- dims[[idx]]
    refvar <- names(args)[idx]
    for (nm in names(args)) {
      if (sizes[nm] > 1 && !identical(dims[[nm]], ref))
        stop("Input variables are incompatible: ", nm, " is of size ", 
             formatDim(dims[[nm]]), " and ", refvar, " is of size ", formatDim(ref))
    }
    Tw <- wetbulb_vec(p, T, rh, psychrometric, icebulb, verbose, lewis)
    if (length(ref) > 1) dim(Tw) <- ref
    return(Tw)
  }
  wetbulb_vec(p, T, rh, psychrometric, icebulb, verbose, lewis)
}
