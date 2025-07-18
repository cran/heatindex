#' @useDynLib heatindex, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @title
#' Relative humidity from wet-bulb temperature
#' @description
#' `rh_from_wetbulb` calculates the relative humidity from the thermodynamic (or psychrometric) wet-bulb (or ice-bulb) temperature using the Rankine-Kirchhoff approximations.
#' @author
#' David M. Romps \email{romps@berkeley.edu}
#' @param p The total air pressuire in Pa. This can be a single number, a vector, a matrix, or an array.
#' @param T The absolute air temperature in Kelvin. This can be a single number, a vector, a matrix, or an array.
#' @param Tw The thermodynamic (or psychrometric) wet-bulb (or ice-bulb) temperature in Kelving. This can be a single number, a vector, a matrix, or an array.
#' @param psychrometric A logical indicating whether to interpret `Tw` as the psychrometric (if `TRUE`) or thermodynamic (if `FALSE`) version. Default is `FALSE`.
#' @param icebulb A logical indicating whether to interpret `Tw` as the ice-bulb (if `TRUE`) or wet-bulb (if `FALSE`) version. Default is `FALSE`.
#' @param verbose A logical indicating whether or not to print warning messages. Default is `TRUE`.
#' @param lewis The Lewis number for moist air. Default is `0.85`.
#' @return Relative humidity in the same shape as `p`, `T`, and `Tw`. The relative humidity is reported with respect to liquid water if `T` is greater than or equal to 273.16 K and with respect to ice if `T` is less than 273.16 K.
#' @examples
#' rh_from_wetbulb(1e5,300,290)
#' rh_from_wetbulb(1e5,301:310,290)
#' rh_from_wetbulb(1e5,301:310,291:300)
#' rh_from_wetbulb(1:10*1e4,301:310,291:300)
#' @references
#' Romps, D. M. (2025). Wet-bulb temperature from pressure, relative humidity, and air temperature. In review.
#' @export
rh_from_wetbulb <- function(p, T, Tw, psychrometric=FALSE, icebulb=FALSE, verbose=TRUE, lewis=0.85) {
  args <- list(p = p, T = T, Tw = Tw)
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
    Tw <- rh_from_wetbulb_vec(p, T, Tw, psychrometric, icebulb, verbose, lewis)
    if (length(ref) > 1) dim(Tw) <- ref
    return(Tw)
  }
  rh_from_wetbulb_vec(p, T, Tw, psychrometric, icebulb, verbose, lewis)
}
