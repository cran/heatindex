#' @useDynLib heatindex, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @title
#' Heat index
#' @description
#' `heatindex` is a simpler and faster version of the heat index that was originally defined in 1979, used by the U.S. National Weather Service, extended to all combinations of temperature and humidity in 2022, and then made simpler and faster in 2025. This simpler and faster version uses a simpler set of physiological equations and a faster computational algorithm without altering the values of the heat index above 300 K (27 C, 80 F) and with only minor changes in the heat index at lower temperatures.
#' @author
#' Yi-Chuan Lu \email{yclu@berkeley.edu} and David M. Romps \email{romps@berkeley.edu}
#' @param T The absolute air temperature in Kelvin. This can be a single number, a vector, a matrix, or an array, but its dimensions must match those of `rh`.
#' @param rh The relative humidity of the air, with values in the range of 0 to 1, with respect to saturation over liquid water for air temperatures over 273.16 K and with respect to saturation over ice for air temperatures lower than 273.16 K. This can be a single number, a vector, a matrix, or an array, but its dimensions must match those of `T`.
#' @return The values of the heat index, in Kelvin, in the same shape as `T` and `rh`.
#' @examples
#' heatindex(300,0.5)
#' heatindex(295:305,0:10/10)
#' @references
#' Steadman, R. G. (1979). The assessment of sultriness. Part I: A temperature-humidity index based on human physiology and clothing science. Journal of Applied Meteorology, 18, 861-873. \doi{10.1175/1520-0450(1979)018<0861:TAOSPI>2.0.CO;2}
#'
#' Lu, X. and Romps, D. M. (2022). Extending the heat index. Journal of Applied Meteorology, 61, 10, 1367--1383. \doi{10.1175/jamc-d-22-0021.1}
#'
#' Lu et al. (2025). Simpler and faster: an improved heat index. In review. For citation details, see \url{https://heatindex.org/docs/citation/}.
#' @export
heatindex <- function(T, rh) {
  args <- list(T = T, rh = rh)
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
    hi <- heatindex_vec(T, rh)
    if (length(ref) > 1) dim(hi) <- ref
    return(hi)
  }
  heatindex_vec(T, rh)
}
