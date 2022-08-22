#' A periodic function
#'
#' Computes the function \eqn{\nu(x) = 1 + cos(2 \pi x)}.
#'
#' @param x an arbitrary numeric vector or scalar.
#'
#' @return The value of the function \eqn{1 + cos(2 \pi x)}.
#'
#' @noRd

nu <- function(x){
##periodic function
	1+cos(2*pi*x)
}

