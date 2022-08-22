#' \eqn{c(t)} at all points
#'
#' Computes c(t) at all arbitrary points in a time vector.
#'
#' @param t a numeric vector of time points at which to evaluate the function.
#' @param terminal.points a numeric vector containing the time points at which c(t) changes.
#' @param ct a numeric vector containing the estimated piecewise constant estimated intensity function c(t). The length of ct should be a whole number power of 2.
#'
#' @return A numeric vector containing the value of c(t) at all time points at which to evaluate it.
#'
#' @noRd

CtAllPoints <- function(t,terminal.points,ct){
##wrapper to compute ct.function at all time points in the time vector
##note that this t is now the time vector and not a single time point
	sapply(X=t,FUN=CtIndividualPoint,terminal.points=terminal.points,ct=ct)
}
