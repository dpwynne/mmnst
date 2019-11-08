#' Additive \theta(t)
#'
#' Calculates the additive intensity function introduced in Ramezan et al. (2014)
#'
#' @param t a numeric vector of time points at which the intensity function is to be calculated
#' @param f a numeric vector containing frequency values
#' @param w0 a numeric vector containing initial phase values
#' @param eta a numeric vector containing eta values (contribution of each periodic component to the intensity function)
#' @param gamma a numeric vector containing gamma values (amplitude of each periodic component in the function)
#' @param terminal.points a numeric vector containing the endpoints of the dyadic partitioning
#' @param ct a numeric vector containing the estimated piecewise constant intensity function c(t). The length of ct should be a whole number power of 2.
#'
#' @return a numeric vector containing the values of the additive intensity function calculated at given time points
#'
#' @export

theta.a <- function(t,f,w0,eta,gamma,terminal.points,ct){
	ct.all.points(t,terminal.points,ct)*(1-sum(eta))+sum(eta*gamma*nu(f*t+w0))
}
