#' Additive intensity function \eqn{\theta(t)}
#'
#' Calculates the additive intensity function \eqn{\theta(t)} introduced in Ramezan *et al.* (2014).
#'
#' @param t a numeric vector of time points at which to evaluate the function.
#' @param f a numeric vector containing frequency values.
#' @param w0 a numeric vector containing initial phase values.
#' @param eta a numeric vector containing \eqn{\eta} values (contribution of each periodic component to the intensity function).
#' @param gamma a numeric vector containing \eqn{\gamma} values (amplitude of each periodic component in the function).
#' @param terminal.points a numeric vector containing the endpoints of the dyadic partitioning.
#' @param ct a numeric vector containing the estimated piecewise constant intensity function \eqn{c(t)}. The length of ct should be a whole number power of 2.
#'
#' @return A numeric vector containing the values of the additive intensity function calculated at given time points.
#'
#' @references Ramezan, R., Marriott, P., and Chenouri, S. (2014), *Statistics in Medicine*, **33**(2), 238-256. doi: 10.1002/sim.5923.
#'
#' @export

ThetaAdditive <- function(t,f,w0,eta,gamma,terminal.points,ct){
	CtAllPoints(t,terminal.points,ct)*(1-sum(eta))+sum(eta*gamma*nu(f*t+w0))
}
