#' Setup likelihood function
#'
#' Create a list of common parameters necessary to calculate likelihood functions for additive and multiplicative models.
#'
#' @param terminal.points a numeric vector containing the endpoints of the dyadic partitioning
#'
#' @return A list of length 5 is returned. The list contains the following scalars:
#'
#' DeltaDi: the delta D value in the likelihood derivation; see Ramezan *et al*. (2014).
#'
#' Di.1: the first Di value in the likelihood derivation; see Ramezan *et al*. (2014).
#'
#' Di.0: the initial Di value in the likelihood derivation; see Ramezan *et al*. (2014).
#'
#' T.data: The length of the data recording window.
#'
#' J: The resolution based on which the piecewise constant intensity function \eqn{c(t)} has been estimated.
#'
#' @references Ramezan, R., Marriott, P., and Chenouri, S. (2014), *Statistics in Medicine*, **33**(2), 238-256. doi: 10.1002/sim.5923.
#'
#' @export

SetupLikelihoods <-function(terminal.points){
  D.i.plus.one <- terminal.points[-c(1)]
  D.i <- terminal.points[-length(terminal.points)]
  DeltaDi <- D.i.plus.one-D.i
	T.data <- max(terminal.points)-min(terminal.points)
	J <- log(length(terminal.points)-1, 2)

	return(list(DeltaDi=DeltaDi,Di.1=D.i.plus.one,Di.0=D.i,T.data=T.data,J=J))
}
