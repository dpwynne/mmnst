#' \eqn{c(t)} at Individual Points
#'
#' Computes the estimated piecewise constant intensity function c(t) at a particular time point t1. The function is defined to be right-continuous.
#'
#' @param t1 a single time value.
#' @param terminal.points a numeric vector containing the time points at which c(t) changes. This vector is usually of length \eqn{2^J+1}.
#' @param ct a numeric vector containing the estimated piecewise constant intensity function c(t).
#'
#' @return A single numeric value of c(t).
#'
#' @export

CtIndividualPoint<-function(t1, terminal.points, ct){
	##if t1 is outside the time frame specified by t.start and t.end, c(t1) is 0
	if( t1 < min(terminal.points) | t1 > max(terminal.points) ) return(0)
  if (t1 == min(terminal.points)) return( ct[1] )
  if (t1 == max(terminal.points)) return(ct[length(ct)])

	##otherwise find the terminal point immediately before t1
	ct.indx <- which(terminal.points[-length(terminal.points)] <= t1  & terminal.points[-1] > t1)
	return(ct[ct.indx])

}


