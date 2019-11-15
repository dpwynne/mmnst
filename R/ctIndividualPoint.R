#' c(t) at Individual Points
#'
#' Computes the estimated piecewise constant intensity function c(t) at a particular time point t0
#'
#' @param t0 a single time value
#' @param terminal.points a numeric vector containing the time points at which c(t) changes
#' @param ct a numeric vector containing the estimated piecewise constant intensity function c(t)
#'
#' @return a single numeric value of c(t)
#'
#' @export

ct.individual.point<-function(t0, terminal.points, ct){
	##if t0 is outside the time frame specified by t.start and t.end, c(t0) is 0
	if( t0 < min(terminal.points) | t0 > max(terminal.points) ) return(0)
	##if t0 is t.start or t.end, c(t0) is the appropriate value
	if (t0 == min(terminal.points)) return( ct[1] )
	if (t0 == max(terminal.points)) return( ct[length(ct)] )
	##otherwise find the terminal point immediately before t0
	ct.indx <- max(which(terminal.points <= t0))
	return(ct[ct.indx])
}
