#' c(t) at Individual Points
#'
#' Computes the estimated piecewise constant intensity function c(t) at a particular time point t
#'
#' @param t a single time value
#' @param terminal.points a numeric vector containing the time points at which c(t) changes
#' @param ct a numeric vector containing the estimated piecewise constant intensity function c(t)
#'
#' @return a single numeric value of c(t)
#'
#' @export
ct.individual.point<-function(t,terminal.points,ct){
	##if t is outside the time frame specified by t.start and t.end, c(t) is 0
	if( t<min(terminal.points) | t>max(terminal.points) ) return(0)
	##if t is t.start or t.end, c(t) is the appropriate value
	if (t==min(terminal.points)) return( ct[1] )
	if (t==max(terminal.points)) return( ct[length(ct)] )
	##otherwise find the terminal point immediately before t
	ct.indx <- max(which(terminal.points <= t))
	return(ct[ct.indx])
}
