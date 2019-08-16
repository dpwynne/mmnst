ct.individual.point<-function(t,terminal.points,ct){
##t is a single time point
##terminal.points a vector of points representing the boundaries of each time partition
##ct is the estimated multiscale intensity function?; must be numeric vector
	##if t is outside the time frame specified by t.start and t.end
	##c(t) is 0
	if( t<min(terminal.points) | t>max(terminal.points) ) return(0)
	##if t is t.start or t.end, c(t) is the appropriate value
	if (t==min(terminal.points)) return( ct[1] )
	if (t==max(terminal.points)) return( ct[length(ct)] )
	##otherwise find the terminal point immediately before t
	ct.indx <- max(which(terminal.points <= t))
	return(ct[ct.indx])
}