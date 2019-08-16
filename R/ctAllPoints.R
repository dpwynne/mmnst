ct.all.points <- function(t,terminal.points,ct){
##wrapper to compute ct.function at all time points in the time vector
##note that this t is now the time vector and not a single time point
	sapply(X=t,FUN=ct.individual.point,terminal.points=terminal.points,ct=ct)
}