theta.m <- function(t,f,w0,eta,gamma,terminal.points,ct){
##t is time vector from time.start to time.end with user-defined resolution
#f is frequencies, w0 is phases
##eta and gamma represent contribution and amplitude of nu to intensity function
	ct.all.points(t,terminal.points,ct)*( (1-sum(eta))+sum(eta*gamma*nu(f*t+w0)) )
}

