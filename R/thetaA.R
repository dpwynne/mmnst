theta.a <- function(t,f,w0,eta,gamma,terminal.points,ct){
	ct.all.points(t,terminal.points,ct)*(1-sum(eta))+sum(eta*gamma*nu(f*t+w0))
}
