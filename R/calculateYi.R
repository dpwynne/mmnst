calculate.yi<-function(spike.train,endpoints.bins,q){
##spike.train is the spike train
##endpoints.bins is the vector of start/endpoints of the bins

N <- length(spike.train)
xi <- numeric(N-1)

for(i in 1:(N-1)){
	##fix to get around "no non-missing arguments to min/max" bug
	indx1 <- min(length(endpoints.bins),which(endpoints.bins>spike.train[i]))
	indx2 <- max(1,which(endpoints.bins<spike.train[i+1]))-1
	sum.term <- (indx2>indx1)*sum(q[indx1:indx2])

	r <- runif(1)
	xi[i] <- sum.term-log(1-r*(1-exp(-q[(indx2+1)])))
}##end for loop

y <- 1-exp(-xi)

return(y)
}
