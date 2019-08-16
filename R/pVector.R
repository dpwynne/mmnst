p.vector<-function(endpoints.bins,theta,resolution){
##input vector of bins to integrate over and intensity function to integrate
##get out the p-values

L <- length(endpoints.bins)
p <- numeric((L-1))
cat(paste('total count =',L,'\n'))

##integrate theta(t) over bin i to i+1 to get p_k
for (i in 1:(L-1)){
	int.factor<-integrate(theta.function,endpoints.bins[i],endpoints.bins[i+1],theta=theta,t.start=endpoints.bins[1],resolution=resolution)$value
	#print(int.factor)
	p[i] <- 1-exp(-int.factor)
	if(i%%5000==0) cat(paste(i,'bins processed\n'))
}##end for loop

q <- -log(1-p) ##return qk - but this is just int.factor, I think?
}
