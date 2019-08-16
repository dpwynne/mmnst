spike.simulation<-function(nruns,Time,pass.arg,ct,type="m",envelope.function=function(t) return(1)){
##returns a list of vectors containing simulated spike times
##nruns is number of simulated runs desired
##Time is vector of at least start and end times
##type is "m" for muliplicative and "a" for additive
##pass.arg is a list of, in order, f, w0, eta, and gamma (parameters)
##ct is the simulated piecewise constant intensity function
##envelope.function is the envelope function

simulated.spikes<-vector("list",length=nruns)

##unpack arguments from pass.arg
f<-pass.arg[[1]]
w0<-pass.arg[[2]]
eta<-pass.arg[[3]]
gamma<-pass.arg[[4]]

t.start<-min(Time)
t.end<-max(Time)
terminal.points <- identify.terminal.points(t.start,t.end,log(length(ct),2))
##get terminal points from ct

##type determines which function to call
if (type == "m"){
spike.fun<-theta.m
}
if (type == "a"){
spike.fun<-theta.a
}

simulated.intensity<-sapply(Time,spike.fun,f=f,w0=w0,eta=eta,gamma=gamma,terminal.points=terminal.points,ct=ct)

envelope<-sapply(Time,envelope.function)

if(min(envelope)<=0){
warning("Envelope function must be strictly positive")
envelope<-envelope+(1-min(envelope))
#stop("Envelope function must be strictly positive")
}

##get constant for rejection sampling
e<-max(simulated.intensity/envelope)

for (i in 1:nruns){
	if (i%%10 == 0){
		cat("Run",i,"\n")
	}else{
      cat("Run",i,"\t")
	}
      x<-numeric(1000*(t.end-t.start))
	##initialize x to avoid memory issue
	##maximum of average 1000 spikes per second - can change this if desired
	spike.counter<-0
	t0<-t.start
      while(t0<t.end){
      	u1<-runif(1)
		lambda.u<-envelope.function(t0)*e
		t0 <- t0-(1/lambda.u)*log(u1)
      	if(t0>t.end) break
      	u2<-runif(1)
		lambda.t <- spike.fun(t0,f,w0,eta,gamma,terminal.points,ct)
      	if(u2<=(lambda.t/lambda.u)){
			spike.counter <- spike.counter+1
			x[spike.counter] <- t0
		}
      } ##end while loop
	x<-x[1:spike.counter] ##only nonzero numbers included
	simulated.spikes[[i]] <- unique(round(x,3))
} ##end for loop
cat("\n") ##for aesthetic purposes
return(simulated.spikes)
}