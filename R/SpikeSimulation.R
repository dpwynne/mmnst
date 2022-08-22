#' Spike train simulation
#'
#' Simulates spike trains from an inhomogeneous Poisson process with a given time-varying intensity function of one of the forms discussed in Ramezan *et al.* (2014).
#'
#' @importFrom stats runif
#'
#' @param nruns a scalar determining the number of spike trains (trials) to be simulated.
#' @param t.start the starting time of the simulated data collection.
#' @param t.end the ending time of the simulated data collection.
#' @param resolution the time resolution of the spike train. The default is 0.001 (1 ms).
#' @param intensity.function any intensity function for the corresponding process from which you want to simulate data.
#' @param pass.arg a list of arguments to pass to the intensity function. In the additive and multiplicative models, this contains four numeric vectors;
#' the first vector determines the frequencies, the second vector determines initial phases, the third vector determines the \eqn{\eta} values and the fourth vector determines the \eqn{\gamma} values.
#' @param envelope.function The envelope function used for acceptance-rejection sampling to simulate the data.
#'
#' @return A list of numeric vectors, each of which contains a simulated spike train.
#'
#' @references Ramezan, R., Marriott, P., and Chenouri, S. (2014), *Statistics in Medicine*, **33**(2), 238-256. doi: 10.1002/sim.5923.
#'
#' @export

SpikeSimulation<-function(nruns, t.start = 0, t.end,
                           resolution = 0.001,
                           intensity.function,
                           pass.arg,
                           envelope.function=function(t) return(1)){
##returns a list of vectors containing simulated spike times
##nruns is number of simulated runs desired
  ## intensity.function is any intensity function for the corresponding process from which you want to simulate data
##resolution: 0.001 = 1 ms assuming time in seconds
  ##pass.arg is a list of arguments to pass to intensity.function
##ct is the simulated piecewise constant intensity function
##envelope.function is the envelope function


  #### BIG NOTE: FIGURE OUT HOW TO USE ... ARGUMENTS TO PASS DIRECTLY TO intensity.function
  # We want to fix this to work with arbitrary intensity function with arbitrary parameters
  ## Use do.call with an argslist

simulated.spikes<-vector("list",length=nruns)

##unpack arguments from pass.arg
f<-pass.arg[[1]]
w0<-pass.arg[[2]]
eta<-pass.arg[[3]]
gamma<-pass.arg[[4]]
ct <- pass.arg[[5]]
# figure out how to use ... notation; then we can pass these arguments directly as f, w0, eta, gamma, ct instead of pass.arg

terminal.points <- IdentifyTerminalPoints(t.start,t.end,log(length(ct),2))
##get terminal points from ct


Time.vector <- seq(t.start, t.end, by = resolution)  # by default 1 ms assuming time in seconds

simulated.intensity<-sapply(Time.vector,intensity.function,f=f,w0=w0,eta=eta,gamma=gamma,terminal.points=terminal.points,ct=ct)

envelope<-sapply(Time.vector,envelope.function)

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
		lambda.t <- intensity.function(t0,f,w0,eta,gamma,terminal.points,ct)
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
