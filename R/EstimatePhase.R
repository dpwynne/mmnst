#' Estimate initial phases of the periodic function \eqn{\nu}
#'
#' Computes the initial phases for spike trains based on the given frequencies.
#'
#' @param spikes a list of spike trains.
#' @param f.hat a numeric vector containing frequency estimates for a particular model.
#'
#' @return A numeric matrix containing the estimates of the initial phases for each frequency and spike train. The rows represent individual spike trains and the columns are individual frequency values.
#'
#' @export

EstimatePhase <- function(spikes,f.hat){
##spikes is the list of spike trains
##f.hat is the identified best frequencies

cat("Estimating Phase\n")

w.hat<-2*pi*f.hat
# w.hat is the angular frequency

ntrials<-length(spikes)

w0.hat <- matrix(NA,ntrials,length(w.hat))
# w0.hat is a matrix where 1 row = 1 spike train and 1 column = 1 frequency
# w0.hat contains the phase offsets for each trial for each frequency

for(itr in 1:ntrials){
	w0.hat.itr <- numeric(length(w.hat))
	spike.train <- spikes[[itr]]  # this goes outside the second for loop
	  for(i in 1:length(w.hat)){
	    sum.sin<- sum(sin(w.hat[i]*spike.train))
	    sum.cos<-sum(cos(w.hat[i]*spike.train))
	    phi.hat<-atan(-sum.sin/sum.cos)
	    if( sin(phi.hat)*sum.sin > cos(phi.hat)*sum.cos ) phi.hat<-phi.hat+pi
	    w0.hat.itr[i] <- phi.hat/(2*pi)
    }
	w0.hat[itr,] <- w0.hat.itr
}

cat("Phase Estimated\n")
return(w0.hat)
}

