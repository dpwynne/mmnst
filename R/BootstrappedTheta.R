#' Bootstrapped Theta
#'
#' Computes a 95\% pointwise confidence interval for \eqn{\theta(t)} using the bootstrap percentile method.
#'
#' @importFrom stats quantile
#'
#' @param theta.t a matrix containing estimates of \eqn{\theta(t)}. Each row represents a time point and each column represents a spike train.
#' Thus, each value in the matrix represents the value of \eqn{\theta(t)} at time point t estimated from a particular spike train.
#' @param B the number of desired bootstrap resamples.
#'
#' @return A matrix containing the lower and upper bounds of the percentile confidence interval at each time point.
#'
#' @noRd

BootstrappedTheta<-function(theta.t,B=1000){

ntrains<-dim(theta.t)[2]
bootstrap.matrix<-matrix(0,dim(theta.t)[1],B)

if(ntrains < 10){
  warning("Fewer than 10 spike trains provided; bootstrap estimates are unreliable.",
          call. = FALSE, immediate. = FALSE)
}

for (i in 1:B){
	sampled.intfns<-sample(1:ntrains,ntrains,replace=TRUE)
	new.sample.matrix<-theta.t[,sampled.intfns, drop = FALSE] # one-train matrix fix
	new.mean<-apply(new.sample.matrix,1,mean)
	bootstrap.matrix[,i]<-new.mean
}##end for loop

bootstrapped.limits<-apply(bootstrap.matrix, 1, quantile, probs=c(0.025,0.975))
if (dim(bootstrapped.limits)[2]>=dim(bootstrapped.limits)[1]) bootstrapped.limits<-t(bootstrapped.limits)
return(bootstrapped.limits)
}
