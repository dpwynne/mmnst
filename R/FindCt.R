#' Find \eqn{c(t)}
#'
#' Compute the piecewise constant estimate of the intensity function, \eqn{c(t)}, for each spike train
#'
#' @param spikes a list of spike trains.
#' @param t.start the starting time of the recording window; the default value is 0.
#' @param t.end the ending time of the recording window.
#' @param lambda a penalty term used in estimating the piecewise constant intensity function c(t). Larger values of lambda result in "smoother" estimates.
#' @param J the maximum size of the tree in the initial dyadic partitioning used to estimate c(t). The final estimate of c(t) will have 2^J values.
#'
#' @return A list of length 2 is returned.
#' The first item in the list is the estimate of \eqn{c(t)} created by averaging the estimates for each spike train.
#' The second item is a matrix in which each row represents the estimate of c(t) for an individual spike train.
#'
#' @export
FindCt<-function(spikes, t.start, t.end, lambda, J){

T.data<-t.end-t.start
val <- floor(2^J)
by.terminal<-T.data/val

terminal.points <- seq(t.start,t.end,by.terminal)
ct.best<-matrix(NA,nrow=length(spikes),ncol=val)

for (i in 1:length(spikes)){
	xi <- spikes[[i]]
	count.points<-numeric(val)
	for (ii in 1:val){
		count.points[ii]<-length(xi[xi>=terminal.points[ii] & xi<terminal.points[ii+1]])
	}
      if (J == 0){
        ct.best[i,] <- count.points
      } else {
        ct.best[i,]<-PoissonRDP(count.points,lambda)
      }
}

##output of PoissonRDP is not exactly c(t), needs to be scaled by multiplying by val/T.data
ct.best<-ct.best/(by.terminal)

## We have 1 c(t) estimate per spike train. This means that the process of merging partitions may not be the same on every spike train.
## To get a single "best" estimate of c(t), we average the c(t) estimates at each time point.
## Because the partitions may merge at different points for different spike trains, averaging the spike trains will "destroy" the apparent merge.
ct.avg<-apply(ct.best,2,mean)
return(list(ct.avg = ct.avg, ct.best = ct.best))
}
