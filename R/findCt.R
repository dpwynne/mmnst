#' Find c(t)
#'
#' Compute the piecewise constant estimate of the intensity function, c(t), for each spike train
#'
#' @param spikes a list of spike trains
#' @param Time a numeric vector containing, at minimum, the start and end points of the spike train recording
#' @param lambda a penalty term used in estimating the piecewise constant intensity function c(t). Larger values of lambda result in "smoother" estimates.
#' @param J the maximum size of the tree in the initial dyadic partitioning used to estimate c(t). The final estimate of c(t) will have 2^J values.
#'
#' @return a list of length 2
#' The first item in the list is the estimate of c(t) created by averaging the estimates for each spike train
#' The second item is a matrix in which each row represents the estimate of c(t) for an individual spike train
#'
#' @export
find.ct<-function(spikes,Time,lambda,J){

time.start<-min(Time)
time.end<-max(Time)
T.data<-time.end-time.start
val <- floor(2^J)
by.terminal<-T.data/val

terminal.points <- seq(time.start,time.end,by.terminal)
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
        ct.best[i,]<-Poisson.RDP(count.points,lambda)
      }
}

##output of Poisson.RDP is not exactly c(t), needs to be scaled by multiplying by val/T.data
ct.best<-ct.best/(by.terminal)

## We have 1 c(t) estimate per spike train. This means that the process of merging partitions may not be the same on every spike train.
## To get a single "best" estimate of c(t), we average the c(t) estimates at each time point.
## Because the partitions may merge at different points for different spike trains, averaging the spike trains will "destroy" the apparent merge.
ct.avg<-apply(ct.best,2,mean)
return(list(ct.avg = ct.avg, ct.best = ct.best))
}
