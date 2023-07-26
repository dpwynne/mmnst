#' Find \eqn{c(t)}
#'
#' Compute the piecewise constant estimate of the intensity function, \eqn{c(t)}, for each spike train
#'
#' @param spikes a list of spike trains.
#' @param t.start the starting time of the recording window; the default value is 0.
#' @param t.end the ending time of the recording window.
#' @param lambda a penalty term used in estimating the piecewise constant intensity function c(t). Larger values of lambda result in "smoother" estimates.
#' @param J the maximum size of the tree in the initial dyadic partitioning used to estimate c(t). The final estimate of c(t) will have 2^J values.
#' @param PSTH if TRUE, will aggregate the spikes across all trials in each and every bin and then estimate c(t) using the post-stimulus-time histogram.
#' If FALSE, will estimate a c(t) for each trial and average the c(t) estimates.
#'
#' @return A list of length 2 is returned.
#' If PSTH=TRUE, the first item in the list (`ct.avg`) is the single estimate of c(t) based on the PSTH (scaled appropriately so that it integrates to the average number of spikes per trial),
#' and the second item in the list (`ct.best`) is a matrix, containing one row per spike train / trial, whose rows each replicate `ct.avg`.
#' If PSTH=FALSE, each row of `ct.best` contains the \eqn{c(t)} estimates for the individual spike trains / trials, and
#' `ct.avg` contains the single estimate of \eqn{c(t)} obtained by averaging the `ct.best` estimates.
#'
#' @export
FindCt<-function(spikes, t.start, t.end, lambda, J, PSTH = FALSE){

T.data<-t.end-t.start
val <- floor(2^J)
by.terminal<-T.data/val

terminal.points <- seq(t.start,t.end,by.terminal)

if (PSTH){

  xi <- sort(unlist(spikes)) # this is a single vector containing all spike times across all trials
  count.points<-numeric(val)

  for (ii in 1:val){
    count.points[ii]<-length(xi[xi>=terminal.points[ii] & xi<terminal.points[ii+1]])
  }

  if (J == 0){
    ct.best.PSTH <- count.points
  } else {
    ct.best.PSTH <-PoissonRDP(count.points,lambda)
  }

  ct.best <- matrix(ct.best.PSTH, nrow = length(spikes), ncol = val, byrow = TRUE)/length(spikes)

} else {

  ct.best<-matrix(NA,nrow=length(spikes),ncol=val)

  if(length(spikes) == 0){
    stop("No spike trains provided; cannot estimate c(t).")
  }

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
}

# This is the new version. It counts the aggregate total of the number of spikes within the
# ct.best is now a matrix that contains the aggregate best ct estimate in each row

##output of PoissonRDP is not exactly c(t), needs to be scaled by multiplying by val/T.data
ct.best<-ct.best/(by.terminal)

## We have 1 c(t) estimate per spike train. This means that the process of merging partitions may not be the same on every spike train.
## To get a single "best" estimate of c(t), we average the c(t) estimates at each time point.
## Because the partitions may merge at different points for different spike trains, averaging the spike trains will "destroy" the apparent merge.
ct.avg<-apply(ct.best,2,mean)


if (length(spikes) == 1){
  warning("c(t) was estimated using only one spike train.", call. = FALSE)
}

return(list(ct.avg = ct.avg, ct.best = ct.best))
}
