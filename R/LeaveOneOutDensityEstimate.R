#' Leave-one-out density function cross-validation
#'
#' Calculates the average of m-1 density function estimates for a given set of time points with the ith density estimate removed.
#' This function implements equation (1) of the supplementary material in Ramezan *et al.* (2014).
#'
#' @param til the time of the lth spike in trial i of a set of spike trains.
#' @param i the index of the trial whose estimated density function is to be removed (in leave-one-out cross-validation).
#' @param points a numerical vector containing the end points of the dyadic partitioning for a given resolution.
#' @param f.hat.minus.i a numeric matrix where the columns represent the m density estimates and the rows represent the time points at which the density estimates are calculated.
#'
#' @return A numeric vector of the average of m-1 density function estimates for a given set of time points with the ith density estimate removed.
#'
#' @references Ramezan, R., Marriott, P., and Chenouri, S. (2014), *Statistics in Medicine*, **33**(2), 238-256. doi: 10.1002/sim.5923.
#'
#' @noRd


LeaveOneOutDensityEstimate <- function(til,i,points,f.hat.minus.i){

if (til < min(points) | til > max(points)) return(0)
if (til==min(points)) return(f.hat.minus.i[i,1])
if (til==max(points)) return(f.hat.minus.i[i,dim(f.hat.minus.i)[2]])
f.minus.i.t.index<-min(which(points>til))-1
return(f.hat.minus.i[i,f.minus.i.t.index])
}
