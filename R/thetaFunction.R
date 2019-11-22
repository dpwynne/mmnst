#' \theta(t)
#'
#' Calculate the value of theta(t) at an individual time point t
#'
#' @param t the time at which to calculate \theta(t). The function \code{q.vector} integrates this function, so t can be of any length (including 1).
#' @param theta a numeric vector containing the estimated intensity function \theta(t).
#' @param t.start the starting time of the recording
#' @param resolution the time difference between successive points in the theta vector
#'
#' @return the value of the intensity function \theta(t) at a specific time point
#'
#' @export

theta.function<-function(t,theta,t.start = 0,resolution) theta[(1+floor((t-t.start)/resolution))]
##pull the value of the intensity function that corresponds to a particular time value
