#' Identify terminal points
#'
#' Computes a vector of times at which the estimated piecewise constant intensity function \eqn{c(t)} changes.
#'
#' @param t.start the starting time of the data collection.
#' @param t.end the ending time of the data collection.
#' @param J the maximum size of the tree in the initial dyadic partitioning used to estimate \eqn{c(t)}.
#'
#' @return A numeric vector of points at which \eqn{c(t)} changes.
#'
#' @export

IdentifyTerminalPoints<-function(t.start,t.end,J){
##t.start and t.end are the starting and ending times of the data collection
##J is identified from N*
terminal.points<-seq(t.start,t.end,length=(2^J+1)) ##vector of terminal points
}
