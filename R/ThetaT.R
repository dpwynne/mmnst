#' Intensity function \eqn{\theta(t)}
#'
#' Calculates the value of \eqn{\theta(t)} at an individual time point \eqn{t}.
#'
#' @param t the time at which to calculate \eqn{\theta(t)}. The function [HaslingerQ()] integrates this function, so t can be of any length (including 1).
#' @param theta a numeric vector containing the estimated intensity function \eqn{\theta(t)}.
#' @param t.start the starting time of the recording.
#' @param resolution the time difference between successive points in the theta vector.
#'
#' @return The value of the intensity function \eqn{\theta(t)} at (a) specific time point(s).
#'
#' @noRd

ThetaT <- function(t,theta,t.start = 0,resolution){
  ##pull the value of the intensity function that corresponds to a particular time value
  theta[(1+floor((t-t.start)/resolution))]
}

