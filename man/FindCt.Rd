% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FindCt.R
\name{FindCt}
\alias{FindCt}
\title{Find \eqn{c(t)}}
\usage{
FindCt(spikes, t.start, t.end, lambda, J, PSTH = FALSE)
}
\arguments{
\item{spikes}{a list of spike trains.}

\item{t.start}{the starting time of the recording window; the default value is 0.}

\item{t.end}{the ending time of the recording window.}

\item{lambda}{a penalty term used in estimating the piecewise constant intensity function c(t). Larger values of lambda result in "smoother" estimates.}

\item{J}{the maximum size of the tree in the initial dyadic partitioning used to estimate c(t). The final estimate of c(t) will have 2^J values.}

\item{PSTH}{if TRUE, will aggregate the spikes across all trials in each and every bin and then estimate c(t) using the post-stimulus-time histogram.
If FALSE, will estimate a c(t) for each trial and average the c(t) estimates.}
}
\value{
A list of length 2 is returned.
If PSTH=TRUE, the first item in the list (\code{ct.avg}) is the single estimate of c(t) based on the PSTH (scaled appropriately so that it integrates to the average number of spikes per trial),
and the second item in the list (\code{ct.best}) is a matrix, containing one row per spike train / trial, whose rows each replicate \code{ct.avg}.
If PSTH=FALSE, each row of \code{ct.best} contains the \eqn{c(t)} estimates for the individual spike trains / trials, and
\code{ct.avg} contains the single estimate of \eqn{c(t)} obtained by averaging the \code{ct.best} estimates.
}
\description{
Compute the piecewise constant estimate of the intensity function, \eqn{c(t)}, for each spike train
}
