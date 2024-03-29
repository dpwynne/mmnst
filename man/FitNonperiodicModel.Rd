% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FitNonperiodicModel.R
\name{FitNonperiodicModel}
\alias{FitNonperiodicModel}
\title{Fit nonperiodic model}
\usage{
FitNonperiodicModel(spikes, setup.pars, terminal.points, ct)
}
\arguments{
\item{spikes}{a list of spike trains.}

\item{setup.pars}{a list of additional parameters for the likelihood function, computed by the \code{\link[=SetupLikelihoods]{SetupLikelihoods()}} function.}

\item{terminal.points}{a numeric vector containing the time points at which c(t) changes.}

\item{ct}{a numeric vector containing the estimated piecewise constant intensity function c(t). The length of ct should be a whole number power of 2.}
}
\value{
A list of length 3 is returned.
The first item in the list is a matrix whose rows each contain the MLEs of eta for a single spike train (this should be 0).
The second item in the list is a matrix whose rows each contain the MLEs of gamma for a single spike train (this should also be 0).
The third item in the list is a matrix whose rows each contain the AIC, AICc, BIC, and log-likelihood for the model for a single spike train.
}
\description{
Fits the nonperiodic multiscale model.
}
\references{
Ramezan, R., Marriott, P., and Chenouri, S. (2014), \emph{Statistics in Medicine}, \strong{33}(2), 238-256. doi: 10.1002/sim.5923.
}
