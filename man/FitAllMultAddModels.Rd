% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FitAllMultAddModels.R
\name{FitAllMultAddModels}
\alias{FitAllMultAddModels}
\title{Fit all multiplicative and additive models}
\usage{
FitAllMultAddModels(
  K,
  spikes,
  f.common.table,
  setup.pars,
  terminal.points,
  ct,
  user.defined.frequencies = NULL
)
}
\arguments{
\item{K}{the number of frequency components in the largest model to fit.}

\item{spikes}{a list of spike trains.}

\item{f.common.table}{a table whose names contain the high-amplitude frequency components as computed by \code{\link[=FindTopFrequencies]{FindTopFrequencies()}}.}

\item{setup.pars}{a list of additional parameters for the likelihood function, computed by \code{\link[=SetupLikelihoods]{SetupLikelihoods()}}.}

\item{terminal.points}{a numeric vector containing the time points at which \eqn{c(t)} changes.}

\item{ct}{a numeric vector containing the estimated piecewise constant intensity function \eqn{c(t)}. The length of \eqn{c(t)} should be a whole number power of 2.}

\item{user.defined.frequencies}{(= NULL) a numeric vector of length at least K containing the frequencies (in decreasing order of importance) to try in the model.
By default this is NULL, and frequencies will be directly read from f.common.table.}
}
\value{
A list of length 3 is returned.
The first item in the list is a list of frequency estimates for each model.
The second item in the list is a list of phase estimates for each model.
The third item in the list is a list of eta/gamma estimates and fit criteria for each model.
}
\description{
Wrapper function to fit all 2K+1 models (multiplicative, additive, and nonperiodic).
}
