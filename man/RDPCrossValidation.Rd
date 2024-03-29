% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RDPCrossValidation.R
\name{RDPCrossValidation}
\alias{RDPCrossValidation}
\title{Recursive dyadic partitioning cross-validation}
\usage{
RDPCrossValidation(
  spikes,
  t.start = 0,
  t.end,
  poss.lambda = seq(0, 10, by = 0.1),
  max.J = 7,
  PSTH = FALSE,
  max.diff = 0.005,
  pct.diff.plot = TRUE,
  print.J.value = TRUE
)
}
\arguments{
\item{spikes}{a list of spike trains.}

\item{t.start}{the starting time of the recording window; the default value is 0.}

\item{t.end}{the ending time of the recording window.}

\item{poss.lambda}{a numeric vector containing a grid of penalty values.}

\item{max.J}{the maximum resolution of the dyadic partitioning used the estimate the piecewise constant intensity function \eqn{c(t)}.}

\item{PSTH}{if TRUE, performs leave-one-train-out cross-validation for the c(t) estimate based on PSTH data.
If FALSE, performs leave-one-spike-out cross-validation for the c(t) estimate from each individual train.}

\item{max.diff}{the maximum allowance for the integrated squared error (ISE) of a smaller model to deviate from the overall minimum ISE.}

\item{pct.diff.plot}{a logical value indicating whether to produce a plot of the percentage difference (above minimum ISE) vs. J.}

\item{print.J.value}{a logical value indicating whether to print off the J value at each step of the cross-validation or not.}
}
\value{
A list of length 3 is returned returned.
The first item in the list is the optimal partition depth as computed by ISE (\eqn{\lambda}).
The second item in the list is the optimal penalty term as corresponding to that partition depth (J).
The third item in the list is a matrix containing the ISE values for all combinations of partition depth and penalty term.
}
\description{
Optimizes a penalized log-likelihood to find the optimal number of partitions for the recursive dyadic partitioning.
}
