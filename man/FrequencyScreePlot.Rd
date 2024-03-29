% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FrequencyScreePlot.R
\name{FrequencyScreePlot}
\alias{FrequencyScreePlot}
\title{Scree plot-like graph of most common frequencies}
\usage{
FrequencyScreePlot(freq.table, spikes = NULL, n = length(spikes))
}
\arguments{
\item{freq.table}{a named numeric vector in which the names are the frequencies and the values are the counts of spike trains
in which those frequencies are among the highest peaks in the periodogram. This vector is typically
output by \code{\link[=FindTopFrequencies]{FindTopFrequencies()}}.}

\item{spikes}{(= NULL) a list of spike trains. Can be omitted if n is given.}

\item{n}{(= length(spikes)) the number of spike trains analyzed by \code{\link[=FindTopFrequencies]{FindTopFrequencies()}}}
}
\value{
the \code{freq.table} argument (invisibly)
}
\description{
Produce a scree plot-like graph to determine number of frequency components
}
