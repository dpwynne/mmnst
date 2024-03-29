% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RasterPlot.R
\name{RasterPlot}
\alias{RasterPlot}
\title{Raster plot}
\usage{
RasterPlot(
  spike.train,
  time.highlight = c(),
  trial.highlight = c(),
  graph.title = NULL
)
}
\arguments{
\item{spike.train}{a list of numeric vectors, each of which contains the spike times from a single trial of the experiment.}

\item{time.highlight}{a numeric vector indicating times (if any) to highlight with dashed vertical lines; for example, the onset/offset times of a stimulus.}

\item{trial.highlight}{a numeric vector indicating trials (if any) to highlight by showing them in a different color.}

\item{graph.title}{the title of the plot}
}
\value{
A ggplot object containing the parameters of the Raster plot.
}
\description{
Produces a raster plot for a single neuron.
}
