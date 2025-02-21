#' Scree plot-like graph of most common frequencies
#'
#' Produce a scree plot-like graph to determine number of frequency components
#'
#' @importFrom rlang .data
#'
#' @param freq.table a named numeric vector in which the names are the frequencies and the values are the counts of spike trains
#' in which those frequencies are among the highest peaks in the periodogram. This vector is typically
#' output by [FindTopFrequencies()].
#' @param spikes (= NULL) a list of spike trains. Can be omitted if n is given.
#' @param n (= length(spikes)) the number of spike trains analyzed by [FindTopFrequencies()]
#' @param frac a vector of values where dashed horizontal lines will be added to the plot. As a single value, this can represent the proportion of spike trains that must exhibit a frequency to be considered significant.
#' @return the `freq.table` argument (invisibly)
#'
#' @export

FrequencyScreePlot <- function(freq.table, spikes = NULL, n = length(spikes), frac = 0.5){
  if(n <= 0){
    stop("Number of spike trains to analyze is 0.")
  }
  if(length(freq.table) <= 0){
    stop("Number of frequencies to analyze is 0.")
  }

  #f.table.sorted <- rev(sort(freq.table))
  f.table.sorted <- sort(freq.table, decreasing = TRUE)
  # decreasing = TRUE keeps same order if the table is already sorted in order
  # rev flips the order of tied frequencies

  freq.df <- data.frame(f.table.sorted/n,
                        seq(1,length(f.table.sorted)))
  names(freq.df) <- c("f", "p", "k")

  plot_basics<-ggplot2::ggplot(freq.df,ggplot2::aes(x = .data$k, y = .data$p))+
   #plot_basics<-ggplot2::ggplot(freq.df,ggplot2::aes(x = k, y = p))+
    ggplot2::theme_classic()+
    ggplot2::theme(panel.grid.major=ggplot2::element_blank(),
                   panel.grid.minor=ggplot2::element_blank(),
                   plot.margin=ggplot2::unit(c(5.08,5.08,5.08,5.08),"mm"),
                   #axis.text=ggplot2::element_text(size=axis.label.size),
                   axis.line.x=ggplot2::element_line(size=0.5),
                   axis.line.y=ggplot2::element_line(size=0.5))+
    ggplot2::scale_x_continuous(limits=c(1,min(10, length(f.table.sorted))) , breaks = 1:10)+ # plot at most the top 10 frequencies
    ggplot2::scale_y_continuous(limits=c(0,1.05))

  plot_labeled<-plot_basics +
    ggplot2::labs(x = "Number of Frequencies", y = "Proportion of Spike Trains") #+
    #ggplot2::theme(axis.title=ggplot2::element_text(size=axis.label.size))

  plot_freq <- suppressWarnings( # suppress warnings because we are cutting off after 10 points
    plot_labeled + ggplot2::geom_point() + ggplot2::geom_line() +
    ggplot2::geom_text(ggplot2::aes(label = paste(.data$f,"")), nudge_y = 0.05) +
    #ggplot2::geom_text(aes(label = paste(as.character(f),"")), nudge_y = 0.05) +
  ggplot2::geom_hline(yintercept = frac, linetype = "dotted")
  )

  suppressWarnings(print(plot_freq))

  invisible(freq.table)
}

