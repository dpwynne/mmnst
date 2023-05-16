#' Goodness of fit pp-plot
#'
#' Runs a goodness-of-fit test for a particular model, i.e., a given intensity function estimate.
#'
#' @importFrom rlang .data
#'
#' @param spikes a list of spike trains.
#' @param theta a numeric matrix in which the *j*th column contains the intensity function estimate for the *j* spike train in the list.
#' @param t.start the starting time of the recording window; the default value is 0.
#' @param t.end the ending time of the recording window. The default value is 10, corresponding to a 10-second recording.
#' @param neuron.name a string containing the name of the neuron being plotted, used only to title the plot.
#' @param resolution a scalar determining the bin width. This is equivalent to Delta in Haslinger *et al.* (2010).
#' @param axis.label.size a scalar determining the font size of the x- and y-axis labels of the plot.
#' @param title.size a scalar determining the font size of the title of the plot (if one exists).
#'
#' @return a list of data frames. Each data frame includes the model and empirical quantiles for one spike train.
#' This data is used to generate the pp-plot outlined in Haslinger *et al.* (2010).
#'
#' @references Haslinger, R., Pipa G., and Brown, E. (2010). Discrete time rescaling theorem: determining goodness of fit for discrete time statistical models of neural spiking. *Neural Computation*. **22**(10):2477-506. doi: 10.1162/NECO_a_00015.
#'
#' @export

GOFPlot <-function(spikes, theta, t.start = 0, t.end = 10,
                   neuron.name = NULL, resolution= (t.end - t.start)/(nrow(theta)-1),
                   axis.label.size = 18, title.size = 24){
  ##spikes is the list of spike trains you want to check goodness of fit for
  ##theta is a matrix in which the rows represent Time and the columns represent spike trains
  ##neuron.name is the name of the neuron being plotted, used only to title the plot if you want it
  ##resolution corresponds to bin width, Delta, in Haslinger (2010)
  ##our default value corresponds to the intensity function length in estimate.theta.t


#library(ggplot2)

  endpoints.bins<-seq(t.start, t.end, by=resolution)

  ##This error should never occur with the default resolution but might with user-defined resolution
  if (length(endpoints.bins)> nrow(theta)){
    stop("The number of bins for GOF plot cannot exceed that of the intensity function; plotting stopped.", call. = FALSE)
  }
  if(length(spikes) == 0){
    stop("No spike trains provided; cannot create a GOF plot.", call. = FALSE)
  }

  spike.gof<-vector("list",length(spikes))
  names(spike.gof)<-paste0("Trial",seq(1,length(spikes)))

  # theta is now the list of individual theta(t) estimates
  for(i in 1:length(spikes)){
    qk <- HaslingerQ(endpoints.bins,theta[,i],resolution)
    transformed.y <- HaslingerYi(spikes[[i]],endpoints.bins,qk)
    transformed.y<-transformed.y[!is.na(transformed.y) & !is.nan(transformed.y)]
    n <- length(transformed.y)
    if (n > 0){
      b <- ((1:n)-0.5)/n ##midpoints of each bin
    }else{
      cat("No numerical arguments found, plotting stopped\n")
      return()
    }
    spike.gof[[i]]<-data.frame(model.quantiles=sort(transformed.y), empirical.quantiles=b)
  }

  n <- max(unlist(lapply(spikes,length)))


  ##setup the plot
  plot.basics<-ggplot2::ggplot(spike.gof[[1]],ggplot2::aes(x = .data$model.quantiles,
                                                           y = .data$empirical.quantiles))+
    ggplot2::theme_classic()+
    ggplot2::theme(panel.grid.major=ggplot2::element_blank(),
                   panel.grid.minor=ggplot2::element_blank(),
                   plot.margin=ggplot2::unit(c(5.08,5.08,5.08,5.08),"mm"),
                   axis.text=ggplot2::element_text(size=axis.label.size),
                   axis.line.x=ggplot2::element_line(size=0.5),
                   axis.line.y=ggplot2::element_line(size=0.5))+
    ggplot2::scale_x_continuous(limits=c(0,1))+
    ggplot2::scale_y_continuous(limits=c(0,1))

  #it looks like it automatically does xlim(0,1) and ylim(0,1) but include just to be sure
  plot.labeled<-plot.basics+ggplot2::labs(x="Empirical CDF",y="Uniform CDF")+ggplot2::theme(axis.title=ggplot2::element_text(size=axis.label.size))

  plot.cdfs<-plot.labeled


  for (i in 1:length(spike.gof)){
    plot.cdfs<-plot.cdfs+ggplot2::geom_line(data=spike.gof[[i]],ggplot2::aes(x=.data$model.quantiles,y=.data$empirical.quantiles),color="gray",linetype=1)
  }

  alpha.corrected <- 0.05/length(spikes)  # bonferroni correction
  ks.critical.value <- NSM3::qKolSmirnLSA(alpha.corrected)

  plot.bands<-plot.cdfs+ggplot2::geom_abline(color="black",linetype=2,size=1,slope=1,intercept=ks.critical.value/sqrt(n))+ggplot2::geom_abline(color="black",linetype=2,size=1,slope=1,intercept=-ks.critical.value/sqrt(n))+ggplot2::geom_abline(color="black",linetype=1,size=1,slope=1,intercept=0)

  if(!is.null(neuron.name)){
    plot.bands <- plot.bands + ggplot2::ggtitle(neuron.name) + ggplot2::theme(plot.title=ggplot2::element_text(size=title.size,face="bold", hjust=0.5))
  }

  print(plot.bands)

  invisible(spike.gof)
}
