#' Raster plot
#'
#' Produces a raster plot for a single neuron.
#'
#' @import ggplot2
#' @importFrom rlang .data
#'
#' @param NeuronNumber the ID for the neuron to be plotted (only used in the title of the raster plot).
#' @param spike.train a list of numeric vectors, each of which contains the spike times from a single trial of the experiment.
#'
#' @return A ggplot object containing the parameters of the Raster plot.
#'
#' @export

RasterPlot <- function(NeuronNumber,spike.train){
##new Raster function for Raster plots in ggplot2

spike.times<-unlist(spike.train)
spike.trials<-sapply(spike.train,length)
spike.numbers<-rep(seq(1,length(spike.trials)),spike.trials)

spike.data<-data.frame(Trial=spike.numbers,Time=spike.times)

tick.separator.x<-(ceiling(max(spike.data$Time))-floor(min(spike.data$Time)))/5
y.high<-50*ceiling(max(spike.data$Trial/50))
y.low<-50*floor(min(spike.data$Trial/50))
y.breaks<-seq(y.low,y.high,length=6)
y.breaks<-y.breaks[which(y.breaks<=max(spike.data$Trial))]

plot.basics<-ggplot(spike.data,aes(x=.data$Time,y=.data$Trial))+
  theme_classic()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),plot.margin=unit(c(5.08,5.08,5.08,5.08),"mm"),axis.text=element_text(size=24),axis.line.x=element_line(size=0.5),axis.line.y=element_line(size=0.5))

plot.ticks<-plot.basics+scale_x_continuous(breaks=seq(floor(min(spike.data$Time)),ceiling(max(spike.data$Time)),by=tick.separator.x))+scale_y_continuous(breaks=y.breaks)

plot.raster<-plot.ticks+geom_point(aes(x=.data$Time,y=.data$Trial),size=0.2)

plot.labeled<-plot.raster+labs(x="Time (sec.)",y="Trial Number")+theme(axis.title=element_text(size=24))

title.graph <- paste0("Neuron: ",as.character(NeuronNumber))
plot.titled<-plot.labeled+ggtitle(title.graph)+theme(plot.title=element_text(size=18,face="bold"))

return(plot.titled)
}
