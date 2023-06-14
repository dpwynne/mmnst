#' Raster plot
#'
#' Produces a raster plot for a single neuron.
#'
#' @import ggplot2
#' @importFrom rlang .data
#'
#' @param spike.train a list of numeric vectors, each of which contains the spike times from a single trial of the experiment.
#' @param time.highlight a numeric vector indicating times (if any) to highlight with dashed vertical lines; for example, the onset/offset times of a stimulus.
#' @param trial.highlight a numeric vector indicating trials (if any) to highlight by showing them in a different color.
#' @param graph.title the title of the plot
#'
#' @return A ggplot object containing the parameters of the Raster plot.
#'
#' @export

RasterPlot <- function(spike.train, time.highlight = c(), trial.highlight = c(), graph.title = NULL){
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

trial.colors <- data.frame(spike.data, Color = rep("plain", nrow(spike.data)))
trial.colors$Color[trial.colors$Trial %in% trial.highlight] <- "highlight"

plot.basics<-ggplot(trial.colors,aes(x=.data$Time,y=.data$Trial)) +
  theme_classic() +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),plot.margin=unit(c(5.08,5.08,5.08,5.08),"mm"),axis.text=element_text(size=24),axis.line.x=element_line(size=0.5),axis.line.y=element_line(size=0.5))

plot.ticks<-plot.basics +
  scale_x_continuous(breaks = seq(floor(min(spike.data$Time)), ceiling(max(spike.data$Time)), by=tick.separator.x)) +
  scale_y_continuous(breaks=y.breaks)

plot.raster<-plot.ticks +
  geom_point(aes(x = .data$Time, y = .data$Trial, color = .data$Color), size=0.2) +
  scale_color_manual(values = c(highlight = "red", plain = "black")) +
  theme(legend.position = "none")

plot.labeled<-plot.raster +
  labs(x="Time (sec.)", y="Trial Number") +
  theme(axis.title=element_text(size=24))

plot.highlighted <- plot.labeled +
  geom_vline(xintercept = time.highlight, size = 1, linetype = "dashed", color = "blue")

plot.titled<-plot.highlighted +
  ggtitle(as.character(graph.title)) +
  theme(plot.title=element_text(size=18,face="bold"))

return(plot.titled)
}
