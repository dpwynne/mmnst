#' Inter-spike interval plot
#'
#' Produces an inter-spike interval (ISI) plot for a single trial.
#'
#' @param NeuronNumber an identifier for the neuron being plotted.
#' @param TrialNumber the trial number whose spike train is being plotted.
#' @param spike.train a list of spike trains.
#'
#' @import ggplot2
#' @importFrom rlang .data
#'
#' @return A ggplot object containing the parameters of the ISI plot.
#'
#' @references Ramezan, R., Marriott, P., and Chenouri, S. (2014), *Statistics in Medicine*, **33**(2), 238-256. doi: 10.1002/sim.5923.
#'
#' @export

ISIPlot <- function(NeuronNumber, TrialNumber, spike.train){

x <- spike.train[[TrialNumber]]
log.diff.x<-log(diff(x)) #base e log, we're just looking to see when points line up
time.fire<-x[2:length(x)]

ISI.data <- data.frame(Time = time.fire, Log.Diff = log.diff.x)

plot.basics<-ggplot(ISI.data,aes(x=.data$Time,y=.data$Log.Diff))+
  theme_classic()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),plot.margin=unit(c(5.08,5.08,5.08,5.08),"mm"),axis.text=element_text(size=24),axis.line.x=element_line(size=0.5),axis.line.y=element_line(size=0.5))
plot.ticks<-plot.basics+scale_x_continuous(breaks=seq(floor(min(ISI.data$Time)),ceiling(max(ISI.data$Time)),length=6))

plot.ISI<-plot.ticks+geom_point(aes(x=.data$Time,y=.data$Log.Diff),size=0.2)
plot.rug<-plot.ISI+geom_rug(sides="b",size=0.5)

plot.labeled<-plot.rug+labs(x="Time (sec.)",y=expression(log(T[i] - T[i-1])))+theme(axis.title=element_text(size=24))

title.graph <- paste0("Neuron: ", as.character(NeuronNumber)," - Trial #",as.character(TrialNumber))
plot.titled<-plot.labeled+ggtitle(title.graph)+theme(plot.title=element_text(size=18,face="bold"))

return(plot.titled)
}

##OLD FUNCTION
#ISIPlot <- function(NeuronNumber,TrialNumber,spikes=spikes){
#  title.graph <- paste("Neuron: ", as.character(NeuronNumber)," - Trial #",as.character(TrialNumber),sep='')
#  x <- spikes[[TrialNumber]]
#  log.diff.x<-log(diff(x))
#  time.fire<-x[2:length(x)]
#  par(mar=c(5, 4, 4, 2) + 1.2)
#  plot(time.fire,log.diff.x,pch=".",cex=2,axes=F,xlab='Time (Sec.)',ylab=expression(log(T[i] - T[i-1])),font.lab=1,cex.lab=2)
#  axis(1,lwd=1,font=1,cex.axis=2)
#  axis(2,lwd=1,font=1,cex.axis=2)
#  rug(time.fire,lwd=1)
#  title(main=title.graph,cex.main=1.5)
#}
