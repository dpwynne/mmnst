##OLD RASTER FUNCTION
#Raster <- function(NeuronNumber,spikes){
#  par(mar=c(5, 4, 4, 2) + 1.2)
#  title.graph <- paste("Neuron: ",as.character(NeuronNumber),sep='')
#  n <- length(spikes)
#  x.max <- ceiling(max(unlist(spikes)))
#  plot(c(0,0),ylim=c(0,n),xlim=c(0,x.max),type='n',axes=FALSE,pch=16,cex.lab=2,ylab='Trial Number',xlab='Time (sec.)')
#  par(new=T)
#  for(trial in 1:n){
#    x<-spikes[[trial]]
#    vert.axis<-rep(trial,length(x))
#    points(vert.axis~x,pch='.',cex=2)
#    if(trial!=n) par(new=T)
#  }
#  axis(1,lwd=1,font=1,cex.axis=2)
#  axis(2,lwd=1,font=1,cex.axis=2)
#  title(main=title.graph,cex.main=1.5)
#}

## Need to delete the library(ggplot2) line and replace ggplot2 functions with ggplot2::function

Raster <- function(NeuronNumber,spike.train=spikes){
##new Raster function for Raster plots in ggplot2

library(ggplot2)

spike.times<-unlist(spike.train)
spike.trials<-sapply(spike.train,length)
spike.numbers<-rep(seq(1,length(spike.trials)),spike.trials)

spike.data<-data.frame(Trial=spike.numbers,Time=spike.times)

tick.separator.x<-(ceiling(max(spike.data$Time))-floor(min(spike.data$Time)))/5
y.high<-50*ceiling(max(spike.data$Trial/50))
y.low<-50*floor(min(spike.data$Trial/50))
y.breaks<-seq(y.low,y.high,length=6)
y.breaks<-y.breaks[which(y.breaks<=max(spike.data$Trial))]

#plot.basics<-ggplot(spike.data,aes(x=Time,y=Trial))+theme_classic()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),plot.margin=unit(c(20.828,10.668,25.908,20.828),"mm"),axis.text=element_text(size=18),axis.line.x=element_line(size=0.5),axis.line.y=element_line(size=0.5))
plot.basics<-ggplot(spike.data,aes(x=Time,y=Trial))+theme_classic()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),plot.margin=unit(c(5.08,5.08,5.08,5.08),"mm"),axis.text=element_text(size=24),axis.line.x=element_line(size=0.5),axis.line.y=element_line(size=0.5))

plot.ticks<-plot.basics+scale_x_continuous(breaks=seq(floor(min(spike.data$Time)),ceiling(max(spike.data$Time)),by=tick.separator.x))+scale_y_continuous(breaks=y.breaks)

plot.raster<-plot.ticks+geom_point(aes(x=Time,y=Trial),size=0.2)

plot.labeled<-plot.raster+labs(x="Time (sec.)",y="Trial Number")+theme(axis.title=element_text(size=24))

title.graph <- paste0("Neuron: ",as.character(NeuronNumber))
plot.titled<-plot.labeled+ggtitle(title.graph)+theme(plot.title=element_text(size=18,face="bold"))

return(plot.titled)
}
