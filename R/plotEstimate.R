plot.estimate<-function(f.hat,w0.hat,eta.hat,gamma.hat,Time,terminal.points,ct,graph.title){
##does the actual plotting of the intensity function estimate
type.of.model<-substr(graph.title,1,1)

if (type.of.model=="M"){
y<-sapply(Time,theta.m,f=f.hat,w0=w0.hat,eta=eta.hat,gamma=gamma.hat,terminal.points=terminal.points,ct=ct)
}else{ ##model is additive
y<-sapply(Time,theta.a,f=f.hat,w0=w0.hat,eta=eta.hat,gamma=gamma.hat,terminal.points=terminal.points,ct=ct)}
y.graphmax<-floor(max(y)/100)*100+100

##do the actual plotting
plot(y~Time,type='s',ylim=c(0,y.graphmax),lty=2,lwd=2,ylab='intensity function',font.lab=1.5,cex.lab=1.5,axes=F,main=paste("Best Model:",graph.title),font.main=1.5)
axis(1,lwd=1,font=1,cex.axis=1.5)
axis(2,lwd=1,font=1,cex.axis=1.5)

return(cbind(Time,y))
}

