plot.theta.t<-function(spikes,f.hat.list,w0.hat.list,K.list,best.models,model.names,Time,terminal.points,ct,intensity.function.length=(1+max(Time)*10^(4-ceiling(log10(diff(range(Time))))))){
#plots average intensity estimate
##f.hat.list, w0.hat.list are lists of f.hat and w0.hat for each model
##K.list is the list of matrices of eta.hat,gama.hat,fit for each model
##intensity.function.length is the number of points in the discretized intensity function
##for 10 seconds or less it defaults to 1 ms; for 10-100 seconds it defaults to 10 ms
##increase the value of intensity.function.length for better resolution

cat("Determining Intensity Function for Best Models\n")

#selected.models<-unique(best.models)
#selected.names<-unique(model.names)
##these two should have the same length

##change to this
selected.models<-best.models
selected.names<-model.names

nmodels<-length(selected.models)

etas<-lapply(K.list,"[[",1)
gamas<-lapply(K.list,"[[",2)

int.estimate<-vector("list",length=nmodels)
names(int.estimate)<-c("AIC","AICc","BIC")

#par(mfrow=c(nmodels,1))

#for (i in 1:nmodels){
for (i in c(2,3)){
##this ensures that only AICc (i=2) and BIC (i=3) run

list.number<-selected.models[i]

theta.t<-matrix(NA,intensity.function.length,length(spikes))

Time.vector<-seq(min(Time),max(Time),length=intensity.function.length)

for(ntrains in 1:length(spikes)){

f.hat<-f.hat.list[[list.number]]
w0.hat<-w0.hat.list[[list.number]][ntrains,]
eta.hat<-etas[[list.number]][ntrains,]
gamma.hat<-gamas[[list.number]][ntrains,]

if (list.number%%2==1){
theta.t[,ntrains]<-y<-sapply(Time.vector,theta.m,f=f.hat,w0=w0.hat,eta=eta.hat,gamma=gamma.hat,terminal.points=terminal.points,ct=ct)
}else{
theta.t[,ntrains]<-y<-sapply(Time.vector,theta.a,f=f.hat,w0=w0.hat,eta=eta.hat,gamma=gamma.hat,terminal.points=terminal.points,ct=ct)
}

}##end ntrains for loop

#return(theta.t)

selection.criteria<-c("AIC","AICc","BIC")

theta.avg<-apply(theta.t,1,mean)


cat(paste("Determining Bootstrap Confidence Interval for",selected.names[i],"Model\n"))
theta.bootstrap<-bootstrapped.theta(theta.t)
theta.matrix<-cbind(Time.vector,theta.bootstrap[,1],theta.avg,theta.bootstrap[,2])
colnames(theta.matrix)<-c("Time","Lower","Average","Upper")
int.estimate[[i]]<-theta.matrix
#y.graphmax<-floor(max(theta.bootstrap)/100)*100+100
y.graphmax<-max(theta.bootstrap)
dev.new()
par(mar=c(5, 4, 4, 2) + 1.2)
plot(theta.avg~Time.vector,type='s',ylim=c(0,y.graphmax),lty=1,lwd=2,ylab='intensity function',xlab='Time',cex.lab=1.5,axes=F,main=paste0("Best Model by ",selection.criteria[i],": ",selected.names[i]),cex.main=1.5,font.lab=1,cex.lab=2)
axis(1,lwd=1,font=1,cex.axis=2)
axis(2,lwd=1,font=1,cex.axis=2)
lines(theta.bootstrap[,1]~Time.vector,type='s',lty=2,lwd=1,col='blue')
lines(theta.bootstrap[,2]~Time.vector,type='s',lty=2,lwd=1,col='blue')

# C.L.T. limits
#l <- theta.avg - 1.96*sd(theta.avg)/sqrt(10)
#u <- theta.avg + 1.96*sd(theta.avg)/sqrt(10)
#lines(u~Time.vector,type='s',lty=2,lwd=1,col='red')
#lines(l~Time.vector,type='s',lty=2,lwd=1,col='red')

}##end i for loop
cat("Intensity Function Plotted\n")
#return(int.estimate)
return(int.estimate[c(2,3)])
}