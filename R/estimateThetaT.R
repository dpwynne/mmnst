estimate.theta.t<-function(spikes,f.hat.list,w0.hat.list,K.list,best.models,Time,terminal.points,ct,intensity.function.length=(1+diff(range((Time)))*10^(3-ceiling(log10(diff(range(Time))))))){
#produces average intensity estimate and bootstrap confidence intervals without plotting
##f.hat.list, w0.hat.list are lists of f.hat and w0.hat for each model
##K.list is the list of matrices of eta.hat,gama.hat,fit for each model
##intensity.function.length is the number of points in the discretized intensity function
##for 10 seconds or less it defaults to 10 ms; for 10-100 seconds it defaults to 100 ms
##increase the value of intensity.function.length for better resolution

cat("Determining Intensity Function for Models\n")

#selected.models<-unique(best.models)
#selected.names<-unique(model.names)
##these two should have the same length

##change to this
selected.models<-best.models$model.numbers
selected.names<-best.models$model.names

model.unique<-unique(best.models)
model.bootstrapped<-numeric(length(model.unique))

nmodels<-length(selected.models)

etas<-lapply(K.list,"[[",1)
gamas<-lapply(K.list,"[[",2)

int.estimate<-vector("list",length=nmodels)
names(int.estimate)<-names(selected.names)

#par(mfrow=c(nmodels,1))

fitted.models<-c()

theta.t.list <- vector("list", length = nmodels)
#names(theta.t.list) <- c("AIC","AICc","BIC")

for (i in 1:nmodels){
#for (i in c(2,3)){
##this ensures that only AICc (i=2) and BIC (i=3) run

list.number<-selected.models[i]

theta.t<-matrix(NA,intensity.function.length,length(spikes))

Time.vector<-seq(min(Time),max(Time),length=intensity.function.length)

if (is.null(dim(ct))){ # if this is a vector of average ct, not a matrix of individual ct
  ct <- matrix(rep(ct, length(spikes)), nrow = length(spikes), byrow = TRUE)
}

for(ntrains in 1:length(spikes)){

f.hat<-f.hat.list[[list.number]]
w0.hat<-w0.hat.list[[list.number]][ntrains,]
eta.hat<-etas[[list.number]][ntrains,]
gamma.hat<-gamas[[list.number]][ntrains,]

if (list.number%%2==1){
theta.t[,ntrains]<-sapply(Time.vector,theta.m,f=f.hat,w0=w0.hat,eta=eta.hat,gamma=gamma.hat,terminal.points=terminal.points,ct=ct[ntrains,])
}else{
theta.t[,ntrains]<-sapply(Time.vector,theta.a,f=f.hat,w0=w0.hat,eta=eta.hat,gamma=gamma.hat,terminal.points=terminal.points,ct=ct[ntrains,])
}

}##end ntrains for loop

theta.t.list[[i]] <- theta.t

theta.avg<-apply(theta.t,1,mean)


if (!(i %in% fitted.models)){
cat(paste("Determining Bootstrap Confidence Interval for",selected.names[i],"Model\n"))
theta.bootstrap<-bootstrapped.theta(theta.t)
theta.matrix<-cbind(Time.vector,theta.bootstrap[,1],theta.avg,theta.bootstrap[,2])
colnames(theta.matrix)<-c("Time","Lower","Average","Upper")
for (j in which(selected.models==selected.models[i])){
int.estimate[[j]]<-theta.matrix
fitted.models<-c(fitted.models,j)
}##end for 
}##end if 
}##end i for loop

names(theta.t.list) <- selected.names

cat("Intensity Function Estimated\n")
#return(int.estimate)
#return(int.estimate[c(2,3)])
#return(list(AIC = int.estimate[[1]],
#	AICc = int.estimate[[2]], 
#	BIC = int.estimate[[3]],
#	individual.thetas = theta.t.list))
return(list(int.estimate = int.estimate,
	individual.thetas = theta.t.list))
}