fit.add.model<-function(spikes,f.hat,w0.hat,setup.pars,terminal.points,ct){
K<-length(f.hat)
J<-setup.pars$J
K.hat <- (K*4)+(2^J)

n.frequencies <- ifelse(f.hat[1] == 0, 0, K)

cat("Fitting Additive Model with",n.frequencies,"Frequencies\n")

add.gama <- add.eta <- matrix(NA,nrow=length(spikes),ncol=K)
fit.matrix<-matrix(NA,nrow=length(spikes),ncol=4)
colnames(fit.matrix) <- c("AIC","AICc","BIC","Log-Likelihood")

##initial guesses for eta and gama
eta.init<-rep(min(0.4,1/K),K)
gama.init<-rep(0.8,K)
init.par<-c(eta.init,gama.init)

eta.names<-paste("Eta",as.character(seq(1:K)))
gama.names<-paste("Gamma",as.character(seq(1:K)))
val.names<-c(eta.names,gama.names)

if (is.null(dim(ct))){ # if this is a vector of average ct, not a matrix of individual ct
  ct <- matrix(rep(ct, length(spikes)), nrow = length(spikes), byrow = TRUE)
}

if (!(0 %in% f.hat)){

for(itr in 1:length(spikes)){
#cat("Spike Train",itr,"\n")

ct.spike.times<-sapply(spikes[[itr]],ct.all.points,terminal.points=terminal.points,ct=ct[itr,])
threshold<-sqrt(sum((init.par)^2))
threshold.counter <- 0
par.new<-init.par

while(threshold>(1e-5) && threshold.counter<=10){
	par.old<-par.new
	optimization<-optim((par.old),AdditiveLogLikelihood.Multiple.f,f.hat=f.hat,w0.hat.itr=w0.hat[itr,],setup.pars=setup.pars,ct=ct[itr,],ct.spike.times=ct.spike.times,individual.spike.train=spikes[[itr]],method="Nelder-Mead", control=list(maxit=2000,fnscale=-1))
  	par.new<-optimization$par
	threshold<-sqrt(sum(par.new-par.old)^2)
	threshold.counter <- threshold.counter+1
  	new.ll<-AdditiveLogLikelihood.Multiple.f(par.new,f.hat=f.hat,w0.hat.itr=w0.hat[itr,],setup.pars=setup.pars,ct=ct[itr,],ct.spike.times=ct.spike.times,individual.spike.train=spikes[[itr]])

	values<-round(c(par.new,threshold,new.ll),6)
	names(values)<-c(val.names,"Convergence Criterion","Log-Likelihood")
#	print(values)
}##end while loop

add.eta[itr,]<-par.new[1:K]
add.gama[itr,]<-par.new[(K+1):(2*K)]
fit.matrix[itr,]<-check.fit(new.ll,K.hat,length(spikes[[itr]]))
}##end for loop
}else{
for (itr in 1:length(spikes)){
ct.spike.times<-sapply(spikes[[itr]],ct.all.points,terminal.points=terminal.points,ct=ct[itr,])
add.eta[itr,]<-rep(0,K)
add.gama[itr,]<-rep(0,K)
DeltaDi <- setup.pars$DeltaDi
new.ll<- sum(log(ct.spike.times+1e-10))-sum(ct*DeltaDi)
fit.matrix[itr,]<-check.fit(new.ll,(2^J),length(spikes[[itr]]))
}##end for loop
}##end if-else loop
cat("Additive Model with",n.frequencies,"Frequencies Fitted\n")
return(list(eta=add.eta,gama=add.gama,fit=fit.matrix))
}