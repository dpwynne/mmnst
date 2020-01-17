################################################
### Compute the Hessian of the loglikelihood ###
################################################
ComputeHessian<-function(spikes,K.list,f.hat.list,w0.hat.list,terminal.points,ct,type,k){
##spikes is a list of spike trains
##K.list is a list of eta and gamma, from ____ function

setup.pars<-SetupLikelihoods(terminal.points)

etas<-lapply(K.list,"[[",1) ##list of eta matrices
gammas<-lapply(K.list,"[[",2) ##list of gamma matrices

type<-substr(type,1,1)

if (type == "m" | type== "M"){
f.hat.actual<-f.hat.list[[2*k-1]]
w0.hat.actual<-w0.hat.list[[2*k-1]]
etas.actual<-etas[[2*k-1]]
gammas.actual<-gammas[[2*k-1]]
ll.function<-MultiplicativeLogLikelihood
}else{
f.hat.actual<-f.hat.list[[2*k]]
w0.hat.actual<-w0.hat.list[[2*k]]
if (type== "a" | type == "A"){
etas.actual<-etas[[2*k]]
gammas.actual<-gammas[[2*k]]
ll.function<-AdditiveLogLikelihood
}else{
stop("Invalid type, please select \"additive\" or \"multiplicative\" \n")
}
}
param.true<-cbind(etas.actual,gammas.actual)

hess <- matrix(0, nrow=2*k, ncol=2*k)
max.itr<-dim(param.true)[1]

for (itr in 1:max.itr){
param.eta.gamma<-param.true[itr,]
spike.train<-spikes[[itr]]
ct.spike.times<-sapply(spike.train,CtAllPoints,terminal.points=terminal.points,ct=ct)
#8/8 Dwight changed the below line
#hess <- hess + hessian(ll.function,param.eta.gamma,f.hat.actual,w0.hat.actual[itr],
hess <- hess + numDeriv::hessian(ll.function,param.eta.gamma,f.hat=f.hat.actual,w0.hat.itr=w0.hat.actual[itr],
setup.pars=setup.pars,ct=ct,individual.spike.train=spikes[[itr]],ct.spike.times=ct.spike.times)
}
Observed.Fisher.Information<- -hess/max.itr
Observed.eigenvalues<-eigen(Observed.Fisher.Information)$value ##comes pre-sorted
Observed.eigenvalues[which(Observed.eigenvalues<0)]<-0
eigen.ratio<-sum(Observed.eigenvalues[1:k])/sum(Observed.eigenvalues)
#eigen.ratio<-Observaed.eigenvalues[1]/Observed.eigenvalues[k+1]
#eigen.ratio<-Observed.eigenvalues[1]/Observed.eigenvalues[2*k]

return(eigen.ratio)
}

