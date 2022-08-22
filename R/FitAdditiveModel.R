#' Fit additive model
#'
#' Fits the additive multiscale model.
#'
#' @importFrom stats optim
#'
#' @param spikes a list of spike trains.
#' @param f.hat a numeric vector containing frequency estimates for a particular model.
#' @param w0.hat a matrix containing phase estimates for each spike train for a particular model.
#' @param setup.pars a list of additional parameters for the likelihood function, computed by the [SetupLikelihoods()] function.
#' @param terminal.points a numeric vector containing the time points at which \eqn{c(t)} changes.
#' @param ct a numeric vector containing the estimated piecewise constant intensity function c(t). The length of ct should be a whole number power of 2.
#'
#' @return A list of length 3.
#' The first item in the list is a matrix whose rows each contain the MLEs of eta for a single spike train.
#' The second item in the list is a matrix whose rows each contain the MLEs of gamma for a single spike train.
#' The third item in the list is a matrix whose rows each contain the AIC, AICc, BIC, and log-likelihood for the model for a single spike train
#'
#' @references Ramezan, R., Marriott, P., and Chenouri, S. (2014), *Statistics in Medicine*, **33**(2), 238-256. doi: 10.1002/sim.5923.
#'
#' @export

FitAdditiveModel<-function(spikes,f.hat,w0.hat,setup.pars,terminal.points,ct){

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

ct.spike.times<-sapply(spikes[[itr]],CtAllPoints,terminal.points=terminal.points,ct=ct[itr,])
threshold<-sqrt(sum((init.par)^2))
threshold.counter <- 0
par.new<-init.par

while(threshold>(1e-5) && threshold.counter<=10){
	par.old<-par.new
	optimization<-optim((par.old),AdditiveLogLikelihood,f.hat=f.hat,w0.hat.itr=w0.hat[itr,],setup.pars=setup.pars,ct=ct[itr,],ct.spike.times=ct.spike.times,individual.spike.train=spikes[[itr]],method="Nelder-Mead", control=list(maxit=2000,fnscale=-1))
  	par.new<-optimization$par
	threshold<-sqrt(sum(par.new-par.old)^2)
	threshold.counter <- threshold.counter+1
  	new.ll<-AdditiveLogLikelihood(par.new,f.hat=f.hat,w0.hat.itr=w0.hat[itr,],setup.pars=setup.pars,ct=ct[itr,],ct.spike.times=ct.spike.times,individual.spike.train=spikes[[itr]])

	values<-round(c(par.new,threshold,new.ll),6)
	names(values)<-c(val.names,"Convergence Criterion","Log-Likelihood")
#	print(values)
}##end while loop

add.eta[itr,]<-par.new[1:K]
add.gama[itr,]<-par.new[(K+1):(2*K)]
fit.matrix[itr,]<-EvaluateModelFit(new.ll,K.hat,length(spikes[[itr]]))
}##end for loop
}else{
for (itr in 1:length(spikes)){
ct.spike.times<-sapply(spikes[[itr]],CtAllPoints,terminal.points=terminal.points,ct=ct[itr,])
add.eta[itr,]<-rep(0,K)
add.gama[itr,]<-rep(0,K)
DeltaDi <- setup.pars$DeltaDi
new.ll<- sum(log(ct.spike.times+1e-10))-sum(ct*DeltaDi)
fit.matrix[itr,]<-EvaluateModelFit(new.ll,(2^J),length(spikes[[itr]]))
}##end for loop
}##end if-else loop
cat("Additive Model with",n.frequencies,"Frequencies Fitted\n")
return(list(eta=add.eta,gama=add.gama,fit=fit.matrix))
}
