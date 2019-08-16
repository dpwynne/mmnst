AdditiveLogLikelihood.Multiple.f<-function(param,f.hat,w0.hat.itr,setup.pars,ct,ct.spike.times,individual.spike.train){
##param is eta.hat and gama.hat - parameters to optimize over
##f.hat is vector of K frequencies
##w0.hat.itr is a vector of K phases
##setup.pars is list of necessary things passed from setup.likelihoods
##ct is c(t), 2^J distinct points coming from Kolaczyk(2003) RDP
##individual.spike.train is a single spike train
##ct.spike.times is c(t) evaluated at each spike time on the individual spike train

##unpack setup.pars
DeltaDi <- setup.pars$DeltaDi
D.i.plus.one<-setup.pars$Di.1
D.i<-setup.pars$Di.0
T.data<-setup.pars$T.data

##unpack param
K<-length(f.hat)
eta.hat<-param[1:K]; gama.hat<-param[(K+1):(2*K)]; 

	# parameter check
	warning.message <- parameter.check(f.hat, w0.hat.itr, eta.hat, gama.hat)
	if(!is.null(warning.message)){
	#	warning.message <- c(warning.message, "NaN returned for log-likelihood.")
	#	cat(paste(warning.message, collapse = "\n"))
		return(NaN)
	}

## fit the actual model
      A1 <- (1-sum(eta.hat))*ct.spike.times
      A2 <- sum(eta.hat*gama.hat)
      A3 <- sapply(individual.spike.train,function(xx){sum(eta.hat*gama.hat*cos(2*pi*(f.hat*xx+w0.hat.itr)))}) 
      A <- sum(log(A1+A2+A3))
      B1<-1-sum(eta.hat)
      B2<-sum(ct*DeltaDi)
      B<-B1*B2
      C<-sum( (eta.hat*gama.hat)*(T.data+(sin(2*pi*(f.hat*T.data+w0.hat.itr))-sin(2*pi*w0.hat.itr))/(2*pi*f.hat))  ) 
      return(A-B-C)
}

