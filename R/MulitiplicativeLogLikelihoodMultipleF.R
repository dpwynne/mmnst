MultiplicativeLogLikelihood.Multiple.f<-function(param,f.hat,w0.hat.itr,setup.pars,ct,ct.spike.times,individual.spike.train){
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
#		warning.message <- c(warning.message, "NaN returned for log-likelihood.")
#		cat(paste(warning.message, collapse = "\n"))
		return(NaN)
	}

## fit the actual model
	constant<-sum(log(ct.spike.times+1e-10))-sum(ct*DeltaDi)
      A<-sum( 
		log(
		1 - sum(eta.hat) + sum(eta.hat*gama.hat) +
		sapply(individual.spike.train, function(t) sum(eta.hat*gama.hat*cos(2*pi*(f.hat*t + w0.hat.itr))))
		)
	) 
      B1<-sum(eta.hat*(gama.hat-1))
      B2<-sum(ct*DeltaDi)
      B<-B1*B2
      C.matrix <- matrix(0, nrow = length(ct), ncol = K)

     for(j in 1:length(ct)){
        for(k in 1:K){
		C.coef <- ct[j]*(eta.hat[k]*gama.hat[k])/(2*pi*f.hat[k])
		C.periodic1 <- sin(2*pi*(f.hat[k]*D.i.plus.one[j]+w0.hat.itr[k]))
		C.periodic2 <- sin(2*pi*(f.hat[k]*D.i[j]+w0.hat.itr[k]))
		C.matrix[j,k] <- C.coef*(C.periodic1 - C.periodic2)
       }
      }
	C <- sum(C.matrix)
   return(constant+A-B-C)
}
