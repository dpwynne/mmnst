#' Additive Model Log Likelihood
#'
#' Computes the log-likelihood of a spike train under the assumptions of the additive model.
#'
#' @param param a numeric vector containing the value of the eta and gamma parameters in the model.
#' For a model containing K frequencies, the first K numbers are eta values and the second K numbers are gamma values.
#' This must be the first argument as we wish to maximize the value of the log-likelihood function over the entire set of eta and gamma parameters.
#' @param f.hat a numeric vector containing frequency estimates for a particular model.
#' @param w0.hat.itr a numeric vector containing phase estimates for a particular model and spike train.
#' @param setup.pars a list of additional parameters for the likelihood function, computed by the [SetupLikelihoods()] function.
#' @param ct a numeric vector containing the estimated piecewise constant intensity function. The length of `ct` should be a whole number power of 2.
#' @param ct.spike.times a numeric vector containing the values of `ct` at the specific times a spike was recorded.
#' @param individual.spike.train a numeric vector containing the spike times for that spike train.
#'
#' @return The value of the log-likelihood function for the additive model.
#'
#' @references Ramezan, R., Marriott, P., and Chenouri, S. (2014), *Statistics in Medicine*, **33**(2), 238-256. doi: 10.1002/sim.5923.
#'
#' @export

AdditiveLogLikelihood<-function(param,f.hat,w0.hat.itr,setup.pars,ct,ct.spike.times,individual.spike.train){

  if (any(f.hat <= 0)){
    stop("Frequencies should all be positive")
  }

##unpack setup.pars
DeltaDi <- setup.pars$DeltaDi
D.i.plus.one<-setup.pars$Di.1
D.i<-setup.pars$Di.0
T.data<-setup.pars$T.data

##unpack param
K<-length(f.hat)
eta.hat<-param[1:K]; gama.hat<-param[(K+1):(2*K)];

	# parameter check
	warning.message <- ParameterFeasibilityCheck(f.hat, w0.hat.itr, eta.hat, gama.hat)
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

