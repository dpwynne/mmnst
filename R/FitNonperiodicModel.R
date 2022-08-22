#' Fit nonperiodic model
#'
#' Fits the nonperiodic multiscale model.
#'
#' @importFrom stats optim
#'
#' @param spikes a list of spike trains.
#' @param setup.pars a list of additional parameters for the likelihood function, computed by the [SetupLikelihoods()] function.
#' @param terminal.points a numeric vector containing the time points at which c(t) changes.
#' @param ct a numeric vector containing the estimated piecewise constant intensity function c(t). The length of ct should be a whole number power of 2.
#'
#' @return A list of length 3 is returned.
#' The first item in the list is a matrix whose rows each contain the MLEs of eta for a single spike train (this should be 0).
#' The second item in the list is a matrix whose rows each contain the MLEs of gamma for a single spike train (this should also be 0).
#' The third item in the list is a matrix whose rows each contain the AIC, AICc, BIC, and log-likelihood for the model for a single spike train.
#'
#' @references Ramezan, R., Marriott, P., and Chenouri, S. (2014), *Statistics in Medicine*, **33**(2), 238-256. doi: 10.1002/sim.5923.
#'
#' @export

FitNonperiodicModel<-function(spikes,setup.pars,terminal.points,ct){

  K<-1
  J<-setup.pars$J
  K.hat <- (K*4)+(2^J)

  cat("Fitting Nonperiodic Model\n")

  mult.gama <- mult.eta <- matrix(NA,nrow=length(spikes),ncol=K)
  fit.matrix<-matrix(NA,nrow=length(spikes),ncol=4)
  colnames(fit.matrix) <- c("AIC","AICc","BIC","Log-Likelihood")

  if (is.null(dim(ct))){ # if this is a vector of average ct, not a matrix of individual ct
    ct <- matrix(rep(ct, length(spikes)), nrow = length(spikes), byrow = TRUE)
  }

  for (itr in 1:length(spikes)){
    ct.spike.times<-sapply(spikes[[itr]],CtAllPoints,terminal.points=terminal.points,ct=ct[itr,])
    mult.eta[itr,]<-rep(0,K)
    mult.gama[itr,]<-rep(0,K)
    DeltaDi <- setup.pars$DeltaDi
    new.ll<- sum(log(ct.spike.times+1e-10))-sum(ct*DeltaDi)
    fit.matrix[itr,]<-EvaluateModelFit(new.ll,(2^J),length(spikes[[itr]]))
  }##end for loop

  cat("Nonperiodic Model Fitted\n")
  return(list(eta=mult.eta,gama=mult.gama,fit=fit.matrix))
}
