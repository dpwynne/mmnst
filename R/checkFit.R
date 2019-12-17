#' Check model fit
#'
#' Computes the AIC, corrected AIC (AICc), and BIC for a model
#'
#' @param l the value of the log-likelihood function evaluated at MLEs of the parameters
#' @param k the number of parameters in the model to be estimated
#' @param n the number of observations used to estimate the parameters
#'
#' @return a numeric vector containing (in order) AIC, AICc, BIC, and maximum log-likelihood

check.fit<-function(l,k,n){
model.AIC<-(-2*l+2*k)
model.AICc<-(model.AIC+(2*(k+1)*(k+2)/(n-k-2)))
model.BIC<-(-2*l+log(n)*k)
return(c(model.AIC,model.AICc,model.BIC,l))
}
