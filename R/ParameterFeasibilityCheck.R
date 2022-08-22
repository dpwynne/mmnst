#' Parameter check function
#'
#' A test function to ensure that appropriate constraints on parameters of the likelihood function are met
#'
#' @param f.hat a numeric vector containing frequency estimates for a particular model.
#' @param w0.hat.itr a numeric vector containing initial phase estimates for every spike train used to fit a particular model.
#' @param eta.hat a numeric vector containing the contribution of each periodic term to the intensity function.
#' @param gama.hat a numeric vector containing the amplitudes of each periodic term in the intensity function.
#' @param f.max a scalar indicating the absolute highest frequency of spiking activity.
#'
#' @return Warning message(s) if the constraints on parameters are violated, otherwise NULL.
#'
#' @noRd

ParameterFeasibilityCheck <- function(f.hat, w0.hat.itr, eta.hat, gama.hat, f.max = 100){
## Note: this parameter check assumes no frequency above f.max = 100 Hz

warning.message <- NULL
if(sum(f.hat<0)>0||sum(f.hat>f.max)>0){
	warning.message <- c(warning.message, "At least one frequency is outside the required range.")
}
    if(sum(eta.hat<0)>0 || sum(eta.hat>1)>0 || sum(eta.hat)>1) {
	warning.message <- c(warning.message, "Constraints on eta hat (from Ramezan et al. 2014) are not met.")
   }
     if(sum(w0.hat.itr<(-0.25))>0 || sum(w0.hat.itr>0.75)>0) {
	warning.message <- c(warning.message, "Constraints on omega (from Ramezan et al. 2014) are not met.")
     }
     if(sum(gama.hat<0)>0){
	warning.message <- c(warning.message, "Constraints on omega (from Ramezan et al. 2014) are not met.")
}

return(warning.message)
}
