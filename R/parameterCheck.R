parameter.check <- function(f.hat, w0.hat.itr, eta.hat, gama.hat, f.max = 100){
##parameter checks for log-likelihood function
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