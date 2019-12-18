#' Fit all multiplicative and additive models
#'
#' Wrapper function to fit all 2K+1 models (multiplicative, additive, and nonperiodic)
#'
#' @param K the number of frequency components in the largest model to fit
#' @param spikes a list of spike trains
#' @param f.common.table a table whose names contain the high-amplitude frequency components as computed by \code{find.top.freqs}
#' @param setup.pars a list of additional parameters for the likelihood function, computed by \code{setup.likelihoods}
#' @param terminal.points a numeric vector containing the time points at which c(t) changes
#' @param ct a numeric vector containing the estimated piecewise constant intensity function c(t). The length of c(t) should be a whole number power of 2.
#' @param user.select whether to allow the user to select the frequencies for the model; must be FALSE for this function to run effectively
#'
#' @return a list of length 3
#' The first item in the list is a list of frequency estimates for each model
#' The second item in the list is a list of phase estimates for each model
#' The third item in the list is a list of eta/gamma estimates and fit criteria for each model
#'
#' @export

fit.all.mult.add.models <- function(K, spikes, f.common.table, setup.pars, terminal.points, ct, user.select = FALSE) {
## This abstracts Step 6 and 7 of the original script to run for each neuron
## For some value K, it fits Multiplicative 1-K, Additive 1-K, and the no periodicity model (i.e. c(t))


  # this is bad, we know EXACTLY the length of f.hat.list, etc. = 2*K + 1
 f.hat.list<-list(c())
 w0.hat.list<-list(c())
 K.list<-list(c())

for(k in 1:K){
   f.hat.list[[2*k]] <- f.hat.list[[2*(k-1)+1]] <- pick.top.K.freqs(f.common.table, k, user.select = user.select)
   w0.hat.list[[2*(k-1)+1]] <-  w0.hat.list[[2*k]] <- phase.estimation(spikes,f.hat.list[[2*k]])
   K.list[[2*(k-1)+1]] <- fit.mult.model(spikes, f.hat.list[[2*(k-1)+1]], w0.hat.list[[2*(k-1)+1]], setup.pars, terminal.points, ct)
   K.list[[2*k]] <- fit.add.model(spikes, f.hat.list[[2*k]], w0.hat.list[[2*k]], setup.pars, terminal.points, ct)
   names(f.hat.list)[c(2*k-1, 2*k)] <- paste(c("Multiplicative", "Additive"), k)
}

   f.hat.list[[2*K+1]] <- 0
   w0.hat.list[[2*K+1]] <-  phase.estimation(spikes,f.hat.list[[2*K+1]])
   K.list[[2*K+1]] <- fit.mult.model(spikes, f.hat.list[[2*K+1]], w0.hat.list[[2*K+1]], setup.pars, terminal.points, ct)

   names(f.hat.list)[2*K+1] <- "Nonperiodic"
   names(w0.hat.list) <- names(f.hat.list)
   names(K.list) <- names(f.hat.list)

return(list(f.hat.list = f.hat.list, w0.hat.list = w0.hat.list, K.list = K.list))
}
