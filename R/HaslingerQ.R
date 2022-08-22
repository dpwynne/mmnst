#' Quantile Vector
#'
#' Calculate the q-values for the goodness-of-fit plot.
#'
#' @importFrom stats integrate
#'
#' @param endpoints.bins a numeric vector containing the endpoints of the bins for the binned spike train.
#' @param theta a numeric vector containing the average of intensity function estimates across trials.
#' @param resolution the time resolution of the intensity function, related to the number of points at which the function is evaluated.
#'
#' @return A numeric vector of q values from Haslinger *et al.* (2010).
#'
#' @references Haslinger, R., Pipa G., and Brown, E. (2010). Discrete time rescaling theorem: determining goodness of fit for discrete time statistical models of neural spiking. *Neural Computation*. **22**(10):2477-506. doi: 10.1162/NECO_a_00015.
#'
#' @noRd


HaslingerQ <- function(endpoints.bins,theta,resolution){
  ##input vector of bins and intensity function to integrate to calculate the q-values

  L <- length(endpoints.bins)
  q <- numeric((L-1))
  cat(paste('total count =',L,'\n'))

  ##integrate theta(t) over bin i to i+1 to get q_k
  for (i in 1:(L-1)){
    q[i]<-integrate(ThetaT,endpoints.bins[i],endpoints.bins[i+1],theta=theta,t.start=endpoints.bins[1],resolution=resolution)$value
    if(i%%5000==0) cat(paste(i,'bins processed\n'))
  }##end for loop

  return(q)
}
