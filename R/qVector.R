#' Quantile Vector
#'
#' Calculate the q-values for the goodness-of-fit plot
#'
#' @importFrom stats integrate
#'
#' @param endpoints.bins a numeric vector containing the endpoints of the bins for the binned spike train
#' @param theta a numeric vector containing the the average of intensity function estimates across trials
#' @param resolution the time resolution of the intensity function, related to the number of points at which the function is evaluated
#'
#' @return a numeric vector of q values from Haslinger (2010)


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
