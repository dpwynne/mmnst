q.vector<-function(endpoints.bins,theta,resolution){
  ##input vector of bins and intensity function to integrate to calculate the q-values

  L <- length(endpoints.bins)
  q <- numeric((L-1))
  cat(paste('total count =',L,'\n'))

  ##integrate theta(t) over bin i to i+1 to get q_k
  for (i in 1:(L-1)){
    q[i]<-integrate(theta.function,endpoints.bins[i],endpoints.bins[i+1],theta=theta,t.start=endpoints.bins[1],resolution=resolution)$value
    if(i%%5000==0) cat(paste(i,'bins processed\n'))
  }##end for loop

  return(q)
}
