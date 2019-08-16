Poisson.RDP<-function(sig,gamma){
# This is a translation of the MATLAB code from Nowak and Kolaczyk (2005)
# adjusting for constant, not polynomial, intensity function; i.e., m = 0 in the original MATLAB code
  
  n <- length(sig)
  J <- log(n,2)
  lam <- gamma*log(n)
  bestFit <- sig + (sig==0)*1e-50
  bestPL <- sig*log(bestFit) - bestFit

 ## outer for loop: start at the smallest division and build up
 for(j in seq(J-1, 0, by = -1)){
 
   ## inner for loop: have to see whether each pair at the current level can be combined
   for(k in 0:(2^j-1)){

     xind <- (2^(J-j)*k + 1):(2^(J-j)*(k+1))
     x <- sig[xind]
     
     bestFit2 <- rep(max(mean(x), 1e-50), length(xind))
     pl0 <- sum(x*log(bestFit2)) - sum(bestFit2)  # don't need matrix multiplication here after all
     pl1 <- sum(bestPL[xind]*2/length(xind)) - lam
     
     bestFit[xind] <- bestFit[xind]*(pl1>pl0) + bestFit2*(pl1<=pl0)
     bestPL[xind] <- max(pl1, pl0)
     }
 }
  # If there are no spikes in a trial, this if clause will take care of it
  # in the form of returning intensity estimate of 0.
  if(sum(sig) < 1){
    bestFit = rep(0,length(bestFit))
  }
  return(bestFit)
}