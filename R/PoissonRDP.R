#' Poisson Recursive Dyadic Partitioning
#'
#'
#' Calculates the piecewise constant intensity function using penalized likelihood recursive dyadic partitioning method in a Poisson model.
#'
#' @param sig a numeric vector determining the raw data, e.g., the heights (counts) of the histograms to be smoothed.
#' @param gamma a scalar determining the penalty factor in the penalized likelihood method.
#'
#' @return A numeric vector of estimated c(t), the smoothed histogram of data.
#'
#' @source This code is based on \href{http://math.bu.edu/people/kolaczyk/software/msglmcode.zip}{MATLAB code} from Kolaczyk and Nowak (2005). Drs. Kolaczyk and Nowak have agreed to allow our translation in the package.
#'
#' @references Kolaczyk, E.D. and Nowak, R.D. (2004). Multiscale likelihood analysis and complexity penalized estimation, \emph{The Annals of Statistics}, \strong{32}(2), 500-527. doi: 10.1214/009053604000000076.
#'
#' Kolaczyk, E.D. and Nowak, R.D. (2005). Multiscale generalized linear models for nonparametric function estimation. \emph{Biometrika},  \strong{92}(1), 119–133. doi: 10.1093/biomet/92.1.119.
#'
#' @export

PoissonRDP<-function(sig,gamma){
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
