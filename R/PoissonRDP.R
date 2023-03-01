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
#' @source This code is based on [MATLAB code](http://math.bu.edu/people/kolaczyk/software/msglmcode.zip) from Kolaczyk and Nowak (2005). Drs. Kolaczyk and Nowak have agreed to allow our translation in the package.
#'
#' @references Kolaczyk, E.D. and Nowak, R.D. (2004). Multiscale likelihood analysis and complexity penalized estimation, *The Annals of Statistics*, **32**(2), 500-527. doi: 10.1214/009053604000000076.
#'
#' Kolaczyk, E.D. and Nowak, R.D. (2005). Multiscale generalized linear models for nonparametric function estimation. *Biometrika*,  **92**(1), 119–133. doi: 10.1093/biomet/92.1.119.
#'
#' @export

PoissonRDP<-function (sig, gamma) {
    # This is a translation and modification (computationally faster version) of the MATLAB code from Nowak and Kolaczyk (2005)
    # adjusting for constant, not polynomial, intensity function; i.e., m = 0 in the original MATLAB code
    # note that the length of the vector sig must be a power of 2

    if (sum(sig) < 1) {
      return(rep(0, length(sig)))
    }

    n <- length(sig)

    if (log2(n)!=round(log2(n))){
      stop("The length of sig must be a power of 2")
    }

    J <- log2(n)
    zeropadding = rep(0, 2^ceiling(J) - n)
    sumx = c(sig, zeropadding)
    n2 = length(sumx)

    lam <- gamma * log(n)
    bestFit <- sig + (sig == 0) * 1e-50
    bestPL <- sig * log(bestFit) - bestFit

    ind_pad = matrix(1:length(sumx), 1)

    for (j in (J - 1):0) {
      n            = 2 ^ (J - j)
      dim(ind_pad) = c(n, length(ind_pad) / n)
      sumx         = .Internal(colSums(sumx, 2L, length(sumx) / 2L, FALSE)) # if this line replaces the next one, the code is faster, but it will mess up R Check for CRAN because of the .Internal function.
      #sumx         = .colSums(sumx, 2L, length(sumx) / 2L, FALSE)
      bestFit2     = sumx / n
      bestFit2[bestFit2 < 1e-50] = 1e-50
      pl0          = sumx * log(bestFit2) - bestFit2 * n
      pl1          = .Internal(colSums(c(bestPL, zeropadding), n, n2 / n, FALSE)) # if this line replaces the next one, the code is faster, but it will mess up R Check for CRAN because of the .Internal function.
      #pl1          = .colSums(c(bestPL, zeropadding), n, n2 / n, FALSE)
      pl1          = pl1 * 2 / n - lam

      for (k in which(pl1[1:2^j] <= pl0[1:2^j])) bestFit[ind_pad[,k]] = bestFit2[k]

      comp = pl1 > pl0
      plmax = pl1 * comp + pl0 * !comp
      bestPL = rep(plmax, rep(n, length(pl0)))[1:length(bestPL)]
    }
    return(bestFit)
  }
