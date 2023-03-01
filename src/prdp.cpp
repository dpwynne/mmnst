#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

/* PoissonRDP<-function (sig, gamma) {
 if (sum(sig) < 1) {return(rep(0, length(sig)))}

 n <- length(sig)
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
 sumx         = .Internal(colSums(sumx, 2L, length(sumx) / 2L, FALSE))
 bestFit2     = sumx / n
 bestFit2[bestFit2 < 1e-50] = 1e-50
 pl0          = sumx * log(bestFit2) - bestFit2 * n
 pl1          = .Internal(colSums(c(bestPL, zeropadding), n, n2 / n, FALSE))
 pl1          = pl1 * 2 / n - lam

 for (k in which(pl1[1:2^j] <= pl0[1:2^j])) bestFit[ind_pad[,k]] = bestFit2[k]

 comp = pl1 > pl0
 plmax = pl1 * comp + pl0 * !comp
 bestPL = rep(plmax, rep(n, length(pl0)))[1:length(bestPL)]
 }
 return(bestFit)
 } */

// [[Rcpp::export]]
arma::vec cpp_PoissonRDP(arma::vec sig, double gamma) {
  int n = sig.n_elem;

  arma::vec bestFit(n, arma::fill::zeros);
  if (sum(sig) < 1) return bestFit; // All zeros

  double J = log2(n);
  int n2 = pow(2, ceil(J));
  arma::vec sumx(n2, arma::fill::zeros);

  sumx.head(n) = sig;

  double lam = gamma * log(n);
  bestFit = arma::clamp(sig, 1e-50, arma::datum::inf);
  arma::vec bestPL(n2, arma::fill::zeros);
  bestPL.head(n) = sig % arma::log(bestFit) - bestFit;

  arma::umat ind_pad = arma::regspace<arma::uvec>(0L, n2 - 1L);

  for (double j = (J) - 1; j >= 0L; j--) { // floor(J)? ceil(J)? int j?
    n = pow(2L, (J) - j); // floor(J)? ceil(J)?
    ind_pad.reshape(n, ind_pad.n_elem / n);
    sumx = arma::sum(arma::reshape(sumx, 2, sumx.n_elem / 2), 0).t(); // colsum of reshaped
    arma::vec bestFit2 = sumx / n;
    bestFit2 = arma::clamp(bestFit2, 1e-50, arma::datum::inf);
    arma::vec pl0 = sumx % arma::log(bestFit2) - bestFit2 * n; // elementwise mult. %
    arma::vec pl1 = arma::sum(arma::reshape(bestPL, n, bestPL.n_elem / n), 0).t(); // colsum of reshaped
    pl1 = pl1 * 2 / n - lam;
    arma::uvec indices = arma::find(pl1.head(pow(2, j)) <= pl0.head(pow(2, j))); // 2^j not guaranteed
    for (unsigned int k : indices) {
      bestFit(ind_pad.col(k)).fill(bestFit2(k));
    }

    arma::vec plmax = arma::max(pl1, pl0); // element-wise max
    arma::vec bestPLupdate = arma::vectorise(arma::repelem(plmax, n, 1));
    bestPL.head(sig.n_elem) = bestPLupdate.head(sig.n_elem);
  }
  return bestFit;
}

/*** R
set.seed(1)
x = sort(rpois(1000, 10))
max(abs(PoissonRDP(x, 1) - as.numeric(cpp_PoissonRDP(x, 1))))

library(microbenchmark)
microbenchmark( PoissonRDP(x, 1),
                cpp_PoissonRDP(x, 1), times = 1)
*/
