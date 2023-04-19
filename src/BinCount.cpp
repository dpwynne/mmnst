#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::uvec Cpp_BinCount_sorted(arma::vec x, arma::vec s) {
  if (s.n_elem <= 1) {
    arma::uvec binned;
    return binned;
  }
  arma::uvec binned = arma::histc(x, s);
  // Armadillo will give the last element of s a singleton-sized bin on its own.
  // Roll this into the second last bin so it looks like (2, 3] instead of (2,3)+{3}
  binned(binned.n_elem - 2) += binned(binned.n_elem - 1);
  // Truncate bins to exclude singleton
  binned = binned.head(binned.n_elem - 1);
  return binned;
}