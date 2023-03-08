#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
IntegerVector fast_first_index_vec_ordered_x (NumericVector x, NumericVector v) {
  IntegerVector y(x.length());
  R_xlen_t last_index = 0;
  for (R_xlen_t i = 0; i < x.length(); i++) {
    y[i] = 0;
    for (R_xlen_t j = last_index; j < v.length(); j++) {
      if (v[j] > x[i]) {
        last_index = j;
        y[i] = j;
        break;
      }
    }
  }
  return y;
}
