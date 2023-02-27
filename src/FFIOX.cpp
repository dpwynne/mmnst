#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
IntegerVector fast_first_index_vec_ordered_x (NumericVector x, NumericVector v) {
  IntegerVector y(x.length());
  int last_index = 0;
  for (unsigned int i = 0; i < x.length(); i++) {
    y[i] = 0;
    for (unsigned int j = last_index; j < v.length(); j++) {
      if (v[j] > x[i]) {
        last_index = j;
        y[i] = j;
        break;
      }
    }
  }
  return y;
}
