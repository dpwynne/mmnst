#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
IntegerVector Cpp_BinCount_sorted(NumericVector x, NumericVector s) {
    IntegerVector counts(s.length() - 1);
    for (R_xlen_t binnum = 0; binnum < s.length(); binnum++){
        for (R_xlen_t pos = 0; pos < x.length(); pos++){
            if (x[pos] >= s[binnum + 1]) break;
            counts[binnum]++;
        }
    }  return counts;
}
