#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
IntegerVector Cpp_BinCount_sorted(NumericVector x, NumericVector s) {
    IntegerVector counts(s.length() - 1);
    unsigned int pos = 0;
    for (unsigned int binnum = 0; binnum < s.length(); binnum++){
        for (pos; pos < x.length(); pos++){
            if (x[pos] >= s[binnum + 1]) break;
            counts[binnum]++;
        }
    }  return counts;
}
