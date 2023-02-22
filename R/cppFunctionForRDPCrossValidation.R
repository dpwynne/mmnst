#This file is just to create the cpp function fast_first_index_vec_ordered_x
# which is used in the RDPCrossValidation function. We must find an appropriate way of including
# this function somewhere in the package. It seems that all files which need compiling should be put
# in the scr folder of the package, but then we should write this as a cpp file.

# Assumes input x is sorted in ascending order so that the inner loop can start
# searching from the previous index. Another sub-OoM improvement on top and
# should scale better for large v.
cppFunction("IntegerVector fast_first_index_vec_ordered_x (NumericVector x, NumericVector v) {
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
}", depends = "Rcpp")
