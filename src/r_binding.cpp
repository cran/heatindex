#include <Rcpp.h>
#include "heatindex.h"

// [[Rcpp::export]]
Rcpp::NumericVector heatindex_vec(Rcpp::NumericVector Ta,
                                  Rcpp::NumericVector RH) {
  size_t n = std::max({Ta.size(), RH.size()});
  if ( (Ta.size() != 1 && Ta.size() != n) ||
       (RH.size() != 1 && RH.size() != n) ) {
    Rcpp::stop("Sizes of Ta and RH do not match");
  }
  Rcpp::NumericVector out(n);
  for (size_t i = 0; i < n; i++) {
    double ta_val = (Ta.size() == 1) ? Ta[0] : Ta[i];
    double rh_val = (RH.size() == 1) ? RH[0] : RH[i];
    out[i] = heatindex(ta_val, rh_val);
  }
  return out;
}
