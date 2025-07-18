#include <Rcpp.h>
#include "heatindex.h"
#include "thermo.h"

// [[Rcpp::export]]
Rcpp::NumericVector heatindex_vec(Rcpp::NumericVector T,
                                  Rcpp::NumericVector rh) {
  size_t n = std::max({T.size(), rh.size()});
  if ( (T .size() != 1 && T .size() != n) ||
       (rh.size() != 1 && rh.size() != n) ) {
    Rcpp::stop("Sizes of T and rh do not match");
  }
  Rcpp::NumericVector out(n);
  for (size_t i = 0; i < n; i++) {
    double  t_val = (T .size() == 1) ? T [0] : T [i];
    double rh_val = (rh.size() == 1) ? rh[0] : rh[i];
    out[i] = heatindex(t_val, rh_val);
  }
  return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector wetbulb_vec(Rcpp::NumericVector p,
                                Rcpp::NumericVector T,
                                Rcpp::NumericVector rh,
                                Rcpp::LogicalVector psychrometric, 
                                Rcpp::LogicalVector icebulb, 
                                Rcpp::LogicalVector verbose,
                                Rcpp::NumericVector Lewis) {
   size_t n = std::max({p.size(), T.size(), rh.size()});
   if ( (p .size() != 1 && p .size() != n) ||
        (T .size() != 1 && T .size() != n) ||
        (rh.size() != 1 && rh.size() != n) ) {
      Rcpp::stop("Sizes of p, T, and rh do not match");
   }
   Rcpp::NumericVector out(n);
   for (size_t i = 0; i < n; i++) {
      double  p_val = (p .size() == 1) ? p [0] : p [i];
      double  t_val = (T .size() == 1) ? T [0] : T [i];
      double rh_val = (rh.size() == 1) ? rh[0] : rh[i];
      out[i] = wetbulb(p_val, t_val, rh_val, psychrometric[0], icebulb[0], verbose[0], Lewis[0]);
   }
   return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector rh_from_wetbulb_vec(Rcpp::NumericVector p,
                                        Rcpp::NumericVector T,
                                        Rcpp::NumericVector Tw,
                                        Rcpp::LogicalVector psychrometric, 
                                        Rcpp::LogicalVector icebulb, 
                                        Rcpp::LogicalVector verbose,
                                        Rcpp::NumericVector Lewis) {
   size_t n = std::max({p.size(), T.size(), Tw.size()});
   if ( (p .size() != 1 && p .size() != n) ||
        (T .size() != 1 && T .size() != n) ||
        (Tw.size() != 1 && Tw.size() != n) ) {
      Rcpp::stop("Sizes of p, T, and Tw do not match");
   }
   Rcpp::NumericVector out(n);
   for (size_t i = 0; i < n; i++) {
      double  p_val = (p .size() == 1) ? p [0] : p [i];
      double  t_val = (T .size() == 1) ? T [0] : T [i];
      double Tw_val = (Tw.size() == 1) ? Tw[0] : Tw[i];
      out[i] = rh_from_wetbulb(p_val, t_val, Tw_val, psychrometric[0], icebulb[0], verbose[0], Lewis[0]);
   }
   return out;
}
