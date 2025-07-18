#include <vector>
#include <iostream>
#include <cmath>
#include <string>
#include <functional>
#include "output.h"
#include "thermo.h"
#include "lambert_w.h"
#include "solve.h"
   
double Le(double T) {
   return E0v + (cvv-cvl)*(T-Ttrip) + rgasv*T;
}
   
double Lm(double T) {
   return E0s + (cvl-cvs)*(T-Ttrip);
}  

// Saturation vapor pressure
// [[Rcpp::export]]
double pvstarl(double T) {
   if (T <= 0.) {
      return 0.;
   } else {
      return ptrip * pow(T/Ttrip,(cpv-cvl)/rgasv) *
         exp( (E0v       - (cvv-cvl)*Ttrip)/rgasv * (1./Ttrip - 1./T) );
   }
}
// [[Rcpp::export]]
double pvstars(double T) {
   if (T <= 0.) {
      return 0.;
   } else {
      return ptrip * pow(T/Ttrip,(cpv-cvs)/rgasv) *
         exp( (E0v + E0s - (cvv-cvs)*Ttrip)/rgasv * (1./Ttrip - 1./T) );
   }
}
// [[Rcpp::export]]
double pvstar(double T) {
   if (T < Ttrip) {
      return pvstars(T);
   } else {
      return pvstarl(T);
   }
}

// The mass fraction of water vapor at saturation with pressure p and temperature T assuming there is no condensate. 
// [[Rcpp::export]]
double qvstarl(double p, double T) {
   if (T <= 0.) {
      return 0.;
   } else if (pvstarl(T) > p) {
      return pvstarl(T)/p;
   } else {
      return 1. / ( rgasv*p/(rgasa*ptrip) * pow(Ttrip/T,(cpv-cvl)/rgasv) *
         exp( - (E0v       - (cvv-cvl)*Ttrip) / rgasv * (1./Ttrip - 1./T) ) - rgasv/rgasa + 1. );
   }
}
// [[Rcpp::export]]
double qvstars(double p, double T) {
   if (T <= 0.) {
      return 0.;
   } else if (pvstars(T) > p) {
      return pvstars(T)/p;
   } else {
      return 1. / ( rgasv*p/(rgasa*ptrip) * pow(Ttrip/T,(cpv-cvs)/rgasv) *
         exp( - (E0v + E0s - (cvv-cvs)*Ttrip) / rgasv * (1./Ttrip - 1./T) ) - rgasv/rgasa + 1. );
   }
}
// [[Rcpp::export]]
double qvstar(double p, double T) {
   if (T < Ttrip) {
      return qvstars(p,T);
   } else {
      return qvstarl(p,T);
   }
}

// The temperature at which air would be saturated with a water-vapor partial pressure of pv.
// [[Rcpp::export]]
double Tstarl(double pv) {
   double a, b, c;
   int branch;
   if (pv <= 0.) {
      return 0.;
   }
   c = (E0v - (cvv-cvl)*Ttrip) / (rgasv*Ttrip);
   b = (cpv-cvl)/rgasv;
   a = ptrip*exp(c);
   branch = -1;
   return Ttrip*c/(b*lambert_w(c/b*pow(pv/a,-1./b),branch,0));
}
// [[Rcpp::export]]
double Tstars(double pv) {
   double c = (E0v + E0s - (cvv-cvs)*Ttrip) / (rgasv*Ttrip);
   double b = (cpv-cvs)/rgasv;
   double a = ptrip*exp(c);
   double l1 = log(c/b) - log(pv/a)/b;
   int branch = 0;
   if (std::isnan(pv)) {
      return std::numeric_limits<double>::quiet_NaN();
   } else if (pv < 0.) {
      return std::numeric_limits<double>::quiet_NaN();
   } else if (pv == 0.) {
      return 0.;
   } else if (l1 < 709.) {
      return Ttrip*c/(b*lambert_w(exp(l1),branch,0));
   } else {
      double l2 = log(l1);
      return Ttrip*c/(b*(l1-l2+l2/l1+l2*(-2+l2)/(2*l1*l1) +
         l2*(6-9*l2+2*l2*l2)/(6*l1*l1*l1) +
         l2*(-12+36*l2-22*l2*l2+3*l2*l2*l2)/(12*l1*l1*l1*l1)));
   }
}
// [[Rcpp::export]]
double Tstar(double pv) {
   if (pv < ptrip) {
      return Tstars(pv);
   } else {
      return Tstarl(pv);
   }
}


// rh is the relative humidity with respect to liquid for T > 273.15 K
//    and with respect to solid for T <= 273.15 K
// psychrometric=TRUE -> Calculate psychrometric wet-bulb temperature
//    instead of thermodynamic, achieved by setting lewistwothirds to its actual value
// icebulb=TRUE -> Calculate ice-bulb temperature instead of wet-bulb
double wetbulb(double p,double T,double rh,bool psychrometric,bool icebulb,bool verbose,double lewis) {
   if (std::isnan(p) || std::isnan(T) || std::isnan(rh)) {
      return std::nan("");
   }
   double lewistwothirds;
   if (psychrometric) { // psychrometric (a.k.a., aspirated, ventilated)
      lewistwothirds = pow(lewis,2./3.);
   } else {
      lewistwothirds = 1.; // thermodynamic
   }
   auto qvstar = icebulb ? +qvstars : +qvstarl;
   auto Tstar  = icebulb ?  +Tstars :  +Tstarl;
   auto L      = icebulb ? std::function<double(double)>(
                              [](double T){ return Le(T)+Lm(T); })
                         : std::function<double(double)>(Le);
   double pv;
   if (T > Ttrip) {
      pv = rh * pvstarl(T);
   } else {
      pv = rh * pvstars(T);
   }
   if (pv > p) {
      if (verbose) {
         WARN("pv = " << pv << " Pa exceeds p = " << p);
      }
      return std::numeric_limits<double>::quiet_NaN();
   }
   double qv = rgasa*pv / (rgasv*p - rgasv*pv + rgasa*pv);
   double cpm = (1.-qv)*cpa + qv*cpv;
   double minTw = Tstar(pv);
   double maxTw = std::min(T,Tstar(p));
   double lewistwothirdscpm = lewistwothirds*cpm;
   auto error = [p,lewistwothirdscpm,T,qv,qvstar,L](double Tw) {
      double qvs = qvstar(p,Tw);
      return lewistwothirdscpm*(Tw-T)*(1.-qvs) + (qvs-qv)*L(Tw);
   };
   double minTwerror = error(minTw);
   double maxTwerror = error(maxTw);
   if (minTwerror == 0.) {
      return minTw;
   } else if (maxTwerror == 0. || minTwerror*maxTwerror > 0.) {
      return maxTw;
   } else {
      return solve_core(error, minTw, maxTw, minTwerror, maxTwerror, 1.e-10, 1000);
   }
}

// Returns the relative humidity with respect to liquid for T > 273.16 K
//    and with respect to solid for T <= 273.16 K
// psychrometric=TRUE -> Interprets Tw to be psychrometric 
//    instead of thermodynamic
// icebulb=TRUE -> Interprets Tw to be the ice-bulb temperature
double rh_from_wetbulb(double p,double T,double Tw,bool psychrometric,bool icebulb,bool verbose,double lewis) {
   if (std::isnan(p) || std::isnan(T) || std::isnan(Tw)) {
      return std::nan("");
   }
   double lewistwothirds;
   if (psychrometric) { // psychrometric (a.k.a., aspirated, ventilated)
      lewistwothirds = pow(lewis,2./3.);
   } else {
      lewistwothirds = 1.; // thermodynamic
   }
   double pvsTw, L, pvsT;
   if (icebulb) {
      pvsTw = pvstars(Tw);
      L = Le(Tw) + Lm(Tw);
   } else {
      pvsTw = pvstarl(Tw);
      L = Le(Tw);
   }
   if (T<Ttrip) {
      pvsT = pvstars(T);
   } else {
      pvsT = pvstarl(T);
   }
   double rh = p/pvsT*(
      (rgasa*L*pvsTw - lewistwothirds*rgasv*cpa*(T-Tw)*(p-pvsTw)) /
      (rgasa*L*p + lewistwothirds*(rgasa*cpv-rgasv*cpa)*(T-Tw)*(p-pvsTw)) );
   if ((rh < 1.e-14) & (rh > -1.e-14)) {
      return 0.;
   } else if (rh <= -1.e-14) {
      if (verbose) {
         char wi = icebulb ? 'i' : 'w';
         WARN("(p,T,T" << wi << ") = (" << p << ", " << T << ", " << Tw << 
            ") cannot be achieved even with zero relative humidity");
      }
      return std::numeric_limits<double>::quiet_NaN();
   } else {
      return rh;
   }
}
