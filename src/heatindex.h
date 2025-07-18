#include "thermo.h"

#ifndef HEATINDEX_H
#define HEATINDEX_H

// Human constants and functions
const double Q        = 180.               ; // W/m^2  , metabolic rate per skin area
const double phi_salt = 0.9                ; // none   , vapor saturation pressure level of saline solution
const double Tc       = 310.               ; // K      , core temeprature
const double Pc       = phi_salt*pvstar(Tc); // Pa     , core vapor pressure
const double Pa0      = 1.6e3              ; // Pa     , reference air vapor pressure in regions III, IV, V, VI, chosen by Steadman
const double hc       = 12.3               ; // W/m^2/K, heat transfer coefficient of the whole body, chosen by Steadman
const double Za       = 60.6/hc            ; // m^2Pa/W, mass transfer resistance of the whole body, chosen by Steadman

double heatindex(
   double T,  // Air temperature, K
   double rh  // Relative humidity, unitless, 0 to 1
);

#endif // HEATINDEX_H
