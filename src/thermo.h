#ifndef THERMO_H
#define THERMO_H

// Thermodynamic constants and functions
const double Ttrip = 273.16      ; // K     , vapor temperature at triple point
const double ptrip = 611.65      ; // Pa    , vapor pressure at triple point
const double E0v   = 2.3740e6    ; // J/kg  , specific internal energy of vapor at the triple point
const double E0s   = 0.3337e6    ; // J/kg  , specific internal energy of solid water at the triple point
const double rgasa = 287.04      ; // J/kg/K, specific gas constant of dry air
const double rgasv = 461.        ; // J/kg/K, specific gas constant of water vapor
const double cva   = 719.        ; // J/kg/K, specific heat capacity of dry air at constant volume
const double cvv   = 1418.       ; // J/kg/K, specific heat capacity of water vapor at constant volume
const double cvl   = 4119.       ; // J/kg/K, specific heat capacity of liquid water at constant volume
const double cvs   = 1861.       ; // J/kg/K, specific heat capacity of solid water at constant volume
const double cpa   = cva + rgasa ; // J/kg/K, specific heat capacity of dry air at constant pressure
const double cpv   = cvv + rgasv ; // J/kg/K, specific heat capacity of water vapor at constant pressure

double pvstar(double T);
double qvstar(double p, double T);
double Tstar(double pv);
double wetbulb(double p,double T,double rh,bool psychrometric=false,bool icebulb=false,bool verbose=false,double lewis=0.85);
double rh_from_wetbulb(double p,double T,double Tw,bool psychrometric=false,bool icebulb=false,bool verbose=false,double lewis=0.85);

#endif // THERMO_H
