#include <vector>
#include <iostream>
#include <cmath>
#include <string>
#include <functional>
#include "heatindex.h"
#include "output.h"

using std::min;

// Thermaodynamic constants and functions
const double Ttrip       = 273.16                  ;// K         , vapor temperature at triple point
const double ptrip       = 611.65                  ;// Pa        , vapor pressure at triple point
const double E0v         = 2.3740e6                ;// J/kg      , specific internal energy of vapor at the triple point
const double E0s         = 0.3337e6                ;// J/kg      , specific internal energy of solid water at the triple point
const double rgasa       = 287.04                  ;// J/kg/K    , specific gas constant of dry air
const double rgasv       = 461.                    ;// J/kg/K    , specific gas constant of water vapor
const double cva         = 719.                    ;// J/kg/K    , specific heat capacity of dry air at constant volume
const double cvv         = 1418.                   ;// J/kg/K    , specific heat capacity of water vapor at constant volume
const double cvl         = 4119.                   ;// J/kg/K    , specific heat capacity of liquid water at constant volume
const double cvs         = 1861.                   ;// J/kg/K    , specific heat capacity of solid water at constant volume
const double cpa         = cva + rgasa             ;// J/kg/K    , specific heat capacity of dry air at constant pressure
const double cpv         = cvv + rgasv             ;// J/kg/K    , specific heat capacity of water vapor at constant pressure

double pvstar(double T){
    if      (T <= 0.) { return 0.;} 
    else if (T<Ttrip) { return ptrip * pow(T/Ttrip,(cpv-cvs)/rgasv) * exp( (E0v + E0s -(cvv-cvs)*Ttrip)/rgasv * (1./Ttrip - 1./T) ); } 
    else {              return ptrip * pow(T/Ttrip,(cpv-cvl)/rgasv) * exp( (E0v       -(cvv-cvl)*Ttrip)/rgasv * (1./Ttrip - 1./T) ); }
}

// Root solver and constants
const int maxIter = 100;
const double tol = 1.e-8; // This will be the precision of the heat index

double solve(const std::function<double(double)>& f, double a, double b, double tol, int maxIter) { // Brent's method
    double fa = f(a), fb = f(b);
    if (fa*fb >= 0){ 
        STOP("Error: root not bracketed.");
    }
    if (fabs(fa) < fabs(fb)){std::swap(a,b); std::swap(fa,fb);}
    double c = a, fc = fa, s = b, d = b - a;
    bool mflag = true;
    
    for (int i = 0; i < maxIter; i++){        
        if (fa != fc && fb != fc){ // inverse quadratic interpolation
              s = a*fb*fc/((fa-fb)*(fa-fc)) + b*fa*fc/((fb-fa)*(fb-fc)) + c*fa*fb/((fc-fa)*(fc-fb));}
        else{ s = b-fb*(b-a)/(fb-fa);} //secant

        if (!( (s>(3*a+b)/4 && s<b) || (s<(3*a+b)/4 && s>b)) || (mflag && fabs(s - b) >= fabs(b - c) / 2) || (!mflag && fabs(s - b) >= fabs(c - d) / 2))
        {s = (a+b)/2.; mflag = true;} // use bisection
        else {mflag = false;}         // accept s
        
        double fs = f(s);
        d = c; c = b; fc = fb;
        if (fa*fs < 0) {b=s; fb=fs;}
        else {a=s; fa=fs;}
        if (fabs(fa)<fabs(fb)){std::swap(a,b); std::swap(fa,fb);}
        if (fabs(b-a) < tol){return b;}
    }
    STOP("Max iterations reached.");
}

// Human constants and functions
const double Q           = 180.                    ;// W/m^2     , metabolic rate per skin area
const double phi_salt    = 0.9                     ;// none      , vapor saturation pressure level of saline solution
const double Tc          = 310.                    ;// K         , core temeprature
const double Pc          = phi_salt*pvstar(Tc)     ;// Pa        , core vapor pressure
const double Pa0         = 1.6e3                   ;// Pa        , reference air vapor pressure in regions III, IV, V, VI, chosen by Steadman
const double hc          = 12.3                    ;// W/m^2/K   , heat transfer coefficient of the whole body, chosen by Steadman
const double Za          = 60.6/hc                 ;// m^2Pa/W   , mass transfer resistance of the whole body, chosen by Steadman

double Qv(double Ta, double Pa){ // respiratory heat loss, W/m^2
    const double p    = 1.013e5      ;// Pa    , atmospheric pressure
    const double eta  = 1.43e-6      ;// kg/J  , "inhaled mass" / "metabolic rate"
    const double L    = 2417405.2    ;// J/kg  , latent heat of vaporization of water 
    return eta * Q *(cpa*(Tc-Ta) + L*rgasa/(p*rgasv) * ( Pc-Pa ) );
}

double Zs(double Rs){ // mass transfer resistance through skin, Pa m^2/W
    return 6.0e8 * Rs*Rs*Rs*Rs*Rs;
}

double Ra(double T1, double T2){ // heat transfer resistance, K m^2/W
    const double epsilon = 0.97     ;
    const double phi_rad = 0.80     ;
    const double sigma   = 5.67e-8  ;
    double hr  = epsilon * phi_rad * sigma* (T1*T1 + T2*T2)*(T1 + T2) ;
    return 1./(hc+hr);
}

bool check_input(double Ta, double RH){
    bool error = false;
    if (Ta < 0.){CERR << "Ta = " << Ta << " K. " << "Air temperature is in Kelvin, and must be positive." << std::endl; error = true;}
    if (RH < 0. || RH > 1.){CERR << "RH = " << RH << ". " << "Relative humidity must be between 0 and 1." << std::endl; error = true;}
    return error;
}

std::vector<double> physiology(double Ta, double RH) {
    if (check_input(Ta, RH)) {STOP("Inputs out of range.");};
    double Pa, Rs, CdTcdt;
    Pa    = RH*pvstar(Ta);
    CdTcdt= Q-Qv(Ta,Pa)-(Tc-Ta)/Ra(Tc,Ta)-(Pc-Pa)/Za;
    Rs    = 0.;
    if (CdTcdt < 0.){
        CdTcdt = 0.;
        auto f = [Ta,Pa](double Ts){return (Ts-Ta)/Ra(Ts,Ta)+min((Pc-Pa)/(Zs((Tc-Ts)/(Q-Qv(Ta,Pa)))+Za),(phi_salt*pvstar(Ts)-Pa)/Za)-(Q-Qv(Ta,Pa));};
        double Ts = solve(f,0.,Tc,tol/100.,maxIter); // the "tol" is divided by 100 here so that the heat index will have precision of "tol" 
        Rs = (Tc-Ts)/(Q-Qv(Ta,Pa));}
    return {Rs,CdTcdt};
}

double heatindex(double Ta, double RH) {
    std::vector<double> physio = physiology(Ta,RH);
    double Rs, CdTcdt;
    Rs    = physio[0];
    CdTcdt= physio[1];
    if (Ta==0.) {return 0.;}
    if (Rs > 0.){
        auto f = [Rs](double Ta){
            double Pa = min(Pa0,pvstar(Ta));
            double Ts = Tc - Rs*(Q-Qv(Ta,Pa));
            double Ps = min((Zs(Rs)*Pa+Za*Pc)/(Zs(Rs)+Za),phi_salt*pvstar(Ts));
            return Q-Qv(Ta,Pa)-(Ts-Ta)/Ra(Ts,Ta)-(Ps-Pa)/Za;};
        return solve(f,0.,345.,tol,maxIter);}
    else {
        auto f = [CdTcdt](double Ta){
            double Pa = min(Pa0,pvstar(Ta));
            return Q-Qv(Ta,Pa)-(Tc-Ta)/Ra(Tc,Ta)-(Pc-Pa)/Za - CdTcdt;};      
        return solve(f,340.,Ta+3500.,tol,maxIter);}
}
