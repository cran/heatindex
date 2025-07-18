#include <vector>
#include <iostream>
#include <cmath>
#include <string>
#include <functional>
#include "heatindex.h"
#include "output.h"
#include "thermo.h"
#include "solve.h"

// Root solver and constants
const int maxIter = 100;
const double tol = 1.e-8; // This will be the precision of the heat index

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

bool check_input(double T, double rh){
    bool error = false;
    if (T < 0.){CERR << "T = " << T << " K. " << "Air temperature is in Kelvin, and must be positive." << std::endl; error = true;}
    if (rh < 0. || rh > 1.){CERR << "rh = " << rh << ". " << "Relative humidity must be between 0 and 1." << std::endl; error = true;}
    return error;
}

std::vector<double> physiology(double T, double rh) {
    if (check_input(T, rh)) {STOP("Inputs out of range.");};
    double Ta, Pa, Rs, CdTcdt;
    Ta    = T;
    Pa    = rh*pvstar(T);
    CdTcdt= Q-Qv(Ta,Pa)-(Tc-Ta)/Ra(Tc,Ta)-(Pc-Pa)/Za;
    Rs    = 0.;
    if (CdTcdt < 0.){
        CdTcdt = 0.;
        auto f = [Ta,Pa](double Ts){return (Ts-Ta)/Ra(Ts,Ta)+std::min((Pc-Pa)/(Zs((Tc-Ts)/(Q-Qv(Ta,Pa)))+Za),(phi_salt*pvstar(Ts)-Pa)/Za)-(Q-Qv(Ta,Pa));};
        double Ts = solve(f,0.,Tc,tol/100.,maxIter); // the "tol" is divided by 100 here so that the heat index will have precision of "tol" 
        Rs = (Tc-Ts)/(Q-Qv(Ta,Pa));}
    return {Rs,CdTcdt};
}

double heatindex(double T, double rh) {
    if (std::isnan(T) || std::isnan(rh)) {
        return std::nan("");
    }
    std::vector<double> physio = physiology(T,rh);
    double Rs, CdTcdt;
    Rs    = physio[0];
    CdTcdt= physio[1];
    if (T==0.) {return 0.;}
    if (Rs > 0.){
        auto f = [Rs](double Ta){
            double Pa = std::min(Pa0,pvstar(Ta));
            double Ts = Tc - Rs*(Q-Qv(Ta,Pa));
            double Ps = std::min((Zs(Rs)*Pa+Za*Pc)/(Zs(Rs)+Za),phi_salt*pvstar(Ts));
            return Q-Qv(Ta,Pa)-(Ts-Ta)/Ra(Ts,Ta)-(Ps-Pa)/Za;};
        return solve(f,0.,345.,tol,maxIter);}
    else {
        auto f = [CdTcdt](double Ta){
            double Pa = std::min(Pa0,pvstar(Ta));
            return Q-Qv(Ta,Pa)-(Tc-Ta)/Ra(Tc,Ta)-(Pc-Pa)/Za - CdTcdt;};      
        return solve(f,340.,T+3500.,tol,maxIter);}
}
