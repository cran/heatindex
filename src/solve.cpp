#include <vector>
#include <iostream>
#include <cmath>
#include <string>
#include <functional>
#include "output.h"
#include "solve.h"

double solve(const std::function<double(double)>& f, double a, double b, double tol, int maxIter) {
    double fa = f(a), fb = f(b);
    if (fa*fb >= 0){ 
        STOP("Error: root not bracketed.");
    }
    else {
        return solve_core(f, a, b, fa, fb, tol, maxIter);
    }
}

double solve_core(const std::function<double(double)>& f, double a, double b, double fa, double fb, double tol, int maxIter) { // Brent's method
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
