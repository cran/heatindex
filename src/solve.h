#ifndef SOLVE_H
#define SOLVE_H

double solve(const std::function<double(double)>& f, double a, double b, double tol, int maxIter);
double solve_core(const std::function<double(double)>& f, double a, double b, double fa, double fb, double tol, int maxIter);

#endif // SOLVE_H
