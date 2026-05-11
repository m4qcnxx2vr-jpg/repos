#ifndef INTEGRATE_H
#define INTEGRATE_H

#include <functional>
#include <cmath>


double integrate(
    const std::function<double(double)>& f,
    double a, double b,
    double acc = 1e-3,
    double eps = 1e-3,
    double f2 = std::numeric_limits<double>::quiet_NaN(),
    double f3 = std::numeric_limits<double>::quiet_NaN()
);

// Error function via integral representation using the adaptive integrator
// Uses three-branch formula for numerical stability:
//   erf(z<0)  = -erf(-z)
//   erf(0≤z≤1)= 2/√π ∫₀ᶻ exp(-x²) dx
//   erf(z>1)  = 1 - 2/√π ∫₀¹ exp(-(z+(1-t)/t)²)/t² dt
double my_erf(double z, double acc = 1e-6, double eps = 1e-6);

#endif 