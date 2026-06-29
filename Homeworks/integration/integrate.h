#ifndef INTEGRATE_H
#define INTEGRATE_H

#include <functional>
#include <cmath>
#include <limits>

double integrate(
    const std::function<double(double)>& f,
    double a,
    double b,
    double acc = 1e-3,
    double eps = 1e-3,
    double f2 = std::numeric_limits<double>::quiet_NaN(),
    double f3 = std::numeric_limits<double>::quiet_NaN()
);

double my_erf(
    double z,
    double acc = 1e-6,
    double eps = 1e-6
);

double integrate_open(
    const std::function<double(double)>& f,
    double a,
    double b,
    double acc = 1e-3,
    double eps = 1e-3
);

double integrate_cc(
    const std::function<double(double)>& f,
    double a,
    double b,
    double acc = 1e-3,
    double eps = 1e-3
);

double integrate_inf(
    const std::function<double(double)>& f,
    double a,
    double b,
    double acc = 1e-6,
    double eps = 1e-6
);

#endif