#include "approx.h"
#include <cmath>
#include <complex>

bool approx(double a, double b, double acc, double eps) {
    if (std::abs(a - b) < acc) return true;
    if (std::abs(a - b) < eps * (std::abs(a) + std::abs(b))) return true;
    return false;
}

bool approx(std::complex<double> a, std::complex<double> b, double acc, double eps) {
    return approx(a.real(), b.real(), acc, eps)
        && approx(a.imag(), b.imag(), acc, eps);
}