#include "SpecialFunctions.h"

#include <cmath>
#include <limits>
#include <vector>

namespace {
    constexpr double PI = 3.14159265358979323846;
}

double erf_approx(double x) {
    if (x < 0.0) {
        return -erf_approx(-x);
    }

    std::vector<double> a {
        0.254829592,
       -0.284496736,
        1.421413741,
       -1.453152027,
        1.061405429
    };

    double t = 1.0 / (1.0 + 0.3275911 * x);

    double sum = t * (
        a[0] + t * (
        a[1] + t * (
        a[2] + t * (
        a[3] + t * a[4]))));

    return 1.0 - sum * std::exp(-x * x);
}

double sgamma(double x) {
    if (x < 0.0) {
        return PI / std::sin(PI * x) / sgamma(1.0 - x);
    }

    if (x < 9.0) {
        return sgamma(x + 1.0) / x;
    }

    double lnsgamma =
        std::log(2.0 * PI) / 2.0
        + (x - 0.5) * std::log(x)
        - x
        + (1.0 / 12.0) / x
        - (1.0 / 360.0) / (x * x * x)
        + (1.0 / 1260.0) / (x * x * x * x * x);

    return std::exp(lnsgamma);
}

double lngamma(double x) {
    if (x <= 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    if (x < 9.0) {
        return lngamma(x + 1.0) - std::log(x);
    }

    return x * std::log(x + 1.0 / (12.0 * x - 1.0 / x / 10.0))
         - x
         + std::log(2.0 * PI / x) / 2.0;
}

double factorial(int n) {
    double result = 1.0;

    for (int i = 2; i <= n; i++) {
        result *= i;
    }

    return result;
}

bool approx(double x, double y, double tolerance) {
    return std::abs(x - y) < tolerance;
}