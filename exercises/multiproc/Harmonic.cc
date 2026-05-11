#include "Harmonic.h"

#include <cmath>

void harm(Data& arg) {
    int a = arg.a;
    int b = arg.b;

    double sum = 0.0;

    for (int i = a; i < b; i++) {
        sum += 1.0 / i;
    }

    arg.sum = sum;
}

double harmonic_asymptotic(int n) {
    const double gamma = 0.5772156649015328606;

    double x = static_cast<double>(n);

    return std::log(x)
         + gamma
         + 1.0 / (2.0 * x)
         - 1.0 / (12.0 * x * x);
}

bool approx(double x, double y, double tolerance) {
    return std::abs(x - y) < tolerance;
}