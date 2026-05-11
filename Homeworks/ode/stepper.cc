#include "stepper.h"

namespace pp {

std::tuple<vector, vector> rkstep12(
    std::function<vector(double, const vector&)> f,
    double x,
    const vector& y,
    double h
) {
    vector k0 = f(x, y);
    vector k1 = f(x + h/2, y + k0 * (h/2));

    vector yh = y + k1 * h;
    vector dy = (k1 - k0) * h;

    return {yh, dy};
}

}