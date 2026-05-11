#include "solver.h"
#include "stepper.h"
#include <cmath>
#include <algorithm>

namespace pp {

std::tuple<std::vector<double>, std::vector<vector>> driver(
    std::function<vector(double, const vector&)> f,
    double a,
    double b,
    const vector& yinit,
    double h,
    double acc,
    double eps
) {
    double x = a;
    vector y = yinit;

    std::vector<double> xlist = {x};
    std::vector<vector> ylist = {y};

    while (true) {
        if (x >= b) return {xlist, ylist};

        if (x + h > b) h = b - x;

        auto [yh, dy] = rkstep12(f, x, y, h);

        double tol = (acc + eps * yh.norm()) * std::sqrt(h / (b - a));
        double err = dy.norm();

        if (err <= tol) {
            x += h;
            y = yh;
            xlist.push_back(x);
            ylist.push_back(y);
        }

        if (err > 0) {
            double factor = std::pow(tol / err, 0.25) * 0.95;
            factor = std::min(factor, 2.0);
            h *= factor;
        } else {
            h *= 2;
        }
    }
}

}