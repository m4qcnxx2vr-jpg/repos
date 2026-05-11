#pragma once
#include "matrix.h"
#include <tuple>
#include <functional>
#include <vector>

namespace pp {

std::tuple<std::vector<double>, std::vector<vector>> driver(
    std::function<vector(double, const vector&)> f,
    double a,
    double b,
    const vector& yinit,
    double h = 0.125,
    double acc = 0.01,
    double eps = 0.01
);

}