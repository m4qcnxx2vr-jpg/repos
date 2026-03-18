#pragma once
#include <vector>
#include <functional>

std::function<double(double)> make_qspline(
    std::vector<double> x,
    std::vector<double> y
);