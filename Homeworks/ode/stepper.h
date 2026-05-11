#pragma once
#include "matrix.h"   // vector is defined in matrix.h from previous homeworks
#include <tuple>
#include <functional>

namespace pp {

std::tuple<vector, vector> rkstep12(
    std::function<vector(double, const vector&)> f,
    double x,
    const vector& y,
    double h
);

}