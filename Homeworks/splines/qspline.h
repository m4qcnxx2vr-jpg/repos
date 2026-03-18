#pragma once 
#include <vector>

struct qspline {
    using vec = std::vector<double>;

    int n;
    vec x, y, b, c;

    qspline(const vec& xs, const vec& ys);

    int binsearch(double z) const;
    double eval(double z) const;
    double deriv(double z) const;
    double integ(double z) const;
};