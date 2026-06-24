#pragma once
#include "matrix.h"

namespace pp {

struct ann {
    int n;
    vector p;

    ann(int n);

    double f(double z) const;
    double df(double z) const;
    double ddf(double z) const;
    double F(double z) const;

    double response(double x) const;
    double response(double x, const vector& params) const;

    double derivative(double x) const;
    double second_derivative(double x) const;
    double antiderivative(double x) const;

    void train(const vector& xs, const vector& ys);
};

}