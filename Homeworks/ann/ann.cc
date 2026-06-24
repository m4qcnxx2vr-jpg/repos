#include "ann.h"
#include "minimize.h"
#include <cmath>
#include <iostream>

namespace pp {

ann::ann(int n) : n(n), p(3 * n) {
    for (int i = 0; i < n; ++i) {
        p[3 * i + 0] = -1.0 + 2.0 * i / (n - 1.0); // a_i
        p[3 * i + 1] = std::log(0.5);               // log(b_i)
        p[3 * i + 2] = 1.0;                         // w_i
    }
}

double ann::f(double z) const {
    return z * std::exp(-z * z);
}

double ann::df(double z) const {
    return (1.0 - 2.0 * z * z) * std::exp(-z * z);
}

double ann::ddf(double z) const {
    return (4.0 * z * z * z - 6.0 * z) * std::exp(-z * z);
}

double ann::F(double z) const {
    return -0.5 * std::exp(-z * z);
}

double ann::response(double x) const {
    return response(x, p);
}

double ann::response(double x, const vector& params) const {
    double sum = 0.0;

    for (int i = 0; i < n; ++i) {
        double a = params[3 * i + 0];
        double b = std::exp(params[3 * i + 1]);
        double w = params[3 * i + 2];

        double z = (x - a) / b;
        sum += w * f(z);
    }

    return sum;
}

double ann::derivative(double x) const {
    double sum = 0.0;

    for (int i = 0; i < n; ++i) {
        double a = p[3 * i + 0];
        double b = std::exp(p[3 * i + 1]);
        double w = p[3 * i + 2];

        double z = (x - a) / b;
        sum += w * df(z) / b;
    }

    return sum;
}

double ann::second_derivative(double x) const {
    double sum = 0.0;

    for (int i = 0; i < n; ++i) {
        double a = p[3 * i + 0];
        double b = std::exp(p[3 * i + 1]);
        double w = p[3 * i + 2];

        double z = (x - a) / b;
        sum += w * ddf(z) / (b * b);
    }

    return sum;
}

double ann::antiderivative(double x) const {
    double sum = 0.0;

    for (int i = 0; i < n; ++i) {
        double a = p[3 * i + 0];
        double b = std::exp(p[3 * i + 1]);
        double w = p[3 * i + 2];

        double z = (x - a) / b;
        sum += w * b * F(z);
    }

    return sum;
}

void ann::train(const vector& xs, const vector& ys) {
    int m = xs.size();

    double xmin = xs[0];
    double xmax = xs[m - 1];
    double width = (xmax - xmin) / n;

    for (int i = 0; i < n; ++i) {
        double a = xmin + (xmax - xmin) * i / (n - 1.0);

        p[3 * i + 0] = a;
        p[3 * i + 1] = std::log(width);
        p[3 * i + 2] = 1.0;
    }

    auto cost = [&](const vector& params) {
        double sum = 0.0;

        for (int k = 0; k < m; ++k) {
            double diff = response(xs[k], params) - ys[k];
            sum += diff * diff;
        }

        return sum;
    };

    int steps = 0;
    p = newton_minimize(cost, p, 1e-4, 200, false, &steps);

    std::cerr << "training steps = " << steps << "\n";
    std::cerr << "final cost = " << cost(p) << "\n";
}

}