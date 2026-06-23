// main.cc -- Homework "Artificial Neural Networks", parts A and B.
//
// Trains a 3-layer ANN (1 input neuron, n hidden neurons, 1 summation
// output neuron) to interpolate a tabulated function g(x), then writes
//   - the trained parameters (a_i, b_i, w_i) to params.txt
//   - a dense table x, g(x), F_p(x), F_p'(x), F_p''(x), antideriv(x)
//     to curve.txt, suitable for plotting (e.g. with gnuplot/python).
//
// Usage:
//   ./main                       (defaults: g(x)=cos(5x-1)exp(-x^2), n=9, [-1,1])
//   ./main n_hidden               (choose number of hidden neurons)
//   ./main n_hidden npoints       (choose number of training points too)

#include "ann.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>

using pp::vector;
using pp::ann;

// The function we train the network to approximate.
static double g(double x) {
    return std::cos(5.0 * x - 1.0) * std::exp(-x * x);
}

int main(int argc, char** argv) {
    int n_hidden = (argc > 1) ? std::atoi(argv[1]) : 9;
    int N        = (argc > 2) ? std::atoi(argv[2]) : 41;
    double xmin = -1.0, xmax = 1.0;

    if (n_hidden < 1 || N < 2) {
        std::cerr << "usage: " << argv[0] << " [n_hidden] [n_points]\n";
        return 1;
    }

    // ---- build the training table {x_k, y_k} ----
    vector xs(N, 0.0), ys(N, 0.0);
    for (int k = 0; k < N; k++) {
        double xk = xmin + (xmax - xmin) * k / (N - 1);
        xs[k] = xk;
        ys[k] = g(xk);
    }

    // ---- build and train the network (part A) ----
    ann net(n_hidden);
    net.init(xmin, xmax, /*seed=*/1);

    int steps = 0;
    double C = net.train(xs, ys, /*acc=*/1e-10, /*max_iter=*/2000, &steps);

    std::cout << "Trained network: " << n_hidden << " hidden neurons, "
              << N << " training points on [" << xmin << ", " << xmax << "]\n";
    std::cout << "  converged in " << steps << " iterations, final cost C(p) = "
              << C << "\n";

    double max_err = 0.0;
    for (int k = 0; k < N; k++)
        max_err = std::max(max_err, std::fabs(net.response(xs[k]) - ys[k]));
    std::cout << "  max |F_p(x_k)-y_k| at training points: " << max_err << "\n";

    // ---- write trained parameters to params.txt ----
    {
        std::ofstream out("params.txt");
        out << "# trained ann parameters: n_hidden = " << n_hidden << "\n";
        out << "# i  a_i  b_i  w_i\n";
        for (int i = 0; i < n_hidden; i++) {
            out << i << "  " << net.p[3*i+0] << "  " << net.p[3*i+1]
                << "  " << net.p[3*i+2] << "\n";
        }
        std::cout << "  wrote params.txt\n";
    }

    // ---- write a dense table for plotting (part A + part B) ----
    {
        std::ofstream out("curve.txt");
        out << "# x  g(x)  F_p(x)  F_p'(x)  F_p''(x)  antideriv(x)\n";
        int M = 400;
        for (int k = 0; k <= M; k++) {
            double x = xmin + (xmax - xmin) * k / M;
            out << x << "  " << g(x) << "  " << net.response(x) << "  "
                << net.derivative(x) << "  " << net.second_derivative(x) << "  "
                << net.antiderivative(x) << "\n";
        }
        std::cout << "  wrote curve.txt (" << (M+1) << " points)\n";
    }

    // ---- a quick interactive-style report: response at a few points ----
    std::cout << "\nSample evaluations:\n";
    std::printf("%8s %12s %12s %12s\n", "x", "F_p(x)", "F_p'(x)", "F_p''(x)");
    for (int k = 0; k <= 8; k++) {
        double x = xmin + (xmax - xmin) * k / 8.0;
        std::printf("%8.3f %12.6f %12.6f %12.6f\n",
                    x, net.response(x), net.derivative(x), net.second_derivative(x));
    }

    return 0;
}
