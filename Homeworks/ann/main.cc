#include "ann.h"
#include <cmath>
#include <fstream>
#include <iostream>

double g(double x) {
    return std::cos(5.0 * x - 1.0) * std::exp(-x * x);
}

double dg(double x) {
    double u = 5.0 * x - 1.0;
    return std::exp(-x * x) *
           (-5.0 * std::sin(u) - 2.0 * x * std::cos(u));
}

double ddg(double x) {
    double u = 5.0 * x - 1.0;
    return std::exp(-x * x) *
           ((4.0 * x * x - 27.0) * std::cos(u)
            + 20.0 * x * std::sin(u));
}

int main() {
    const int ndata = 40;
    const int nneurons = 4;

    const double xmin = -1.0;
    const double xmax = 1.0;

    pp::vector xs(ndata);
    pp::vector ys(ndata);

    for (int i = 0; i < ndata; ++i) {
        xs[i] = xmin + (xmax - xmin) * i / (ndata - 1.0);
        ys[i] = g(xs[i]);
    }

    pp::ann network(nneurons);
    network.train(xs, ys);

    std::ofstream data("ann.tsv");

    data << "# x "
         << "g "
         << "ann "
         << "dg "
         << "ann_dg "
         << "ddg "
         << "ann_ddg "
         << "ann_antiderivative\n";

    double A0 = network.antiderivative(xmin);

    for (int i = 0; i <= 400; ++i) {
        double x = xmin + (xmax - xmin) * i / 400.0;

        data << x << " "
             << g(x) << " "
             << network.response(x) << " "
             << dg(x) << " "
             << network.derivative(x) << " "
             << ddg(x) << " "
             << network.second_derivative(x) << " "
             << network.antiderivative(x) - A0
             << "\n";
    }

    data.close();

    std::cout << "ANN exercise A+B\n";
    std::cout << "neurons = " << nneurons << "\n";
    std::cout << "training points = " << ndata << "\n";
    std::cout << "data written to ann.tsv\n";
    std::cout << "plots written to ann_fit.svg, ann_derivative.svg, ann_second_derivative.svg, ann_antiderivative.svg\n";

    return 0;
}