#include "ann.h"
#include <cmath>
#include <cstdio>
#include <random>

using pp::vector;
using pp::ann;

double g(double x) {
    return std::cos(5.0 * x - 1.0) * std::exp(-x * x);
}
double dg(double x) {
    // d/dx [cos(5x-1) exp(-x^2)] = -5 sin(5x-1) exp(-x^2) - 2x cos(5x-1) exp(-x^2)
    return -5.0 * std::sin(5.0 * x - 1.0) * std::exp(-x * x)
           - 2.0 * x * std::cos(5.0 * x - 1.0) * std::exp(-x * x);
}

int main() {
    // ---------- Part A: train network to approximate g(x) ----------
    int N = 41;
    vector xs(N, 0.0), ys(N, 0.0);
    for (int k = 0; k < N; k++) {
        double xk = -1.0 + 2.0 * k / (N - 1);
        xs[k] = xk;
        ys[k] = g(xk);
    }

    int n_hidden = 9;
    ann net(n_hidden);
    net.init(-1.0, 1.0, /*seed=*/42);

    int steps = 0;
    double C0 = 0.0;
    for (int k = 0; k < N; k++) { double r = net.response(xs[k]) - ys[k]; C0 += r*r; }

    double Cfinal = net.train(xs, ys, 1e-10, 500, &steps);

    std::printf("Part A: n_hidden=%d, points=%d\n", n_hidden, N);
    std::printf("  initial cost = %.6e\n", C0);
    std::printf("  final cost   = %.6e  (after %d iterations)\n", Cfinal, steps);

    std::printf("\n  x        g(x)        F_p(x)      |error|\n");
    for (int k = 0; k < 9; k++) {
        double xk = -1.0 + 2.0 * k / 8.0;
        double yk = g(xk);
        double Fk = net.response(xk);
        std::printf("  %6.3f  %10.6f  %10.6f  %10.3e\n", xk, yk, Fk, std::fabs(Fk - yk));
    }

    // out-of-sample check (between training points)
    double max_err = 0.0;
    for (int k = 0; k < 200; k++) {
        double xk = -1.0 + 2.0 * k / 199.0;
        max_err = std::max(max_err, std::fabs(net.response(xk) - g(xk)));
    }
    std::printf("\n  max |F_p(x)-g(x)| on dense grid [-1,1]: %.4e\n", max_err);

    // ---------- Part B: derivatives and antiderivative ----------
    std::printf("\nPart B: derivatives at sample points\n");
    std::printf("  x        F'(x)       g'(x)       |error|\n");
    double max_derr = 0.0;
    for (int k = 0; k < 9; k++) {
        double xk = -1.0 + 2.0 * k / 8.0;
        double Fp = net.derivative(xk);
        double gp = dg(xk);
        std::printf("  %6.3f  %10.6f  %10.6f  %10.3e\n", xk, Fp, gp, std::fabs(Fp - gp));
        max_derr = std::max(max_derr, std::fabs(Fp - gp));
    }
    std::printf("  max derivative error (9-pt sample): %.4e\n", max_derr);

    // Check antiderivative via numerical integration consistency:
    // antiderivative(b) - antiderivative(a) should equal the definite
    // integral of F_p over [a,b]. We verify against a Simpson estimate
    // of the integral of F_p itself (self-consistency of analytic formula).
    double a_pt = -0.5, b_pt = 0.5;
    double exact_FTC = net.antiderivative(b_pt) - net.antiderivative(a_pt);

    int M = 2000; // Simpson's rule, must be even
    double h = (b_pt - a_pt) / M;
    double simpson = net.response(a_pt) + net.response(b_pt);
    for (int i = 1; i < M; i++) {
        double xi = a_pt + i * h;
        simpson += (i % 2 == 0 ? 2.0 : 4.0) * net.response(xi);
    }
    simpson *= h / 3.0;

    std::printf("\n  Antiderivative check on [-0.5, 0.5]:\n");
    std::printf("    F_antideriv(b)-F_antideriv(a) = %.8f\n", exact_FTC);
    std::printf("    Simpson integral of F_p(x)dx  = %.8f\n", simpson);
    std::printf("    difference                    = %.3e\n", std::fabs(exact_FTC - simpson));

    return 0;
}
