#include "minimize.h"
#include "qr.h"
#include <cmath>
#include <algorithm>

namespace pp {

vector numerical_gradient(scalar_function phi, vector x, bool central) {
    int n = x.size();
    vector g(n);

    if (!central) {
        double phix = phi(x);
        double scale = std::pow(2.0, -26.0); // eps^(1/2)

        for (int i = 0; i < n; ++i) {
            double dxi = (1.0 + std::abs(x[i])) * scale;
            x[i] += dxi;
            g[i] = (phi(x) - phix) / dxi;
            x[i] -= dxi;
        }
    } else {
        double scale = std::pow(2.0, -52.0 / 3.0); // eps^(1/3)

        for (int i = 0; i < n; ++i) {
            double dxi = (1.0 + std::abs(x[i])) * scale;

            x[i] += dxi;
            double fplus = phi(x);
            x[i] -= 2.0 * dxi;
            double fminus = phi(x);
            x[i] += dxi;

            g[i] = (fplus - fminus) / (2.0 * dxi);
        }
    }

    return g;
}

matrix numerical_hessian(scalar_function phi, vector x, bool central) {
    int n = x.size();
    matrix H(n, n);

    if (!central) {
        double scale = std::pow(2.0, -13.0); // eps^(1/4)
        vector gx = numerical_gradient(phi, x, false);

        for (int j = 0; j < n; ++j) {
            double dxj = (1.0 + std::abs(x[j])) * scale;
            x[j] += dxj;
            vector dg = numerical_gradient(phi, x, false) - gx;
            x[j] -= dxj;

            for (int i = 0; i < n; ++i)
                H(i, j) = dg[i] / dxj;
        }
    } else {
        double scale = std::pow(2.0, -52.0 / 3.0); // eps^(1/3)
        double phix = phi(x);

        // Diagonal terms: central second-difference, one extra eval each.
        for (int i = 0; i < n; ++i) {
            double dxi = (1.0 + std::abs(x[i])) * scale;

            x[i] += dxi;
            double fplus = phi(x);
            x[i] -= 2.0 * dxi;
            double fminus = phi(x);
            x[i] += dxi;

            H(i, i) = (fplus - 2.0 * phix + fminus) / (dxi * dxi);
        }

        // Off-diagonal terms: mixed central differences.
        for (int j = 0; j < n; ++j) {
            double dxj = (1.0 + std::abs(x[j])) * scale;

            for (int i = j + 1; i < n; ++i) {
                double dxi = (1.0 + std::abs(x[i])) * scale;

                x[i] += dxi; x[j] += dxj;
                double fpp = phi(x);
                x[j] -= 2.0 * dxj;
                double fpm = phi(x);
                x[i] -= 2.0 * dxi;
                double fmm = phi(x);
                x[j] += 2.0 * dxj;
                double fmp = phi(x);
                x[i] += dxi; x[j] -= dxj;

                double hij = (fpp - fpm - fmp + fmm) / (4.0 * dxi * dxj);
                H(i, j) = hij;
                H(j, i) = hij;
            }
        }
    }

    return H;
}

vector newton_minimize(
    scalar_function phi,
    vector x,
    double acc,
    int max_iter,
    bool central,
    int* steps
) {
    int n = x.size();
    int iter = 0;
    double mu = 1e-6; // Levenberg regularization, adapted as we go

    for (; iter < max_iter; ++iter) {
        vector g = numerical_gradient(phi, x, central);
        if (g.norm() < acc) break;

        double phix = phi(x);
        bool improved = false;
        vector xnew = x;

        // Try increasingly strong regularization until a downhill step is found.
        for (int tries = 0; tries < 30; ++tries) {
            matrix H = numerical_hessian(phi, x, central);
            for (int i = 0; i < n; ++i) H(i, i) += mu;

            vector dx;
            try {
                qr QRH(H);
                dx = QRH.solve(-g);
            } catch (...) {
                mu *= 10.0;
                continue;
            }

            double lambda = 1.0;
            while (lambda >= 1.0 / 1024.0) {
                vector z = x + lambda * dx;
                if (phi(z) < phix) {
                    xnew = z;
                    improved = true;
                    break;
                }
                lambda /= 2.0;
            }

            if (improved) {
                mu = std::max(mu / 5.0, 1e-10);
                break;
            }

            mu *= 10.0;
        }

        if (!improved) {
            // Could not find a downhill step even with strong regularization;
            // take the smallest allowed step anyway and keep going.
            matrix H = numerical_hessian(phi, x, central);
            for (int i = 0; i < n; ++i) H(i, i) += mu;
            qr QRH(H);
            vector dx = QRH.solve(-g);
            xnew = x + (1.0 / 1024.0) * dx;
        }

        x = xnew;
    }

    if (steps) *steps = iter;
    return x;
}

double rosenbrock(const vector& x) {
    double X = x[0];
    double Y = x[1];
    return (1.0 - X) * (1.0 - X) + 100.0 * (Y - X * X) * (Y - X * X);
}

double himmelblau(const vector& x) {
    double X = x[0];
    double Y = x[1];
    double a = X * X + Y - 11.0;
    double b = X + Y * Y - 7.0;
    return a * a + b * b;
}

double breit_wigner(double E, double m, double Gamma, double A) {
    return A / ((E - m) * (E - m) + Gamma * Gamma / 4.0);
}

}
