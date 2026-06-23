#include "roots.h"
#include "qr.h"
#include <cmath>
#include <algorithm>
#include <vector>

namespace pp {

static vector default_dx(const vector& x) {
    vector dx(x.size());
    double scale = std::pow(2.0, -26);

    for (int i = 0; i < x.size(); ++i) {
        dx[i] = std::max(std::abs(x[i]), 1.0) * scale;
    }

    return dx;
}

static void numerical_jacobian(
    vector_function f,
    vector x,
    const vector& fx,
    matrix& J
) {
    int n = x.size();

    vector dx = default_dx(x);

    for (int j = 0; j < n; ++j) {
        x[j] += dx[j];
        vector df = f(x) - fx;

        for (int i = 0; i < n; ++i) {
            J(i, j) = df[i] / dx[j];
        }

        x[j] -= dx[j];
    }
}

static bool step_too_small(const vector& Dx, const vector& x) {
    vector dx = default_dx(x);

    for (int i = 0; i < x.size(); ++i) {
        if (std::abs(Dx[i]) > dx[i]) return false;
    }

    return true;
}

vector newton(
    vector_function f,
    vector x,
    double acc,
    double alpha_min,
    int max_iter
) {
    vector fx = f(x);

    matrix J(x.size(), x.size());

    for (int iter = 0; iter < max_iter; ++iter) {
        if (fx.norm() < acc) break;

        numerical_jacobian(f, x, fx, J);

        qr QRJ(J);
        vector Dx = QRJ.solve(-fx);

        if (step_too_small(Dx, x)) break;

        double alpha = 1.0;
        vector z;
        vector fz;

        while (true) {
            z = x + alpha * Dx;
            fz = f(z);

            if (fz.norm() < fx.norm()) break;
            if (alpha < alpha_min) break;

            alpha /= 2.0;
        }

        x = z;
        fx = fz;
    }

    return x;
}
vector newton_quadratic(
    vector_function f,
    vector x,
    double acc,
    double alpha_min,
    int max_iter
) {
    vector fx = f(x);

    matrix J(x.size(), x.size());

    for (int iter = 0; iter < max_iter; ++iter) {
        if (fx.norm() < acc) break;

        numerical_jacobian(f, x, fx, J);

        qr QRJ(J);
        vector Dx = QRJ.solve(-fx);

        if (step_too_small(Dx, x)) break;

        double lambda = 1.0;
        double phi0 = fx.norm();

        vector z = x + lambda * Dx;
        vector fz = f(z);
        double phi_lambda = fz.norm();

        while (
            phi_lambda >= (1.0 - lambda / 2.0) * phi0
            && lambda > alpha_min
        ) {
            double denominator =
                2.0 * (phi_lambda - phi0 + lambda * phi0 / 2.0);

            double lambda_new;

            if (std::abs(denominator) < 1e-14) {
                lambda_new = lambda / 2.0;
            } else {
                lambda_new = -lambda * lambda * phi0 / denominator;
            }

            if (lambda_new < lambda / 4.0 || lambda_new > lambda / 2.0) {
                lambda_new = lambda / 2.0;
            }

            lambda = lambda_new;

            z = x + lambda * Dx;
            fz = f(z);
            phi_lambda = fz.norm();
        }

        x = z;
        fx = fz;
    }

    return x;
}
vector gradient_rosenbrock(const vector& x) {
    double X = x[0];
    double Y = x[1];

    vector g(2);

    g[0] = -2.0 * (1.0 - X) - 400.0 * X * (Y - X * X);
    g[1] = 200.0 * (Y - X * X);

    return g;
}

vector gradient_himmelblau(const vector& x) {
    double X = x[0];
    double Y = x[1];

    vector g(2);

    g[0] = 4.0 * X * (X * X + Y - 11.0)
         + 2.0 * (X + Y * Y - 7.0);

    g[1] = 2.0 * (X * X + Y - 11.0)
         + 4.0 * Y * (X + Y * Y - 7.0);

    return g;
}

// ---------- small RK12 ODE solver for hydrogen ----------

static vector rkstep12(
    std::function<vector(double, const vector&)> F,
    double r,
    const vector& y,
    double h,
    vector& err
) {
    vector k0 = F(r, y);
    vector yt = y + h * k0;
    vector k1 = F(r + h, yt);

    vector yh(y.size());
    err.resize(y.size());

    for (int i = 0; i < y.size(); ++i) {
        yh[i] = y[i] + h * 0.5 * (k0[i] + k1[i]);
        err[i] = yh[i] - yt[i];
    }

    return yh;
}

static vector ode_driver(
    std::function<vector(double, const vector&)> F,
    double a,
    vector ya,
    double b,
    double h,
    double acc,
    double eps
) {
    double r = a;
    vector y = ya;

    while (r < b) {
        if (r + h > b) h = b - r;

        vector err(y.size());
        vector yh = rkstep12(F, r, y, h, err);

        double tol = (acc + eps * yh.norm()) * std::sqrt(h / (b - a));
        double err_norm = err.norm();

        if (err_norm <= tol) {
            r += h;
            y = yh;
        }

        if (err_norm > 0) {
            h *= std::min(std::pow(tol / err_norm, 0.25) * 0.95, 2.0);
        } else {
            h *= 2.0;
        }
    }

    return y;
}

double hydrogen_mismatch(double E, double rmin, double rmax, double acc, double eps) {
    auto F = [E](double r, const vector& y) {
        vector dydr(2);

        dydr[0] = y[1];
        dydr[1] = -2.0 * (E + 1.0 / r) * y[0];

        return dydr;
    };

    vector y0(2);
    y0[0] = rmin - rmin * rmin;
    y0[1] = 1.0 - 2.0 * rmin;

    vector yfinal = ode_driver(F, rmin, y0, rmax, 0.01, acc, eps);

    return yfinal[0];
}

std::vector<hydrogen_point> hydrogen_wavefunction(
    double E, double rmin, double rmax, double acc, double eps
) {
    auto F = [E](double r, const vector& y) {
        vector dydr(2);
        dydr[0] = y[1];
        dydr[1] = -2.0 * (E + 1.0 / r) * y[0];
        return dydr;
    };

    vector y0(2);
    y0[0] = rmin - rmin * rmin;
    y0[1] = 1.0 - 2.0 * rmin;

    std::vector<hydrogen_point> points;
    points.push_back({rmin, y0[0], y0[1]});

    double r = rmin;
    vector y = y0;
    double h = 0.01;
    double a = rmin, b = rmax;

    while (r < b) {
        if (r + h > b) h = b - r;

        vector err(y.size());
        vector yh = rkstep12(F, r, y, h, err);

        double tol = (acc + eps * yh.norm()) * std::sqrt(h / (b - a));
        double err_norm = err.norm();

        if (err_norm <= tol) {
            r += h;
            y = yh;
            points.push_back({r, y[0], y[1]});
        }

        if (err_norm > 0) {
            h *= std::min(std::pow(tol / err_norm, 0.25) * 0.95, 2.0);
        } else {
            h *= 2.0;
        }
    }

    return points;
}

double find_hydrogen_energy(double rmin, double rmax, double acc, double eps) {
    auto M = [=](const vector& x) {
        vector y(1);
        y[0] = hydrogen_mismatch(x[0], rmin, rmax, acc, eps);
        return y;
    };

    vector start(1);
    start[0] = -0.7;

    vector root = newton(M, start, 1e-6, 1e-3, 100);

    return root[0];
}

} // namespace pp