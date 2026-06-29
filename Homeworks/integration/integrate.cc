#include "integrate.h"
#include <cmath>
#include <functional>
#include <limits>

// ──────────────────────────────────────────────────────────────────────────────
double integrate(
    const std::function<double(double)>& f,
    double a, double b,
    double acc,
    double eps,
    double f2,
    double f3)
{
    double h = b - a;

    // First call: no cached points yet
    if (std::isnan(f2)) {
        f2 = f(a + 2.0*h/6.0);
        f3 = f(a + 4.0*h/6.0);
    }

    double f1 = f(a +   h/6.0);
    double f4 = f(a + 5.0*h/6.0);

    // Higher-order quadrature rule
    double Q = (2.0*f1 + f2 + f3 + 2.0*f4) / 6.0 * h;
    // Embedded lower-order rule (error estimate)
    double q = (     f1 + f2 + f3 +      f4) / 4.0 * h;

    double err = std::abs(Q - q);
    double tol = acc + eps * std::abs(Q);

    if (err < tol) {
        return Q;
    } else {
        double mid = (a + b) / 2.0;
        if (mid == a || mid == b) {
            return Q;
        }
        double acc2 = acc / std::sqrt(2.0);  // distribute accuracy budget
        // Left half reuses f1, f2; right half reuses f3, f4
        return integrate(f, a,   mid, acc2, eps, f1, f2)
             + integrate(f, mid, b,   acc2, eps, f3, f4);
    }
}

// ──────────────────────────────────────────────────────────────────────────────
// Error function via integral representation
// ──────────────────────────────────────────────────────────────────────────────
double my_erf(double z, double acc, double eps)
{
    const double two_over_sqrtpi = 2.0 / std::sqrt(M_PI);

    if (z < 0.0) {
        return -my_erf(-z, acc, eps);
    } else if (z <= 1.0) {
        auto integrand = [](double x) { return std::exp(-x*x); };
        return two_over_sqrtpi * integrate(integrand, 0.0, z, acc, eps);
    } else {
        // erf(z) = 1 - 2/√π ∫₀¹ exp(-(z+(1-t)/t)²) / t² dt
        // Substitution maps the tail [z,∞) onto (0,1], avoiding the infinite limit.
        auto integrand = [z](double t) {
            double u = z + (1.0 - t) / t;
            return std::exp(-u*u) / (t*t);
        };
        return 1.0 - two_over_sqrtpi * integrate(integrand, 0.0, 1.0, acc, eps);
    }
}

// ──────────────────────────────────────────────────────────────────────────────
// Part B: Clenshaw–Curtis variable transformation
double integrate_cc(
    const std::function<double(double)>& f,
    double a, double b,
    double acc,
    double eps)
{
    double half = (b - a) / 2.0;
    auto transformed = [&f, a, b, half](double theta) {
        double x;
        if (theta < M_PI / 2.0) {
            // near theta=0 -> x near b; cos(theta) is well-conditioned here,
            // but compute x as b - half*(1-cos(theta)) to stay precise near x=b
            double s = std::sin(theta / 2.0);
            x = b - half * 2.0 * s * s;
        } else {
            // near theta=pi -> x near a; compute x as a + half*(1+cos(theta))
            double c = std::cos(theta / 2.0);
            x = a + half * 2.0 * c * c;
        }
        return f(x) * std::sin(theta) * half;
    };

    return integrate(transformed, 0.0, M_PI, acc, eps);
}

// ──────────────────────────────────────────────────────────────────────────────
// Part B: infinite/semi-infinite limits, built on integrate_cc
double integrate_inf(
    const std::function<double(double)>& f,
    double a, double b,
    double acc,
    double eps)
{
    bool a_inf = std::isinf(a);
    bool b_inf = std::isinf(b);

    if (a_inf && b_inf) {
        // x = t/(1-t^2), t in (-1, 1)
        auto transformed = [&f](double t) {
            double t2 = t * t;
            double denom = 1.0 - t2;
            double x = t / denom;
            double dxdt = (1.0 + t2) / (denom * denom);
            return f(x) * dxdt;
        };
        return integrate_cc(transformed, -1.0, 1.0, acc, eps);
    } else if (b_inf) {
        // [a, +inf): x = a + t/(1-t), t in [0,1)
        auto transformed = [&f, a](double t) {
            double denom = 1.0 - t;
            double x = a + t / denom;
            double dxdt = 1.0 / (denom * denom);
            return f(x) * dxdt;
        };
        return integrate_cc(transformed, 0.0, 1.0, acc, eps);
    } else if (a_inf) {
        // (-inf, b]: x = b - t/(1-t), t in [0,1)
        auto transformed = [&f, b](double t) {
            double denom = 1.0 - t;
            double x = b - t / denom;
            double dxdt = 1.0 / (denom * denom);
            return f(x) * dxdt;
        };
        return integrate_cc(transformed, 0.0, 1.0, acc, eps);
    } else {
        // Both finite: no transformation needed
        return integrate_cc(f, a, b, acc, eps);
    }
}