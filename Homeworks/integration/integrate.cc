#include "integrate.h"
#include <cmath>
#include <functional>
#include <limits>

// ──────────────────────────────────────────────────────────────────────────────
// Recursive open 4-point adaptive integrator
//
// Higher-order rule (degree-4 precision, 4 points):
//   Q = (2*f1 + f2 + f3 + 2*f4) / 6 * (b-a)
//
// Embedded lower-order rule (degree-2 precision, same 4 points):
//   q = (f1 + f2 + f3 + f4) / 4 * (b-a)
//
// Points at x = a + h/6, a + 2h/6, a + 4h/6, a + 5h/6  (open: avoid endpoints)
// On subdivision, f2 becomes f1 of the right child, f3 becomes f4 of the left
// child — so only 2 new evaluations are needed per recursive split instead of 4.
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
        // erf(z) = 2/√π ∫₀ᶻ exp(-x²) dx
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