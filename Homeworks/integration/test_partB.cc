#include "integrate.h"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>
#include <limits>

struct CountedF {
    std::function<double(double)> f;
    mutable int ncalls = 0;
    double operator()(double x) const { ++ncalls; return f(x); }
};

void compare(const std::string& label,
             const std::function<double(double)>& f,
             double a, double b, double exact,
             double acc = 1e-6, double eps = 1e-6)
{
    CountedF cf1{f};
    double r1 = integrate(std::ref(cf1), a, b, acc, eps);

    CountedF cf2{f};
    double r2 = integrate_cc(std::ref(cf2), a, b, acc, eps);

    std::cout << std::left << std::setw(28) << label << "\n";
    std::cout << "    plain:  result=" << std::setprecision(10) << r1
              << "  |err|=" << std::scientific << std::setprecision(2) << std::abs(r1-exact)
              << "  ncalls=" << std::fixed << cf1.ncalls << "\n";
    std::cout << "    CC:     result=" << std::setprecision(10) << r2
              << "  |err|=" << std::scientific << std::setprecision(2) << std::abs(r2-exact)
              << "  ncalls=" << std::fixed << cf2.ncalls << "\n\n";
}

void test_inf(const std::string& label,
               const std::function<double(double)>& f,
               double a, double b, double exact,
               double acc = 1e-6, double eps = 1e-6)
{
    CountedF cf{f};
    double r = integrate_inf(std::ref(cf), a, b, acc, eps);
    std::cout << std::left << std::setw(28) << label
              << "  result=" << std::setprecision(10) << r
              << "  exact=" << exact
              << "  |err|=" << std::scientific << std::setprecision(2) << std::abs(r-exact)
              << "  ncalls=" << std::fixed << cf.ncalls << "\n";
}

int main() {
    const double inf = std::numeric_limits<double>::infinity();

    std::cout << std::string(90,'=') << "\n";
    std::cout << "Part B: plain integrator vs Clenshaw-Curtis on endpoint singularities\n";
    std::cout << std::string(90,'=') << "\n\n";

    compare("∫0^1 1/sqrt(x) dx = 2",
            [](double x){ return 1.0/std::sqrt(x); }, 0.0, 1.0, 2.0);

    compare("∫0^1 ln(x)/sqrt(x) dx = -4",
            [](double x){ return std::log(x)/std::sqrt(x); }, 0.0, 1.0, -4.0);

    compare("∫0^1 1/sqrt(1-x) dx = 2 (other endpoint)",
            [](double x){ return 1.0/std::sqrt(1.0-x); }, 0.0, 1.0, 2.0);

    compare("∫0^1 sqrt(x) dx = 2/3 (smooth, sanity check)",
            [](double x){ return std::sqrt(x); }, 0.0, 1.0, 2.0/3.0);

    std::cout << std::string(90,'=') << "\n";
    std::cout << "Part B: infinite-limit integrals\n";
    std::cout << std::string(90,'=') << "\n\n";

    // Gaussian: ∫_{-inf}^{inf} exp(-x^2) dx = sqrt(pi)
    test_inf("∫_-inf^inf exp(-x²) dx = √π",
              [](double x){ return std::exp(-x*x); }, -inf, inf, std::sqrt(M_PI));

    // ∫_0^inf exp(-x) dx = 1
    test_inf("∫0^inf exp(-x) dx = 1",
              [](double x){ return std::exp(-x); }, 0.0, inf, 1.0);

    // ∫_0^inf 1/(1+x^2) dx = pi/2
    test_inf("∫0^inf 1/(1+x²) dx = π/2",
              [](double x){ return 1.0/(1.0+x*x); }, 0.0, inf, M_PI/2.0);

    // ∫_-inf^0 exp(x) dx = 1
    test_inf("∫-inf^0 exp(x) dx = 1",
              [](double x){ return std::exp(x); }, -inf, 0.0, 1.0);

    // ∫_-inf^inf 1/(1+x^2) dx = pi
    test_inf("∫-inf^inf 1/(1+x²) dx = π",
              [](double x){ return 1.0/(1.0+x*x); }, -inf, inf, M_PI);

    return 0;
}
