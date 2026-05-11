#include "integrate.h"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

// ──────────────────────────────────────────────────────────────────────────────
// Helper: count function evaluations via a wrapper
// ──────────────────────────────────────────────────────────────────────────────
struct CountedF {
    std::function<double(double)> f;
    mutable int ncalls = 0;
    double operator()(double x) const { ++ncalls; return f(x); }
};

// ──────────────────────────────────────────────────────────────────────────────
// Print a test result row
// ──────────────────────────────────────────────────────────────────────────────
void test(const std::string& label,
          const std::function<double(double)>& f,
          double a, double b,
          double exact,
          double acc = 1e-4, double eps = 1e-4)
{
    CountedF cf{f};
    double result = integrate(cf, a, b, acc, eps);
    double err    = std::abs(result - exact);
    bool   ok     = err < acc + eps * std::abs(exact);

    std::cout << std::left << std::setw(34) << label
              << "  result=" << std::setw(12) << std::setprecision(8) << result
              << "  exact="  << std::setw(12) << exact
              << "  |err|="  << std::scientific << std::setprecision(2) << err
              << "  ncalls=" << std::setw(5)   << cf.ncalls
              << "  " << (ok ? "OK" : "FAIL")
              << std::fixed << "\n";
}

int main()
{
    std::cout << std::string(100, '=') << "\n";
    std::cout << "  Recursive Open 4-Point Adaptive Integrator — Test Suite\n";
    std::cout << std::string(100, '=') << "\n\n";

    // ── Basic integral tests ──────────────────────────────────────────────────
    std::cout << "── Standard test integrals (acc=eps=1e-4) ──────────────────────\n";

    test("∫₀¹ √x dx  = 2/3",
         [](double x){ return std::sqrt(x); },
         0.0, 1.0, 2.0/3.0);

    test("∫₀¹ 1/√x dx  = 2",
         [](double x){ return 1.0/std::sqrt(x); },
         0.0, 1.0, 2.0);

    // ∫₀¹ √(1-x²) dx = π/4  (quarter-circle area; π/2 would be the full semicircle ∫₋₁¹)
    test("∫₀¹ √(1-x²) dx  = π/4",
         [](double x){ return std::sqrt(1.0 - x*x); },
         0.0, 1.0, M_PI/4.0);

    test("∫₀¹ ln(x)/√x dx  = -4",
         [](double x){ return std::log(x)/std::sqrt(x); },
         0.0, 1.0, -4.0);

    // A couple of extra interesting ones
    test("∫₀^π sin(x) dx  = 2",
         [](double x){ return std::sin(x); },
         0.0, M_PI, 2.0);

    test("∫₁^e ln(x) dx  = 1",
         [](double x){ return std::log(x); },
         1.0, M_E, 1.0);

    // ── erf(1) accuracy study ─────────────────────────────────────────────────
    std::cout << "\n── erf(1) accuracy vs acc (eps=0) ──────────────────────────────\n";
    const double erf1_exact = 0.84270079294971486934;
    std::cout << std::setprecision(14)
              << "  Reference erf(1) = " << erf1_exact << "\n\n";

    std::ofstream acc_file("erf_accuracy.dat");
    acc_file << "# acc  |erf(1)-exact|\n";

    std::cout << std::left
              << std::setw(12) << "acc"
              << std::setw(22) << "my_erf(1)"
              << std::setw(22) << "|error|"
              << "\n";
    std::cout << std::string(56, '-') << "\n";

    for (double acc = 0.1; acc > 1e-10; acc *= 0.1) {
        double val = my_erf(1.0, acc, 0.0);
        double err = std::abs(val - erf1_exact);
        std::cout << std::setw(12) << std::scientific << std::setprecision(1) << acc
                  << std::setw(22) << std::fixed      << std::setprecision(14) << val
                  << std::setw(22) << std::scientific << std::setprecision(4)  << err
                  << "\n";
        acc_file << acc << " " << err << "\n";
    }
    acc_file.close();

    // ── erf plot data ─────────────────────────────────────────────────────────
    std::cout << "\n── Writing erf plot data ────────────────────────────────────────\n";
    std::ofstream erf_file("erf_data.dat");
    erf_file << "# z  my_erf(z)  std::erf(z)  |diff|\n";

    for (int i = -40; i <= 40; ++i) {
        double z   = i * 0.1;
        double my  = my_erf(z, 1e-6, 1e-6);
        double ref = std::erf(z);
        erf_file << z << " " << my << " " << ref
                 << " " << std::abs(my - ref) << "\n";
    }
    erf_file.close();
    std::cout << "  erf_data.dat written (z from -4 to 4, step 0.1)\n";

    std::cout << "\n" << std::string(100, '=') << "\n";
    std::cout << "  Done.  Data files: erf_data.dat, erf_accuracy.dat\n";
    std::cout << std::string(100, '=') << "\n";
    return 0;
}