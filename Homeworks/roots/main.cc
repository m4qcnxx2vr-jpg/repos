#include "roots.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <fstream>

int main() {
    std::cout << std::setprecision(12);

    std::cout << "============================================================\n";
    std::cout << "A. Newton's method with numerical Jacobian and line-search\n";
    std::cout << "============================================================\n\n";

    auto f1 = [](const pp::vector& x) {
        pp::vector y(1);
        y[0] = x[0] * x[0] - 2.0;
        return y;
    };

    pp::vector x1(1);
    x1[0] = 1.0;

    pp::vector r1 = pp::newton(f1, x1);

    std::cout << "Test 1: root of x^2 - 2 = 0\n";
    std::cout << "numerical root = " << r1[0] << "\n";
    std::cout << "exact root     = " << std::sqrt(2.0) << "\n";
    std::cout << "error          = " << std::abs(r1[0] - std::sqrt(2.0)) << "\n\n";

    pp::vector xr(2);
    xr[0] = 0.5;
    xr[1] = 0.5;

    pp::vector rr = pp::newton(pp::gradient_rosenbrock, xr, 1e-6, 1e-3,1000);

    std::cout << "Rosenbrock valley function\n";
    std::cout << "minimum found  = (" << rr[0] << ", " << rr[1] << ")\n";
    std::cout << "expected       = (1, 1)\n\n";

    std::cout << "Himmelblau function\n";
    std::cout << "minima found:\n";

    std::vector<pp::vector> starts;

    pp::vector h1(2);
    h1[0] = 3.0;
    h1[1] = 2.0;
    starts.push_back(h1);

    pp::vector h2(2);
    h2[0] = -3.0;
    h2[1] = 3.0;
    starts.push_back(h2);

    pp::vector h3(2);
    h3[0] = -4.0;
    h3[1] = -3.0;
    starts.push_back(h3);

    pp::vector h4(2);
    h4[0] = 4.0;
    h4[1] = -2.0;
    starts.push_back(h4);

    for (const pp::vector& s : starts) {
        pp::vector rh = pp::newton(pp::gradient_himmelblau, s);
        std::cout << "(" << rh[0] << ", " << rh[1] << ")\n";
    }

    std::cout << "\n";
    std::cout << "Known Himmelblau minima approximately:\n";
    std::cout << "(3, 2)\n";
    std::cout << "(-2.805118, 3.131312)\n";
    std::cout << "(-3.779310, -3.283186)\n";
    std::cout << "(3.584428, -1.848126)\n\n";


    std::cout << "============================================================\n";
    std::cout << "B. Hydrogen atom bound state by shooting method\n";
    std::cout << "============================================================\n\n";

    double rmin = 1e-3;
    double rmax = 8.0;
    double acc = 1e-6;
    double eps = 1e-6;

    double E0 = pp::find_hydrogen_energy(rmin, rmax, acc, eps);

    std::cout << "Parameters:\n";
    std::cout << "rmin = " << rmin << "\n";
    std::cout << "rmax = " << rmax << "\n";
    std::cout << "acc  = " << acc << "\n";
    std::cout << "eps  = " << eps << "\n\n";

    std::cout << "Lowest energy root:\n";
    std::cout << "numerical E0 = " << E0 << "\n";
    std::cout << "exact E0     = -0.5\n";
    std::cout << "error        = " << std::abs(E0 + 0.5) << "\n\n";

    // Dump the wavefunction at E0 for plotting against the exact solution
    {
        auto pts = pp::hydrogen_wavefunction(E0, rmin, rmax, acc, eps);
        std::ofstream out("wavefunction.txt");
        out << std::setprecision(12);
        out << "# r  f_numerical(r)  f_exact(r)=r*exp(-r)\n";
        for (const auto& p : pts) {
            double f_exact = p.r * std::exp(-p.r);
            out << p.r << " " << p.f << " " << f_exact << "\n";
        }
    }

    std::cout << "Convergence study (writing data files for plotting):\n\n";

    // --- vs rmax (rmin, acc, eps fixed) ---
    {
        std::ofstream out("conv_rmax.txt");
        out << std::setprecision(12);
        out << "# rmax  E0  error=|E0+0.5|\n";
        std::cout << "vs rmax:\n";
        std::vector<double> rmax_vals = {2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 15.0};
        for (double rmx : rmax_vals) {
            double e = pp::find_hydrogen_energy(rmin, rmx, acc, eps);
            double err = std::abs(e + 0.5);
            out << rmx << " " << e << " " << err << "\n";
            std::cout << "  rmax = " << rmx << "  E0 = " << e << "  err = " << err << "\n";
        }
        std::cout << "\n";
    }

    // --- vs rmin (rmax, acc, eps fixed) ---
    {
        std::ofstream out("conv_rmin.txt");
        out << std::setprecision(12);
        out << "# rmin  E0  error=|E0+0.5|\n";
        std::cout << "vs rmin:\n";
        std::vector<double> rmin_vals = {1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6};
        for (double rmn : rmin_vals) {
            double e = pp::find_hydrogen_energy(rmn, rmax, acc, eps);
            double err = std::abs(e + 0.5);
            out << rmn << " " << e << " " << err << "\n";
            std::cout << "  rmin = " << rmn << "  E0 = " << e << "  err = " << err << "\n";
        }
        std::cout << "\n";
    }

    // --- vs acc (rmin, rmax, eps fixed) ---
    {
        std::ofstream out("conv_acc.txt");
        out << std::setprecision(12);
        out << "# acc  E0  error=|E0+0.5|\n";
        std::cout << "vs acc:\n";
        std::vector<double> acc_vals = {1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8};
        for (double a : acc_vals) {
            double e = pp::find_hydrogen_energy(rmin, rmax, a, eps);
            double err = std::abs(e + 0.5);
            out << a << " " << e << " " << err << "\n";
            std::cout << "  acc = " << a << "  E0 = " << e << "  err = " << err << "\n";
        }
        std::cout << "\n";
    }

    // --- vs eps (rmin, rmax, acc fixed) ---
    {
        std::ofstream out("conv_eps.txt");
        out << std::setprecision(12);
        out << "# eps  E0  error=|E0+0.5|\n";
        std::cout << "vs eps:\n";
        std::vector<double> eps_vals = {1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8};
        for (double e_ : eps_vals) {
            double e = pp::find_hydrogen_energy(rmin, rmax, acc, e_);
            double err = std::abs(e + 0.5);
            out << e_ << " " << e << " " << err << "\n";
            std::cout << "  eps = " << e_ << "  E0 = " << e << "  err = " << err << "\n";
        }
        std::cout << "\n";
    }


    std::cout << "============================================================\n";
    std::cout << "C. Quadratic interpolation line-search / optimization\n";
    std::cout << "============================================================\n\n";

    std::cout << "Testing quadratic interpolation line-search on Rosenbrock:\n";

    pp::vector xc(2);
    xc[0] = 0.5;
    xc[1] = 0.5;

    pp::vector rc = pp::newton_quadratic(
        pp::gradient_rosenbrock,
        xc,
        1e-6,
        1.0 / 128.0,
        1000
    );

    std::cout << "minimum found  = (" << rc[0] << ", " << rc[1] << ")\n";
    std::cout << "expected       = (1, 1)\n\n";

    std::cout << "Part C implementation:\n";
    std::cout << "1. The Jacobian matrix is allocated once before the Newton loop.\n";
    std::cout << "2. The line-search uses quadratic interpolation instead of simple halving.\n";

    return 0;
}