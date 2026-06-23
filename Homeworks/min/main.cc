#include "minimize.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>

int main() {
    std::cout << std::setprecision(10);

    std::cout << "============================================================\n";
    std::cout << "A. Newton minimization, numerical gradient/Hessian, forward diff\n";
    std::cout << "============================================================\n\n";

    {
        pp::vector x(2);
        x[0] = -1.0; x[1] = 1.5;
        int steps = 0;
        pp::vector xmin = pp::newton_minimize(pp::rosenbrock, x, 1e-3, 1000, false, &steps);

        std::cout << "Rosenbrock function:\n";
        std::cout << "minimum found = (" << xmin[0] << ", " << xmin[1] << ")\n";
        std::cout << "expected      = (1, 1)\n";
        std::cout << "phi(min)      = " << pp::rosenbrock(xmin) << "\n";
        std::cout << "steps taken   = " << steps << "\n\n";
    }

    {
        pp::vector x(2);
        x[0] = 1.0; x[1] = 1.0;
        int steps = 0;
        pp::vector xmin = pp::newton_minimize(pp::himmelblau, x, 1e-3, 1000, false, &steps);

        std::cout << "Himmelblau function (single start):\n";
        std::cout << "minimum found = (" << xmin[0] << ", " << xmin[1] << ")\n";
        std::cout << "phi(min)      = " << pp::himmelblau(xmin) << "\n";
        std::cout << "steps taken   = " << steps << "\n\n";

        std::cout << "Himmelblau minima from several starting points:\n";
        std::vector<pp::vector> starts;
        for (auto p : std::vector<std::pair<double,double>>{
                {3.0, 2.0}, {-3.0, 3.0}, {-4.0, -3.0}, {4.0, -2.0}}) {
            pp::vector s(2);
            s[0] = p.first; s[1] = p.second;
            starts.push_back(s);
        }
        for (auto& s : starts) {
            int st = 0;
            pp::vector m = pp::newton_minimize(pp::himmelblau, s, 1e-3, 1000, false, &st);
            std::cout << "  (" << m[0] << ", " << m[1] << ")  steps = " << st << "\n";
        }
        std::cout << "\n";
    }

    std::cout << "============================================================\n";
    std::cout << "B. Higgs boson: Breit-Wigner fit\n";
    std::cout << "============================================================\n\n";

    std::vector<double> energy, signal, error;
    {
        double e, s, d;
        while (std::cin >> e >> s >> d) {
            energy.push_back(e);
            signal.push_back(s);
            error.push_back(d);
        }
    }

    if (energy.empty()) {
        std::cout << "No Higgs data read from stdin -- run as:\n";
        std::cout << "  ./main < higgs.data.txt\n\n";
    } else {
        auto deviation = [&](const pp::vector& p) {
            double m = p[0], Gamma = p[1], A = p[2];
            double D = 0.0;
            for (std::size_t i = 0; i < energy.size(); ++i) {
                double F = pp::breit_wigner(energy[i], m, Gamma, A);
                double r = (F - signal[i]) / error[i];
                D += r * r;
            }
            return D;
        };

        pp::vector p0(3);
        p0[0] = 125.0; // m
        p0[1] = 4.0;   // Gamma
        p0[2] = 10.0;  // A

        int steps = 0;
        pp::vector p = pp::newton_minimize(deviation, p0, 1e-3, 1000, false, &steps);

        double m = p[0], Gamma = std::abs(p[1]), A = p[2];

        std::cout << "Fit results:\n";
        std::cout << "mass  m     = " << m << " GeV\n";
        std::cout << "width Gamma = " << Gamma << " GeV (upper limit on physical width)\n";
        std::cout << "scale A     = " << A << "\n";
        std::cout << "D(min)      = " << deviation(p) << "\n";
        std::cout << "steps taken = " << steps << "\n\n";

        std::ofstream out("higgs_fit.txt");
        out << std::setprecision(10);
        out << "# E  signal  error  fit\n";
        for (std::size_t i = 0; i < energy.size(); ++i) {
            out << energy[i] << " " << signal[i] << " " << error[i] << " "
                << pp::breit_wigner(energy[i], m, Gamma, A) << "\n";
        }

        std::ofstream curve("higgs_curve.txt");
        curve << std::setprecision(10);
        curve << "# E  fit_curve\n";
        double emin = energy.front(), emax = energy.back();
        int N = 400;
        for (int i = 0; i <= N; ++i) {
            double E = emin + (emax - emin) * i / N;
            curve << E << " " << pp::breit_wigner(E, m, Gamma, A) << "\n";
        }
    }

    std::cout << "============================================================\n";
    std::cout << "C. Central vs forward finite differences\n";
    std::cout << "============================================================\n\n";

    {
        pp::vector x(2);
        x[0] = -1.0; x[1] = 1.5;

        int steps_fwd = 0;
        pp::vector xf = pp::newton_minimize(pp::rosenbrock, x, 1e-3, 1000, false, &steps_fwd);

        int steps_ctr = 0;
        pp::vector xc = pp::newton_minimize(pp::rosenbrock, x, 1e-3, 1000, true, &steps_ctr);

        std::cout << "Rosenbrock function:\n";
        std::cout << "forward diff : min = (" << xf[0] << ", " << xf[1]
                  << ")  steps = " << steps_fwd
                  << "  phi = " << pp::rosenbrock(xf) << "\n";
        std::cout << "central diff : min = (" << xc[0] << ", " << xc[1]
                  << ")  steps = " << steps_ctr
                  << "  phi = " << pp::rosenbrock(xc) << "\n\n";
    }

    {
        pp::vector x(2);
        x[0] = 1.0; x[1] = 1.0;

        int steps_fwd = 0;
        pp::vector xf = pp::newton_minimize(pp::himmelblau, x, 1e-3, 1000, false, &steps_fwd);

        int steps_ctr = 0;
        pp::vector xc = pp::newton_minimize(pp::himmelblau, x, 1e-3, 1000, true, &steps_ctr);

        std::cout << "Himmelblau function:\n";
        std::cout << "forward diff : min = (" << xf[0] << ", " << xf[1]
                  << ")  steps = " << steps_fwd
                  << "  phi = " << pp::himmelblau(xf) << "\n";
        std::cout << "central diff : min = (" << xc[0] << ", " << xc[1]
                  << ")  steps = " << steps_ctr
                  << "  phi = " << pp::himmelblau(xc) << "\n\n";
    }

  

    return 0;
}
