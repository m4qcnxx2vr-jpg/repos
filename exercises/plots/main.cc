#include "SpecialFunctions.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

int main() {
    std::ofstream erf_data("erf.data");
    std::ofstream erf_tab("erf_tab.data");

    std::ofstream gamma_data("gamma.data");
    std::ofstream gamma_tab("gamma_tab.data");

    std::ofstream lngamma_data("lngamma.data");
    std::ofstream lngamma_tab("lngamma_tab.data");

    if (!erf_data || !erf_tab || !gamma_data || !gamma_tab || !lngamma_data || !lngamma_tab) {
        std::cerr << "Error: could not open one or more data files.\n";
        return 1;
    }

    erf_data << std::setprecision(16);
    erf_tab << std::setprecision(16);
    gamma_data << std::setprecision(16);
    gamma_tab << std::setprecision(16);
    lngamma_data << std::setprecision(16);
    lngamma_tab << std::setprecision(16);

    // Error function curve
    for (double x = -3.0; x <= 3.0; x += 0.02) {
        erf_data << x << " " << erf_approx(x) << "\n";
    }

    // Tabulated error-function values
    std::vector<double> erf_points {
        -2.0,
        -1.0,
        -0.5,
         0.0,
         0.5,
         1.0,
         2.0
    };

    for (double x : erf_points) {
        erf_tab << x << " " << std::erf(x) << "\n";
    }

    // Gamma function curve
    // The gamma function has poles at 0, -1, -2, -3, ...
    // We skip points close to the poles so gnuplot does not draw ugly connecting lines.
    for (double x = 0.05; x <= 6.0; x += 0.01) {
    double y = sgamma(x);

    if (std::isfinite(y) && std::abs(y) < 100.0) {
        gamma_data << x << " " << y << "\n";
    } else {
        gamma_data << "\n";
    }
}

    // Gamma tabulated values:
    // Gamma(n) = (n - 1)!
    for (int n = 1; n <= 7; n++) {
        double x = n;
        double y = factorial(n - 1);

        gamma_tab << x << " " << y << "\n";
    }

    // Log-gamma curve
    for (double x = 0.05; x <= 10.0; x += 0.01) {
        lngamma_data << x << " " << lngamma(x) << "\n";
    }

    // Log-gamma tabulated values:
    // ln(Gamma(n)) = ln((n - 1)!)
    for (int n = 1; n <= 10; n++) {
        double x = n;
        double y = std::log(factorial(n - 1));

        lngamma_tab << x << " " << y << "\n";
    }

    std::cout << std::setprecision(16);

    std::cout << "Error function checks\n";
    for (double x : erf_points) {
        double calculated = erf_approx(x);
        double expected = std::erf(x);

        std::cout << "x = " << x
                  << ", erf_approx(x) = " << calculated
                  << ", std::erf(x) = " << expected
                  << ", approx = " << std::boolalpha
                  << approx(calculated, expected, 1e-6)
                  << "\n";
    }

    std::cout << "\nGamma function factorial checks\n";
    for (int n = 1; n <= 7; n++) {
        double calculated = sgamma(n);
        double expected = factorial(n - 1);

        std::cout << "Gamma(" << n << ") = " << calculated
                  << ", expected " << (n - 1) << "! = " << expected
                  << ", approx = " << std::boolalpha
                  << approx(calculated, expected, 1e-8)
                  << "\n";
    }

    std::cout << "\nLog-gamma function checks\n";
    for (int n = 1; n <= 10; n++) {
        double calculated = lngamma(n);
        double expected = std::log(factorial(n - 1));

        std::cout << "ln Gamma(" << n << ") = " << calculated
                  << ", expected ln(" << (n - 1) << "!) = " << expected
                  << ", approx = " << std::boolalpha
                  << approx(calculated, expected, 1e-8)
                  << "\n";
    }

    return 0;
}