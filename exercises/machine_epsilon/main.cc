#include "epsilon.h"
#include "compare.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>

// This is AI generated, i made the other files myself. I gave a rough description of what i wanted and it generated the code.
int main() {
    std::cout << std::scientific;
    std::cout << std::setprecision(20);

    std::cout << "Machine epsilon\n";
    std::cout << "===============\n\n";

    float f_eps = machine_epsilon_float();
    double d_eps = machine_epsilon_double();
    long double l_eps = machine_epsilon_long_double();

    std::cout << "Calculated epsilons:\n";
    std::cout << "float       eps = " << f_eps << "\n";
    std::cout << "double      eps = " << d_eps << "\n";
    std::cout << "long double eps = " << l_eps << "\n\n";

    std::cout << "System epsilons from std::numeric_limits:\n";
    std::cout << "float       eps = " << std::numeric_limits<float>::epsilon() << "\n";
    std::cout << "double      eps = " << std::numeric_limits<double>::epsilon() << "\n";
    std::cout << "long double eps = " << std::numeric_limits<long double>::epsilon() << "\n\n";

    std::cout << "Theoretical checks:\n";
    std::cout << "2^(-23) = " << std::pow(2.0, -23) << "   for float\n";
    std::cout << "2^(-52) = " << std::pow(2.0, -52) << "   for double\n\n";

    std::cout << "Compare calculated values with system values using approx:\n";
    std::cout << "float  correct? " 
              << (approx(f_eps, std::numeric_limits<float>::epsilon()) ? "true" : "false") 
              << "\n";

    std::cout << "double correct? " 
              << (approx(d_eps, std::numeric_limits<double>::epsilon()) ? "true" : "false") 
              << "\n\n";

    std::cout << "Floating-point order of addition\n";
    std::cout << "================================\n\n";

    double epsilon = std::pow(2.0, -52);
    double tiny = epsilon / 2.0;

    double a = 1.0 + tiny + tiny;
    double b = tiny + tiny + 1.0;

    std::cout << "a == b ? " << (a == b ? "true" : "false") << "\n";
    std::cout << "a > 1  ? " << (a > 1.0 ? "true" : "false") << "\n";
    std::cout << "b > 1  ? " << (b > 1.0 ? "true" : "false") << "\n\n";

    std::cout << std::fixed << std::setprecision(17);
    std::cout << "       tiny = " << tiny << "\n";
    std::cout << "1+tiny+tiny = " << a << "\n";
    std::cout << "tiny+tiny+1 = " << b << "\n\n";

    std::cout << "Comparing doubles\n";
    std::cout << "=================\n\n";

    double d1 = 0.1 + 0.1 + 0.1 + 0.1 + 0.1 + 0.1 + 0.1 + 0.1;
    double d2 = 8.0 * 0.1;

    std::cout << "d1 == d2 ? " << (d1 == d2 ? "true" : "false") << "\n";
    std::cout << "approx(d1, d2) ? " << (approx(d1, d2) ? "true" : "false") << "\n\n";

    std::cout << "d1 = " << d1 << "\n";
    std::cout << "d2 = " << d2 << "\n";

    return 0;
}