#include "mc.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <string>

void line()
{
    std::cout << "------------------------------------------------------------\n";
}

void print_result(
    const std::string& name,
    double value,
    double error,
    double exact
)
{
    std::cout << std::left << std::setw(18) << name
              << " estimate = " << std::setw(16) << value
              << " estimated error = " << std::setw(16) << error
              << " actual error = " << std::abs(value - exact)
              << "\n";
}

int main()
{
    std::cout << std::setprecision(12);

    const double pi = std::numbers::pi;

    pp::vector square_a = {-1.0, -1.0};
    pp::vector square_b = { 1.0,  1.0};

    pp::Func unit_circle = [](pp::vector x) {
        double r2 = x[0] * x[0] + x[1] * x[1];
        return r2 <= 1.0 ? 1.0 : 0.0;
    };

    line();
    std::cout << "Task A: Plain Monte Carlo integration\n";
    line();

    std::cout << "\nA1: Area of the unit circle\n";
    std::cout << "Exact value: pi = " << pi << "\n\n";

    for (int N : {1000, 3000, 10000, 30000, 100000}) {
        auto [value, error] = pp::plainmc(unit_circle, square_a, square_b, N);

        std::cout << "N = " << std::setw(8) << N
                  << " estimate = " << std::setw(16) << value
                  << " estimated error = " << std::setw(16) << error
                  << " actual error = " << std::abs(value - pi)
                  << "\n";
    }

    std::ofstream circle_data("circle_errors.txt");
    circle_data << "# N estimated_error actual_error actual_error_times_sqrtN\n";

    for (int N = 10; N <= 10000; N *= 2) {
        auto [value, error] = pp::plainmc(unit_circle, square_a, square_b, N);
        double actual_error = std::abs(value - pi);

        circle_data << N << " "
                    << error << " "
                    << actual_error << " "
                    << actual_error * std::sqrt(static_cast<double>(N))
                    << "\n";
    }

    std::cout << "\nA2: Error scaling\n";
    std::cout << "See plot A2.png for a plot of the estimated and actual errors as a function of N.\n";

    pp::vector ellipsoid_a = {-1.0, -2.0, -3.0};
    pp::vector ellipsoid_b = { 1.0,  2.0,  3.0};

    pp::Func ellipsoid = [](pp::vector x) {
        double test =
            x[0] * x[0] / 1.0 +
            x[1] * x[1] / 4.0 +
            x[2] * x[2] / 9.0;

        return test <= 1.0 ? 1.0 : 0.0;
    };

    double exact_ellipsoid = 4.0 * pi * 1.0 * 2.0 * 3.0 / 3.0;

    auto [ell_value, ell_error] = pp::plainmc(
        ellipsoid,
        ellipsoid_a,
        ellipsoid_b,
        200000
    );

    std::cout << "\nA3: Volume of ellipsoid with axes a=1, b=2, c=3\n";
    std::cout << "Exact value: " << exact_ellipsoid << "\n";
    print_result("Plain MC", ell_value, ell_error, exact_ellipsoid);

    line();
    std::cout << "Task B: Quasi-random Monte Carlo\n";
    line();

    pp::vector cube_a = {0.0, 0.0, 0.0};
    pp::vector cube_b = {pi, pi, pi};

    pp::Func singular = [pi](pp::vector x) {
        double denominator =
            1.0 - std::cos(x[0]) * std::cos(x[1]) * std::cos(x[2]);

        return 1.0 / (pi * pi * pi * denominator);
    };

    double exact_singular =
        std::pow(std::tgamma(0.25), 4.0) / (4.0 * std::pow(pi, 3.0));

    std::cout << "\nB1: Difficult singular integral\n";
    std::cout << "Exact value: " << exact_singular << "\n\n";

    int N_singular = 200000;

    auto [lcg_value, lcg_error] = pp::lcgmc(
        singular,
        cube_a,
        cube_b,
        N_singular
    );

    auto [std_value, std_error] = pp::plainmc(
        singular,
        cube_a,
        cube_b,
        N_singular
    );

    auto [q_value, q_error] = pp::quasimc(
        singular,
        cube_a,
        cube_b,
        N_singular
    );

    print_result("LCG MC", lcg_value, lcg_error, exact_singular);
    print_result("std::mt19937 MC", std_value, std_error, exact_singular);
    print_result("Halton quasi MC", q_value, q_error, exact_singular);

    std::cout << "\nB2: Comment\n";
    std::cout << "The quasi-random error estimate is obtained by comparing two Halton estimates.\n";
    std::cout << "The LCG and std::mt19937 estimates use pseudo-random numbers.\n";

    line();
    std::cout << "Task C: Recursive stratified sampling\n";
    line();

    std::cout << "\nC1: Stratified sampling test on the unit circle\n";
    std::cout << "Exact value: pi = " << pi << "\n\n";

    auto [strat_value, strat_error] = pp::stratifiedmc(
        unit_circle,
        square_a,
        square_b,
        100000
    );

    print_result("Stratified MC", strat_value, strat_error, pi);

    std::cout << "\nDone.\n";

    return 0;
}