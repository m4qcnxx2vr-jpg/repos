#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <numbers>
#include "sfuns.h"

bool approx(double a, double b, double acc = 1e-6, double eps = 1e-6) {
	return std::abs(a - b) < acc + eps * std::abs(b);
}

int main() {
	const double pi = std::numbers::pi;
	const double e = std::numbers::e;

	const std::complex<double> i(0, 1);

	std::cout << std::setprecision(15);

	std::cout << "Simple math with <cmath>\n";
	std::cout << "sqrt(2) = " << std::sqrt(2.0) << "\n";
	std::cout << "21.0 / 5.0 = " << 21.0 / 5.0 << "\n";
	std::cout << "exp(pi) = " << std::exp(pi) << "\n";
	std::cout << "pow(pi, e) = " << std::pow(pi, e) << "\n";

	std::cout << "\nComplex math with <complex>\n";
	std::cout << "exp(i) = " << std::exp(i) << "\n";
	std::cout << "pi * i = " << pi * i << "\n";
	std::cout << "e * i = " << e * i << "\n";
	std::cout << "pow(i, i) = " << std::pow(i, i) << "\n";
	std::cout << "log(i) = " << std::log(i) << "\n";

	std::cout << "\nGamma function\n";
	std::cout << "x\tfgamma(x)\texact\t\tcheck\n";

	double factorial = 1.0;

	for (int x = 1; x <= 10; ++x) {
		if (x > 1) {
			factorial *= x - 1;
		}

		double calculated = sfuns::fgamma(x);
		double exact = factorial;

		std::cout << x << "\t"
		          << calculated << "\t"
		          << exact << "\t";

		if (approx(calculated, exact)) {
			std::cout << "OK\n";
		} else {
			std::cout << "not OK\n";
		}
	}

	std::cout << "\nLog-gamma function\n";
	std::cout << "x\tlngamma(x)\tlog(exact)\tcheck\n";

	factorial = 1.0;

	for (int x = 1; x <= 10; ++x) {
		if (x > 1) {
			factorial *= x - 1;
		}

		double calculated = sfuns::lngamma(x);
		double exact = std::log(factorial);

		std::cout << x << "\t"
		          << calculated << "\t"
		          << exact << "\t";

		if (approx(calculated, exact)) {
			std::cout << "OK\n";
		} else {
			std::cout << "not OK\n";
		}
	}

	return 0;
}