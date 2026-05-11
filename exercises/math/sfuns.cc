#include <cmath>
#include <limits>
#include <numbers>
#include "sfuns.h"

namespace sfuns {

constexpr double PI = std::numbers::pi;
constexpr double NaN = std::numeric_limits<double>::quiet_NaN();

double fgamma(double x) {
	if (x < 0) {
		return PI / std::sin(PI * x) / fgamma(1 - x);
	}

	if (x < 9) {
		return fgamma(x + 1) / x;
	}

	double lnfgamma =
		x * std::log(x + 1 / (12 * x - 1 / (10 * x)))
		- x
		+ std::log(2 * PI / x) / 2;

	return std::exp(lnfgamma);
}

double lngamma(double x) {
	if (x <= 0) {
		return NaN;
	}

	if (x < 9) {
		return lngamma(x + 1) - std::log(x);
	}

	double result =
		x * std::log(x + 1 / (12 * x - 1 / (10 * x)))
		- x
		+ std::log(2 * PI / x) / 2;

	return result;
}

}