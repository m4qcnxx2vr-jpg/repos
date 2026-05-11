#include "compare.h"

#include <algorithm>
#include <cmath>

bool approx(double a, double b, double acc, double eps) {
    double diff = std::abs(a - b);

    if (diff <= acc) {
        return true;
    }

    double max_ab = std::max(std::abs(a), std::abs(b));

    return diff <= eps * max_ab;
}