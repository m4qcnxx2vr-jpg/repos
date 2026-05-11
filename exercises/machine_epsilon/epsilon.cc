#include "epsilon.h"

float machine_epsilon_float() {
    float f = 1.0f;

    while (1.0f + f != 1.0f) {
        f /= 2.0f;
    }

    return 2.0f * f;
}

double machine_epsilon_double() {
    double d = 1.0;

    while (1.0 + d != 1.0) {
        d /= 2.0;
    }

    return 2.0 * d;
}

long double machine_epsilon_long_double() {
    long double l = 1.0L;

    while (1.0L + l != 1.0L) {
        l /= 2.0L;
    }

    return 2.0L * l;
}