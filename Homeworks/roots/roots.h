#pragma once
#include "matrix.h"
#include <functional>
#include <vector>

namespace pp {

using vector_function = std::function<vector(const vector&)>;

vector newton(
    vector_function f,
    vector x,
    double acc = 1e-6,
    double alpha_min = 1e-3,
    int max_iter = 100
);

vector newton_quadratic(
    vector_function f,
    vector x,
    double acc = 1e-6,
    double alpha_min = 1.0 / 128.0,
    int max_iter = 100
);

vector gradient_rosenbrock(const vector& x);
vector gradient_himmelblau(const vector& x);

double hydrogen_mismatch(double E, double rmin, double rmax, double acc, double eps);
double find_hydrogen_energy(double rmin, double rmax, double acc, double eps);

// Integrates the radial equation at energy E and records (r, f(r), f'(r))
// at every accepted step, so the wavefunction can be plotted.
struct hydrogen_point { double r, f, fp; };
std::vector<hydrogen_point> hydrogen_wavefunction(
    double E, double rmin, double rmax, double acc, double eps
);

}