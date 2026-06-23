#pragma once
#include "matrix.h"
#include <functional>

namespace pp {

using scalar_function = std::function<double(const vector&)>;

// Numerical gradient and Hessian of a scalar function phi at point x.
// central=false uses forward differences (dx ~ eps^(1/2) for gradient,
// eps^(1/4) for Hessian); central=true uses central differences
// (dx ~ eps^(1/3) for both, as derived in the assignment).
vector numerical_gradient(scalar_function phi, vector x, bool central = false);
matrix numerical_hessian(scalar_function phi, vector x, bool central = false);

// Newton's method for minimization with numerical gradient/Hessian,
// Levenberg regularization, and backtracking line-search.
// If steps != nullptr, the number of Newton iterations taken is stored there.
vector newton_minimize(
    scalar_function phi,
    vector x,
    double acc = 1e-3,
    int max_iter = 1000,
    bool central = false,
    int* steps = nullptr
);

double rosenbrock(const vector& x);
double himmelblau(const vector& x);

// Breit-Wigner resonance lineshape: F(E|m,Gamma,A) = A / [(E-m)^2 + Gamma^2/4]
double breit_wigner(double E, double m, double Gamma, double A);

}
