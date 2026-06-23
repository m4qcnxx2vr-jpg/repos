#pragma once
#include "matrix.h"
#include <functional>

namespace pp {

// ------------------------------------------------------------
// ann: a simple 3-layer artificial neural network
//
//   input (identity) -> n hidden neurons -> output (summation)
//
//   hidden neuron i:  y_i = w_i * f( (x - a_i) / b_i )
//   network:          F_p(x) = sum_i y_i
//
// Parameters are packed as p = [a_0,b_0,w_0, a_1,b_1,w_1, ...]
// (3*n numbers total), so the whole network is one point in R^{3n}.
//
// The activation f is supplied together with its first and second
// derivatives and its antiderivative, so that the network can return
// not just F_p(x) but also F_p'(x), F_p''(x) and an antiderivative of
// F_p (part B). All of these are analytic if f is analytic, since the
// network response is just a linear combination of (rescaled) copies
// of f.
//
// Training (part A) minimizes
//   C(p) = sum_k ( F_p(x_k) - y_k )^2
// using a Levenberg-Marquardt / Gauss-Newton scheme built on top of
// the analytic Jacobian of the residuals r_k(p) = F_p(x_k) - y_k.
// We use the *analytic* gradient/Jacobian rather than the numerical
// one from minimize.h, because the cost surface in (a_i,b_i) is
// highly non-convex and numerical differentiation is too noisy to
// reliably drive convergence (this is explicitly warned about in the
// assignment text).
// ------------------------------------------------------------
struct ann {
    int n;                              // number of hidden neurons
    std::function<double(double)> f;    // activation function
    std::function<double(double)> df;   // f'
    std::function<double(double)> d2f;  // f''
    std::function<double(double)> If;   // antiderivative of f (any fixed constant of integration)
    vector p;                           // parameters, length 3n: (a_i,b_i,w_i)

    // Build a network with n hidden neurons using the given activation
    // and its derivatives/antiderivative. Defaults to the Gaussian
    // wavelet f(x) = x*exp(-x^2), which is analytic and well suited
    // to interpolation (part A/B recommend it).
    explicit ann(int n_);
    ann(int n_,
        std::function<double(double)> f_,
        std::function<double(double)> df_,
        std::function<double(double)> d2f_,
        std::function<double(double)> If_);

    // initialize parameters spread over [xmin,xmax] with small random weights
    void init(double xmin, double xmax, unsigned seed = 1);

    // --- network evaluation -------------------------------------------------
    double response(double x) const;       // F_p(x)
    double derivative(double x) const;      // F_p'(x)
    double second_derivative(double x) const; // F_p''(x)
    double antiderivative(double x) const;  // an antiderivative of F_p, value at x

    // gradient of F_p(x) with respect to the parameters p (length 3n)
    vector grad_params(double x) const;

    // --- training (part A) --------------------------------------------------
    // Train the network so that F_p(x_k) ~= y_k for the given table.
    // Returns the final cost C(p). If steps != nullptr, stores iteration count.
    double train(const vector& x, const vector& y,
                 double acc = 1e-8, int max_iter = 2000, int* steps = nullptr);

private:
    double cost(const vector& x, const vector& y) const;
};

// Some convenient activation functions (Gaussian wavelet is the default).
double gaussian_wavelet(double x);      // f(x)  = x*exp(-x^2)
double gaussian_wavelet_d1(double x);   // f'(x) = (1-2x^2)*exp(-x^2)
double gaussian_wavelet_d2(double x);   // f''(x)= (4x^3-6x)*exp(-x^2)
double gaussian_wavelet_I(double x);    // int f = -0.5*exp(-x^2)

} // namespace pp
