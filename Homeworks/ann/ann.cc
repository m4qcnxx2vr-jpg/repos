#include "ann.h"
#include "qr.h"
#include <cmath>
#include <random>
#include <algorithm>

namespace pp {

// ---------------- default activation: Gaussian wavelet ----------------
// f(x)  = x * exp(-x^2)
// f'(x) = (1 - 2x^2) * exp(-x^2)
// f''(x)= (4x^3 - 6x) * exp(-x^2)
// I(x)  = -0.5 * exp(-x^2)      (antiderivative, constant of integration = 0)
double gaussian_wavelet(double x) {
    return x * std::exp(-x * x);
}
double gaussian_wavelet_d1(double x) {
    return (1.0 - 2.0 * x * x) * std::exp(-x * x);
}
double gaussian_wavelet_d2(double x) {
    return (4.0 * x * x * x - 6.0 * x) * std::exp(-x * x);
}
double gaussian_wavelet_I(double x) {
    return -0.5 * std::exp(-x * x);
}

// ---------------- constructors ----------------
ann::ann(int n_)
    : n(n_),
      f(gaussian_wavelet),
      df(gaussian_wavelet_d1),
      d2f(gaussian_wavelet_d2),
      If(gaussian_wavelet_I),
      p(3 * n_, 0.0)
{}

ann::ann(int n_,
         std::function<double(double)> f_,
         std::function<double(double)> df_,
         std::function<double(double)> d2f_,
         std::function<double(double)> If_)
    : n(n_), f(f_), df(df_), d2f(d2f_), If(If_), p(3 * n_, 0.0)
{}

void ann::init(double xmin, double xmax, unsigned seed) {
    std::mt19937 rng(seed);
    double span = std::max(xmax - xmin, 1e-3);

    // Spread centers evenly over the domain (with a touch of jitter)
    // rather than fully at random: gives every neuron a sensible,
    // non-overlapping starting region to specialize in.
    std::uniform_real_distribution<double> jitter(-0.4, 0.4);
    std::uniform_real_distribution<double> uw(-0.5, 0.5);

    double b0 = span / n; // matches the b_target scale used by train()

    for (int i = 0; i < n; i++) {
        double center = xmin + span * (i + 0.5) / n;
        p[3 * i + 0] = center + jitter(rng) * (span / n); // a_i
        p[3 * i + 1] = b0;                                // b_i
        p[3 * i + 2] = uw(rng);                            // w_i
    }
}

// ---------------- evaluation ----------------
double ann::response(double x) const {
    double F = 0.0;
    for (int i = 0; i < n; i++) {
        double a = p[3 * i + 0], b = p[3 * i + 1], w = p[3 * i + 2];
        double u = (x - a) / b;
        F += w * f(u);
    }
    return F;
}

double ann::derivative(double x) const {
    double F = 0.0;
    for (int i = 0; i < n; i++) {
        double a = p[3 * i + 0], b = p[3 * i + 1], w = p[3 * i + 2];
        double u = (x - a) / b;
        F += w * df(u) / b;
    }
    return F;
}

double ann::second_derivative(double x) const {
    double F = 0.0;
    for (int i = 0; i < n; i++) {
        double a = p[3 * i + 0], b = p[3 * i + 1], w = p[3 * i + 2];
        double u = (x - a) / b;
        F += w * d2f(u) / (b * b);
    }
    return F;
}

double ann::antiderivative(double x) const {
    // d/dx [ w*b*If(u) ] = w*b*If'(u)*(1/b) = w*f(u). Good.
    double F = 0.0;
    for (int i = 0; i < n; i++) {
        double a = p[3 * i + 0], b = p[3 * i + 1], w = p[3 * i + 2];
        double u = (x - a) / b;
        F += w * b * If(u);
    }
    return F;
}

vector ann::grad_params(double x) const {
    vector g(3 * n, 0.0);
    for (int i = 0; i < n; i++) {
        double a = p[3 * i + 0], b = p[3 * i + 1], w = p[3 * i + 2];
        double u = (x - a) / b;
        double fu = f(u);
        double dfu = df(u);

        g[3 * i + 0] = -w * dfu / b;          // dF/da_i
        g[3 * i + 1] = -w * dfu * u / b;      // dF/db_i
        g[3 * i + 2] = fu;                    // dF/dw_i
    }
    return g;
}

double ann::cost(const vector& x, const vector& y) const {
    double C = 0.0;
    for (int k = 0; k < x.size(); k++) {
        double r = response(x[k]) - y[k];
        C += r * r;
    }
    return C;
}

// ---------------- training: Levenberg-Marquardt on the residuals -------
// r_k(p) = F_p(x_k) - y_k,  J_{k,j} = d r_k / d p_j = d F_p(x_k) / d p_j
//
// LM step solves  (J^T J + mu*diag(J^T J)) dp = -J^T r
// We build J^T J + mu*diag and -J^T r directly and solve with the qr
// class (this matrix is m x m, symmetric positive definite for mu>0).
//
// See the regularization comment below for why two extra penalty terms
// are added to the plain sum-of-squares cost.
double ann::train(const vector& x, const vector& y,
                   double acc, int max_iter, int* steps) {
    int N = x.size();
    int m = 3 * n;

    double xmin_ = x[0], xmax_ = x[0];
    for (int k = 1; k < N; k++) {
        xmin_ = std::min(xmin_, x[k]);
        xmax_ = std::max(xmax_, x[k]);
    }
    double xspan = std::max(xmax_ - xmin_, 1.0);

    // Width regularization: a neuron is penalized smoothly as its width
    // b_i strays from a "sensible" scale comparable to the spacing
    // between neurons across the domain. Without this, the optimizer
    // can collapse a neuron onto a single point (b_i -> 0, large
    // compensating w_i -- a spike that fits one data point and ruins
    // everything in between) or spread it out far beyond the domain.
    // Penalizing log(b_i / b_target)^2 discourages both extremes
    // symmetrically in log-width.
    //
    // Weight regularization: penalizing w_i directly (a standard L2 /
    // ridge term) is also needed. Width control alone can be evaded by
    // a *pair* of narrow, opposite-sign-weighted neurons sitting next
    // to each other: individually each width is only moderately small
    // and escapes a strong width penalty, but their large cancelling
    // weights still reconstruct a sharp, non-smooth spike between two
    // data points. Penalizing |w_i| directly removes the incentive for
    // such large cancellations regardless of what the widths are doing.
    double b_target = xspan / (2.0 * n);
    double lambda_b = 0.01; // strength of the log-width penalty
    double lambda_w = 0.001; // strength of the weight-magnitude penalty
    // scale for the weight penalty: typical |y| values, so the penalty
    // is comparable in size to a typical data residual
    double yscale = 0.0;
    for (int k = 0; k < N; k++) yscale = std::max(yscale, std::fabs(y[k]));
    yscale = std::max(yscale, 1e-3);

    int Ntot = N + 2 * n; // data residuals + width-reg + weight-reg residuals

    double mu = 1e-3;
    double C = cost(x, y);
    for (int i = 0; i < n; i++) {
        double rho_b = std::sqrt(lambda_b) * std::log(p[3*i+1] / b_target);
        double rho_w = std::sqrt(lambda_w) * p[3*i+2] / yscale;
        C += rho_b * rho_b + rho_w * rho_w;
    }
    int iter = 0;

    for (; iter < max_iter; iter++) {
        matrix J(Ntot, m, 0.0);
        vector r(Ntot, 0.0);

        for (int k = 0; k < N; k++) {
            double xk = x[k];
            r[k] = response(xk) - y[k];
            vector gk = grad_params(xk);
            for (int j = 0; j < m; j++) J(k, j) = gk[j];
        }
        // width-regularization rows: rho_i = sqrt(lambda_b)*log(b_i/b_target)
        for (int i = 0; i < n; i++) {
            double bi = p[3 * i + 1];
            r[N + i] = std::sqrt(lambda_b) * std::log(bi / b_target);
            J(N + i, 3 * i + 1) = std::sqrt(lambda_b) / bi;
        }
        // weight-regularization rows: rho_i = sqrt(lambda_w)*w_i/yscale
        for (int i = 0; i < n; i++) {
            r[N + n + i] = std::sqrt(lambda_w) * p[3 * i + 2] / yscale;
            J(N + n + i, 3 * i + 2) = std::sqrt(lambda_w) / yscale;
        }

        vector JTr(m, 0.0);
        for (int j = 0; j < m; j++) {
            double s = 0.0;
            for (int k = 0; k < Ntot; k++) s += J(k, j) * r[k];
            JTr[j] = s;
        }

        if (JTr.norm() < acc) break;

        matrix A(m, m, 0.0);
        for (int i = 0; i < m; i++) {
            for (int j = i; j < m; j++) {
                double s = 0.0;
                for (int k = 0; k < Ntot; k++) s += J(k, i) * J(k, j);
                A(i, j) = s;
                A(j, i) = s;
            }
        }

        auto total_cost = [&](const vector& pp_) {
            ann trial = *this;
            trial.p = pp_;
            double c = trial.cost(x, y);
            for (int i = 0; i < n; i++) {
                double bi = pp_[3 * i + 1];
                if (bi <= 0.0) return 1e300; // never accept non-positive width
                double rho_b = std::sqrt(lambda_b) * std::log(bi / b_target);
                double rho_w = std::sqrt(lambda_w) * pp_[3 * i + 2] / yscale;
                c += rho_b * rho_b + rho_w * rho_w;
            }
            return c;
        };

        bool improved = false;
        vector pnew = p;

        for (int tries = 0; tries < 50 && !improved; tries++) {
            matrix M = A;
            for (int i = 0; i < m; i++) M(i, i) += mu * (A(i, i) + 1e-9);

            vector dp;
            try {
                qr QRM(M);
                dp = QRM.solve(-1.0 * JTr);
            } catch (...) {
                mu *= 10.0;
                continue;
            }

            vector ptry = p;
            ptry += dp;

            double Cnew = total_cost(ptry);
            if (Cnew < C) {
                pnew = ptry;
                C = Cnew;
                improved = true;
                mu = std::max(mu / 7.0, 1e-12);
                break;
            }
            mu *= 10.0;
        }

        if (!improved) break; // stuck: stop iterating
        p = pnew;
    }

    if (steps) *steps = iter;
    return cost(x, y); // report the actual data-fit cost, not the regularized one
}

} // namespace pp
